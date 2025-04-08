import numpy as np
import networkx as nx
import scanpy as sc
import scipy.sparse as sp
from scipy.spatial import Delaunay
from sklearn.neighbors import kneighbors_graph, radius_neighbors_graph
import scipy.sparse as sp
from scipy.sparse import csr_matrix
import community as community_louvain
import igraph as ig
import leidenalg
from sklearn.cluster import KMeans
from scipy.linalg import eigh

def generate_spatial_graph(adata, method="knn", k=6, radius=50, obsp_key="spatial_connectivities"):
    """
    Generate a spatial graph from an AnnData object and store it in `adata.obsp`.

    Parameters:
    - adata (AnnData): The AnnData object containing spatial transcriptomics data.
    - method (str): Graph construction method. Options: "knn", "delaunay", "radius".
    - k (int): Number of neighbors for kNN graph (only used if method="knn").
    - radius (float): Radius cutoff for radius-based graph (only used if method="radius").
    - obsp_key (str): Key to store the adjacency matrix in `adata.obsp`.

    Returns:
    - Updated AnnData object with the spatial graph stored in `adata.obsp[obsp_key]`.
    """
    if 'spatial' not in adata.obsm:
        raise ValueError("AnnData object must contain 'spatial' coordinates in adata.obsm['spatial'].")

    coords = adata.obsm['spatial']

    if coords.shape[1] != 2:
        raise ValueError("Spatial coordinates must have shape (n, 2) for x, y spatial locations.")

    num_cells = coords.shape[0]

    if method == "knn":
        # k-Nearest Neighbors Graph
        graph_matrix = kneighbors_graph(coords, k, mode="connectivity", include_self=False)

    elif method == "delaunay":
        # Delaunay Triangulation Graph
        tri = Delaunay(coords)
        edges = np.vstack([tri.simplices[:, [0, 1]],
                           tri.simplices[:, [1, 2]],
                           tri.simplices[:, [2, 0]]])
        graph_matrix = np.zeros((num_cells, num_cells))
        for edge in edges:
            graph_matrix[edge[0], edge[1]] = 1
            graph_matrix[edge[1], edge[0]] = 1  # Undirected graph
        graph_matrix = csr_matrix(graph_matrix)  # Convert to sparse matrix

    elif method == "radius":
        # Radius-Based Graph
        graph_matrix = radius_neighbors_graph(coords, radius, mode="connectivity", include_self=False)

    else:
        raise ValueError("Invalid method. Choose from 'knn', 'delaunay', or 'radius'.")

    # Store the spatial graph in AnnData's obsp slot
    adata.obsp[obsp_key] = graph_matrix

    return adata

def load_spatial_graph(adata):
    """
    Load the spatial graph from a file or memory.
    Ensure it is structured as an undirected weighted graph.
    """
    # load file as network
    spatial_graph = nx.from_numpy_array(adata.obsp['spatial_connectivities'], create_using=nx.Graph())
    # Check node/edge integrity; ensures node indices = row indices of adata.obsp
    assert not nx.is_directed(spatial_graph), "Graph is directed! Convert to undirected using nx.Graph()."
    # ensure edges weighted by spatial connectivity
    adj_matrix = sp.coo_matrix(adata.obsp['spatial_connectivities'])
    for i, j, w in zip(adj_matrix.row, adj_matrix.col, adj_matrix.data):
        spatial_graph[i][j]['weight'] = w
    # ensure graph is weighted
    for u, v, data in spatial_graph.edges(data=True):
        assert 'weight' in data, f"Edge ({u}, {v}) is missing weight!"

    return spatial_graph

def cluster_spatial_graph(graph, method='leiden'):

    """
    Perform graph-based clustering using the specified method.
    Supported methods: 'louvain', 'leiden', 'spectral'.
    """

    if method == "louvain":
        cluster_map = community_louvain.best_partition(graph)
    elif method == "leiden":
        # Convert networkx to igraph
        node_list = list(graph.nodes)
        node_index = {node: idx for idx, node in enumerate(node_list)}
        ig_graph = ig.Graph()
        ig_graph.add_vertices(len(node_list))
        ig_graph.add_edges([(node_index[u], node_index[v]) for u, v in graph.edges])

        # Perform Leiden clustering
        partition = leidenalg.find_partition(ig_graph, leidenalg.ModularityVertexPartition)

        # Map back cluster assignments
        cluster_map = {node_list[i]: cluster_id for i, cluster_id in enumerate(partition)}
    elif method == "spectral":
        # compute Laplacian matrix
        A = nx.adjacency_matrix(graph).toarray()
        D = np.diag(A.sum(axis=1))
        L = D - A
        # get eigenvectors
        _, eigvecs = eigh(L, subset_by_index=[1, 3]) # filler for now: ask Shu what he would want num_clusters to be
        # apply k-means clustering
        kmeans = KMeans(n_clusters=3, random_state=42, n_init=10)
        labels = kmeans.fit_predict(eigvecs)
        cluster_map = {node: labels[i] for i, node in enumerate(graph.nodes())}
    else:
        raise ValueError("Unsupported clustering method")

    # Assign cluster labels to graph nodes
    for node, cluster in cluster_map.items():
        graph.nodes[node]["cluster"] = cluster

    return graph
