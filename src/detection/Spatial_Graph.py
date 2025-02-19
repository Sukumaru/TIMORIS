import numpy as np
import networkx as nx
import scanpy as sc
from scipy.spatial import Delaunay
from sklearn.neighbors import kneighbors_graph, radius_neighbors_graph
from scipy.sparse import csr_matrix

def generate_spatial_graph(adata, method="knn", k=6, radius=50, obsp_key="spatial_graph"):
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
