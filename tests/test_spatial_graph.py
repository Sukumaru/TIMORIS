import pytest
import numpy as np
import anndata
from timoris.spatial_graph import generate_spatial_graph, load_spatial_graph, cluster_spatial_graph
import networkx as nx
import scanpy as sc
import scipy.sparse as sp
from tempfile import NamedTemporaryFile


@pytest.fixture
def synthetic_adata():
    """Generate a synthetic AnnData object with spatial coordinates for testing."""
    np.random.seed(42)
    coords = np.random.rand(50, 2) * 100  # 50 random cells
    adata = anndata.AnnData(X=np.random.rand(50, 10))  # 50 cells, 10 gene features
    adata.obsm['spatial'] = coords
    return adata

def test_knn_graph(synthetic_adata):
    adata = generate_spatial_graph(synthetic_adata, method="knn", k=6, obsp_key="knn_graph")
    assert "knn_graph" in adata.obsp
    assert adata.obsp["knn_graph"].shape == (50, 50)

def test_delaunay_graph(synthetic_adata):
    adata = generate_spatial_graph(synthetic_adata, method="delaunay", obsp_key="delaunay_graph")
    assert "delaunay_graph" in adata.obsp
    assert adata.obsp["delaunay_graph"].shape == (50, 50)

def test_radius_graph(synthetic_adata):
    adata = generate_spatial_graph(synthetic_adata, method="radius", radius=20, obsp_key="radius_graph")
    assert "radius_graph" in adata.obsp
    assert adata.obsp["radius_graph"].shape == (50, 50)

def test_load_spatial_graph():
    adj_matrix = sp.coo_matrix(([1, 2, 3], ([0, 1, 2], [1, 2, 0])), shape=(3,3))

    # dummy AnnData
    adata = sc.AnnData(np.zeros((3, 3)))
    adata.obsp['spatial_connectivities'] = adj_matrix.toarray()

    # temp h5ad file
    with NamedTemporaryFile(suffix='.h5ad') as temp_file:
        adata.write(temp_file.name)
        graph = load_spatial_graph(temp_file.name)
        # assertions
        assert isinstance(graph, nx.Graph), 'Returned object is not a NetworkX graph.'
        assert not nx.is_directed(graph), 'Graph should be undirected'
        assert len(graph.nodes) == 3, 'Graph should have 3 nodes.'
        assert len(graph.edges) == 3, 'Graph should have 3 edges.'
        # check expected edge weights
        expected_weights = {(0, 1): 1, (1, 2): 2, (2, 0): 3}
        for u, v, data in graph.edges(data=True):
            assert 'weight' in data, f"Edge ({u}, {v}) has no weight."
            assert data['weight'] == expected_weights.get((u, v)) or expected_weights.get((v, u)), f"Unexpected weight {data['weight']} for edge ({u}, {v})."

@pytest.mark.parametrize('method', ['louvain', 'leiden', 'spectral'])
def test_cluster_spatial_graph(method):
    test_graph = nx.karate_club_graph()
    cluster_spatial_graph(test_graph)
    clusters = nx.get_node_attributes(test_graph, 'cluster')
    assert len(clusters) == len(test_graph.nodes())
    # ensure at least 2 clusters
    unique_clusters = set(clusters.values())
    assert len(unique_clusters) >= 2


if __name__ == "__main__":
    pytest.main()
