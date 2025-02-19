import pytest
import numpy as np
import anndata
from detection.spatial_graph import generate_spatial_graph

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

if __name__ == "__main__":
    pytest.main()
