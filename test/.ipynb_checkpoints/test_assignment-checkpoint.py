import sys
import os
import numpy as np
import pytest
from anndata import AnnData

# Ensure timoris/ is importable
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../timoris")))

from timoris import assign_tubule_regions_by_rays

@pytest.fixture
def synthetic_ring_tubule():
    """
    Generate a ring of cells around a centroid to simulate a tubule.
    """
    n_cells = 100
    radius = 100
    angle = np.linspace(0, 2 * np.pi, n_cells, endpoint=False)

    # Create ring coordinates
    x = 500 + radius * np.cos(angle)
    y = 500 + radius * np.sin(angle)
    coords = np.vstack([x, y]).T

    # Create AnnData object with empty gene matrix, just for spatial testing
    adata = AnnData(np.zeros((n_cells, 10)))
    adata.obsm["spatial"] = coords

    # Define a centroid at the center of the ring
    centroid = np.array([[500, 500]])

    return adata, centroid

def test_assign_tubule_regions_by_rays(synthetic_ring_tubule):
    adata, centroid = synthetic_ring_tubule
    adata = assign_tubule_regions_by_rays(adata, centroid, n_directions=72, max_distance=150)

    assert "tubule_id" in adata.obs.columns
    labels = adata.obs["tubule_id"].astype(int).values

    # At least half of the ring cells should be captured
    assert (labels != -1).sum() > len(labels) // 2

    # All assigned cells should have the same tubule ID
    assert len(set(labels[labels != -1])) == 1
