import numpy as np
from src.detection.spatial_analysis import compute_nearest_neighbor_distances

def test_compute_nearest_neighbor_distances():
    coords = np.array([[0, 0], [1, 1], [2, 2]])
    distances = compute_nearest_neighbor_distances(coords)
    assert len(distances) == len(coords),
    assert all(d > 0 for d in distances), 
