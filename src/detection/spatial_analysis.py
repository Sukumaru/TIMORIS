import numpy as np
from scipy.spatial import KDTree

def compute_nearest_neighbor_distances(coords):
    """
    Compute the nearest neighbor distances for a set of spatial coordinates.

    Parameters:
    coords (array-like): A list of (x, y) coordinates.

    Returns:
    np.ndarray: An array of nearest neighbor distances.
    """
    if len(coords) < 2:
        raise ValueError("At least two coordinates are required.")
    
    tree = KDTree(coords)
    distances, _ = tree.query(coords, k=2)  # k=2 to get the first nearest neighbor (excluding self)
    return distances[:, 1]  # Return the second column, which contains the nearest neighbor distances
