import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN

def sample_empty_centroids(adata, grid_size=200, min_distance=100, density_threshold=1e-5):
    """
    Sample potential tubule centroids by scanning for empty areas in a grid.

    Parameters:
    - adata: AnnData object with adata.obsm["spatial"]
    - grid_size: spacing between grid points (in spatial units)
    - min_distance: minimum distance from nearest cell
    - density_threshold: minimum spatial density required to count as "empty"

    Returns:
    - centroids: ndarray of shape (n, 2), coordinates of empty centroids
    """
    coords = adata.obsm["spatial"]
    density = adata.obs.get("spatial_density", None)
    if density is None:
        raise ValueError("Run `estimate_spatial_density()` first.")

    # Create a grid over spatial space
    x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
    y_min, y_max = coords[:, 1].min(), coords[:, 1].max()
    xx, yy = np.meshgrid(
        np.arange(x_min, x_max, grid_size),
        np.arange(y_min, y_max, grid_size)
    )
    grid_points = np.c_[xx.ravel(), yy.ravel()]

    # Use nearest neighbor distance to check isolation
    nn = NearestNeighbors(n_neighbors=1)
    nn.fit(coords)
    dist, _ = nn.kneighbors(grid_points)

    # Filter grid points that are far from any cell
    keep = dist[:, 0] > min_distance

    return grid_points[keep]

def cluster_centroids(centroids, eps=300, min_samples=2):
    """
    Cluster nearby centroid candidates to collapse them into single tubule centers.

    Parameters:
    - centroids: np.ndarray of shape (n, 2)
    - eps: maximum distance between points to be considered in the same cluster
    - min_samples: minimum points to form a cluster

    Returns:
    - cluster_centers: list of [x, y] points (one per detected tubule)
    """
    if len(centroids) == 0:
        return np.array([])

    db = DBSCAN(eps=eps, min_samples=min_samples)
    labels = db.fit_predict(centroids)

    cluster_centers = []
    for label in set(labels):
        if label == -1:
            continue  # skip noise
        points = centroids[labels == label]
        mean_point = points.mean(axis=0)
        cluster_centers.append(mean_point)

    return np.array(cluster_centers)
