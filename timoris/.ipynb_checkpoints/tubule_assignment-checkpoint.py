import numpy as np
from collections import defaultdict
from sklearn.neighbors import KDTree
import pandas as pd 

def assign_tubule_regions_by_rays(adata, centroids, n_directions=72, max_distance=500):
    """
    Assign tubule regions by casting rays from centroids in 360 degrees.

    Parameters:
    - adata: AnnData with adata.obsm["spatial"]
    - centroids: ndarray of shape (n_centroids, 2)
    - n_directions: number of angles to sweep (e.g., 72 → every 5°)
    - max_distance: how far each ray can search

    Returns:
    - AnnData with adata.obs["tubule_id"] assigned
    """
    coords = adata.obsm["spatial"]
    tree = KDTree(coords)
    angle_step = 2 * np.pi / n_directions

    # Map: cell index → set of centroid IDs that reached it
    cell_hits = defaultdict(list)

    for tid, center in enumerate(centroids):
        for i in range(n_directions):
            theta = i * angle_step
            direction = np.array([np.cos(theta), np.sin(theta)])
            for dist in np.linspace(10, max_distance, num=100):
                probe = center + dist * direction
                dist_to_cell, idx = tree.query([probe], k=1)
                if dist_to_cell[0][0] < 20:  # within range of a cell
                    cell_hits[idx[0][0]].append(tid)
                    break  # stop after first hit

    # Resolve overlaps: assign each cell to closest centroid
    cell_to_tubule = -1 * np.ones(adata.n_obs, dtype=int)

    for cell_idx, centroid_ids in cell_hits.items():
        if len(centroid_ids) == 1:
            cell_to_tubule[cell_idx] = centroid_ids[0]
        else:
            # Break tie by choosing closest centroid
            dists = [np.linalg.norm(coords[cell_idx] - centroids[cid]) for cid in centroid_ids]
            closest_id = centroid_ids[np.argmin(dists)]
            cell_to_tubule[cell_idx] = closest_id

    adata.obs["tubule_id"] = pd.Categorical(cell_to_tubule)
    return adata