from sklearn.neighbors import kneighbors_graph
from scipy.sparse import csr_matrix
import networkx as nx
import numpy as np
from scipy.spatial import ConvexHull
import pandas as pd

def refine_tubules(adata, k=6, min_neighbor_agreement=0.5, max_fill_distance=40, use_hulls=True):
    """
    Refine ray-assigned tubule boundaries using spatial smoothing, graph consistency,
    and (optional) convex hull-based filling.

    Parameters:
    - adata: AnnData with `tubule_id` and `obsm["spatial"]`
    - k: neighbors for graph smoothing
    - min_neighbor_agreement: threshold for removing outliers
    - max_fill_distance: max distance for filling in missed cells near tubule boundaries
    - use_hulls: if True, expand tubules using convex hulls

    Returns:
    - adata: updated AnnData with refined `tubule_id`
    """

    coords = adata.obsm["spatial"]
    tubule_ids = adata.obs["tubule_id"].astype("object").copy()

    # 1. Remove outliers
    graph = kneighbors_graph(coords, k, mode="connectivity", include_self=False)
    G = nx.from_scipy_sparse_array(graph)

    for node in G.nodes:
        tid = tubule_ids.iloc[node]
        if tid == -1:
            continue
        neighbors = list(G.neighbors(node))
        neighbor_ids = tubule_ids.iloc[neighbors]
        match = (neighbor_ids == tid).sum()
        if match / k < min_neighbor_agreement:
            tubule_ids.iloc[node] = -1

    # 2. Fill nearby cells based on neighborhood
    new_ids = tubule_ids.copy()
    for node in G.nodes:
        if new_ids.iloc[node] != -1:
            continue
        neighbors = list(G.neighbors(node))
        valid = tubule_ids.iloc[neighbors][tubule_ids.iloc[neighbors] != -1]
        if len(valid) == 0:
            continue
        most_common = valid.value_counts().idxmax()
        if (valid == most_common).sum() >= 3:
            new_ids.iloc[node] = most_common

    tubule_ids = new_ids.copy()

    # 3. Optional: Convex hull fill-in
    if use_hulls:
        for tid in set(tubule_ids):
            if tid == -1:
                continue
            mask = (tubule_ids == tid)
            if mask.sum() < 3:
                continue
            cluster_coords = coords[mask.values]
            try:
                hull = ConvexHull(cluster_coords)
                polygon = cluster_coords[hull.vertices]

                from matplotlib.path import Path
                poly_path = Path(polygon)

                for i in range(adata.n_obs):
                    if tubule_ids.iloc[i] != -1:
                        continue
                    if poly_path.contains_point(coords[i]):
                        dist = np.min(np.linalg.norm(coords[i] - polygon, axis=1))
                        if dist <= max_fill_distance:
                            tubule_ids.iloc[i] = tid
            except:
                continue

    adata.obs["tubule_id"] = pd.Categorical(tubule_ids)
    return adata