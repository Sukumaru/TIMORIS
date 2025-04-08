# Work Breakdown Structure (WBS) for TIMORIS

## 1. First Basic Output: Tubule-Like Structure Detection

### 1.1 Goal

Detect tubule-like structures in spatial transcriptomics data using spatial geometry and graph-based refinement — without requiring predefined marker genes.

---

### 1.2 Steps to Compute Output

#### 1.2.1 Data Preparation ✅
- **Dataset selection**: Small-scale spatial transcriptomics dataset for testing (e.g., mouse testis).
- **Format verification**: Ensure data is in `.h5ad` (AnnData) format.
- **Preprocessing**:
  - Normalize gene expression values.
  - Filter low-quality or low-expression cells.
  - Extract spatial coordinates per cell.

#### 1.2.2 Feature Extraction ✅
- **Compute spatial statistics**:
  - Estimate local density via kernel density estimation.
  - Compute nearest-neighbor distances.
  - (Optional) Apply spatial autocorrelation metrics (e.g., Moran’s I).

#### 1.2.3 Structure Detection ⚙️
- **Centroid identification**:
  - Sample centroids from spatial density valleys.
  - Cluster centroids using DBSCAN to ensure one per tubule.
- **Tubule assignment**:
  - Use 360° radial ray-casting to assign boundary cells.
  - Resolve overlaps via distance-based conflict resolution.
- **Refinement**:
  - Remove spatial outliers.
  - Expand or smooth boundaries using kNN graph and convex hulls.

#### 1.2.4 Output Generation ⚙️
- **Visualization**:
  - Overlay detected tubules on spatial plots.
  - Color-code by `tubule_id` or confidence.
- **Export**:
  - Save annotated `.h5ad` with `obs["tubule_id"]` and metadata.
  - Optionally export boundary outlines or convex hull polygons.

---

## 2. Additional Tasks

### 1.2.5 Optimization of Detection Algorithms
- **Performance testing**:
  - Benchmark clustering and ray assignment on various datasets.
  - ✅ Deliverable: Performance summary report.
- **Parameter tuning**:
  - Optimize settings for centroid clustering, ray density, and graph refinement.
  - ✅ Deliverable: Documented best parameters + validation set metrics.

### 1.2.6 Integration with Visualization Tools
- **Interactive visualization**:
  - Build support for viewing tubules with Napari or Scanpy.
  - 🎯 Completion: Users can inspect and compare tubules visually.
- **Static plot export**:
  - Generate image overlays for use in figures and validation reports.

### 1.2.7 Validation with Real Experimental Data
- **Dataset selection**:
  - Use internal MERFISH testis dataset (or public equivalent).
- **Expert comparison**:
  - Collaborate with domain expert(s) for manual annotations.
  - Deliverable: Agreement metrics (precision, recall, IoU).

---