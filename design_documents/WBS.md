# Work Breakdown Structure (WBS) for TIMORIS

## 1. First Basic Output: Tubule-Like Structure Detection

### 1.1 Goal
The first basic output of TIMORIS will be the identification of tubule-like structures from spatial transcriptomics data. This involves detecting spatial patterns based on gene expression and spatial coordinates.

### 1.2 Steps to Compute Output
#### 1.2.1 Data Preparation
- **Identify dataset**: Select a small-scale spatial transcriptomics dataset for initial testing.
- **Format verification**: Ensure the dataset is in `.h5ad` format.
- **Data preprocessing**:
  - Normalize gene expression values.
  - Filter low-quality or low-expression cells.
  - Extract spatial coordinates of each cell.

#### 1.2.2 Feature Extraction
- **Identify key marker genes**:
  - Use known tubule marker genes (if available) or perform gene selection based on clustering.
- **Compute spatial statistics**:
  - Calculate neighborhood density and nearest-neighbor distances.
  - Apply spatial autocorrelation metrics (e.g., Moran’s I) to detect local structure.

#### 1.2.3 Structure Detection
- **Apply clustering algorithms**:
  - Use DBSCAN, k-means, or graph-based clustering to detect spatially structured regions.
- **Shape analysis**:
  - Compute elongation and curvature metrics to determine tubule-like morphology.
- **Validate detected structures**:
  - Compare detected clusters with known tubule locations (if annotated data is available).
  - Compute performance metrics such as precision and recall.

#### 1.2.4 Output Generation
- **Generate visualization**:
  - Overlay detected tubules on spatial plots.
  - Color-code clusters based on spatial properties.
- **Export results**:
  - Annotated `.h5ad` file with tubule labels.

---
## Additional Tasks
### 1.2.5: Optimization of Detection Algorithms
- **Performance Testing**
  - Benchmark clustering algorithms on test datasets.
  - Deliverable: Performance report summarizing results.
  - Completion Criteria: Performance improvements documented and integrated.

- **Parameter Tuning**
  - Test different parameter settings for clustering.
  - Optimize for highest precision and recall.
  - Deliverable: Optimized parameter settings.
  - Completion Criteria: Selected parameters improve accuracy on validation set.

### 1.2.6: Integration with Visualization Tools
- **Develop Visualization**
  - Implement visualization in **Napari** or **Scanpy**.
  - Deliverable: Interactive visualization module.
  - Completion Criteria: Users can inspect detected tubules interactively.

- **Export Visualization Outputs**
  - Generate static plots for reports.
  - Deliverable: Collection of annotated images.
  - Completion Criteria: Output files successfully generated and documented.

### 1.2.7: Validation with Real Experimental Data
- **Select Experimental Dataset**
  - Choose a publicly available dataset with known tubule structures.
  - Deliverable: Dataset documentation and metadata.
  - Completion Criteria: Dataset is selected and loaded into the pipeline.

- **Compare with Expert Annotations**
  - Collaborate with domain experts to manually annotate structures.
  - Deliverable: Validation report comparing results.
  - Completion Criteria: Agreement metrics (e.g., precision, recall) computed and analyzed.

---
