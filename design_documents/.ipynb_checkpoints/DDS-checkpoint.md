# Design Document Specification (DDS)

## 1. Goals and Milestones

### Goals
- Develop TIMORIS, a tool for detecting tubule-like structures in spatial transcriptomics data using geometric and graph-based reasoning.
- Implement efficient, biologically grounded algorithms for tubule detection (e.g., spatial density mapping, centroid-based region growing).
- Provide intuitive spatial visualizations and diagnostic plots.
- Ensure compatibility with single-cell spatial datasets (e.g., MERFISH, seqFISH, Visium).
- Enable downstream marker discovery and morphological analysis.

### Milestones
- **Milestone 1**: Identify and preprocess spatial transcriptomics datasets. **(✅ complete)**
- **Milestone 2**: Develop and optimize spatial centroid + ray-based tubule detection algorithm. **(✅ functional, tuning in progress)**
- **Milestone 3**: Implement adaptive refinement and boundary smoothing using graph theory. **(✅ complete)**
- **Milestone 4**: Build and validate spatial visualizations and overlays. **(🟡 partial, ongoing)**
- **Milestone 5**: Package and document TIMORIS for general use. **(🟢 75% complete)**

---

## 2. Modules and Their Dependencies

### **Module 1: Data Preprocessing**
- **Scope**: Load and normalize spatial transcriptomics datasets.
- **Contents**:
  - `preprocess_data()` – filtering, log-normalization, HVG selection
  - `estimate_spatial_density()` – KDE-based local density scoring
- **Dependencies**: Downstream modules rely on this for cleaned inputs and spatial structure.

---

### **Module 2: Centroid Detection**
- **Scope**: Identify candidate tubule centers (non-cell spatial voids).
- **Contents**:
  - `sample_empty_centroids()` – grid-based sampling of low-density space
  - `cluster_centroids()` – collapse nearby candidates into single centroid per tubule
- **Dependencies**: Uses spatial density from Module 1, provides seeds for tubule region assignment.

---

### **Module 3: Spatial Graph and Region Assignment**
- **Scope**: Build spatial graphs and assign cells to tubules via geometry.
- **Contents**:
  - `generate_spatial_graph()` / `load_spatial_graph()` – construct kNN spatial graph
  - `assign_tubule_regions_by_rays()` – 360° ray-casting from each centroid to identify boundary cells
- **Dependencies**: Centroid list (Module 2), spatial coordinates (Module 1)

---

### **Module 4: Tubule Refinement**
- **Scope**: Refine tubule boundaries using connectivity, density, and hull logic.
- **Contents**:
  - `refine_tubules()` – removes outliers, fills boundary gaps, applies convex hull expansion
- **Dependencies**: Works on `tubule_id` assignments from Module 3

---

### **Module 5: Visualization and Reporting**
- **Scope**: Generate summary plots and support evaluation/diagnostics.
- **Contents** (to be expanded):
  - Spatial overlays of `tubule_id`
  - Radial profiles or convex hulls (future)
- **Dependencies**: All prior modules

---

## 3. Future Enhancements

- Add tubule-specific marker gene discovery and radial maturity scoring
- Export boundary outlines (e.g., convex hull polygons) for figure-ready overlays
- Integrate automatic maturity gradient extraction per tubule
- GUI support via napari or Shiny (optional)

