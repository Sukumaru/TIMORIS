# Design Document Specification (DDS)

## 1. Goals and Milestones

### Goals
- Develop TIMORIS, a tool for detecting tubule-like structures in spatial transcriptomics data.
- Implement efficient algorithms for spatial pattern recognition.
- Enable visualization and validation of detected tubules.
- Ensure scalability from small test datasets to full-scale datasets.

### Milestones
- **Milestone 1**: Identify and preprocess spatial transcriptomics datasets. **(done)**
- **Milestone 2**: Develop and optimize tubule detection algorithms. **(done 75%)**
- **Milestone 3**: Implement spatial visualization tools. 
- **Milestone 4**: Validate detection accuracy against expert annotations.
- **Milestone 5**: Release a documented, scalable pipeline. **(done 50%)**

---

## 2. Modules and Their Dependencies

### **Module 1: Data Preprocessing**
- **Scope**: Loads, formats, and normalizes spatial transcriptomics datasets.
- **Content**:
  - File parsers for `.h5ad`
  - Data filtering (low-quality cells removal)
  - Normalization and transformation functions
- **Dependencies**: Required by Feature Extraction and Detection modules.

### **Module 2: Feature Extraction**
- **Scope**: Identifies key marker genes and computes spatial statistics.
- **Content**:
  - Marker gene selection
  - Spatial clustering analysis
  - Nearest-neighbor computations
- **Dependencies**: Relies on preprocessed data from Module 1, feeds into Detection.

### **Module 3: Tubule Structure Detection**
- **Scope**: Identifies and segments tubule-like structures.
- **Content**:
  - Clustering algorithms (DBSCAN, k-means, graph-based)
  - Shape and morphology analysis
  - Tubule structure annotation
- **Dependencies**: Uses features from Module 2, provides output for Visualization.

### **Module 4: Visualization**
- **Scope**: Displays spatial plots and annotated images.
- **Content**:
  - Heatmaps of detected tubules
  - Spatial overlay of tubule structures
  - Visualization (Scanpy)
- **Dependencies**: Uses detected tubule structures from Module 3.

### **Module 5: Validation and Reporting**
- **Scope**: Compares results with expert annotations and outputs summary statistics.
- **Content**:
  - Benchmarking metrics (precision, recall, spatial clustering scores)
  - Performance evaluation
- **Dependencies**: Uses outputs from Detection and Visualization
