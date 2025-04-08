The data will look like:
- Spatial Transcriptomics data (with paired gene expression and spatial location data)

Repositories where this data is available:
## 1. Mouse Testis Spatial Transcriptomics (Example Data)

- **Description**: Spatial transcriptomics dataset capturing germ cell organization and tubule structures in adult mouse testis using high-resolution MERFISH.
- **Source**: Jun Z. Li Lab, University of Michigan
- **Link**: [Spatial Viewer](https://junzlilab-shiny.med.umich.edu/spatialViewer/)
- **Format**: MERFISH spatial transcriptomics with 491 gene panel; single-cell resolution.

### Key metrics were:

- Number of cells per sample: ~10,000–20,000  
- Genes measured: 491 (targeted panel)  
- Imaging system: Vizgen MERSCOPE & Nikon confocal  
- Tubule structures segmented: Yes (custom ray-based geometry + refinement)  
- Spatial resolution: Subcellular  
- Metadata includes: Cell coordinates, tubule assignments, gene expression matrix

### Justification for this Dataset:

This dataset is an ideal choice to test the **TIMORIS** tool because it contains **spatially resolved transcriptomics data** from **mouse testis** — a tissue where **tubule structures** (e.g., seminiferous tubules) are crucial for development and function. Tubule-like structures are readily identifiable in testis tissue, making it an excellent biological context for testing TIMORIS. The high-resolution single-cell data also allows for fine-grained analysis of tubule morphology and gene expression.

### Biological Question:

Using the **TIMORIS** tool, we aim to answer the following question:

- **How can spatially resolved gene expression patterns reveal the structural and functional organization of seminiferous tubules in mouse testis tissue?**

Specifically, we will:
- Detect tubule-like regions in the tissue based purely on **spatial coordinates and density**, without relying on specific gene expression markers.

### Expected Results:

- **Tubule regions** will be detected based on their **spatial geometry** (i.e., centroids, boundary cells).
- **Visualizations** should show tubule-like clusters, with each tubule region identified and labeled.
- **Expected answer**: The tool will assign **a unique `tubule_id`** to each tubule-like structure, and the **spatial layout** of these regions will align with known anatomical features of the seminiferous tubules.

## 2. Human Prostate Cancer Spatial Transcriptomics 
- **Description**: Highlights spatial gene expression in normal and cancerous prostate glandular tubule structures.
- **Source**: [10x Genomics](https://www.10xgenomics.com/)
- **Link**: [Prostate Cancer Visium Dataset](https://www.10xgenomics.com/datasets/human-prostate-cancer-adjacent-normal-section-with-if-staining-ffpe-1-standard)
- **Format**: Spatial transcriptomics, annotated histology images.
### Key metrics were:

- Spots detected under tissue (Loupe Manual Alignment): 3,460
- Median genes per spot: 4,614
- Median UMI Counts per spot: 11,444

## 3. Human Breast Cancer Spatial Transcriptomics
- **Description**: Spatial gene expression in Human Breast Cancer
- **Source**: [Broad Institute - Single Cell Portal](https://singlecell.broadinstitute.org/)
- **Link**: [Human Breast Cancer](https://singlecell.broadinstitute.org/](https://singlecell.broadinstitute.org/single_cell/study/SCP1256/visium-demo-study#study-summary))
- **Alternative Link**: https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Breast_Cancer_Block_A_Section_1
- **Format**: Spatial transcriptomics.

### Key cell metrics were:

- Spots detected under tissue - 3,798
- Median UMI counts per spot - 20,762
- Median genes per spot - 6,026
