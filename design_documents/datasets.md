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
