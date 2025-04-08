# TIMORIS  
**Tubule Identification and Mapping for Optimized Recognition in Spatial transcriptomics**

## Inspiration

Spatial transcriptomics has revolutionized the ability to study gene expression in the spatial context of tissues, uncovering structural and functional relationships. However, detecting complex tissue architectures, such as tubule-like structures, remains a challenge. TIMORIS was developed to address this gap by combining computational power with biological insights, allowing researchers to better analyze and interpret spatial transcriptomics data.

- The need to detect and characterize specialized tissue structures in kidney, liver, germ, and other organ systems where tubules play a crucial role.
- The advancements in spatially resolved transcriptomics technologies that provide unprecedented resolution but lack standardized tools for spatial pattern detection.
- The growing interest in integrating geometric modeling and spatial algorithms with biological data to reveal new insights into tissue organization and disease mechanisms.

### Understanding Tubule Structures

Tubules are **elongated, tube-like anatomical structures** found in various biological systems. They play critical roles in **filtration, secretion, absorption, and transport** within tissues. 

### 🔬 **Examples of Tubule Structures in Human Biology**

| **Tubule Type**           | **Function**                                    | **Example Tissue**       |
|---------------------------|--------------------------------------------------|--------------------------|
| **Renal Tubules**         | Filtration and absorption of waste and nutrients | **Kidneys**              |
| **Mammary Ducts**         | Transport of milk from glands to nipple         | **Breast Tissue**        |
| **Seminiferous Tubules**  | Sperm production and maturation                 | **Testes**               |
| **Bile Ducts**            | Transport of bile for digestion                 | **Liver & Gallbladder**  |

---

## Description

TIMORIS is a Python package for detecting tubule-like structures in high-resolution spatial transcriptomics data, preferably at single-cell resolution. It combines **geometric inference**, **ray-based spatial assignment**, and **graph-based refinement** to produce accurate segmentations of tubular tissue regions — even in the absence of gene expression markers.

---

## Getting Started

### Input

- **File Format**: `.h5ad` (AnnData format)
- **Source**: Data from platforms such as 10x Genomics Visium, Vizgen MERFISH, or custom tissue datasets

#### Required Data Structure

- Gene expression matrix (rows = genes, columns = cells/spots)
- Spatial coordinates per cell/spot (`adata.obsm["spatial"]`)
- Optional: metadata or pseudotime/maturity scores

### Output

- **File Format**: `.h5ad` with updated `obs["tubule_id"]`
- **Includes**:
  - Detected tubule region labels
  - Spatial density estimates
  - Annotated centroid positions
  - Optional maturity scores and gene enrichment metrics
  - Overlay plots for result visualization

---

## Features

- ### 🧭 Ray-Based Tubule Assignment  
  Sweeps 360° radially from candidate centroids to assign boundary cells, simulating geometric tracing of tubular cavities.

- ### 🧠 Adaptive Centroid Detection  
  Uses spatial density valleys and clustering to detect non-cell centroid positions in tubule-like cavities.

- ### 🪢 Graph-Based Refinement  
  Refines coarse boundaries using k-nearest neighbor graphs, outlier filtering, and convex hull-based edge completion.

- ### 🧬 Marker Discovery (Optional)  
  Enables post hoc identification of marker genes for each tubule structure via differential expression analysis.

- ### 🖼 Visual Overlay  
  Easily visualizes detected tubules over spatial coordinates or histology backgrounds using `scanpy` or `matplotlib`.

---

## Tutorial

A Jupyter tutorial is available at `tutorial/tutorial.ipynb` demonstrating the full tubule detection pipeline on example data.

---