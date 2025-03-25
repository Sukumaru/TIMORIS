# TIMORIS
**Tubule Identification and Mapping for Optimized Recognition in Spatial transcriptomics**
## Inspiration

Spatial transcriptomics has revolutionized the ability to study gene expression in the spatial context of tissues, uncovering structural and functional relationships. However, detecting complex tissue architectures, such as tubule-like structures, remains a challenge. TIMORIS was developed to address this gap by combining computational power with biological insights, allowing researchers to better analyze and interpret spatial transcriptomics data.


The need to detect and characterize specialized tissue structures in kidney, liver, germ and other organ systems where tubules play a crucial role.

The advancements in spatially resolved transcriptomics technologies that provide unprecedented resolution but lack standardized tools for spatial pattern detection.

The growing interest in integrating machine learning and statistical modeling with biological data to reveal new insights into tissue organization and disease mechanisms.

### Understanding Tubule Structures

Tubules are **elongated, tube-like anatomical structures** found in various biological systems. They play critical roles in **filtration, secretion, absorption, and transport** within tissues. 

### 🔬 **Examples of Tubule Structures in Human Biology**

| **Tubule Type**           | **Function**                                    | **Example Tissue**       |
|---------------------------|------------------------------------------------|--------------------------|
| **Renal Tubules**         | Filtration and absorption of waste and nutrients | **Kidneys**              |
| **Mammary Ducts**         | Transport of milk from glands to nipple         | **Breast Tissue**        |
| **Seminiferous Tubules**  | Sperm production and maturation                 | **Testes**               |
| **Bile Ducts**            | Transport of bile for digestion                 | **Liver & Gallbladder**  |

## Description

This repository contains a tool developed as part of the BIOINF 576 course to detect tubule-like structures in tissue using spatial transcriptomics data, preferably at single-cell resolution. The tool integrates spatial pattern detection algorithms and transcriptomics data analysis to identify and characterize tubule-like formations within tissue samples.


## Getting Started

### Input

File Format: .h5ad (AnnData format used for spatial transcriptomics)

Source: Data obtained from public repositories such as 10x Genomics Visium, SpatialDB, or custom experimental datasets.

#### Structure:

- Gene expression matrix (rows = genes, columns = cells/spots)

- Spatial coordinates of each cell/spot

- Metadata containing sample details

### Output

File Format: .h5ad with additional annotations

Contains:

- Detected tubule-like structures with labeled clusters

- Spatial statistics including neighborhood density and curvature scores

- Visualization overlays for detected tubule regions

## Features

- ### Spatial Pattern Detection: Leverages advanced spatial statistics and machine learning techniques to detect elongated tubule-like structures.

- Single-Cell Resolution: Optimized for high-resolution spatial transcriptomics data, ensuring that fine-grained details of tubule-like structures are preserved and accurately detected.

- Visualization: Provides intuitive visual outputs to overlay detected structures on tissue images and expression maps.


## Tutorial

Tutorial is available at `tutorial/tutorial.pynb`
