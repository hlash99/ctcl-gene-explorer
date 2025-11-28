# CTCL Single-Cell Gene Explorer 🧬

**A Clinical Decision Support Tool for Cutaneous T-Cell Lymphoma**

[Link to Live App](https://ctcl-gene-explorer-h8jfajpu28z7ndpimyg8kr.streamlit.app/)

## Overview
Differentiating advanced Cutaneous T-Cell Lymphoma (CTCL) from benign inflammatory conditions like Eczema is a significant clinical challenge. This dashboard visualizes single-cell RNA sequencing (scRNA-seq) data from **40,000+ cells** (GSE128531) to identify potential diagnostic markers.

## Features
- **Single-Cell Atlas:** Interactive UMAP visualization of 40k cells across Tumor, Normal, and Eczema samples.
- **Differential Expression:** Real-time comparison of gene expression between malignant and benign T-cells.
- **Automated Insights:** Logic layer that interprets expression trends (e.g., "Significantly Higher in Tumor").
- **Optimized Performance:** Uses sparse matrix compression to serve ~6GB of genomic data in a <100MB lightweight app.

## Tech Stack
- **Python** (Pandas, NumPy, Scipy)
- **Bioinformatics:** Scanpy (Single-cell analysis suite)
- **Visualization:** Streamlit, Matplotlib, Seaborn
- **Data:** NCBI GEO (Accession GSE128531)

## How to Run Locally
1. Clone the repo:
   ```bash
   git clone [https://github.com/YourUsername/ctcl-gene-explorer.git](https://github.com/YourUsername/ctcl-gene-explorer.git)
