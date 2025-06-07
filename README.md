# Tumor Deconvolution on LUAD Datasets from TCGA

This repository contains the code and data used for tumor deconvolution focusing on mutant EGFR cells. The analysis is primarily based on LUAD datasets from the TCGA project. This was analyzed for the Hannun Lab, Stony Brook Medicine.

## Repository Structure

This https://zenodo.org/records/10042128 houses the bulk mRNA sequencing data required for deconvolution. The reference for the data is sourced from iHLCA but has been further processed to suit our specific use case. Notably, mutations in EGFR and KRAS were identified. Within EGFR, the exon 19 deletions (del745-750, 45%) and exon 21 substitutions (L858R, 40-45%) are of particular significance, as referenced from the GDC portal.

### `scripts/`

Within this directory, you will find scripts categorized into three primary functions:

-   **Preprocess**: Preparation of data for deconvolution.
-   **Deconvolution**: Actual tumor deconvolution operations.
-   **Postprocess**: Processes following the deconvolution to finalize results.

## Usage

1.  Ensure that you have all necessary prerequisites installed.
2.  Clone this repository.
3.  Navigate to the `scripts/` directory.
4.  Run the preprocessing scripts followed by the deconvolution and then the post-processing scripts.

> *PS - Change the paths according to your folder location.*
