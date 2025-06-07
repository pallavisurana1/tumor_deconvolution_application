# 🧬 Tumor Deconvolution of LUAD TCGA Datasets

This repository contains code and data used for tumor deconvolution analysis focused on EGFR-mutant cells, using LUAD (Lung Adenocarcinoma) datasets from the TCGA project. The project was conducted for the **Hannun Lab** at **Stony Brook Medicine**.

---

## 📁 Repository Overview

The analysis utilizes bulk mRNA sequencing data available on Zenodo:

🔗 [TCGA-LUAD bulk data (Zenodo)](https://zenodo.org/records/10042128)

The raw data originates from iHLCA and has been processed to support tumor deconvolution workflows. Mutational profiling includes detection of common driver mutations, specifically in **EGFR** and **KRAS**.

**EGFR mutation highlights**:
- Exon 19 deletions (e.g., del745-750) — ~45%
- Exon 21 substitutions (e.g., L858R) — ~40–45%

Data source and annotations are based on the [GDC Portal](https://portal.gdc.cancer.gov/).

---

## 📂 Project Structure

```
scripts/
├── preprocess/      # Scripts for preparing input data
├── deconvolution/   # Scripts performing the deconvolution
└── postprocess/     # Scripts for downstream analysis and visualization
```

Each folder contains self-contained scripts organized by function.

---

## ⚙️ How to Use

1. **Install Dependencies**: Ensure required packages are installed (R or Python, depending on the scripts).
2. **Clone the Repository**:
   ```bash
   git clone https://github.com/your-username/luad-deconvolution.git
   cd luad-deconvolution/scripts
   ```
3. **Run Analysis Pipeline**:
   - Execute preprocessing scripts.
   - Run deconvolution.
   - Follow up with postprocessing scripts for final results.

> ⚠️ **Note:** Update file paths in the scripts based on your local directory structure.

---

## 📬 Contact

For questions or collaboration inquiries, please reach out to:

**Pallavi Surana**  
📧 pallavi.surana@stonybrook.edu
