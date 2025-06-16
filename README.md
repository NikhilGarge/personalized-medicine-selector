# Personalized Medicine Framework for Type 2 Diabetes Treatment Optimization
Developed a machine learning pipeline to predict optimal diabetes treatments using simulated genomic data. Implemented data preprocessing, feature selection, model training, and evaluation in R and Python. Skills: R, Python, machine learning, genomics

## Overview

This project demonstrates a computational framework for **personalized treatment selection** that integrates genomic risk profiling with machine learning to optimize treatment selection for Type 2 Diabetes. By combining simulated genome-wide association study (GWAS) data and model-based recursive partitioning (MOB), we provide a framework for predicting personalized treatment responses based on genetic markers.

---

## Key Features  

### 1. **Simulated Genomic Data Generation**  
   - Generates synthetic SNP datasets mimicking real-world genetic diversity.  
   - Models gene-treatment interactions to identify treatment-responsive variants[1][4].

### 2. **SNP Prioritization**  
   - Identifies top 100 statistically significant SNPs influencing treatment efficacy through gene-treatment interaction mapping[4].

### 3. **Personalized Treatment Modeling**  
   - Implements **ensemble MOB trees** using the `model4you` R library to:  
   - Stratify patients into subgroups with differential treatment responses[2][5].  
   - Estimate personalized treatment effects via bootstrapped subsampling[2][5].

### 4. **Clinical Decision Support**  
   - Outputs treatment recommendations tailored to individual genetic profiles.  
   - Provides interpretable models for clinical translation[1][5].

---

## Installation  

### Requirements  
- R (≥ 4.0)  
- Python (≥ 3.8, optional for integration)

### R Packages  
install.packages(c("model4you", "party", "caret", "dplyr", "tidyr"))

### Clone Repository  
git clone https://github.com/NikhilGarge/personalized-medicine-selector.git
cd personalized-medicine-selector

## Project Structure

- `Treatment_response_Code_v5.R`  
  Main R script for data simulation, model training, and evaluation.

- `README.md`  
  Project documentation.

- (Optional) `notebooks/`  
  Jupyter or R Markdown notebooks for exploratory analysis.

- `results/`  
  Output files, model summaries, and performance metrics.

---

## How to Run

1. **Clone the repository:**
git clone https://github.com/NikhilGarge/personalized-medicine-selector.git
cd personalized-medicine-selector

2. **Install dependencies:**  
- R packages: `model4you`, `party`, `caret`, `dplyr`, etc.
- (Optional) Python: `rpy2` for R-Python integration.

3. **Run the main script:**
Execute the main R script with default parameters:  

source("Treatment_response_Code_v5.R")

results <- run_full_workflow(
n_patients = 10000,
n_snps = 1000,
n_top_snps = 100,
treatment_effects = c(T1 = 0, T2 = -0.3, T3 = -0.6, T4 = -0.9)
)

4. **Key Outputs:**
- train_data / test_data: training and test sets.
- best_treatment_pred: Optimal treatment selected per patient on test data.
- eval_results: Model evaluation using R-Squared metric on training and test sets. This helps in fine tuning the parameters of model to avoid underfitting/overfitting of models

---

## Applications  
- Demonstrates how **genomic data** can bridge the gap between research and clinical practice.  
- Provides methodology for:  
  - Precision treatment selection 
  - Drug repurposing based on genetic subgroups[1][2][5].

---

## References  
1. Personalized medicine projects using genomic data and MOB models[1][4][5].  
2. Treatment effect analysis with `model4you` and ensemble methods[2][5].
