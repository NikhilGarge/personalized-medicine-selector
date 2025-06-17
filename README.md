# Personalized Medicine Framework for Type 2 Diabetes Treatment Optimization

## Overview

This project demonstrates an advanced computational framework for personalized treatment selection that combines genomic risk profiling with both traditional machine learning and deep learning approaches. The framework optimizes Type 2 Diabetes treatment selection through:

**Multi-Method Architecture**:
- **Machine Learning Pipeline (R)**
  - Model-based recursive partitioning (MOB) for idnetifying gene-treatment interactions
  - GWAS-inspired data simulation with SNP-treatment effect modeling
  - Feature selection using pmforest for genomic marker identification
  - Treatment effect estimation via `model4you` library

- **Deep Learning Pipeline (TensorFlow/Keras)**
  - Deep neural networks for clinical outcome prediction (HbA1c levels)
  - Embedding layers for genetic variant representation
  - Multi-treatment interaction modeling
  - Ensemble prediction of optimal treatment regimens

**Implementation Highlights**:
- Hybrid R/Python architecture for comparative analysis
- Genomic data simulation with realistic SNP-treatment interactions
- MOB-based patient stratification combined with neural network predictions
- Cross-validate performance metrics for both approaches

**Key Components**:
1. Synthetic dataset generation (10k patients, 2000 SNPs)
2. Dimensionality reduction for SNP selection (in R).
3. No dimensionality reduction needed in Deep Learning model for predicting clinical outcomes for different treatments
4. Comparative model evaluation framework
5. Treatment recommendation engine with explainable AI components

Skills: R, Python, TensorFlow, Keras, machine learning, deep learning, genomics, treatment effect analysis, data simulation

---

## Installation  

### Requirements  
- Python (version ≥ 3.8) with packages: tensorflow, keras, numpy, pandas, scikit-learn
- R (≥ 4.0) with packages: model4you, pmforest, data.table, caret
- Git

## Steps

### Clone Repository  
`git clone https://github.com/NikhilGarge/personalized-medicine-selector.git/`  
`cd personalized-medicine-selector/`

### Install Python Modules (from terminal)
`pip install tensorflow keras numpy pandas scikit-learn`

### Install R Packages (from R console)  
`install.packages(c("model4you", "party", "caret", "dplyr", "tidyr"))`

## Project Structure

- `Personalized_Treatment_genomics_T2D.ipynb`  
  Python script for data simulation, deep learning model training, and evaluation. The goal of this code is to predict optimal treatment (tailored to individual's geneitic profile).
  
- `Genomic_Treatment_Selection_T2D.R`  
  R script for data simulation, model training, and evaluation. The goal of this code is Personalized Treatment effect estimation via `model4you` library.

- `README.md`  
  Project documentation.

- `results/`  
  Personalized Treatment effect estimation, model performance metrics, Optimal treatment prediction using deep learning

---

## How to Run

1. **Machine learning model (R) - Personalized treatment effect estimation:**  
  Train and evaluate the MOB-based model:
  `source("Genomic_Treatment_Selection_T2D.R")/`  
  `results <- run_full_workflow(
  n_patients = 10000,
  n_snps = 2000,
  n_top_snps = 20,
  treatment_effects = c(T1 = 0, T2 = -0.3, T3 = -0.6, T4 = -0.9),
  ntree = 300
  )`

2. **Deep Learning Model (Python) - Personalized treatment Prediction:**  
  Train and evaluate the TensorFlow model:
  
  `python Personalized_Treatment_genomics_T2D.py/`

3. **Key Outputs:**
- Model accuracy metrics - MAE and R-Squared - for Deep learning (python) and MOB(R) models on train and test sets. 
- Personalized/Optimal treatment prediction using Deep Learning model - high prediction accuracy.
- Personalized treatment effect estimations using `model4you` - highly explainable model.

---

## Applications  
- Demonstrates how **genomic data** can bridge the gap between research and clinical practice.  
- Provides methodology for:  
  - Precision treatment selection 
  - Drug repurposing based on genetic subgroups.

---

