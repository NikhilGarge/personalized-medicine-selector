#-----------------------------------------------------------
# Title: Personalized Treatment Selection Workflow
# Description: Simulates SNP data, performs gene-treatment interaction analysis,
#              and applies machine learning for treatment selection in T2D.
# Author: Nikhil Garge
# Date: 2025-06-16
# Version: 1.0
# Contact: nikhil.garge@gmail.com
# Dependencies: model4you, party, caret, dplyr, etc.
#-----------------------------------------------------------

#' Simulate SNP Matrix
#'
#' @param n_patients Number of patients
#' @param n_snps Number of SNPs
#' @param maf_min Minimum MAF
#' @param maf_max Maximum MAF
#' @return SNP matrix (patients x SNPs)
#' @export
simulate_snp_matrix <- function(n_patients, n_snps, maf_min = 0.1, maf_max = 0.5) {
  maf_vec <- runif(n_snps, maf_min, maf_max)
  SNP_matrix <- sapply(maf_vec, function(maf) rbinom(n_patients, 2, maf))
  colnames(SNP_matrix) <- paste0("SNP", 1:n_snps)
  return(SNP_matrix)
}

#' Simulate Patient Data
#'
#' @param n_patients Number of patients
#' @param n_snps Number of SNPs
#' @param SNP_matrix SNP matrix
#' @param n_causal_snps Number of causal SNPs
#' @param n_interacting_snps Number of SNPs with treatment interaction
#' @param n_epistatic_pairs Number of epistatic SNP pairs
#' @param treatment_effects Named vector of treatment effects (e.g., c(T1=0, T2=-0.3, ...))
#' @return Data frame with patient info and SNPs
#' @export
simulate_patient_data <- function(
    n_patients,
    n_snps,
    SNP_matrix,
    n_causal_snps = 20,
    n_interacting_snps = 10,
    n_epistatic_pairs = 5,
    treatment_effects = c(T1=0, T2=-0.3, T3=-0.6, T4=-0.9)
) {
  set.seed(123)
  
  # Treatments and assignment
  treatments <- names(treatment_effects)
  n_treatments <- length(treatments)
  treatment_assignment <- sample(treatments, n_patients, replace = TRUE)
  
  # Select causal SNPs for main effects and interactions
  causal_indices <- sample(1:n_snps, n_causal_snps)
  causal_treatments <- sample(treatments, n_causal_snps, replace = TRUE)
  causal_effects <- rnorm(n_causal_snps, mean = -0.6, sd = 0.8)
  
  # Typical baseline HbA1c for T2D
  baseline_hba1c <- 8
  hba1c <- rep(baseline_hba1c, n_patients) + rnorm(n_patients, 0, 0.5)
  
  # Add SNP-treatment interactions
  for (i in 1:n_causal_snps) {
    idx <- which(treatment_assignment == causal_treatments[i])
    hba1c[idx] <- hba1c[idx] + SNP_matrix[idx, causal_indices[i]] * causal_effects[i]
  }
  
  # Add main treatment effects
  for (t in treatments) {
    idx <- which(treatment_assignment == t)
    hba1c[idx] <- hba1c[idx] + treatment_effects[t]
  }
  
  # Add epistatic effects (SNP-SNP interaction)
  if (n_epistatic_pairs > 0) {
    epistatic_indices <- matrix(sample(1:n_snps, 2 * n_epistatic_pairs, replace = FALSE), ncol = 2)
    for (j in 1:n_epistatic_pairs) {
      snp1 <- epistatic_indices[j, 1]
      snp2 <- epistatic_indices[j, 2]
      epi_effect <- rnorm(1, mean = -0.4, sd = 0.15)
      hba1c <- hba1c + (SNP_matrix[, snp1] * SNP_matrix[, snp2]) * epi_effect
    }
  }
  
  # Add random noise
  hba1c <- hba1c + rnorm(n_patients, 0, 0.3)
  
  # Assemble data frame
  geno_data <- data.frame(
    PatientID = 1:n_patients,
    Treatment = treatment_assignment,
    HbA1c = hba1c,
    SNP_matrix
  )
  
  geno_data$Treatment <- as.factor(geno_data$Treatment)
  
  return(geno_data)
}


#' Calculate Interaction P-values for SNPs
#'
#' @param geno_data Data frame with phenotype, treatment, and SNPs
#' @param snp_cols Vector of SNP column names
#' @return Named vector of minimum interaction p-values for each SNP
#' @export
calculate_interaction_pvals <- function(geno_data, snp_cols) {
  sapply(snp_cols, function(snp) {
    model <- lm(HbA1c ~ Treatment * geno_data[[snp]], data = geno_data)
    coef_summary <- summary(model)$coefficients
    interaction_rows <- grep(":", rownames(coef_summary), value = TRUE)
    if (length(interaction_rows) > 0) {
      min(coef_summary[interaction_rows, "Pr(>|t|)"])
    } else {
      NA
    }
  })
}

#' Minor Allele Frequency Filter
#'
#' @param geno SNP vector
#' @return TRUE if MAF >= 0.05, FALSE otherwise
#' @export
maf_filter <- function(geno) {
  p <- mean(geno, na.rm = TRUE) / 2
  maf <- min(p, 1 - p)
  return(maf >= 0.05)
}

#' Filter SNPs by MAF
#'
#' @param geno_data Data frame with SNPs
#' @param snp_cols SNP column names
#' @return Filtered SNP column names
#' @export
filter_snps_by_maf <- function(geno_data, snp_cols) {
  maf_pass <- sapply(snp_cols, function(snp) maf_filter(geno_data[[snp]]))
  snp_cols[maf_pass]
}

#' Prune SNPs for Linkage Disequilibrium
#'
#' @param geno_data Data frame with SNPs
#' @param snp_cols_maf SNP column names after MAF filtering
#' @param threshold LD threshold (R^2)
#' @return Pruned SNP column names
#' @export
prune_ld <- function(geno_data, snp_cols_maf, threshold = 0.8) {
  library(Matrix)
  geno_mat <- as.matrix(geno_data[, snp_cols_maf])
  cor_mat <- cor(geno_mat, use = "pairwise.complete.obs")
  to_remove <- c()
  for (i in 1:(ncol(cor_mat) - 1)) {
    for (j in (i + 1):ncol(cor_mat)) {
      if (abs(cor_mat[i, j])^2 > threshold && !(colnames(cor_mat)[j] %in% to_remove)) {
        to_remove <- c(to_remove, colnames(cor_mat)[j])
      }
    }
  }
  setdiff(snp_cols_maf, to_remove)
}

#' Personalized Treatment Prediction using Model4You
#'
#' @param geno_data Data frame with phenotype, treatment, and top SNPs
#' @param top_snps Vector of top SNP column names
#' @param ntree Number of trees for the forest
#' @return Personalized model object
#' @export
personalized_treatment_forest <- function(geno_data, top_snps, ntree = 10) {
  library(model4you)
  base_model <- lm(HbA1c ~ Treatment, data = geno_data[, c("HbA1c", "Treatment", top_snps)])
  mob_forest <- pmforest(
    base_model,
    data = geno_data[, c("HbA1c", "Treatment", top_snps)],
    control = ctree_control(minbucket = 100),
    ntree = ntree
  )
  #personalized_models <- pm(mob_forest)
  return(mob_forest)
}

#' Evaluate Personalized Model on Test Set
#'
#' For each patient in the test set, fit a personalized model using the model4you forest,
#' extract treatment effect estimates, p-values, and confidence intervals, and recommend the best treatment.
#'
#' @param mob_forest A pmforest object (from model4you)
#' @param train_data Train set data frame
#' @param test_data Test set data frame
#' @return Data frame with predicted effects, CIs, and recommended treatments
#' @export
evaluate_personalized_model <- function(mob_forest, train_data, test_data) {
  R2_train = 1 - (sum((predict(mob_forest, train_data)[1][,1] - train_data$HbA1c)**2) / sum((mean(train_data$HbA1c) - train_data$HbA1c)**2))
  R2_test = 1 - (sum((predict(mob_forest, test_data)[1][,1] - test_data$HbA1c)**2) / sum((mean(test_data$HbA1c) - test_data$HbA1c)**2))
  data.frame(
    R2_train = R2_train,
    R2_test = R2_test
  )
}

#' Computer personalized treatment for each patient in Test Set
#'
#' For each patient in the test set, fit a personalized model using the model4you forest,
#' extract treatment effect estimates, predictions, and recommend the best treatment.
#'
#' @param forest A pmforest object (from model4you)
#' @param test_data Test set data frame
#' @param top_snps Vector of top SNP column names
#' @return Data frame with predicted effects, CIs, and recommended treatments
#' @export
Predict_best_treatment <- function(forest, test_data, top_snps) {
  library(model4you)
  results <- lapply(1:nrow(test_data), function(i) {
    # Prepare new patient data as a data.frame
    new_patient <- test_data[i, c("HbA1c", "Treatment", top_snps), drop = FALSE]
    # Fit personalized model for this patient
    pm <- pmodel(forest, newdata = new_patient)
    coefs <- head(pm)
    # Find the best Treatment (lowest predicted HbA1c)
    # Assuming treatment is coded as factor with T1 as reference
    pred_T1 <- coefs[1, "(Intercept)"]
    pred_T2 <- coefs[1, "(Intercept)"] + coefs[1, "TreatmentT2"]
    pred_T3 <- coefs[1, "(Intercept)"] + coefs[1, "TreatmentT3"]
    pred_T4 <- coefs[1, "(Intercept)"] + coefs[1, "TreatmentT4"]
    
    best_trt <- paste("T", which.min(c(pred_T1, pred_T2, pred_T3, pred_T4)), sep="")
    # Extract info for each treatment
    out <- data.frame(
      patient = test_data$PatientID[i],
      pred_T1 = pred_T1,
      pred_T2 = pred_T2,
      pred_T3 = pred_T3,
      pred_T4 = pred_T4,
      best_treatment = best_trt,
      best_trt_prediction = min(c(pred_T1, pred_T2, pred_T3, pred_T4)),
      stringsAsFactors = FALSE
    )
    out
  })
  # Combine all results
  do.call(rbind, results)
}

#' Benchmark: Standard-of-Care and Random Assignment
#'
#' @param test_data Test set data frame
#' @param treatment_effects Named vector of treatment effects
#' @return Data frame with standard and random assignment outcomes
#' @export
benchmark_treatment_selection <- function(test_data, treatment_effects) {
  set.seed(123)
  n <- nrow(test_data)
  # Standard-of-care: always assign T1 (or best on average)
  soc_trt <- names(which.min(treatment_effects))
  
  soc_outcomes <- test_data$HbA1c -
    as.numeric(treatment_effects[as.character(test_data$Treatment)]) +
    as.numeric(treatment_effects[soc_trt])
  
  # Random assignment
  rand_trt <- sample(names(treatment_effects), n, replace = TRUE)
  
  rand_outcomes <- test_data$HbA1c -
    as.numeric(treatment_effects[as.character(test_data$Treatment)]) +
    as.numeric(treatment_effects[rand_trt])
  
  data.frame(
    patient = test_data$PatientID,
    soc_treatment = soc_trt,
    soc_outcome = soc_outcomes,
    rand_treatment = rand_trt,
    rand_outcome = rand_outcomes
  )
}

#' Full Workflow with Validation and Benchmarking
#'
#' @param n_patients Number of patients
#' @param n_snps Number of SNPs
#' @param n_top_snps Number of top SNPs to select
#' @param train_frac Fraction of data for training
#' @param ... Additional parameters for simulation
#' @return List of results
#' @export
run_full_workflow <- function(
    n_patients = 1000,
    n_snps = 200,
    n_top_snps = 20,
    train_frac = 0.8,
    treatment_effects = c(T1=0, T2=-0.3, T3=-0.6, T4=-0.9),
    n_causal_snps = 20,
    n_interacting_snps = 10,
    n_epistatic_pairs = 5,
    ntree = 100
) {
  set.seed(123)
  SNP_matrix <- simulate_snp_matrix(n_patients, n_snps)
  geno_data <- simulate_patient_data(
    n_patients, n_snps, SNP_matrix,
    n_causal_snps = n_causal_snps,
    n_interacting_snps = n_interacting_snps,
    n_epistatic_pairs = n_epistatic_pairs,
    treatment_effects = treatment_effects
  )
  snp_cols <- grep("^SNP", names(geno_data), value = TRUE)
  interaction_pvals <- calculate_interaction_pvals(geno_data, snp_cols)
  interaction_pvals_sig <- names(which(interaction_pvals < 0.05)) # assuming alpha 0.05 for statistical significance
  
  if (length(interaction_pvals_sig) == 0) {
    stop("No SNPs with gene-treatment interaction p-value < 0.05 found. Please relax alpha > 0.05 or note that no analysis is possible because no significant gene-treatment interaction was identified.")
  }
  
  snp_cols_maf <- filter_snps_by_maf(geno_data, interaction_pvals_sig)
  
  if (length(snp_cols_maf) == 0) {
    stop("No SNPs identfiied with signficant gene-treatment interaction and reasonable MAF.")
  }
  
  interaction_pvals_maf <- interaction_pvals[snp_cols_maf]
  
  snp_cols_ld <- prune_ld(geno_data, snp_cols_maf)
  interaction_pvals_ld <- interaction_pvals_maf[snp_cols_ld]
  top_snps <- names(sort(interaction_pvals_ld))[1:min(n_top_snps, length(interaction_pvals_ld))]
  
  # Split data
  train_idx <- sample(1:n_patients, size = floor(train_frac * n_patients))
  train_data <- geno_data[train_idx, ]
  test_data <- geno_data[-train_idx, ]
  
  # Fit personalized model on training set
  mob_forest <- personalized_treatment_forest(train_data, top_snps, ntree = ntree)
  
  # Evaluate the forest model on test set
  eval_results <- evaluate_personalized_model(mob_forest, train_data, test_data)
  
  # Benchmark
  benchmark_results <- benchmark_treatment_selection(test_data, treatment_effects)
  
  #Personalized best treatment predictions
  best_treatment_pred <- Predict_best_treatment(mob_forest, test_data, top_snps)
  
  return(list(
    train_data = train_data,
    test_data = test_data,
    top_snps = top_snps,
    eval_results = eval_results,
    benchmark_results = benchmark_results,
    best_treatment_pred = best_treatment_pred
  ))
}

# Example usage:
results <- run_full_workflow(n_patients = 1000, n_snps = 200, n_top_snps = 20, ntree = 100)

head(results$eval_results) # how model performs on train and test sets
head(results$benchmark_results) 
head(results$best_treatment_pred)

# Compare overall clinical response when patients are assigned:
# standard care (best known treatment) vs. 
# random treatment assignment vs. 
# best treatment based on personalized model
mean(results$benchmark_results$soc_outcome)
mean(results$benchmark_results$rand_outcome)
mean(results$best_treatment_pred$best_trt_prediction)

