# Load required libraries
library(dplyr)
library(data.table)
library(randomForestSRC)

#' Random Forest Variable Selection for Survival Analysis
#'
#' @description
#' Performs Random Forest variable selection on survival data with gene 
#' expression values. The function integrates survival information with 
#' gene expression data, removes incomplete cases, and uses variable hunting 
#' (vh) method to select the most important genes for survival prediction.
#'
#' @param s.data A data frame containing survival information with three columns:
#'   \describe{
#'     \item{s_id}{Character vector of sample identifiers}
#'     \item{event}{Binary numeric (0 or 1) indicating event occurrence}
#'     \item{months}{Numeric vector (>0) representing survival time in months}
#'   }
#' @param value.data A data frame containing gene expression values with columns:
#'   \describe{
#'     \item{s_id}{Character vector of sample identifiers (must match s.data)}
#'     \item{gene1, gene2, ..., geneN}{Numeric vectors of gene expression values}
#'   }
#' @param outPATH Character string specifying output directory path. 
#'   Default is current working directory from getwd().
#'
#' @return A list object (vh.sur.cox) containing:
#'   \describe{
#'     \item{topvars}{Character vector of top selected gene names}
#'     \item{modelsize}{Integer indicating number of genes selected}
#'     \item{varselect}{Data frame with variable selection statistics}
#'     \item{RF.gene.tsv / can.term.df}{Data frame with term, depth, and rel.freq}
#'     \item{select_genes}{Character vector of selected gene names}
#'     \item{value.data.s}{Merged data frame of survival data and selected genes}
#'   }
#'
#' @export
RF.run <- function(s.data, value.data, outPATH = getwd()) { 
    # Integrate survival and expression data
    sur.sam <- s.data %>% 
        dplyr::left_join(value.data, by = "s_id")
    
    # Check for duplicated patient IDs
    if (sum(duplicated(sur.sam$s_id)) != 0) {
        stop("Duplicated patient s_id detected. Please check your data.")
    }
    
    # Remove incomplete cases
    sur.comp <- sur.sam[complete.cases(sur.sam), ]
    
    # Extract value matrix without s_id
    sam.df <- sur.comp %>% 
        dplyr::select(2:length(.))
    
    # First fit the random forest model
    fit.f.genes <- rfsrc(
        Surv(months, event) ~ ., 
        data = sam.df,
        ntree = 100,
        importance = TRUE
    )
    
    # Random Forest variable selection using variable hunting method
    vh.sur.cox <- var.select(
        fit.f.genes, 
        method = "vh", 
        nstep = 10, 
        nrep = 50
    )
    
    # Extract top variables and model size
    topvars <- vh.sur.cox$topvars
    modelsize <- vh.sur.cox$modelsize
    selected_var <- data.frame(vh.sur.cox$varselect)
    
    # Create summary data frame of selected genes
    vh.sur.cox$RF.gene.tsv <- vh.sur.cox$can.term.df <- data.frame(
        term = topvars,
        depth = selected_var$depth[seq_len(modelsize)],
        rel.freq = selected_var$rel.freq[seq_len(modelsize)]
    )
    
    # Filter gene names present in value.data
    vh.sur.cox$select_genes <- names(value.data)[
        names(value.data) %in% topvars
    ]
    
    # Create filtered dataset with selected genes only
    value.data.s <- value.data %>% 
        dplyr::select(s_id, all_of(vh.sur.cox$select_genes))
    vh.sur.cox$value.data.s <- merge(s.data, value.data.s)
    
    # Write results to file
    fwrite(
        vh.sur.cox$can.term.df, 
        file = file.path(outPATH, "RF.gene.tsv"), 
        sep = '\t', 
        col.names = TRUE, 
        row.names = FALSE
    )
    
    return(vh.sur.cox)
}

# ============================================================================
# A Simulated Dataset as the Test Run
# ============================================================================

set.seed(12345)  # For reproducibility

# Parameters
n_samples <- 100
n_genes <- 200

# Generate sample IDs
sample_ids <- sprintf("SAMPLE-%03d", seq_len(n_samples))

# Generate survival data
# Event: 0 (censored) or 1 (event occurred)
# Months: survival time (positive numeric)
survival_data <- data.frame(
    s_id = sample_ids,
    event = rbinom(n_samples, size = 1, prob = 0.4),  # 40% event rate
    months = runif(n_samples, min = 0.5, max = 60)    # 0.5 to 60 months
)

# Generate gene expression data
# Create gene names: GENE001, GENE002, ..., GENE200
gene_names <- sprintf("GENE%03d", seq_len(n_genes))

# Generate expression matrix (normalized values around 0)
# Using normal distribution with mean=0, sd=1
expression_matrix <- matrix(
    rnorm(n_samples * n_genes, mean = 0, sd = 1),
    nrow = n_samples,
    ncol = n_genes
)

# Convert to data frame and add sample IDs
expression_data <- as.data.frame(expression_matrix)
colnames(expression_data) <- gene_names
expression_data <- cbind(
    data.frame(s_id = sample_ids),
    expression_data
)

# ============================================================================
# Display sample data
# ============================================================================

cat("=== Survival Data (first 10 rows) ===\n")
print(head(survival_data, 10))

cat("\n=== Expression Data (first 6 rows, first 6 genes) ===\n")
print(survival_data[1:6, ])
print(expression_data[1:6, 1:7])

cat("\n=== Data Dimensions ===\n")
cat("Survival data:", nrow(survival_data), "samples\n")
cat("Expression data:", nrow(expression_data), "samples x", 
    ncol(expression_data) - 1, "genes\n")

cat("\n=== Event Distribution ===\n")
print(table(survival_data$event))

cat("\n=== Survival Time Summary ===\n")
print(summary(survival_data$months))

# ============================================================================
# Run RF.run function
# ============================================================================

cat("\n\n=== Running RF.run function ===\n")

# Create output directory if it doesn't exist
output_dir <- file.path(getwd(), "RF_output")
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# Run the function
results <- RF.run(
    s.data = survival_data,
    value.data = expression_data,
    outPATH = output_dir
)

# Display results
cat("\n=== RF.run Results ===\n")
cat("Number of genes selected:", results$modelsize, "\n")
cat("Selected genes:", paste(results$select_genes, collapse = ", "), "\n")

cat("\n=== Top Selected Genes (first 10) ===\n")
print(head(results$can.term.df, 10))

cat("\n=== Output file created ===\n")
cat("File location:", file.path(output_dir, "RF.gene.tsv"), "\n")

cat("\n=== Filtered Dataset Dimensions ===\n")
cat("Final dataset:", nrow(results$value.data.s), "samples x", 
    ncol(results$value.data.s) - 3, "selected genes\n")

cat("\n✓ Function completed successfully!\n")