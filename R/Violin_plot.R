#' Violin Plot
#'
#' @description
#' Creates a violin plot with overlaid jitter points and boxplot to visualize
#' the distribution of different status groups (Low vs High).
#' Combines violin plot density, individual data points, and summary statistics
#' in a single comprehensive visualization.
#'
#' @param RS.signature A data frame containing risk score information with columns:
#'   \describe{
#'     \item{S_ID}{Character vector of sample identifiers}
#'     \item{value}{Numeric vector of risk score values}
#'     \item{status}{Binary numeric (0 or 1) indicating risk status. 
#'                   0 = Low Risk, 1 = High Risk}
#'   }
#' @param dataset.name Character string specifying the name of the dataset for
#'   plot title. Default is NULL.
#' @param low.color Character string for low risk group color. Default is "#3498db" (blue).
#' @param high.color Character string for high risk group color. Default is "#e74c3c" (red).
#' @param jitter.width Numeric value for jitter width. Default is 0.1.
#' @param jitter.alpha Numeric value (0-1) for jitter point transparency. Default is 0.3.
#' @param jitter.size Numeric value for jitter point size. Default is 1.
#' @param violin.alpha Numeric value (0-1) for violin plot transparency. Default is 0.7.
#'
#' @return A ggplot2 object displaying the risk score distribution violin plot.
#'
#' @details
#' The function creates a layered visualization with three components:
#' 1. Violin plot showing the density distribution of risk scores
#' 2. Jittered points showing individual sample values
#' 3. Box plot overlay showing median, quartiles, and outliers
#' 
#' The risk status is automatically converted from binary (0/1) to labeled
#' factors ("Low Risk"/"High Risk") for better interpretability.
#'
#' @note
#' - Requires packages: dplyr, ggplot2
#' - Status must be binary (0 or 1) where 0 = Low Risk, 1 = High Risk
#' - The function assumes S_ID, value, and status columns exist in input data
#' - Box plot width is fixed at 0.1 to avoid obscuring violin plot
#'
#' @examples
#' \dontrun{
#' # Prepare risk score data
#' risk_data <- data.frame(
#'     S_ID = paste0("SAMPLE", seq_len(100)),
#'     value = c(rnorm(50, mean = 2, sd = 0.5),
#'               rnorm(50, mean = 4, sd = 0.8)),
#'     status = rep(c(0, 1), each = 50)
#' )
#' 
#' # Create default plot
#' plot_risk_violin(risk_data, dataset.name = "TCGA-BRCA")
#' 
#' # Customize colors and transparency
#' plot_risk_violin(
#'     RS.signature = risk_data,
#'     dataset.name = "Example Violin",
#'     low.color = "#2ecc71",
#'     high.color = "#e67e22",
#'     jitter.alpha = 0.5,
#'     violin.alpha = 0.5
#' )
#' }
#'
#' @export
plot_risk_violin <- function(RS.signature, 
                             dataset.name = NULL,
                             low.color = "#3498db",
                             high.color = "#e74c3c",
                             jitter.width = 0.1,
                             jitter.alpha = 0.3,
                             jitter.size = 1,
                             violin.alpha = 0.7) {
    
    require(dplyr)
    require(ggplot2)
    
    # Validate input
    if (!all(c("S_ID", "value", "status") %in% colnames(RS.signature))) {
        stop("RS.signature must contain columns: S_ID, value, status")
    }
    
    if (!all(RS.signature$status %in% c(0, 1))) {
        stop("status must be binary (0 or 1)")
    }
    
    # Prepare data
    df <- RS.signature %>% 
        dplyr::select(S_ID, value, status)
    
    # Convert status to factor with labels
    df$status <- factor(
        df$status, 
        levels = c(0, 1), 
        labels = c("Low Risk", "High Risk")
    )
    
    # Create title
    if (is.null(dataset.name)) {
        plot_title <- "Value Distribution by Risk Status"
    } else {
        plot_title <- paste0("Value Distribution by Risk Status (", 
                             dataset.name, ")")
    }
    
    # Create violin plot
    p <- ggplot(df, aes(x = status, y = value, fill = status)) +
        geom_violin(trim = FALSE, alpha = violin.alpha) +
        geom_jitter(width = jitter.width, 
                    alpha = jitter.alpha, 
                    size = jitter.size) +
        geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
        scale_fill_manual(
            values = c("Low Risk" = low.color, "High Risk" = high.color)
        ) +
        labs(
            title = plot_title,
            x = "Risk Status",
            y = "Risk Score",
            fill = "Risk Status"
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10),
            legend.position = "right"
        )
    
    return(p)
}


# ============================================================================
# Simulation Example and Testing
# ============================================================================

#' Generate Simulated Risk Score Data
#'
#' @description
#' Creates a simulated dataset with risk scores for testing plot_risk_violin().
#'
#' @param n.samples Total number of samples. Default is 100.
#' @param mean.low Mean risk score for low risk group. Default is 2.
#' @param sd.low Standard deviation for low risk group. Default is 0.5.
#' @param mean.high Mean risk score for high risk group. Default is 4.
#' @param sd.high Standard deviation for high risk group. Default is 0.8.
#' @param prop.high Proportion of high risk samples (0-1). Default is 0.5.
#'
#' @return A data frame with S_ID, value, and status columns.
#'
#' @export
simulate_risk_data <- function(n.samples = 100,
                               mean.low = 2,
                               sd.low = 0.5,
                               mean.high = 4,
                               sd.high = 0.8,
                               prop.high = 0.5) {
    
    set.seed(123)  # For reproducibility
    
    n.high <- round(n.samples * prop.high)
    n.low <- n.samples - n.high
    
    risk_data <- data.frame(
        S_ID = sprintf("SAMPLE-%03d", seq_len(n.samples)),
        value = c(
            rnorm(n.low, mean = mean.low, sd = sd.low),
            rnorm(n.high, mean = mean.high, sd = sd.high)
        ),
        status = rep(c(0, 1), c(n.low, n.high))
    )
    
    return(risk_data)
}


# ============================================================================
# Example Usage
# ============================================================================

cat("=== Generating Simulated Risk Score Data ===\n")

# Create simulated data with 150 samples
risk_signature <- simulate_risk_data(
    n.samples = 150,
    mean.low = 1.8,
    sd.low = 0.6,
    mean.high = 4.2,
    sd.high = 1.0,
    prop.high = 0.4
)

cat("Data dimensions:", nrow(risk_signature), "samples\n")
cat("Risk status distribution:\n")
print(table(risk_signature$status))

cat("\nRisk score summary by status:\n")
print(
    risk_signature %>%
        group_by(status) %>%
        summarise(
            n = n(),
            mean = mean(value),
            sd = sd(value),
            min = min(value),
            max = max(value)
        )
)

cat("\n=== Creating Violin Plots ===\n\n")

# Example 1: Default plot
cat("Example 1: Default settings\n")
p1 <- plot_risk_violin(
    RS.signature = risk_signature,
    dataset.name = "Example 1"
)
print(p1)

cat("\n")

# Example 2: Custom colors
cat("Example 2: Custom colors (green/orange)\n")
p2 <- plot_risk_violin(
    RS.signature = risk_signature,
    dataset.name = "Example 2",
    low.color = "#2ecc71",
    high.color = "#e67e22"
)
print(p2)

cat("\n")

# Example 3: More visible points
cat("Example 3: Increased point visibility\n")
p3 <- plot_risk_violin(
    RS.signature = risk_signature,
    dataset.name = "Example 3",
    jitter.alpha = 0.6,
    jitter.size = 1.5,
    violin.alpha = 0.5
)
print(p3)

cat("\n✓ All plots created successfully!\n")