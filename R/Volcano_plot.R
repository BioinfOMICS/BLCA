#' Volcano Plot
#'
#' @description
#' Creates a volcano plot displays coefficient estimates
#' versus statistical significance, with optional text labels for significant
#' associations. Coefficients are inverted so that High > Low associations
#' appear on the right side.
#'
#' @param res.df A data frame containing association analysis results with columns:
#'   \describe{
#'     \item{TIL}{Character vector of names}
#'     \item{coef or Estimate}{Numeric vector of coefficient estimates}
#'     \item{p_val or Pr(>|z|)}{Numeric vector of p-values}
#'   }
#' @param p.threshold Numeric value for p-value significance threshold. 
#'   Default is 0.05.
#' @param coef.threshold Numeric value for coefficient significance threshold. 
#'   Default is 0 (any non-zero coefficient).
#' @param add.labels Logical indicating whether to add text labels for 
#'   significant points. Default is TRUE.
#' @param down.color Character string for "High > Low" color. 
#'   Default is "#B2182B" (red).
#' @param up.color Character string for "Low > High" color. 
#'   Default is "#2166AC" (blue).
#' @param ns.color Character string for non-significant points color. 
#'   Default is "grey70".
#' @param point.size Numeric value for point size. Default is 3.
#' @param point.alpha Numeric value (0-1) for point transparency. Default is 0.7.
#' @param label.size Numeric value for label text size. Default is 3.
#' @param max.overlaps Maximum number of overlapping labels. Default is 30.
#' @param plot.title Character string for plot title. 
#'   Default is "Volcano Plot: TIL Association Analysis".
#'
#' @return A ggplot2 object displaying the volcano plot.
#'
#' @details
#' The function performs the following steps:
#' 1. Standardizes column names (handles both 'coef'/'Estimate' and 'p_val'/'Pr(>|z|)')
#' 2. Inverts coefficient signs so High > Low appears on right (positive)
#' 3. Filters data to include only significant associations (p < p.threshold)
#' 4. Calculates -log10(p-value) for y-axis
#' 5. Classifies points as "Up" (Low > High), "Down" (High > Low), or "NS"
#' 6. Creates volcano plot with threshold lines
#' 7. Optionally adds text labels using ggrepel
#'
#' @note
#' - Requires packages: dplyr, ggplot2, ggrepel (if add.labels = TRUE)
#' - Coefficients are automatically inverted (multiplied by -1)
#' - Only points with p < p.threshold are included in the plot
#' - "Down" (red) represents High risk > Low risk associations
#' - "Up" (blue) represents Low risk > High risk associations
#'
#' @examples
#' \dontrun{
#' # Prepare association results
#' til_results <- data.frame(
#'     TIL = c("CD8_T", "CD4_T", "NK", "Macrophage", "B_cell"),
#'     Estimate = c(0.5, -0.3, 0.8, -0.2, 0.1),
#'     `Pr(>|z|)` = c(0.001, 0.02, 0.0001, 0.1, 0.5),
#'     check.names = FALSE
#' )
#' 
#' # Create default volcano plot
#' plot_til_volcano(til_results)
#' 
#' # Custom thresholds without labels
#' plot_til_volcano(
#'     res.df = til_results,
#'     p.threshold = 0.01,
#'     coef.threshold = 0.3,
#'     add.labels = FALSE
#' )
#' 
#' # Custom colors
#' plot_til_volcano(
#'     res.df = til_results,
#'     down.color = "#E63946",
#'     up.color = "#1D3557",
#'     point.size = 4
#' )
#' }
#'
#' @export
plot_til_volcano <- function(res.df,
                             p.threshold = 0.05,
                             coef.threshold = 0,
                             add.labels = TRUE,
                             down.color = "#B2182B",
                             up.color = "#2166AC",
                             ns.color = "grey70",
                             point.size = 3,
                             point.alpha = 0.7,
                             label.size = 3,
                             max.overlaps = 30,
                             plot.title = "Volcano Plot: TIL Association Analysis") {
    
    require(dplyr)
    require(ggplot2)
    if (add.labels) {
        require(ggrepel)
    }
    
    # Standardize column names
    if ("Estimate" %in% colnames(res.df)) {
        colnames(res.df)[colnames(res.df) == "Estimate"] <- "coef"
    }
    if ("Pr(>|z|)" %in% colnames(res.df)) {
        colnames(res.df)[colnames(res.df) == "Pr(>|z|)"] <- "p_val"
    }
    
    # Validate required columns
    required_cols <- c("TIL", "coef", "p_val")
    if (!all(required_cols %in% colnames(res.df))) {
        stop(paste("res.df must contain columns:", 
                   paste(required_cols, collapse = ", ")))
    }
    
    # Convert to data frame if list
    if (is.list(res.df$TIL)) {
        res.df$TIL <- unlist(res.df$TIL)
    }
    if (is.list(res.df$coef)) {
        res.df$coef <- unlist(res.df$coef)
    }
    if (is.list(res.df$p_val)) {
        res.df$p_val <- unlist(res.df$p_val)
    }
    
    # Invert coefficient so High > Low is on right side
    res.df$coef <- -res.df$coef
    
    # Filter for significant results
    til_data <- res.df %>% 
        dplyr::filter(p_val < p.threshold)
    
    if (nrow(til_data) == 0) {
        warning(paste("No significant associations found with p <", p.threshold))
        return(NULL)
    }
    
    # Calculate -log10(p-value)
    til_data$neg_log10_p <- -log10(til_data$p_val)
    
    # Classify significance
    til_data$significance <- ifelse(
        til_data$p_val < p.threshold & abs(til_data$coef) > coef.threshold,
        ifelse(til_data$coef > 0, "Up", "Down"),
        "NS"
    )
    
    # Mark significant points for labeling
    til_data$is_significant <- til_data$p_val < p.threshold
    
    # Create base volcano plot
    p <- ggplot(til_data, aes(x = coef, y = neg_log10_p)) +
        geom_point(
            aes(color = significance),
            size = point.size,
            alpha = point.alpha
        ) +
        scale_color_manual(
            values = c(
                "Down" = down.color, 
                "NS" = ns.color, 
                "Up" = up.color
            ),
            labels = c(
                "Down" = "High > Low", 
                "NS" = "Not Sig", 
                "Up" = "Low > High"
            )
        ) +
        geom_hline(
            yintercept = -log10(p.threshold),
            linetype = "dashed",
            color = "grey30"
        ) +
        geom_vline(
            xintercept = c(-coef.threshold, coef.threshold),
            linetype = "dashed",
            color = "grey30"
        ) +
        labs(
            title = plot.title,
            x = "Coefficient",
            y = "-log10(P-value)",
            color = "Significance"
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            legend.position = "right"
        )
    
    # Add text labels if requested
    if (add.labels) {
        p <- p +
            ggrepel::geom_text_repel(
                data = subset(til_data, is_significant == TRUE),
                aes(label = TIL),
                size = label.size,
                max.overlaps = max.overlaps,
                box.padding = 0.5,
                min.segment.length = 0
            )
    }
    
    return(p)
}


# ============================================================================
# Simulation and Testing
# ============================================================================

#' Generate Simulated TIL Association Data
#'
#' @description
#' Creates simulated TIL association analysis results for testing.
#'
#' @param n.til Number of TIL cell types. Default is 20.
#' @param prop.sig Proportion of significant associations. Default is 0.4.
#' @param seed Random seed for reproducibility. Default is 123.
#'
#' @return A data frame with TIL, Estimate, and Pr(>|z|) columns.
#'
#' @export
simulate_til_data <- function(n.til = 20, 
                              prop.sig = 0.4, 
                              seed = 123) {
    
    set.seed(seed)
    
    n.sig <- round(n.til * prop.sig)
    n.nonsig <- n.til - n.sig
    
    # Generate TIL names
    til_names <- c(
        "CD8_T_cell", "CD4_T_cell", "NK_cell", "B_cell", 
        "Macrophage_M1", "Macrophage_M2", "Dendritic_cell",
        "Neutrophil", "Monocyte", "Treg", "Th1", "Th2", 
        "Plasma_cell", "Mast_cell", "Eosinophil",
        paste0("TIL_", seq_len(n.til - 15))
    )[seq_len(n.til)]
    
    # Generate coefficients
    coef_sig <- rnorm(n.sig, mean = 0, sd = 1)
    coef_nonsig <- rnorm(n.nonsig, mean = 0, sd = 0.3)
    
    # Generate p-values
    p_sig <- runif(n.sig, min = 0.0001, max = 0.04)
    p_nonsig <- runif(n.nonsig, min = 0.05, max = 0.9)
    
    # Combine data
    til_results <- data.frame(
        TIL = til_names,
        Estimate = c(coef_sig, coef_nonsig),
        `Pr(>|z|)` = c(p_sig, p_nonsig),
        check.names = FALSE
    )
    
    # Shuffle rows
    til_results <- til_results[sample(seq_len(nrow(til_results))), ]
    rownames(til_results) <- NULL
    
    return(til_results)
}


# ============================================================================
# Example Usage
# ============================================================================

cat("=== Generating Simulated TIL Association Data ===\n")

# Create simulated data
til_results <- simulate_til_data(n.til = 25, prop.sig = 0.5)

cat("Total TIL types:", nrow(til_results), "\n")
cat("Significant associations (p < 0.05):", 
    sum(til_results$`Pr(>|z|)` < 0.05), "\n\n")

cat("=== Top 10 TIL Associations ===\n")
print(head(til_results[order(til_results$`Pr(>|z|)`), ], 10))

cat("\n\n=== Creating Volcano Plots ===\n\n")

# Example 1: Default plot with labels
cat("Example 1: Default volcano plot with labels\n")
p1 <- plot_til_volcano(til_results)
print(p1)

cat("\n")

# Example 2: Stricter thresholds
cat("Example 2: Stricter thresholds (p < 0.01, |coef| > 0.5)\n")
p2 <- plot_til_volcano(
    res.df = til_results,
    p.threshold = 0.01,
    coef.threshold = 0.5,
    plot.title = "Volcano Plot: Strict Threshold"
)
print(p2)

cat("\n")

# Example 3: Without labels, custom colors
cat("Example 3: No labels with custom colors\n")
p3 <- plot_til_volcano(
    res.df = til_results,
    add.labels = FALSE,
    down.color = "#E63946",
    up.color = "#1D3557",
    ns.color = "grey80",
    point.size = 4,
    plot.title = "Volcano Plot: Custom Styling"
)
print(p3)

cat("\n✓ All volcano plots created successfully!\n")