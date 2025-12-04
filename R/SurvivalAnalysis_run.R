#' Univariate Cox Proportional Hazards Survival Analysis
#'
#' @description
#' Performs univariate Cox proportional hazards analysis on survival data with
#' multiple factor levels. Creates Kaplan-Meier survival curves, cumulative
#' hazard plots, and diagnostic plots. Calculates hazard ratios with confidence
#' intervals comparing each factor level against specified reference groups.
#'
#' @param surv.table A data frame containing survival information with columns:
#'   \describe{
#'     \item{S_ID}{Character vector of sample identifiers}
#'     \item{time}{Numeric vector representing survival time (in months)}
#'     \item{event}{Binary numeric (0 or 1) indicating event occurrence}
#'     \item{factor}{Character/factor vector indicating group membership}
#'   }
#' @param factor.table A data frame containing factor comparison pairs with columns:
#'   \describe{
#'     \item{HR.factor}{Character vector of hazard ratio factor levels}
#'     \item{ref.factor}{Character vector of reference factor levels}
#'   }
#' @param main Character string for plot title. Default is NULL.
#' @param plot.format Character string specifying plot format ('pdf' or 'json').
#'   Default is 'pdf'. Affects censor point size in plots.
#'
#' @return A list containing:
#'   \describe{
#'     \item{coxph.full}{List of Cox proportional hazards model objects}
#'     \item{KM}{List of Kaplan-Meier survival plot objects}
#'     \item{CH}{List of cumulative hazard plot objects}
#'     \item{ZPH}{List of proportional hazards test plots (or NULL if NA)}
#'     \item{DIAG}{List of diagnostic plots (Schoenfeld residuals)}
#'     \item{stat.tab}{Data frame with statistical results including HR, p-values, 
#'                     confidence intervals, concordance index, AIC, BIC}
#'   }
#'
#' @details
#' The function performs the following steps:
#' 1. Validates input data structure and column names
#' 2. Determines reference groups based on frequency in factor.table
#' 3. Fits Cox proportional hazards models for each reference group
#' 4. Calculates hazard ratios, confidence intervals, and p-values
#' 5. Tests proportional hazards assumption using Schoenfeld residuals
#' 6. Creates Kaplan-Meier survival curves with risk tables
#' 7. Generates cumulative hazard plots with event tables
#' 8. Produces diagnostic plots for model assessment
#'
#' @note
#' - Requires packages: tidyverse, data.table, survival, survminer, broom
#' - All factors in surv.table must exist in factor.table
#' - Function handles multiple reference groups automatically
#' - If Cox model returns NA, ZPH plots will be NULL with warning note
#'
#' @examples
#' \dontrun{
#' # Prepare survival data
#' surv_data <- data.frame(
#'     S_ID = paste0("PAT", seq_len(100)),
#'     time = runif(100, 1, 60),
#'     event = rbinom(100, 1, 0.4),
#'     factor = sample(c("Low", "High"), 100, replace = TRUE)
#' )
#' 
#' # Define factor comparisons
#' factor_data <- data.frame(
#'     HR.factor = c("High"),
#'     ref.factor = c("Low")
#' )
#' 
#' # Run analysis
#' results <- plot_unicox_surv(
#'     surv.table = surv_data,
#'     factor.table = factor_data,
#'     main = "Overall Survival",
#'     plot.format = "pdf"
#' )
#' 
#' # View statistical results
#' print(results$stat.tab)
#' 
#' # Display Kaplan-Meier plot
#' print(results$KM$KM)
#' }
#'
#' @export
plot_unicox_surv <- function(surv.table, 
                             factor.table, 
                             main = NULL, 
                             plot.format = "pdf") {
    
    require(tidyverse)
    require(data.table)
    require(survival)
    require(survminer)
    require(broom)
    
    # Input validation
    if (ncol(surv.table) < 4) {
        stop("The surv.table must contain at least 4 columns.")
    }
    if (ncol(factor.table) < 2) {
        stop("The factor.table must contain at least 2 columns.")
    }
    if (colnames(surv.table)[1] != "S_ID") {
        stop("The first column name of surv.table must be 'S_ID'.")
    }
    if (colnames(surv.table)[2] != "time") {
        stop("The second column name of surv.table must be 'time'.")
    }
    if (colnames(surv.table)[3] != "event") {
        stop("The third column name of surv.table must be 'event'.")
    }
    if (colnames(surv.table)[4] != "factor") {
        stop("The fourth column name of surv.table must be 'factor'.")
    }
    if (colnames(factor.table)[1] != "HR.factor") {
        stop("The first column name of factor.table must be 'HR.factor'.")
    }
    if (colnames(factor.table)[2] != "ref.factor") {
        stop("The second column name of factor.table must be 'ref.factor'.")
    }
    
    # Check factor consistency
    unique_factors <- unique(c(factor.table$HR.factor, factor.table$ref.factor))
    if (length(unique(surv.table$factor)) != length(unique_factors)) {
        stop("Please check your 'factor' in surv.table are all in factor.table.")
    }
    if (length(setdiff(unique(factor.table$HR.factor), unique(surv.table$factor))) > 0) {
        stop("Please check your 'HR.factor' in factor.table is in your surv.table!")
    }
    if (length(setdiff(unique(factor.table$ref.factor), unique(surv.table$factor))) > 0) {
        stop("Please check your 'ref.factor' in factor.table is in your surv.table!")
    }
    
    # Determine reference groups by frequency
    ref.count <- factor.table %>% 
        group_by(ref.factor) %>% 
        summarise(N = n()) %>% 
        arrange(desc(N))
    
    reference <- unique(factor.table$ref.factor)[
        order(match(unique(factor.table$ref.factor), ref.count$ref.factor))
    ]
    
    # Perform Cox analysis for each reference group
    surv.result <- NULL
    for (i in seq_along(reference)) {
        
        REF <- reference[i]
        HRF <- setdiff(unique(surv.table$factor), REF)
        fac.level <- c(REF, HRF)
        
        nSample <- surv.table %>% 
            group_by(factor) %>% 
            summarise(N = n()) %>% 
            ungroup() %>% 
            arrange(match(factor, fac.level))
        
        surv.table$factor <- factor(surv.table$factor, levels = fac.level)
        
        # Fit Cox proportional hazards model
        scox <- coxph(Surv(time, event == 1) ~ factor, data = surv.table)
        assign(paste0("cox_", i), scox)
        
        # Summarize model statistics
        scox.tidy <- tidy(scox, exponentiate = TRUE, conf.int = TRUE) %>% 
            mutate(HR.factor = str_remove(string = term, pattern = "factor"), 
                   ref.factor = REF) %>% 
            dplyr::select(HR.factor, ref.factor, 
                          "HR" = "estimate", 
                          "cox.pval" = "p.value", 
                          "conf.low.95" = "conf.low", 
                          "conf.high.95" = "conf.high")
        
        scox.glance <- glance(scox)
        
        # Combine results
        surv.res <- scox.tidy %>% 
            mutate(c.index = scox.glance$concordance, 
                   log.rank.pval = scox.glance$p.value.sc, 
                   wald.pval = scox.glance$p.value.wald, 
                   likelihood.pval = scox.glance$p.value.log, 
                   r.square = scox.glance$r.squared, 
                   AIC = scox.glance$AIC, 
                   BIC = scox.glance$BIC, 
                   num.HR = nSample$N[-1], 
                   num.ref = nSample$N[1]) %>% 
            dplyr::select(HR.factor, ref.factor, log.rank.pval, HR, c.index, 
                          num.HR, num.ref, cox.pval, conf.low.95, conf.high.95, 
                          wald.pval, likelihood.pval, r.square, AIC, BIC)
        
        surv.result <- rbindlist(
            list(surv.result, surv.res), 
            use.names = TRUE, 
            fill = TRUE
        )
    }
    
    # Merge with factor table
    stat.res <- factor.table %>% 
        left_join(surv.result, by = c("HR.factor", "ref.factor")) %>% 
        mutate(note = NA)
    
    # Diagnostic plots
    DIAG <- ggcoxdiagnostics(cox_1, type = "schoenfeld")
    
    # Proportional hazards test
    if (all(!is.na(tidy(cox_1, exponentiate = TRUE, conf.int = TRUE)$estimate))) {
        cox.zph.fit <- cox.zph(cox_1)
        ZPH <- ggcoxzph(
            cox.zph.fit, 
            resid = TRUE, 
            se = TRUE,
            ylab = "Beta(t) for factor", 
            xlab = "Months"
        )
    } else {
        ZPH <- NULL
        stat.res$note <- "coxph() returns NA, so cannot plot the ZPH graph."
    }
    
    # Create survival curves
    surv.table$factor <- as.character(surv.table$factor)
    mfit <- surv_fit(Surv(time, event == 1) ~ factor, data = surv.table)
    
    # Prepare plot labels
    legend.labs <- str_remove(names(mfit$strata), pattern = "factor=")
    main.HR <- paste(
        paste0("\nHR of ", stat.res$HR.factor, "/", stat.res$ref.factor, 
               " = ", format(stat.res$HR, digits = 3)), 
        collapse = ""
    )
    
    # Set censor size based on plot format
    censor_size <- ifelse(plot.format == "pdf", 4, 1.5)
    
    # Kaplan-Meier survival plot
    KM <- ggsurvplot(
        mfit,  
        color = "strata", 
        palette = c("#EA4335", "#4285F4", "#43A047", "#FDD835", 
                    "#8E24AA", "#6D4C41", "#00897B"), 
        linetype = 1, 
        conf.int = FALSE,
        censor = TRUE, 
        censor.shape = "+", 
        censor.size = censor_size,
        pval = FALSE, 
        pval.size = 5,
        surv.median.line = "none",
        risk.table = TRUE, 
        tables.height = 0.25, 
        tables.y.text = FALSE, 
        tables.col = "strata",
        ncensor.plot = FALSE, 
        ncensor.plot.height = 0.25,
        xlab = "Months", 
        ylab = "Survival Probability",
        title = paste0(main, "\nLog-Rank P value = ", 
                       format(unique(stat.res$log.rank.pval), digits = 3), 
                       main.HR),
        legend = "right", 
        legend.title = "", 
        legend.labs = legend.labs
    )
    
    # Cumulative hazard plot
    CH <- ggsurvplot(
        mfit, 
        fun = "cumhaz",
        color = "strata", 
        palette = c("#EA4335", "#4285F4", "#43A047", "#FDD835", 
                    "#8E24AA", "#6D4C41", "#00897B"), 
        linetype = 1, 
        conf.int = FALSE,
        censor = TRUE, 
        censor.shape = "+", 
        censor.size = censor_size,
        pval = FALSE, 
        pval.size = 5,
        surv.median.line = "none",
        cumevents = TRUE, 
        tables.height = 0.25, 
        tables.y.text = FALSE, 
        tables.col = "strata",
        ncensor.plot = FALSE, 
        ncensor.plot.height = 0.25,
        xlab = "Months", 
        ylab = "Cumulative Hazard",
        title = paste0(main, "\nLog-Rank P value = ", 
                       format(unique(stat.res$log.rank.pval), digits = 3), 
                       main.HR),
        legend = "right", 
        legend.title = "", 
        legend.labs = legend.labs
    )
    
    # Return results
    return(list(
        coxph.full = mget(ls(pattern = "^cox_"), envir = environment()),
        KM = mget(ls(pattern = "^KM"), envir = environment()),
        CH = mget(ls(pattern = "^CH"), envir = environment()),
        ZPH = mget(ls(pattern = "^ZPH"), envir = environment()),
        DIAG = mget(ls(pattern = "^DIAG"), envir = environment()),
        stat.tab = stat.res
    ))
}