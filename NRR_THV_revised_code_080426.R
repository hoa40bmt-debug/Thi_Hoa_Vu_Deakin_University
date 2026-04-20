####################################################################################
#
# Nature-related risk perception of smallholder farmers and implications for global agricultural supply chain sustainability
#    By: Hoa Vu
#    Created: July 2025
#    Email: h.vu@deakin.edu.au
#    Accompanied data: 
#                       Clean_data_V2.csv - Survey data
#                       RS_variable_matrix_selected_v2.csv - a CSV table with matrix of selected predictor variable
#
#    Description: This script statistically analyses Tea Farmers’ Perceived Risk And Vulnerability To Nature-Related Risks In Northern Vietnam.
#                 This script includes 3 parts:
#                     Part I: Preparation - loading R packages, functions and survey data
#                     Part II: Tea farmers’ perspectives on ecosystem services, nature-related risk experiences, and future risk perceptions
#                     Part III: Factors influencing farmers’ nature-related risk perception
#
####################################################################################

rm(list = ls())
gc()

####################################################################################
# PART I: PREPARATION ----
####################################################################################

# 1. Working directory and data loading ----
setwd("C:/Users/s223118632/OneDrive - Deakin University/PhD/Paper 1/PP1_Data/Final_code2")

# Survey data
data <- read.csv("Clean_data_riskscore_V2.csv")

# Output folder for all revised results
output_dir <- "C:/Users/s223118632/OneDrive - Deakin University/PhD/Paper 1/Results/FDR_Check"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Loading required packages ----
library(tidyverse)
library(ggplot2)
library(likert)
library(ggpubr)
library(RColorBrewer)
library(ggridges)
library(corrplot)
library(MASS)
library(rsample)
library(rlang)
library(parsnip)
library(recipes)
library(yardstick)
library(ggh4x)
library(tidyr)
library(forcats)
library(patchwork)
library(moments)

# 3. Customised functions ----

## 3.1. convert_to_factor ----
convert_to_factor <- function(data, cols, levels, labels = NULL) {
  existing_cols <- cols[cols %in% colnames(data)]
  missing_cols <- setdiff(cols, existing_cols)
  
  if (length(missing_cols) > 0) {
    warning("Missing columns not converted: ", paste(missing_cols, collapse = ", "))
  }
  
  for (col in existing_cols) {
    if (is.null(labels)) {
      data[[col]] <- factor(data[[col]], levels = levels)
    } else {
      data[[col]] <- factor(data[[col]], levels = levels, labels = labels)
    }
  }
  return(data)
}

## 3.2. recipe_simple ----
recipe_simple <- function(dataset, fml) {
  recipe(fml, data = dataset) %>%
    step_string2factor(all_nominal(), -all_outcomes()) %>%
    prep(data = dataset)
}

## 3.3. risk_model ----
# Revised for journal submission:
# - full dataset only
# - backward stepwise selection
# - no train/test split
risk_model <- function(numeric_data, response, predictor, title, stepwise = TRUE, k = 3.5){
  print(paste("Running model for:", title))
  print("Creating model data")
  
  model_data <- numeric_data %>%
    dplyr::select(all_of(c(response, predictor$Predictor))) %>%
    na.omit()
  
  fml <- as.formula(paste(response, "~ ."))
  lm_model <- lm(fml, data = model_data)
  
  if(stepwise){
    print("Variable selection using backward stepwise AIC")
    
    step_model <- stepAIC(lm_model, direction = "backward", k = k, trace = FALSE)
    step_summary <- summary(step_model)
    step_coefs <- step_summary$coefficients
    
    selected_vars <- rownames(step_coefs)[rownames(step_coefs) != "(Intercept)"]
    
    if(length(selected_vars) == 0) {
      warning(paste("No predictors retained for", response))
      return(NULL)
    }
    
    print(paste("Number of retained predictors:", length(selected_vars)))
    print(selected_vars)
    
    model_data <- model_data[, c(response, selected_vars), drop = FALSE]
    fml <- as.formula(paste(response, "~", paste(selected_vars, collapse = " + ")))
  }
  
  print("Running linear regression model on full dataset")
  
  recipe_prepped <- recipe_simple(dataset = model_data, fml)
  model_baked <- bake(recipe_prepped, new_data = model_data)
  
  glm_fit_obj <- linear_reg() %>%
    set_engine("lm") %>%
    fit(fml, data = model_baked)
  
  print("Producing outputs")
  
  ml_summary <- summary(glm_fit_obj$fit)
  
  # Coefficient table
  coefs <- as.data.frame(ml_summary$coefficients)
  coefs$RS <- title
  coefs$vars <- rownames(coefs)
  
  vars_order <- predictor$order[match(coefs$vars[coefs$vars != "(Intercept)"], predictor$Predictor)]
  coefs$order <- c(0, vars_order)
  
  vars_name <- predictor$Predictor_name[match(coefs$vars[coefs$vars != "(Intercept)"], predictor$Predictor)]
  coefs$vars_name <- c("Intercept", vars_name)
  
  coefs <- coefs[, c(5, 6, 7, 8, 1, 2, 3, 4)]
  names(coefs)[5:8] <- c("Estimate", "Std. Error", "t value", "p")
  
  # Model performance
  glm_perf <- data.frame(
    RS = title,
    n = nrow(model_baked),
    r_squared = round(ml_summary$r.squared, 3),
    adj_r_squared = round(ml_summary$adj.r.squared, 3)
  )
  
  return(list(
    glm_fit = coefs,
    glm_perf = glm_perf
  ))
}

## 3.4. risk_model_batch ----
risk_model_batch <- function(numeric_data, predictor_tbl, responses, titles, stepwise = TRUE){
  glm_fit_list <- list()
  glm_perf_list <- list()
  
  for(i in 1:length(responses)){
    response <- responses[i]
    print(response)
    
    predictor <- predictor_tbl[predictor_tbl[[response]] == 1, c("order", "Predictor", "Predictor_name")] %>%
      na.omit()
    
    title <- titles[i]
    
    model_result <- risk_model(
      numeric_data = numeric_data,
      response = response,
      predictor = predictor,
      title = title,
      stepwise = stepwise,
      k = 3.5
    )
    
    if(!is.null(model_result)){
      glm_fit_list[[title]] <- model_result[["glm_fit"]]
      glm_perf_list[[title]] <- model_result[["glm_perf"]]
    }
  }
  
  return(list(
    glm_fit_list = glm_fit_list,
    glm_perf_list = glm_perf_list
  ))
}

## 3.5. risk_model_outputs ----
risk_model_outputs <- function(risk_models, titles){
  
  if (length(risk_models[["glm_fit_list"]]) == 0) {
    stop("No models were successfully fitted. Please check the regression step.")
  }
  
  glm_fits <- bind_rows(risk_models[["glm_fit_list"]])
  model_perf <- bind_rows(risk_models[["glm_perf_list"]])
  
  glm_fits_fdr <- glm_fits %>%
    group_by(RS) %>%
    group_modify(~{
      tmp <- .x
      
      tmp$p_fdr <- NA_real_
      non_intercept_idx <- which(tmp$vars != "(Intercept)")
      
      if(length(non_intercept_idx) > 0){
        tmp$p_fdr[non_intercept_idx] <- p.adjust(tmp$p[non_intercept_idx], method = "fdr")
      }
      
      tmp$Significant_raw <- ifelse(tmp$vars == "(Intercept)", NA, tmp$p < 0.05)
      tmp$Significant_fdr <- ifelse(tmp$vars == "(Intercept)", NA, tmp$p_fdr < 0.05)
      
      tmp
    }) %>%
    ungroup()
  
  df_sig <- glm_fits_fdr %>%
    filter(vars != "(Intercept)", Significant_fdr == TRUE) %>%
    mutate(
      vars_wrapped = str_wrap(vars_name, width = 40),
      RS_wrapped = factor(RS, levels = titles),
      order = as.numeric(order),
      var_factor = fct_reorder(vars_wrapped, order, .desc = TRUE),
      label = sprintf("%.2f", Estimate)
    )
  
  if(nrow(df_sig) > 0){
    
    limit_sig <- max(abs(df_sig$Estimate), na.rm = TRUE)
    limit_sig <- ceiling(limit_sig)
    
    all_x <- levels(df_sig$RS_wrapped)
    all_y <- levels(df_sig$var_factor)
    
    bg_grid <- expand.grid(
      RS_wrapped = all_x,
      var_factor = all_y,
      stringsAsFactors = FALSE
    )
    
    bg_grid$RS_wrapped <- factor(bg_grid$RS_wrapped, levels = all_x)
    bg_grid$var_factor <- factor(bg_grid$var_factor, levels = all_y)
    
    p_heat_sig <- ggplot() +
      geom_tile(
        data = bg_grid,
        aes(x = RS_wrapped, y = var_factor),
        fill = "white",
        colour = "grey92",
        linewidth = 0.35,
        width = 1,
        height = 1
      ) +
      geom_tile(
        data = df_sig,
        aes(x = RS_wrapped, y = var_factor, fill = Estimate),
        colour = "grey92",
        linewidth = 0.35,
        width = 1,
        height = 1
      ) +
      geom_text(
        data = df_sig,
        aes(x = RS_wrapped, y = var_factor, label = label),
        size = 2.8,
        colour = "black"
      ) +
      scale_fill_steps2(
        low = "#D83016",
        mid = "white",
        high = "#5AB4AC",
        midpoint = 0,
        name = "Coefficient",
        n.breaks = 6,
        limits = c(-limit_sig, limit_sig)
      ) +
      scale_x_discrete(expand = c(0, 0), position = "bottom") +
      scale_y_discrete(expand = c(0, 0)) +
      coord_fixed(ratio = 0.5) +
      labs(
        x = "Perceived Nature-related Risks",
        y = "Explanatory Variables"
      ) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(color = "black", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 12, face = "bold"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5)
      )
    
  } else {
    p_heat_sig <- ggplot() +
      annotate("text", x = 1, y = 1, label = "No predictors passed FDR", size = 6) +
      theme_void()
  }
  
  main_text_table <- glm_fits_fdr %>%
    filter(vars != "(Intercept)", Significant_fdr == TRUE) %>%
    dplyr::select(RS, vars_name, Estimate, `Std. Error`, `t value`, p, p_fdr) %>%
    arrange(RS, p_fdr)
  
  supplementary_table <- glm_fits_fdr %>%
    mutate(
      Significant_raw = case_when(
        is.na(Significant_raw) ~ "",
        Significant_raw ~ "Yes",
        TRUE ~ "No"
      ),
      Significant_fdr = case_when(
        is.na(Significant_fdr) ~ "",
        Significant_fdr ~ "Yes",
        TRUE ~ "No"
      )
    ) %>%
    dplyr::select(RS, vars_name, Estimate, `Std. Error`, `t value`, p, p_fdr,
                  Significant_raw, Significant_fdr) %>%
    arrange(RS, desc(vars_name == "Intercept"), p)
  
  model_perf_table <- model_perf %>%
    dplyr::select(RS, n, r_squared, adj_r_squared) %>%
    arrange(RS)
  
  return(list(
    heatmap_sig_fdr = p_heat_sig,
    glm_fit = glm_fits,
    glm_fit_fdr = glm_fits_fdr,
    main_text_table = main_text_table,
    supplementary_table = supplementary_table,
    model_performance = model_perf_table
  ))
}

####################################################################################
# PART II: TEA FARMERS’ PERSPECTIVES ----
####################################################################################

# 1. Prepare data for Likerts ----
likert_cols_3 <- c("PAST_Soil_Condition", "PAST_Water_availability", "PAST_Water_quality", 
                   "PAST_Temperature", "PAST_Hot_days","PAST_Cold.days",
                   "PAST_Pests",  "PAST_Weeds")

likert_cols_5_experience <- c("RE_Soil_condition","RE_Water_availability","RE_Water_quality", 
                              "RE_Change.in.temp", "RE_Hot_days", "RE_Cold.days","RE_Change_in_rainfall", "RE_Landslides",
                              "RE_Pest_and_disease","RE_Weeds")

likert_cols_5_likelihood <- c("RL_Soil.condition","RL_Water_availability", "RL_Water_quality",
                              "RL_Change.in.temp", "RL_Hot_days", "RL_Cold_days","RL_Change_in_rainfall",  "RL_Landslides",
                              "RL_Pest_and_disease", "RL_Weeds")

likert_cols_5_impact <- c("RI_Soil_condition","RI_Provisioing_Water_availability","RI_Water.quality",
                          "RI_Change.in.temp",  "RI_Hot.days","RI_Cold.days","RI_Change_in_rainfall",  "RI_Landslides",
                          "RI_Pest_and_disease",  "RI_Weeds")

changes_3levels <- convert_to_factor(data, likert_cols_3, c(1, 2, 3), c("Decreased", "No change", "Increased"))
changes_3levels <- na.omit(changes_3levels[, likert_cols_3])

experience_5levels <- convert_to_factor(data, likert_cols_5_experience, c(0,1, 2, 3, 4, 5), 
                                        c("NA","No impact", "Low impact", "Moderate impact", "High impact", "Very high impact"))
experience_5levels <- experience_5levels[, likert_cols_5_experience]
experience_5levels[experience_5levels=="NA"] <- NA
experience_5levels[] <- lapply(experience_5levels, function(x) {
  x <- factor(x, levels = c("No impact", "Low impact", "Moderate impact", "High impact", "Very high impact"), ordered = TRUE)
  return(x)
})

likelihood_5levels <- convert_to_factor(data, likert_cols_5_likelihood, c(1, 2, 3, 4, 5), 
                                        c("No chance", "Low chance", "Moderate chance", "High chance", "Very high chance"))
likelihood_5levels <- na.omit(likelihood_5levels[, likert_cols_5_likelihood])

impact_5levels <- convert_to_factor(data, likert_cols_5_impact, c(1, 2, 3, 4, 5), 
                                    c("No impact", "Low impact", "Moderate impact", "High impact", "Very high impact"))
impact_5levels <- na.omit(impact_5levels[, likert_cols_5_impact])

# 2. Figure: Change perspective----
names(changes_3levels) <- c("Soil health", "Water availability", "Water quality", 
                            "Average temperature ", "Consecutive hot days","Consecutive cold days",
                            "Varieties of pest and disease", "Varieties of weeds")

change_positive_increase <- changes_3levels[, 1:3]
change_negative_increase <- changes_3levels[, 4:8]

change_benefit <- plot(likert(change_positive_increase), group.order = names(change_positive_increase), 
                       colors = c("#D83016","grey90","#5AB4AC"),
                       text.size = 3.5,text.color = "black", centered = TRUE, 
                       legend.position = "right", legend = "Response") +
  ggtitle("(a) Perceived Change in Ecosystem Services Beneficial to Tea Production") +
  theme(axis.text = element_text(size=11, color="black"), legend.text = element_text(size=11))

change_stress <- plot(likert(change_negative_increase), group.order = names(change_negative_increase), 
                      colors = c("#5AB4AC","grey90", "#D83016"),
                      text.size = 3.5,text.color = "black", centered = TRUE, 
                      legend.position = "right", legend = "Response") +
  ggtitle("(b) Perceived Change in Ecosystem Stressors Hindering Tea Production") +
  theme(axis.text = element_text(size=11, color="black"), legend.text = element_text(size=11))

ggarrange(change_benefit, change_stress, nrow=2, align = "v", heights = c(3.8,5))

# 3. Figure: Risk perception ----
var_names <- c("Soil degradation","Water scarcity","Water quality decrease", 
               "Changes of temperature", "Consecutive hot days", "Consecutive cold days",
               "Changes of rainfall pattern",  "Landslides",
               "Pest and disease occurrence","Weeds occurrence")

colors <- c("#CCCCCC", "#FCF0DA", "#FDCB8C", "#FC8E59", "#D83016")

names(experience_5levels) <- var_names
RE_impact <- plot(likert(experience_5levels),
                  colors = colors, 
                  center = 3, plot.percent.neutral = TRUE,
                  text.size = 3.5, text.color = "black", centered = TRUE, 
                  wrap = 50, legend.position = "right", legend = "a, Experienced risk impact") +
  theme(axis.text = element_text(size = 11, color = "black"), 
        legend.text = element_text(size = 11), 
        legend.title = element_text(face = "bold"))

names(likelihood_5levels) <- var_names
FR_likelihood <- plot(likert(likelihood_5levels),
                      colors = colors,
                      text.size = 3.5, text.color = "black", centered = TRUE,
                      legend.position = "right", legend = "b, Future risk likelihood") +
  theme(axis.text = element_text(size = 11, color = "black"), 
        legend.text = element_text(size = 11), 
        legend.title = element_text(face = "bold"))

names(impact_5levels) <- var_names
FR_impact <- plot(likert(impact_5levels),
                  colors = colors,
                  text.size = 3.5, text.color = "black", centered = TRUE,
                  legend.position = "right", legend = "c, Future risk impact") +
  theme(axis.text = element_text(size = 11, color = "black"), 
        legend.text = element_text(size = 11), 
        legend.title = element_text(face = "bold"))

ggarrange(RE_impact, FR_likelihood, FR_impact, 
          ncol = 1, nrow = 3, 
          align = "v", 
          heights = c(4.5, 4.5, 4.5))

####################################################################################
# PART III: FACTORS INFLUENCING RISK PERCEPTION ----
####################################################################################

# 1. Computing risk scores ----
data$RI_overall <- rowMeans(data[, likert_cols_5_impact], na.rm = TRUE)
data$RL_overall <- rowMeans(data[, likert_cols_5_likelihood], na.rm = TRUE)

risk_score <- data[, likert_cols_5_likelihood] * data[, likert_cols_5_impact]
risk_score$Overall <- rowMeans(risk_score, na.rm = TRUE)
names(risk_score) <- c("RS_Soil_condition","RS_Water_availability","RS_Water_quality", 
                       "RS_Change.in.temp", "RS_Hot_days", "RS_Cold.days","RS_Change_in_rainfall", "RS_Landslides",
                       "RS_Pest_and_disease","RS_Weeds", "RS_overall")

data <- cbind(data, risk_score)
####################################################################################
# CHECK DISTRIBUTION OF RISK SCORES (RS variables) ----
####################################################################################

# List of RS variables (11 risks)
rs_vars <- c("RS_Change.in.temp","RS_Change_in_rainfall", "RS_Cold.days", "RS_Hot_days",
             "RS_Pest_and_disease","RS_Weeds", "RS_Landslides",
             "RS_Soil_condition","RS_Water_availability","RS_Water_quality", 
             "RS_overall")

cat("\n================ DISTRIBUTION CHECK FOR RISK SCORES ================\n")

for (var in rs_vars) {
  
  cat("\n--------------------------------------------------\n")
  cat("Variable:", var, "\n")
  
  x <- data[[var]]
  x <- x[!is.na(x)]
  
  # Summary stats
  cat("\nSummary:\n")
  print(summary(x))
  
  # Skewness & kurtosis
  cat("\nSkewness:", round(skewness(x), 3), "\n")
  cat("Kurtosis:", round(kurtosis(x), 3), "\n")
  
  # Shapiro-Wilk test (sample if large n)
  if(length(x) > 5000){
    set.seed(123)
    x_sample <- sample(x, 5000)
    shapiro <- shapiro.test(x_sample)
    cat("\nShapiro-Wilk (sample n=5000): p-value =", shapiro$p.value, "\n")
  } else if(length(x) >= 3 & length(x) <= 5000){
    shapiro <- shapiro.test(x)
    cat("\nShapiro-Wilk: p-value =", shapiro$p.value, "\n")
  } else {
    cat("\nShapiro-Wilk: Not enough observations\n")
  }
}

cat("\n================ END OF DISTRIBUTION CHECK ================\n")
####################################################################################
# APPENDIX TABLE: DISTRIBUTION OF RISK SCORES ----
####################################################################################

library(dplyr)

appendix_table <- data.frame()

for (var in rs_vars) {
  
  x <- data[[var]]
  x <- x[!is.na(x)]
  
  # Shapiro test (with sampling rule)
  if(length(x) > 5000){
    set.seed(123)
    x_sample <- sample(x, 5000)
    shapiro_p <- shapiro.test(x_sample)$p.value
  } else if(length(x) >= 3){
    shapiro_p <- shapiro.test(x)$p.value
  } else {
    shapiro_p <- NA
  }
  
  temp <- data.frame(
    Variable = var,
    Min = min(x),
    Median = median(x),
    Mean = round(mean(x), 2),
    Max = max(x),
    Skewness = round(moments::skewness(x), 3),
    Kurtosis = round(moments::kurtosis(x), 3),
    Shapiro_p = signif(shapiro_p, 3)
  )
  
  appendix_table <- rbind(appendix_table, temp)
}

# View in console
print(appendix_table)

# Save
write.csv(appendix_table,
          file.path(output_dir, "Appendix_Distribution_RiskScores.csv"),
          row.names = FALSE)

# 2. Predictor matrix & correlations between X variables ----
predictor_tbl <- read.csv("RS_variable_matrix_selected_v4_Stepwise.csv")
predictor_tbl$order <- row.names(predictor_tbl)

responses <- c("RS_Change.in.temp","RS_Change_in_rainfall", "RS_Cold.days", "RS_Hot_days",
               "RS_Pest_and_disease","RS_Weeds", "RS_Landslides",
               "RS_Soil_condition","RS_Water_availability","RS_Water_quality", 
               "RS_overall")

titles <- c("Changes of Temperature", "Changes of Rainfall","Consecutive Cold Days","Consecutive Hot Days", 
            "Pest and Disease Occurrence","Weeds Occurrence","Landslides",
            "Soil Degradation","Water Scarcity","Water Quality Decrease",
            "Overall nature-related risks")

all_needed <- unique(c(responses, predictor_tbl$Predictor))
missing_cols <- setdiff(all_needed, names(data))
if (length(missing_cols) > 0) {
  warning("These variables are listed but not found in `data`: ",
          paste(missing_cols, collapse = ", "))
}

numeric_data <- data %>%
  dplyr::select(all_of(intersect(all_needed, names(data)))) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(),
                              ~ suppressWarnings(as.numeric(as.character(.)))))

numeric_data <- numeric_data %>% dplyr::select(where(~ !all(is.na(.))))

pred_cols <- intersect(predictor_tbl$Predictor, names(numeric_data))
if (length(pred_cols) < 2) {
  warning("Not enough numeric predictor columns found to compute a correlation matrix.")
} else {
  X_cors <- stats::cor(numeric_data[, pred_cols, drop = FALSE],
                       use = "pairwise.complete.obs", method = "pearson")
  corrplot(X_cors, method = "square", type = "full",
           tl.col = "black", mar = c(0, 0, 0, 0), tl.srt = 45)
}

# 3. Regression models ----
RS_models  <- risk_model_batch(numeric_data, predictor_tbl, responses, titles, stepwise = TRUE)
RS_results <- risk_model_outputs(RS_models, titles)

# Show only remaining heatmap
RS_results[["heatmap_sig_fdr"]]

# Extract tables
glm_fit_raw         <- RS_results[["glm_fit"]]
glm_fit_fdr         <- RS_results[["glm_fit_fdr"]]
main_text_table     <- RS_results[["main_text_table"]]
supplementary_table <- RS_results[["supplementary_table"]]
model_performance   <- RS_results[["model_performance"]]

# Save tables
write.csv(glm_fit_raw,
          file.path(output_dir, "glm_model_fits_full_data_raw.csv"),
          row.names = FALSE)

write.csv(glm_fit_fdr,
          file.path(output_dir, "glm_model_fits_full_data_fdr.csv"),
          row.names = FALSE)

write.csv(main_text_table,
          file.path(output_dir, "Table_main_text_significant_findings_FDR.csv"),
          row.names = FALSE)

write.csv(supplementary_table,
          file.path(output_dir, "Table_S3_full_regression_coefficients_FDR.csv"),
          row.names = FALSE)

write.csv(model_performance,
          file.path(output_dir, "Table_S4_model_performance_R2.csv"),
          row.names = FALSE)

# =========================
# SAVE FIGURES
# =========================

save_png <- function(plot_obj, file_path, width_mm, height_mm, dpi = 600) {
  png(file_path, width = width_mm, height = height_mm, units = "mm", res = dpi)
  print(plot_obj)
  dev.off()
}

save_pdf <- function(plot_obj, file_path, width_in, height_in) {
  pdf(file_path, width = width_in, height = height_in)
  print(plot_obj)
  dev.off()
}

W_TWO_COL <- 180
H_CHANGE  <- 220
H_RISK    <- 270
H_HEATMAP <- 260
W_HEATMAP <- 260

if (!exists("fig_change")) {
  fig_change <- ggpubr::ggarrange(change_benefit, change_stress,
                                  nrow = 2, align = "v", heights = c(3.8, 5))
}
if (!exists("fig_risk")) {
  fig_risk <- ggpubr::ggarrange(RE_impact, FR_likelihood, FR_impact,
                                ncol = 1, nrow = 3, align = "v",
                                heights = c(4.5, 4.5, 4.5))
}

p_heat_sig <- if (exists("RS_results")) RS_results[["heatmap_sig_fdr"]] else NULL

save_png(fig_change,
         file.path(output_dir, "Figure_changes.png"),
         W_TWO_COL, H_CHANGE)

save_png(fig_risk,
         file.path(output_dir, "Figure_risk_perception.png"),
         W_TWO_COL, H_RISK)

if (!is.null(p_heat_sig)) {
  save_png(p_heat_sig,
           file.path(output_dir, "Figure_regression_heatmap_FDR_significant_only.png"),
           W_HEATMAP, H_HEATMAP)
  
  save_pdf(p_heat_sig,
           file.path(output_dir, "Figure_regression_heatmap_FDR_significant_only.pdf"),
           14, 10)
}

cat("All revised results saved to:\n", output_dir, "\n")
print(list.files(output_dir))