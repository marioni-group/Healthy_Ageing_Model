install.packages("ggplot2")
install.packages("tidyr")
install.packages("performance")
install.packages("patchwork") 
install.packages("qqplotr")
install.packages("showtext")
library(ggplot2)
library(tidyr)
library(performance)
library(patchwork)
library(qqplotr)
library(showtext)

#extract HA scores     
HA_scores <- lavPredict(HA_mod)

#convert this lavaan matrix to HA data frame
HA_scores_df <- as.data.frame(HA_scores)#read in clocks data
clocks_all_data <- read.csv("LBC_clock_output_LBC36.csv")

#make correlation matrix of outcomes
HA_outcomes <- c("cognitive", "physical", "biological", "subjective")
HA_cor <- cor(HA_scores_df[, HA_outcomes], use = "complete.obs")
rownames(HA_cor)[rownames(HA_cor) == "cognitive"] <- "Cognitive"
colnames(HA_cor)[colnames(HA_cor) == "cognitive"] <- "Cognitive"
rownames(HA_cor)[rownames(HA_cor) == "physical"] <- "Physical"
colnames(HA_cor)[colnames(HA_cor) == "physical"] <- "Physical"
rownames(HA_cor)[rownames(HA_cor) == "biological"] <- "Biological"
colnames(HA_cor)[colnames(HA_cor) == "biological"] <- "Biological"
rownames(HA_cor)[rownames(HA_cor) == "subjective"] <- "QOL"
colnames(HA_cor)[colnames(HA_cor) == "subjective"] <- "QOL"

# Plot using the same style as your lavInspect-based example
corrplot(HA_cor,
         type = "lower",
         order = "hclust",
         tl.col = "black",
         tl.cex = 1.7,
         addCoef.col = "black",
         method = "color",
         number.cex = 1.7,
         number.font = 1)

#add lbc36no to HA df (same order so directy corresponds)
HA_scores_df$lbc36no <- my_data_or_r$lbc36no

#create new data frame filtered for LBC36 and wave 1
clocks_lbc36_wave1 <- clocks_all_data %>%
  filter(cohort == "LBC36", WAVE == "1")

#rename columns to match my_data
colnames(clocks_lbc36_wave1)[colnames(clocks_lbc36_wave1) == "ID"] <- "lbc36no"

#merge HA and clocks dfs
scores_clocks_merged <- merge(HA_scores_df, clocks_lbc36_wave1, by = "lbc36no")

#outlier removal for clocks
clocks <- c("DNAmPhenoAge", "DNAmGrimAge2BasedOnRealAge")
for(i in clocks){
  hist(scores_clocks_merged[[i]],
       xlab = paste("X-Axis Label for", i),
       main = paste("Histogram of", i))
}
for (var_name in clocks) {
  cat("Testing for variable:", var_name, "\n")
  check_normality(scores_clocks_merged[[var_name]])
}
# Testing for variable: DNAmPhenoAge 
# Shapiro-Wilk Test p-value: 4.788071e-07 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: DNAmGrimAge2BasedOnRealAge 
# Shapiro-Wilk Test p-value: 1.339386e-09 
# The variable is NOT normally distributed (Shapiro-Wilk Test).

# Initialize a list to store summaries
new_summary <- list()
# Loop through the new clock variables
for (var_name in clocks) {
  var_data <- na.omit(scores_clocks_merged[[var_name]])
  new_summary[[var_name]] <- data.frame(
    Variable = var_name,
    Mean = mean(var_data),
    SD = sd(var_data),
    Min = min(var_data),
    Max = max(var_data),
    Range = max(var_data) - min(var_data),
    n = length(var_data)
  )
}
# Combine into a single data frame
new_clock_stats <- do.call(rbind, new_summary)
# Append to the existing cont_var_stats data frame
cont_var_stats <- rbind(cont_var_stats, new_clock_stats)

# Create a new workbook for summary stats
stats_wb <- createWorkbook()
# Add each summary stats data frame as a sheet
addWorksheet(stats_wb, "parametric stats")
writeData(stats_wb, "parametric stats", cont_var_stats)

addWorksheet(stats_wb, "nonparametric stats")
writeData(stats_wb, "nonparametric stats", quantile_stats)

addWorksheet(stats_wb, "sex stats")
writeData(stats_wb, "sex stats", cat_var_stats)

# Save the new workbook with a distinct name
saveWorkbook(stats_wb, "summary_statistics.xlsx", overwrite = TRUE)

#check data dimensions pre-outlier removal 
dim(scores_clocks_merged)
#[1] 895  12

#check number of NAs pre-outlier removal  
table(is.na(scores_clocks_merged))
#FALSE 
#10740

scores_clocks_merged <- recode_IQR_outliers_to_na(scores_clocks_merged, clocks, multiplier = 3)

#check data dimensions post_outlier removal 
dim(scores_clocks_merged)
#[1] 895  12

#check number of NAs post-outlier removal 
table(is.na(scores_clocks_merged))
#FALSE  TRUE 
#10739     1 

#loop through linear models for outcomes with clocks, make results data frames
outcomes <- c("cognitive", "physical", "biological", "subjective", "HA")

#DNAmPhenoAge loop
pheno_results_list <- list() #initialise model results list
pheno_results <- data.frame(
  Outcome = character(),
  Predictor = character(),
  Estimate = numeric(),
  LowerCI = numeric(),
  UpperCI = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)  #initialize model results data frame

for (outcome in outcomes) {
  print(paste("Processing outcome:", outcome))  # Debugging output
  
  #formula_p <- as.formula(paste(outcome, "~ scale(Age) + as.factor(Sex) + scale(DNAmPhenoAge)")) #create formula
  formula_p <- as.formula(paste("scale(", outcome, ") ~ scale(DNAmPhenoAge)")) #create formula
  model_p <- lm(formula_p, data = scores_clocks_merged) #fit model
  
  #print diagnostics
  clockp_model_title <- paste("Model:", outcome, "~ DNAmPhenoAge")
  clockp_diagnostics <- check_model(model_p)
  clockp_diagnostics_plots <- plot(clockp_diagnostics)
  print(wrap_plots(clockp_diagnostics_plots) + plot_annotation(title = clockp_model_title))
  
  pheno_results_list[[outcome]] <- summary(model_p) #store full model summary
  
  coef_table_p <- summary(model_p)$coefficients  # Get coefficient and p value matrix
  # Check if DNAmPhenoAge exists in the model (avoid errors if it's missing)
  if (!("scale(DNAmPhenoAge)" %in% rownames(coef_table_p))) {
    print(paste("Skipping", outcome, "because DNAmPhenoAge is missing"))
    next  # Skip to next outcome if DNAmPhenoAge was removed
  }
  betas <- coef_table_p["scale(DNAmPhenoAge)", "Estimate"]  # Extract estimates
  p_values <- coef_table_p["scale(DNAmPhenoAge)", "Pr(>|t|)"]  # Extract p-values
  
  confint_table_p <- confint(model_p) #Get confidence interval matrix
  lower_ci <- confint_table_p["scale(DNAmPhenoAge)", "2.5 %"]
  upper_ci <- confint_table_p["scale(DNAmPhenoAge)", "97.5 %"]
  
  # Create a temporary data frame for this outcome
  temp_df <- data.frame(
    Outcome = outcome,
    Predictor = "DNAmPhenoAge",
    Estimate = betas,
    LowerCI = lower_ci,
    UpperCI = upper_ci,
    P_Value = p_values
  )
  
  # Append results to the final data frame
  pheno_results <- rbind(pheno_results, temp_df)
}

for (outcome in names(pheno_results_list)) {
  cat("\nLinear model for", outcome, ":\n")
  print(pheno_results_list[[outcome]])
  cat("\n", strrep("=", 80), "\n")
}

#DNAmGrimAge2BasedOnRealAge loop
grim_results_list <- list() #initialise model results list
grim_results <- data.frame(
  Outcome = character(),
  Predictor = character(),
  Estimate = numeric(),
  LowerCI = numeric(),
  UpperCI = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)  #initialize model results data frame

for (outcome in outcomes) {
  print(paste("Processing outcome:", outcome))  # Debugging output
  
  #formula_g <- as.formula(paste(outcome, "~ scale(Age) + as.factor(Sex) + scale(DNAmGrimAge2BasedOnRealAge)")) #create formula
  formula_g <- as.formula(paste("scale(", outcome, ") ~ scale(DNAmGrimAge2BasedOnRealAge)"))
  model_g <- lm(formula_g, data = scores_clocks_merged) #fit model
  
  #print diagnostics
  clockg_model_title <- paste("Model:", outcome, "~ DNAmGrimAge2")
  clockg_diagnostics <- check_model(model_g)
  clockg_diagnostics_plots <- plot(clockg_diagnostics)
  print(wrap_plots(clockg_diagnostics_plots) + plot_annotation(title = clockg_model_title))
  
  grim_results_list[[outcome]] <- summary(model_g) #store full model summary
  
  coef_table_g <- summary(model_g)$coefficients  # Get coefficient and p value matrix
  # Check if DNAmGrimAge exists in the model (avoid errors if it's missing)
  if (!("scale(DNAmGrimAge2BasedOnRealAge)" %in% rownames(coef_table_g))) {
    print(paste("Skipping", outcome, "because DNAmGrimAge2BasedOnRealAge is missing"))
    next  # Skip to next outcome if DNAmGrimAge2BasedOnRealAge was removed
  }
  betas <- coef_table_g["scale(DNAmGrimAge2BasedOnRealAge)", "Estimate"]  # Extract estimates
  p_values <- coef_table_g["scale(DNAmGrimAge2BasedOnRealAge)", "Pr(>|t|)"]  # Extract p-values
  
  confint_table_g <- confint(model_g) #Get confidence interval matrix
  lower_ci <- confint_table_g["scale(DNAmGrimAge2BasedOnRealAge)", "2.5 %"]
  upper_ci <- confint_table_g["scale(DNAmGrimAge2BasedOnRealAge)", "97.5 %"]
  
  # Create a temporary data frame for this outcome
  temp_df <- data.frame(
    Outcome = outcome,
    Predictor = "DNAmGrimAge2BasedOnRealAge",
    Estimate = betas,
    LowerCI = lower_ci,
    UpperCI = upper_ci,
    P_Value = p_values
  )
  
  # Append results to the final data frame
  grim_results <- rbind(grim_results, temp_df)
}

for (outcome in names(grim_results_list)) {
  cat("\nLinear model for", outcome, ":\n")
  print(grim_results_list[[outcome]])
  cat("\n", strrep("=", 80), "\n")
}

#Combine data frames
clocks_results <- rbind(pheno_results, grim_results)
clocks_results$Outcome <- gsub("HA", "Healthy Ageing", clocks_results$Outcome)
clocks_results$Outcome <- gsub("subjective", "QOL", clocks_results$Outcome)
clocks_results$Outcome <- factor(clocks_results$Outcome, levels = c("Healthy Ageing", "QOL", "biological", "physical", "cognitive"))

#false discovery rate correction
clocks_results$FDR_Pvalue <- p.adjust(clocks_results$P_Value, method = 'BH')

#CI range calculation
clocks_results$CI_range <- clocks_results$UpperCI - clocks_results$LowerCI

# Add the Times New Roman font
font_add(family = "Times New Roman", 
         regular = "/System/Library/Fonts/Supplemental/Times New Roman.ttf")

# Activate showtext to render text properly
showtext_auto()

clocks_results$Outcome <- factor(clocks_results$Outcome,
                                 levels = c("cognitive", "physical", "biological", "QOL", "Healthy Ageing"),
                                 labels = c("Cognitive", "Physical", "Biological", "QOL", "HA"))

clocks_results$Outcome <- factor(clocks_results$Outcome, levels = rev(levels(factor(clocks_results$Outcome))))

#Graph of models
ggplot(clocks_results, aes(x = Outcome, y = Estimate, color = Predictor)) +  
  # Scatter points
  geom_point(size = 3, position = position_dodge(width = 0.5)) +  
  # Error bars
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, position = position_dodge(width = 0.5)) +
  # Horizontal reference line
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
  # Labels and formatting
  labs(title = "Associations Between Epigenetic Age Acceleration and Healthy Ageing",
       x = "Outcome",
       y = "Standardised Beta [95% CI]",
       color = "Predictor") +
  scale_color_manual(values = c("DNAmPhenoAge" = "blue", "DNAmGrimAge2BasedOnRealAge" = "red"),
                    labels = c("DNAmGrimAge2", "DNAmPhenoAge")) +  # Change legend label here
  theme_minimal(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(size = 16, face = "bold", family = "Times New Roman"),
    axis.text.y = element_text(size = 12, family = "Times New Roman"),
    axis.title.y = element_text(size = 14, family = "Times New Roman"),
    axis.text.x = element_text(size = 12, family = "Times New Roman"),
    axis.title.x = element_text(size = 14, family = "Times New Roman"),
    legend.text = element_text(size = 12, family = "Times New Roman"),
    legend.title = element_text(size = 14, family = "Times New Roman")
  ) +
  coord_flip()  # Flip axes for readability

#Phenotypic correlates
# Merge HA_scores_df with selected columns from my_data using 'lbc36no' as the key
scores_phenotypes_merged <- merge(
  HA_scores_df, 
  my_data_or[, c("lbc36no", "packyears_w1", "alcunitwk_w1", "bmi_w1", "depgrp_w1", "yrsedu_w1", "ageyears_w1", "sex")], 
  by = "lbc36no", 
  all = TRUE
)

#define phenoptypes
phenotypes <- c("packyears_w1", "alcunitwk_w1", "bmi_w1", "depgrp_w1", "yrsedu_w1")
# Define a mapping for predictor names
predictor_labels <- c(
  "yrsedu_w1" = "Years Full-Time Education",
  "packyears_w1" = "Smoking Exposure (packyears)",
  "depgrp_w1" = "SIMD Group",
  "bmi_w1" = "BMI",
  "alcunitwk_w1" = "Alcohol Consumption (units per week)"
)

# Create an empty dataframe to store results
phenotypes_results <- data.frame(
  Outcome = character(),
  Predictor = character(),
  Estimate = numeric(),
  P_Value = numeric(),
  LowerCI = numeric(),
  UpperCI = numeric(),
  #FDR_Pvalue = numeric(), 
  stringsAsFactors = FALSE
)  #initialize model results data frame

# Loop through each combination of outcome and phenotype
for (outcome in outcomes) {
  for (phenotype in phenotypes) {
    #Create a formula 
    #formula_pheno <- as.formula(paste("scale(", outcome, ") ~ scale(ageyears_w1) + as.factor(sex) + scale(", phenotype,")"))
    formula_pheno <- as.formula(paste("scale(", outcome, ") ~ scale(", phenotype,")"))
    #Fit the linear model
    model_pheno <- lm(formula_pheno, data = scores_phenotypes_merged)
    
    # Model title
    model_title_pheno <- paste("Model:", outcome, "~", phenotype)
    
    # Fix: wrap the check_model plots before adding a title
    pheno_diagnostics <- check_model(model_pheno)
    pheno_diagnostics_plots <- plot(pheno_diagnostics)
    print(wrap_plots(pheno_diagnostics_plots) + plot_annotation(title = model_title_pheno))

    # Extract beta coefficient and p-value for the outcome predictor
    beta <- coef(model_pheno)[2]  # Second coefficient is for the predictor
    p_value <- summary(model_pheno)$coefficients[2, 4]  # Extract p-value
    # Calculate confidence intervals
    ci <- confint(model_pheno)[2, ]  # Second row for the predictor
    # Store results
    new_row <- data.frame(
      Outcome = outcome, 
      Predictor = phenotype, 
      Estimate = beta, 
      P_Value = p_value,
      LowerCI = ci[1],
      UpperCI = ci[2],
      stringsAsFactors = FALSE
    )
    phenotypes_results <- rbind(phenotypes_results, new_row)
    }
  }
phenotypes_results$FDR_Pvalue <- p.adjust(phenotypes_results$P_Value, method = 'BH')
phenotypes_results$CI_range <- phenotypes_results$UpperCI - phenotypes_results$LowerCI

# Sort phenotypes_results by Predictor
phenotypes_results <- phenotypes_results[order(phenotypes_results$Predictor), ]

#creating new names for the outcomes
outcome_labels <- c("cognitive" = "Cognitive", "physical" = "Physical", 
                    "biological" = "Biological", "subjective" = "QOL", "HA" = "HA")

# Plot using ggplot2
ggplot(phenotypes_results, aes(x = Estimate, y = Predictor, color = Predictor)) +  # color by Predictor
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = LowerCI, xmax = UpperCI), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  facet_wrap(~ Outcome, nrow = 3, ncol = 2, labeller = labeller(Outcome = outcome_labels)) +
  scale_x_continuous(limits = c(-0.4, 0.6)) +
  scale_color_brewer(palette = "Dark2") +  # You can change this to "Set1", etc.
  theme_minimal(base_family = "Times New Roman") +  # Use Times New Roman globally
  theme(
    panel.spacing = unit(1.5, "lines"),
    plot.title = element_text(size = 14, face = "bold", family = "Times New Roman"),
    axis.title = element_text(family = "Times New Roman"),
    axis.text = element_text(family = "Times New Roman"),
    strip.text = element_text(family = "Times New Roman"),  # facet labels
    legend.position = "none"
  ) +
  labs(
    title = "Associations Between Disease Risk Factors and Healthy Ageing", 
    x = "Standardised Beta [95% CI]", 
    y = "Predictor",
    color = "Predictor"
  ) +
  scale_y_discrete(labels = predictor_labels)

library(openxlsx)

# Create a new workbook
LM_wb <- createWorkbook()

# Add each data frame as a sheet
addWorksheet(LM_wb, "clocks results")
writeData(LM_wb, "clocks results", clocks_results)

addWorksheet(LM_wb, "lifestyles results")
writeData(LM_wb, "lifestyles results", phenotypes_results)

# Save the workbook with new name
saveWorkbook(LM_wb, "Linear_models_results.xlsx", overwrite = TRUE)