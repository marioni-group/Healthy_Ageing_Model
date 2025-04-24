#load required libraries
library(haven)

#read in data
my_data_og <- read_sav("LBC1936_HealthyAgeing_PhenotypicAndOmics_ES_22JAN2024.sav")
height_data <- read.csv("lbc36_height.csv")
my_data <- merge(my_data_og, height_data, by = "lbc36no") #merge height into my_data

#create required columns with averages of variables
my_data$ageyears_w1 <- (my_data$agedays_w1/365.25)
my_data$yearssmoked_w1 <- with(my_data, ifelse(
  smokcat_w1 == 2, ageyears_w1 - smokagestart_w1, 
  ifelse(smokcat_w1 == 1, smokagestop_w1 - smokagestart_w1, 
         0)
))
my_data$packyears_w1 <- (my_data$smoknumcigs_w1/20)*my_data$yearssmoked_w1
my_data$sysbp_w1 <- (my_data$sbp1sit_w1 + my_data$sbp2sit_w1 + my_data$sbp3sit_w1)/3
my_data$grip_w1 <- (my_data$griprh_w1 + my_data$griplh_w1)/2
my_data$nlratio_w1 <- my_data$bld_neut_w1 / my_data$bld_lymph_w1

#define lists of variables of interest
var_list <- c("digback_w1", "matreas_w1", "blkdes_w1", 
              "digsym_w1", "symsear_w1", "nart_w1", 
              "lmtotal_w1", "vpatotal_w1", "ittotal_w1",
              "crtmean_w1", "spantot_w1", "vftot_w1", "wtar_w1",
              "grip_w1", "sixmwk_w1", "sysbp_w1", "fvc_w1", 
              "fev_w1", "adl_w1","bld_creat_w1", "bld_hba1c_w1", 
              "bld_hdlchol_w1", "bld_hdlrat_w1", "bld_triglyc_w1", "bld_crprot_w1", 
              "bld_tsh_w1", "bld_fibrin_w1", "nlratio_w1", 
              "Telomere_length_bp_w1","WHOQOLDomain1_w1", "WHOQOLDomain2_w1", 
              "WHOQOLDomain3_w1", "WHOQOLDomain4_w1", "packyears_w1", 
              "alcunitwk_w1", "bmi_w1", "depgrp_w1", "yrsedu_w1")
input_vars <- c("digback_w1", "matreas_w1", "blkdes_w1", 
                "digsym_w1", "symsear_w1", "nart_w1", 
                "lmtotal_w1", "vpatotal_w1", "ittotal_w1",
                "crtmean_w1", "spantot_w1", "vftot_w1", "wtar_w1",
                "grip_w1", "sixmwk_w1", "sysbp_w1", "fvc_w1", 
                "fev_w1", "adl_w1","bld_creat_w1", "bld_hba1c_w1", 
                "bld_hdlchol_w1", "bld_hdlrat_w1", "bld_triglyc_w1", "bld_crprot_w1", 
                "bld_tsh_w1", "bld_fibrin_w1", "nlratio_w1", 
                "Telomere_length_bp_w1","WHOQOLDomain1_w1", "WHOQOLDomain2_w1", 
                "WHOQOLDomain3_w1", "WHOQOLDomain4_w1")
corr_vars <- c("packyears_w1", "alcunitwk_w1", "bmi_w1", "depgrp_w1", "yrsedu_w1")
height_irr <- c("digback_w1", "matreas_w1", "blkdes_w1", 
                "digsym_w1", "symsear_w1", "nart_w1", 
                "lmtotal_w1", "vpatotal_w1", "ittotal_w1",
                "crtmean_w1", "spantot_w1", "vftot_w1", "wtar_w1", "sysbp_w1", 
                "adl_w1","bld_creat_w1", "bld_hba1c_w1", "bld_hdlchol_w1", 
                "bld_hdlrat_w1", "bld_triglyc_w1", "bld_crprot_w1", 
                "bld_tsh_w1", "bld_fibrin_w1", "nlratio_w1", 
                "Telomere_length_bp_w1","WHOQOLDomain1_w1", "WHOQOLDomain2_w1", 
                "WHOQOLDomain3_w1", "WHOQOLDomain4_w1")
height_imp <- c("grip_w1", "sixmwk_w1", "fvc_w1", "fev_w1")

#create histograms to check continuous vs categorical
for(i in var_list){
  hist(my_data[[i]],
       xlab = paste("X-Axis Label for", i),
       main = paste("Histogram of", i))
}
hist(my_data$ageyears_w1)

#defining categorical and continuous variables
cat_var <- c("adl_w1", "WHOQOLDomain1_w1", "WHOQOLDomain2_w1", 
              "WHOQOLDomain3_w1", "WHOQOLDomain4_w1", "depgrp_w1", "yrsedu_w1")
cont_var <- c("digback_w1", "matreas_w1", "blkdes_w1", 
              "digsym_w1", "symsear_w1", "nart_w1", 
              "lmtotal_w1", "vpatotal_w1", "ittotal_w1",
              "crtmean_w1", "spantot_w1", "vftot_w1", "wtar_w1",
              "grip_w1", "sixmwk_w1", "sysbp_w1", "fvc_w1", 
              "fev_w1","bld_creat_w1", "bld_hba1c_w1", 
              "bld_hdlchol_w1", "bld_hdlrat_w1", "bld_triglyc_w1", "bld_crprot_w1", 
              "bld_tsh_w1", "bld_fibrin_w1", "nlratio_w1", 
              "Telomere_length_bp_w1", "packyears_w1", 
              "alcunitwk_w1", "bmi_w1")


#function for testing variables for normality 
check_normality <- function(variable) {
  shapiro_test <- shapiro.test(variable)
  cat("Shapiro-Wilk Test p-value:", shapiro_test$p.value, "\n")
  
  if (shapiro_test$p.value > 0.05) {
    cat("The variable is normally distributed (Shapiro-Wilk Test).\n")
  } else {
    cat("The variable is NOT normally distributed (Shapiro-Wilk Test).\n")
  }
}

#test all variables of interest for normality 
for (var_name in cont_var) {
  cat("Testing for variable:", var_name, "\n")
  check_normality(my_data[[var_name]])
}

#outputs of normality testing
# Testing for variable: digback_w1 
# Shapiro-Wilk Test p-value: 1.397242e-16 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: matreas_w1 
# Shapiro-Wilk Test p-value: 2.307999e-14 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: blkdes_w1 
# Shapiro-Wilk Test p-value: 6.5972e-08 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: digsym_w1 
# Shapiro-Wilk Test p-value: 0.1492106 
# The variable is normally distributed (Shapiro-Wilk Test).
# Testing for variable: symsear_w1 
# Shapiro-Wilk Test p-value: 1.741786e-05 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: nart_w1 
# Shapiro-Wilk Test p-value: 9.180844e-14 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: lmtotal_w1 
# Shapiro-Wilk Test p-value: 1.453825e-07 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: vpatotal_w1 
# Shapiro-Wilk Test p-value: 1.046082e-17 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: ittotal_w1 
# Shapiro-Wilk Test p-value: 6.194166e-14 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: crtmean_w1 
# Shapiro-Wilk Test p-value: 3.423895e-17 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: spantot_w1 
# Shapiro-Wilk Test p-value: 1.603061e-08 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: vftot_w1 
# Shapiro-Wilk Test p-value: 1.468935e-05 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: wtar_w1 
# Shapiro-Wilk Test p-value: 2.834507e-23 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: grip_w1 
# Shapiro-Wilk Test p-value: 4.299079e-12 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: sixmwk_w1 
# Shapiro-Wilk Test p-value: 1.455349e-35 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: sysbp_w1 
# Shapiro-Wilk Test p-value: 4.111928e-08 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: fvc_w1 
# Shapiro-Wilk Test p-value: 4.316695e-08 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: fev_w1 
# Shapiro-Wilk Test p-value: 0.002038109 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: bld_creat_w1 
# Shapiro-Wilk Test p-value: 2.956476e-17 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: bld_hba1c_w1 
# Shapiro-Wilk Test p-value: 3.201573e-37 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: bld_hdlchol_w1 
# Shapiro-Wilk Test p-value: 1.892259e-17 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: bld_hdlrat_w1 
# Shapiro-Wilk Test p-value: 1.874787e-14 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: bld_triglyc_w1 
# Shapiro-Wilk Test p-value: 1.659574e-37 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: bld_crprot_w1 
# Shapiro-Wilk Test p-value: 6.354838e-45 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: bld_tsh_w1 
# Shapiro-Wilk Test p-value: 6.89838e-44 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: bld_fibrin_w1 
# Shapiro-Wilk Test p-value: 5.908021e-12 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: nlratio_w1 
# Shapiro-Wilk Test p-value: 3.39266e-30 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: Telomere_length_bp_w1 
# Shapiro-Wilk Test p-value: 6.304872e-18 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: packyears_w1 
# Shapiro-Wilk Test p-value: 2.445626e-22 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: alcunitwk_w1 
# Shapiro-Wilk Test p-value: 1.76649e-39 
# The variable is NOT normally distributed (Shapiro-Wilk Test).
# Testing for variable: bmi_w1 
# Shapiro-Wilk Test p-value: 1.657698e-15 
# The variable is NOT normally distributed (Shapiro-Wilk Test).

#define list of non-normal distributed variables
non_norm <- c("digback_w1", "matreas_w1", "blkdes_w1", 
               "symsear_w1", "nart_w1", 
              "lmtotal_w1", "vpatotal_w1", "ittotal_w1",
              "crtmean_w1", "spantot_w1", "vftot_w1", "wtar_w1",
              "grip_w1", "sixmwk_w1", "sysbp_w1", "fvc_w1", 
              "fev_w1","bld_creat_w1", "bld_hba1c_w1", 
              "bld_hdlchol_w1", "bld_hdlrat_w1", "bld_triglyc_w1", "bld_crprot_w1", 
              "bld_tsh_w1", "bld_fibrin_w1", "nlratio_w1", 
              "Telomere_length_bp_w1", "packyears_w1", 
              "alcunitwk_w1", "bmi_w1")

#define list of normally distributed variables
norm <- c("digsym_w1")

#function for recoding outliers to NA by IQR if non-normally distributed
recode_IQR_outliers_to_na <- function(data, var_names, multiplier = 3) {
  for (var_name in var_names) {
    # Extract the variable data
    variable_data <- data[[var_name]]
    # Calculate quartiles and IQR
    Q1 <- quantile(variable_data, 0.25, na.rm = TRUE)
    Q3 <- quantile(variable_data, 0.75, na.rm = TRUE)
    IQR_value <- IQR(variable_data, na.rm = TRUE)
    # Compute outlier thresholds
    outlier_threshold_lower <- Q1 - multiplier * IQR_value
    outlier_threshold_upper <- Q3 + multiplier * IQR_value
    # Identify outliers and recode them to NA
    outlier_indices <- which(variable_data < outlier_threshold_lower |
                               variable_data > outlier_threshold_upper)
    # Recode outliers to NA
    variable_data[outlier_indices] <- NA
    data[[var_name]] <- variable_data
    # Print the count of outliers per variable
    cat("\nOutliers for", var_name, ":")
    if (length(outlier_indices) > 0) {
      cat(length(outlier_indices), "outliers recoded to NA.\n")
    } else {
      cat("No outliers identified.\n")
    }
  }
  return(data)
}

#function for recoding outliers to NA by SD if normally distributed
recode_sd_outliers_to_na <- function(data, var_names, sd_threshold = 3) {
  for (var_name in var_names) {
    # Extract the variable data
    variable_data <- data[[var_name]]
    # Calculate the mean and standard deviation
    mean_value <- mean(variable_data, na.rm = TRUE)
    sd_value <- sd(variable_data, na.rm = TRUE)
    # Define the outlier thresholds (mean ± 3 * standard deviation)
    lower_threshold <- mean_value - sd_threshold * sd_value
    upper_threshold <- mean_value + sd_threshold * sd_value
    # Identify outliers and replace them with NA
    outlier_indices <- which(variable_data < lower_threshold |
                               variable_data > upper_threshold)
    # Recode outliers to NA
    variable_data[outlier_indices] <- NA
    data[[var_name]] <- variable_data
    # Print the count of outliers per variable
    cat("\nOutliers for", var_name, ":")
    if (length(outlier_indices) > 0) {
      cat(length(outlier_indices), "outliers recoded to NA.\n")
    } else {
      cat("No outliers identified.\n")
    }
  }
  return(data)
}

#check data dimensions pre-outlier removal 
dim(my_data)
#[1] 1091  332

#check number of NAs pre-outlier removal  
table(is.na(my_data))
# FALSE   TRUE 
#327214  34998

#applying IQR removal function to appropriate variables
my_data <- recode_IQR_outliers_to_na(my_data, non_norm, multiplier = 3)

#applying SD removal function to appropriate variables 
my_data <- recode_sd_outliers_to_na(my_data, norm, sd_threshold = 3)

#check data dimensions post_outlier removal 
dim(my_data)
#[1] 1091  332

#check number of NAs post-outlier removal 
table(is.na(my_data))
# FALSE   TRUE 
#327035  35177 

#save out dataset with outliers removed 
write.csv(my_data, file = "LBC1936_HealthyAgeing_PhenotypicAndOmics_ES_22JAN2024_outliers_removed.csv", row.names = FALSE)
my_data_or <- read.csv("LBC1936_HealthyAgeing_PhenotypicAndOmics_ES_22JAN2024_outliers_removed.csv")

#create data frames with descriptive statistics of data being used
#Create an empty list to store summary statistics
#Continuous variables
cont_summary <- list()
# Variables to summarize
vars_to_summarise <- c("ageyears_w1", var_list)
# Loop through each variable and compute summary statistics
for (var_name in vars_to_summarise) {
  # Remove NAs and compute statistics
  var_data <- na.omit(my_data_or[[var_name]])  # Remove NAs before calculations
  
  cont_summary[[var_name]] <- data.frame(
    Variable = var_name,
    Mean = mean(var_data, na.rm = TRUE),
    SD = sd(var_data, na.rm = TRUE),
    Min = min(var_data, na.rm = TRUE),
    Max = max(var_data, na.rm = TRUE),
    Range = max(var_data, na.rm = TRUE) - min(var_data, na.rm = TRUE),
    n = length(var_data)  # Number of non-missing values
  )
}

# Combine all statistics into a single dataframe
cont_var_stats <- do.call(rbind, cont_summary)

# Print the summary statistics dataframe
#print(cont_var_stats)

#Categorical variables
# Function to summarize categorical variables
categorical_summary <- function(var_name, data) {
  var_data <- as.factor(data[[var_name]])  # Convert to factor if not already
  freq_table <- table(var_data, useNA = "ifany")  # Count occurrences, including NAs
  prop_table <- prop.table(freq_table) * 100  # Convert to percentages
  
  # Create summary dataframe
  summary_df <- data.frame(
    Variable = var_name,
    Category = names(freq_table),
    Count = as.numeric(freq_table),
    Percentage = round(as.numeric(prop_table), 2)
  )
  
  return(summary_df)
}

# Generate summary for the 'sex' column
cat_var_stats <- categorical_summary("sex", my_data_or)

# Filter the dataframe to keep only the row where Category == "1"
cat_var_stats <- cat_var_stats[cat_var_stats$Category == "1", ]

# Rename the Category "1" to "Male"
cat_var_stats$Category <- "Male"

# Print the summary statistics
#print(cat_var_stats)

#Quantile summary stats for variables that don't look normally distributed based on histograms
# Define variables for which you want median, Q1, and Q3
quantile_vars <- c("matreas_w1", "nart_w1", "vpatotal_w1", "wtar_w1", "adl_w1", 
                   "bld_crprot_w1", "WHOQOLDomain1_w1", "WHOQOLDomain2_w1", 
                   "WHOQOLDomain3_w1", "packyears_w1", "alcunitwk_w1", "depgrp_w1",
                   "yrsedu_w1")

# Create an empty list to store summary statistics
quantile_summary <- list()

# Loop through each variable and compute statistics
for (var_name in quantile_vars) {
  var_data <- na.omit(my_data_or[[var_name]])  # Remove NAs before calculations
  
  quantile_summary[[var_name]] <- data.frame(
    Variable = var_name,
    Median = median(var_data, na.rm = TRUE),
    Q1 = quantile(var_data, 0.25, na.rm = TRUE),
    Q3 = quantile(var_data, 0.75, na.rm = TRUE),
    n = length(var_data)  # Number of non-missing values
  )
}

# Combine into a single dataframe
quantile_stats <- do.call(rbind, quantile_summary)

# Print the new summary table
#print(quantile_stats)