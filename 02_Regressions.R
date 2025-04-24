#read in outliers removed dataset
my_data_or <- read.csv("LBC1936_HealthyAgeing_PhenotypicAndOmics_ES_22JAN2024_outliers_removed.csv")
my_data_or_r <- my_data_or  

#regressing age and sex out of most variables
for (var_name in height_irr) {
  #fit linear model with age (agedays_w1) as the predictor
  model <- lm(formula = as.formula(paste0(var_name, "~ agedays_w1 + sex")), data = my_data_or, na.action = na.exclude)
  #get the residuals from the model
  residuals <- scale(resid(model))
  #print a message for each variable being processed
  #cat("Residuals (age and sex regressed out) for", var_name, ":\n")
  #print the residuals for the current variable
  #print(residuals)
  #store the residuals in a new column named 'residuals_<var_name>' to avoid overwriting
  my_data_or_r[,var_name] <- residuals
}

#regressing age, sex and height out of variables where height is important
for (var_name in height_imp) {
  #fit linear model with age (agedays_w1) as the predictor
  model <- lm(formula = as.formula(paste0(var_name, "~ agedays_w1 + sex + height_w1")), data = my_data_or, na.action = na.exclude)
  #get the residuals from the model
  residuals <- scale(resid(model))
  #print a message for each variable being processed
  #cat("Residuals (age and sex regressed out) for", var_name, ":\n")
  #print the residuals for the current variable
  #print(residuals)
  #store the residuals in a new column named 'residuals_<var_name>' to avoid overwriting
  my_data_or_r[,var_name] <- residuals
}

#save out dataset with outliers removed, residualised 
write.csv(my_data_or_r, file = "LBC1936_HealthyAgeing_PhenotypicAndOmics_ES_22JAN2024_outliers_removed_residualised.csv", row.names = FALSE)