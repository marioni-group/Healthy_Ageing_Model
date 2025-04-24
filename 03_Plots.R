#read in outliers removed and residualised dataset
my_data_or_r <- read.csv("LBC1936_HealthyAgeing_PhenotypicAndOmics_ES_22JAN2024_outliers_removed_residualised.csv")

#creating HA character variables
cog_tests_list <- c("digback_w1", "matreas_w1", "blkdes_w1", 
                    "digsym_w1", "symsear_w1", "nart_w1", 
                    "lmtotal_w1", "vpatotal_w1", "ittotal_w1",
                    "crtmean_w1", "spantot_w1", "vftot_w1", "wtar_w1")
phys_tests_list <- c("grip_w1", "sixmwk_w1", "sysbp_w1", "fvc_w1", 
                     "fev_w1", "adl_w1")
bio_tests_list <- c("bld_creat_w1", "bld_hba1c_w1", "bld_hdlchol_w1", 
                    "bld_hdlrat_w1", "bld_triglyc_w1", "bld_crprot_w1", 
                    "bld_tsh_w1", "bld_fibrin_w1", 
                    "nlratio_w1", "Telomere_length_bp_w1")
subj_tests_list <- c("WHOQOLDomain1_w1", "WHOQOLDomain2_w1", 
                     "WHOQOLDomain3_w1", "WHOQOLDomain4_w1")

#function to make pdf bivariate plots
biv_plot <- function(data, var_names, group_name) {
  #open pdf to save plots
  pdf(paste0(group_name, "_plot.pdf"))
  par(mfrow=c(3,2))
  #loop through variable names
  for (i in seq_len(length(var_names) - 1)) {  #seq_len to avoid off-by-one errors
    for (j in (i + 1):length(var_names)) {  #second loop for pairs
      #access the columns by name
      var_name <- var_names[i]
      var_name2 <- var_names[j]
      #extract numeric data
      x <- data[[var_name]]
      y <- data[[var_name2]]
      #plot the numeric data
      plot(x, y,
             xlab = var_name, 
             ylab = var_name2, 
             main = paste(var_name, "vs", var_name2))
      #add LOESS line
      loess_fit <- loess(y ~ x, data = data)
      loess_line <- predict(loess_fit, newdata = data.frame(x = sort(x)))
      lines(sort(x), loess_line, col = "red", lwd = 2)
    }
  }
  dev.off()
}

#function to make pdf density plots
density_plot <- function(data, var_names, group_name) {
  #open PDF to save density plots
  pdf(paste0(group_name, "_density_plot.pdf"))
  par(mfrow=c(3,2))
  #loop through variable names
  for (var_name in var_names) {
    #extract numeric data
    x <- data[[var_name]]
    #check if the variable is numeric
    if (is.numeric(x)) {
      #plot density
      plot(density(x, na.rm = TRUE), 
           main = paste("Density of", var_name),
           xlab = var_name,
           ylab = "Density",
           col = "blue",
           lwd = 2)
    }
  }
  dev.off()
}

#making bivariate plots for each HA group and saving as pdfs
biv_plot(my_data_or_r, cog_tests_list, "cog_tests")
biv_plot(my_data_or_r, phys_tests_list, "phys_tests")
biv_plot(my_data_or_r, bio_tests_list, "bio_tests")
biv_plot(my_data_or_r, subj_tests_list, "subj_tests")

#making density plots for each HA group and saving as pdfs
density_plot(my_data_or_r, cog_tests_list, "cog_tests")
density_plot(my_data_or_r, phys_tests_list, "phys_tests")
density_plot(my_data_or_r, bio_tests_list, "bio_tests")
density_plot(my_data_or_r, subj_tests_list, "subj_tests")