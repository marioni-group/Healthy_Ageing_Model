install.packages("lavaan", dependencies = TRUE)
install.packages("semPlot")
install.packages("dplyr")
install.packages("corrplot")
library(lavaan)
library(semPlot)
library(dplyr)
library(corrplot)

#latent variables
c <- '
  cognitive =~ digback_w1 + vpatotal_w1 + lmtotal_w1 + matreas_w1 + blkdes_w1 
  + spantot_w1 + digsym_w1 + symsear_w1 + ittotal_w1 + crtmean_w1 + 
  vftot_w1 + wtar_w1 + nart_w1
  # within-domain covariances
  blkdes_w1 ~~ matreas_w1 #Visuospatial domain
  blkdes_w1 ~~ spantot_w1
  matreas_w1 ~~ spantot_w1
  nart_w1 ~~ vftot_w1 #Verbal ability domain
  wtar_w1 ~~ vftot_w1
  wtar_w1 ~~ nart_w1
  lmtotal_w1 ~~ vpatotal_w1 #Memory domain
  lmtotal_w1 ~~ digback_w1
  vpatotal_w1 ~~ digback_w1
  ittotal_w1 ~~ digsym_w1 #Processing speed domain
  ittotal_w1 ~~ symsear_w1
  ittotal_w1 ~~ crtmean_w1
  digsym_w1 ~~ symsear_w1
  digsym_w1 ~~ crtmean_w1
  symsear_w1 ~~ crtmean_w1
'

p <- '
  physical =~ grip_w1 + sixmwk_w1 + sysbp_w1 + fvc_w1 + fev_w1 + adl_w1
  fvc_w1 ~~ fev_w1 #Lung function
'

b <-  '
  biological =~ bld_creat_w1 + bld_hba1c_w1 + bld_hdlchol_w1 + 
  bld_hdlrat_w1 + bld_triglyc_w1 + bld_tsh_w1 + bld_crprot_w1 + bld_fibrin_w1 +
  nlratio_w1 + Telomere_length_bp_w1
  bld_hdlrat_w1 ~~ bld_triglyc_w1 #Lipids
  bld_hdlchol_w1 ~~ bld_hdlrat_w1
  bld_hdlchol_w1 ~~ bld_triglyc_w1
  bld_crprot_w1 ~~ bld_fibrin_w1 #Inflammation
  nlratio_w1 ~~ bld_fibrin_w1
  nlratio_w1 ~~ bld_crprot_w1
'

s <-  '
  subjective =~ WHOQOLDomain1_w1 + WHOQOLDomain2_w1 + WHOQOLDomain3_w1 +
  WHOQOLDomain4_w1
'

# List of models and their names
models <- list(
  "Cognitive" = c,
  "Physical" = p,
  "Biological" = b,
  "Subjective" = s
)

#Fit each model, generate SEM plot, print fit measures and summary and generate correlation matrix
cog_mod <- sem(model = c, data = my_data_or_r, missing = "ml.x")

semPaths(cog_mod,
         what = "std",         # Show standardized estimate
         layout = "tree2",        # Arrange nodes in a hierarchical tree structure
         edge.label.cex = 0.8,    # Adjust size of edge labels
         curveAdjacent = TRUE,   # Curve edges for better clarity
         title = TRUE,           # Show title
         #nCharNodes = 6,          # Truncate variable names if needed
         label.cex = 1.6,           # Adjust size of node labels
         fade = FALSE,            #stops edge fading based on loading
         residuals = FALSE,      #removes circle arrows
         intercepts = FALSE,     #removes triangles
         nodeLabels = c("digit\nspan\nbackwards", "verbal\npaired\nassociates", "logical\nmemory", 
                        "matrix\nreasoning","block\ndesign", "spatial\nspan", "digit\nsymbol", 
                        "symbol\nsearch", "inspection\ntime", "reaction\ntime", "verbal\nfluency", 
                        "WTAR","NART", "Cognitive")
)
print(fitMeasures(cog_mod, c("cfi", "tli", "rmsea", "srmr"))
)
#cfi   tli   rmsea  srmr 
#0.964 0.944 0.060 0.040 

print(summary(cog_mod, standardized = TRUE))

#make data frame summary
cog_mod_summary_df <- parameterEstimates(cog_mod)
cog_mod_summary_df$fdr_pvalue <- p.adjust(cog_mod_summary_df$pvalue, method = "fdr")
print(cog_mod_summary_df)

#correlation matrix
cog_cor <- lavInspect(cog_mod, "sampstat")$cov
cog_cor <- cov2cor(cog_cor[cog_tests_list, cog_tests_list])
corrplot(cog_cor, type = "lower", order = "hclust", tl.col = "black", addCoef.col = 'black', 
         method = "color", number.cex = 0.7, number.font = 1)

#pca and scree plot
cog_scores_df <- my_data_or_r[, c(60:64, 67, 72:74, 76:79)]
cog_scores_df <- cog_scores_df[complete.cases(cog_scores_df), ]
cog_pca <- prcomp(cog_scores_df, center = TRUE, scale. = TRUE)
cog_explained_variance <- summary(cog_pca)$importance[2,]
screeplot(cog_pca, main = "Scree Plot of Cog PCA", type = "lines", col = "blue")
title(xlab = "Principal Component")
print(cog_explained_variance)

bio_mod <- sem(model = b, data = my_data_or_r, missing = "ml.x")
semPaths(bio_mod,
         what = "std",         # Show standardized estimate
         layout = "tree2",        # Arrange nodes in a hierarchical tree structure
         edge.label.cex = 0.8,    # Adjust size of edge labels
         curveAdjacent = TRUE,   # Curve edges for better clarity
         title = TRUE,           # Show title
         #nCharNodes = 6,          # Truncate variable names if needed
         label.cex = 1.6,           # Adjust size of node labels
         fade = FALSE,            #stops edge fading based on loading
         residuals = FALSE,      #removes circle arrows
         intercepts = FALSE,     #removes triangles
         nodeLabels = c("creatinine", "hba1c", "hdlc", "hdl\nratio", 
                        "triglycerides", "tsh", "c-reactive\nprotein", "fibrinogen",
                        "neutrophil/\nleukocyte\nratio", "telomere\nlength\nbasepairs", "Biological")
)
print(fitMeasures(bio_mod, c("cfi", "tli", "rmsea", "srmr"))
)
#cfi   tli   rmsea  srmr 
#0.960 0.938 0.042 0.032

print(summary(bio_mod, standardized = TRUE))
#make data frame summary
bio_mod_summary_df <- parameterEstimates(bio_mod)
bio_mod_summary_df$fdr_pvalue <- p.adjust(bio_mod_summary_df$pvalue, method = "fdr")
print(bio_mod_summary_df)

bio_cor <- lavInspect(bio_mod, "sampstat")$cov
bio_cor <- cov2cor(bio_cor[bio_tests_list, bio_tests_list])
corrplot(bio_cor, type = "lower", order = "hclust", tl.col = "black", addCoef.col = 'black', 
         method = "color", number.cex = 0.7, number.font = 1)

bio_scores_df <- my_data_or_r[, c(132, 139, 145:150, 324, 328)]
bio_scores_df <- bio_scores_df[complete.cases(bio_scores_df), ]
bio_pca <- prcomp(bio_scores_df, center = TRUE, scale. = TRUE)
bio_explained_variance <- summary(bio_pca)$importance[2,]
screeplot(bio_pca, main = "Scree Plot of Bio PCA", type = "lines", col = "blue")
title(xlab = "Principal Component")
print(bio_explained_variance)

phys_mod <- sem(model = p, data = my_data_or_r, missing = "ml.x")
semPaths(phys_mod,
         what = "std",         # Show standardized estimate
         layout = "tree2",        # Arrange nodes in a hierarchical tree structure
         edge.label.cex = 0.8,    # Adjust size of edge labels
         curveAdjacent = TRUE,   # Curve edges for better clarity
         title = TRUE,           # Show title
         #nCharNodes = 6,          # Truncate variable names if needed
         label.cex = 1.6,           # Adjust size of node labels
         fade = FALSE,            #stops edge fading based on loading
         residuals = FALSE,      #removes circle arrows
         intercepts = FALSE,     #removes triangles
         nodeLabels = c("grip\nstrength", "six\nmetre\nwalk", "systolic\nblood\npressure",
                        "fvc", "fev", "adl\n(disability)", "Physical")
)
print(fitMeasures(phys_mod, c("cfi", "tli", "rmsea", "srmr"))
)
#cfi   tli   rmsea  srmr 
#0.988 0.977 0.046 0.031

print(summary(phys_mod, standardized = TRUE))

#make data frame summary
phys_mod_summary_df <- parameterEstimates(phys_mod)
phys_mod_summary_df$fdr_pvalue <- p.adjust(phys_mod_summary_df$pvalue, method = "fdr")
print(phys_mod_summary_df)

phys_cor <- lavInspect(phys_mod, "sampstat")$cov
phys_cor <- cov2cor(phys_cor[phys_tests_list, phys_tests_list])
corrplot(phys_cor, type = "lower", order = "hclust", tl.col = "black", addCoef.col = 'black', 
         method = "color", number.cex = 1.3, number.font = 1)

phys_scores_df <- my_data_or_r[, c(91, 95, 108:109, 325:326)]
phys_scores_df <- phys_scores_df[complete.cases(phys_scores_df), ]
phys_pca <- prcomp(phys_scores_df, center = TRUE, scale. = TRUE)
phys_explained_variance <- summary(phys_pca)$importance[2,]
screeplot(phys_pca, main = "Scree Plot of Phys PCA", type = "lines", col = "blue")
title(xlab = "Principal Component")
print(phys_explained_variance)

subj_mod <- sem(model = s, data = my_data_or_r, missing = "ml.x")
semPaths(subj_mod,
         what = "std",         # Show standardized estimate
         layout = "tree2",        # Arrange nodes in a hierarchical tree structure
         edge.label.cex = 0.8,    # Adjust size of edge labels
         curveAdjacent = TRUE,   # Curve edges for better clarity
         title = TRUE,           # Show title
         #nCharNodes = 6,          # Truncate variable names if needed
         label.cex = 1.6,           # Adjust size of node labels
         fade = FALSE,            #stops edge fading based on loading
         residuals = FALSE,      #removes circle arrows
         intercepts = FALSE,     #removes triangles
         nodeLabels = c("WHOQOL\nDomain1\nphysical", "WHOQOL\nDomain2\npsychological",
                        "WHOQOL\nDomain3\nsocial", "WHOQOL\nDomain4\nenvironment", "QOL")
)
print(fitMeasures(subj_mod, c("cfi", "tli", "rmsea", "srmr"))
)
#  cfi   tli rmsea  srmr 
#0.962 0.885 0.147 0.032

print(summary(subj_mod, standardized = TRUE))
#make data frame summary
subj_mod_summary_df <- parameterEstimates(subj_mod)
subj_mod_summary_df$fdr_pvalue <- p.adjust(subj_mod_summary_df$pvalue, method = "fdr")
print(subj_mod_summary_df)

subj_cor <- lavInspect(subj_mod, "sampstat")$cov
subj_cor <- cov2cor(subj_cor[subj_tests_list, subj_tests_list])
corrplot(subj_cor, type = "lower", order = "hclust", tl.col = "black", addCoef.col = 'black', 
         method = "color", number.cex = 1.7, number.font = 1)

subj_scores_df <- my_data_or_r[, c(257:260)]
subj_scores_df <- subj_scores_df[complete.cases(subj_scores_df), ]
subj_pca <- prcomp(subj_scores_df, center = TRUE, scale. = TRUE)
subj_explained_variance <- summary(subj_pca)$importance[2,]
screeplot(subj_pca, main = "Scree Plot of QOL PCA", type = "lines", col = "blue")
title(xlab = "Principal Component")
print(subj_explained_variance)

#Generate HA latent variable and SEM plot
HA_model <- '
  HA =~ cognitive + physical + biological + subjective
'
HA_mod <- sem(c(c, p, b, s, HA_model), data = my_data_or_r, missing = "ml.x")
semPaths(HA_mod,
         what = "std",         # Show standardized estimate
         layout = "tree2",        # Arrange nodes in a hierarchical tree structure
         edge.label.cex = 0.8,    # Adjust size of edge labels
         curveAdjacent = TRUE,   # Curve edges for better clarity
         title = TRUE,           # Show title
         #nCharNodes = 6,          # Truncate variable names if needed
         label.cex = 1.6,           # Adjust size of node labels
         fade = FALSE,            #stops edge fading based on loading
         residuals = FALSE,      #removes circle arrows
         intercepts = FALSE,     #removes triangles
         structural = TRUE,      #removes first layer of loadings
         nodeLabels = c("Cognitive", "Physical", "Biological", "QOL", 
                        "Healthy Ageing")
)
summary(HA_mod, standardized = TRUE)
print(fitMeasures(HA_mod, c("cfi", "tli", "rmsea", "srmr")))
#cfi   tli   rmsea  srmr 
#0.913 0.903 0.041 0.049

#make data frame summary
HA_mod_summary_df <- parameterEstimates(HA_mod)
HA_model <- '
  HA =~ cognitive + physical + biological + subjective
'
HA_mod <- sem(c(c, p, b, s, HA_model), data = my_data_or_r, missing = "ml.x")
semPaths(HA_mod,
         what = "std",         # Show standardized estimate
         layout = "tree2",        # Arrange nodes in a hierarchical tree structure
         edge.label.cex = 0.8,    # Adjust size of edge labels
         curveAdjacent = TRUE,   # Curve edges for better clarity
         title = TRUE,           # Show title
         #nCharNodes = 6,          # Truncate variable names if needed
         label.cex = 1.6,           # Adjust size of node labels
         fade = FALSE,            #stops edge fading based on loading
         residuals = FALSE,      #removes circle arrows
         intercepts = FALSE,     #removes triangles
         structural = TRUE,      #removes first layer of loadings
         nodeLabels = c("Cognitive", "Physical", "Biological", "QOL", 
                        "Healthy Ageing")
)
summary(HA_mod, standardized = TRUE)
print(fitMeasures(HA_mod, c("cfi", "tli", "rmsea", "srmr")))
#cfi   tli   rmsea  srmr 
#0.917 0.907 0.041 0.047

#make data frame summary
HA_mod_summary_df <- parameterEstimates(HA_mod)
HA_mod_summary_df$fdr_pvalue <- p.adjust(HA_mod_summary_df$pvalue, method = "fdr")
print(HA_mod_summary_df)

#saving out to excel
install.packages("openxlsx")
library(openxlsx)

# Create a new workbook
SEM_wb <- createWorkbook()

# Add each data frame as a sheet
addWorksheet(SEM_wb, "cog mod summary")
writeData(SEM_wb, "cog mod summary", cog_mod_summary_df)

addWorksheet(SEM_wb, "phys mod summary")
writeData(SEM_wb, "phys mod summary", phys_mod_summary_df)

addWorksheet(SEM_wb, "bio mod summary")
writeData(SEM_wb, "bio mod summary", bio_mod_summary_df)

addWorksheet(SEM_wb, "QOL mod summary")
writeData(SEM_wb, "QOL mod summary", subj_mod_summary_df)

addWorksheet(SEM_wb, "HA mod summary")
writeData(SEM_wb, "HA mod summary", HA_mod_summary_df)

# Save the workbook
saveWorkbook(SEM_wb, "SEM_model_summaries.xlsx", overwrite = TRUE)