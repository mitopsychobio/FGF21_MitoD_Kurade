####### MiSBE FGF21 Analyses-----
## Prepared by Mangesh Kurade & Alex Behnke

## Interpolation for FGF21 from ELISA plates absorbance readings-----

## This script demonstrates the workflow used to interpolate protein concentration
# from ELISA absorbance readings using data from Plate 8 as an example.

# NOTE:
#   - For each ELISA plate, samples are run in duplicates across two sub-plates,
#     labeled as Plate8A and Plate8B (duplicate measurements).
#   - This code processes both sub-plates simultaneously and later combines their
#     results to obtain average protein concentrations.
#
#   - The same code skeleton was applied to all 12 ELISA plates.
#     Each plate's analysis was executed separately with adjustments to:
#       * File paths for input/output.
#       * Platemap sheet numbers (if applicable).

##########################################################################################################################################################################################################################

# Load required libraries

library(tidyverse)
library(janitor)
library(dplyr)
library(tidyr)
library(stringr)
install.packages("matrixStats")
library(matrixStats)
install.packages("readxl")
packageVersion("dplyr")
install.packages("writexl")
library(writexl)


# ==============================================================================
# PLATEMAP
# ==============================================================================

# Read the platemap for Plate 8. (Adjust the sheet number for other plates as needed)

fgf_8_platemap <- readxl::read_xlsx("/Users/MangeshK/Desktop/Spectramax/FGF_Platemaps_All.xlsx" , sheet = 8)%>%
  pivot_longer(cols = -...1,names_to = "column", values_to = "sampleName")%>%
  rename(row = ...1) %>%
  mutate(across('row', str_replace, 'row', '')) %>%
  mutate(across('column', str_replace, 'C', '')) %>%
  unite(well, row, column, sep = "", remove = F)

# ==============================================================================
# READ AND PROCESS DATA
# ==============================================================================


# Read Plate 8 raw data file (which contains data for both sub-plates: A & B)

fgf_8AB <- readxl::read_xlsx("/Users/MangeshK/Desktop/Spectramax/Plate 8/Plate 8AB.xlsx")

# separate 450 and 540

fgf_8A_450 <- fgf_8AB[3:10, 3:14]
fgf_8A_540 <- fgf_8AB[3:10, 16:27]
fgf_8B_450 <- fgf_8AB[15:22, 3:14]
fgf_8B_540 <- fgf_8AB[15:22, 16:27]

# Rename columns for consistency

names(fgf_8A_450)[1:12] <- as.character(1:12)
names(fgf_8A_540)[1:12] <- as.character(1:12)
names(fgf_8B_450)[1:12] <- as.character(1:12)
names(fgf_8B_540)[1:12] <- as.character(1:12)

# Process data for Plate8A (450 and 540 readings)

fgf_8A_450_final <- fgf_8A_450 %>%
  mutate_all(as.numeric) %>%
  mutate(row = LETTERS[1:8]) %>%
  select(row, everything()) %>%
  pivot_longer(cols = !row, names_to = "column", values_to = "OD_450") %>%
  unite(well, row, column, sep = "", remove = FALSE) %>%
  full_join(fgf_8_platemap,.)


fgf_8A_540_final <- fgf_8A_540 %>%
  mutate_all(as.numeric) %>%
  mutate(row = LETTERS[1:8]) %>%
  select(row, everything()) %>%
  pivot_longer(cols = !row, names_to = "column", values_to = "OD_540") %>%
  unite(well, row, column, sep = "", remove = FALSE) %>%
  full_join(fgf_8_platemap,.)

# Process data for Plate8B (450 and 540 readings)

fgf_8B_450_final <- fgf_8B_450 %>%
  mutate_all(as.numeric) %>%
  mutate(row = LETTERS[1:8]) %>%
  select(row, everything()) %>%
  pivot_longer(cols = !row, names_to = "column", values_to = "OD_450") %>%
  unite(well, row, column, sep = "", remove = FALSE) %>%
  full_join(fgf_8_platemap,.)


fgf_8B_540_final <- fgf_8B_540 %>%
  mutate_all(as.numeric) %>%
  mutate(row = LETTERS[1:8]) %>%
  select(row, everything()) %>%
  pivot_longer(cols = !row, names_to = "column", values_to = "OD_540") %>%
  unite(well, row, column, sep = "", remove = FALSE) %>%
  full_join(fgf_8_platemap,.)


# ==============================================================================
# BACKGROUND & BLANK CORRECTION
# ==============================================================================

# Background correction: subtract 540 nm absorbance from 450 nm absorbance

fgf_8A_bgcorr <- full_join(fgf_8A_450_final, fgf_8A_540_final) %>%
  mutate(OD_bgcorr = OD_450 - OD_540)
fgf_8B_bgcorr <- full_join(fgf_8B_450_final, fgf_8B_540_final) %>%
  mutate(OD_bgcorr = OD_450 - OD_540)

# Blank correction: subtract the blank reading from each well

fgf_8A_blank <- fgf_8A_bgcorr$OD_bgcorr[fgf_8A_bgcorr$sampleName == "BLANK"]
fgf_8A_blankcorr <- fgf_8A_bgcorr %>%
  mutate(OD_blankcorr = OD_bgcorr - fgf_8A_blank) %>%
  select(-OD_450, -OD_540, -OD_bgcorr) %>%
  rename (OD = OD_blankcorr) %>%
  mutate(plate = "A")

fgf_8B_blank <- fgf_8B_bgcorr$OD_bgcorr[fgf_8B_bgcorr$sampleName == "BLANK"]
fgf_8B_blankcorr <- fgf_8B_bgcorr %>%
  mutate(OD_blankcorr = OD_bgcorr - fgf_8B_blank) %>%
  select(-OD_450, -OD_540, -OD_bgcorr) %>%
  rename (OD = OD_blankcorr) %>%
  mutate(plate = "B")

# ==============================================================================
# LINEAR MODEL FITTING & INTERPOLATION
# ==============================================================================

# Combine data from Plate8A and Plate8B for duplicate measurements

fgf_8_all <- rbind(fgf_8A_blankcorr, fgf_8B_blankcorr) %>%
  mutate(OD = case_when(OD >= 0 ~ OD, TRUE ~ NA_real_)) %>%
  mutate(log_OD = ifelse(is.na(OD) | OD <= 0, NA_real_, log(OD)))

# Define the standard samples and their known concentrations

std <- fgf_8_all %>%
  filter(sampleName %in% c("STD 1", "STD 3", "STD 5", "STD 7")) %>%
  mutate(conc = case_when(
    sampleName == "STD 1" ~ 2000,
    sampleName == "STD 3" ~ 500,
    sampleName == "STD 5" ~ 125,
    sampleName == "STD 7" ~ 31.3
  )) %>%
  mutate(log_conc = log(conc))

# Fit a linear model on the logarithmic values of OD and concentration
model <- lm(log_conc ~ log_OD, data = std)

# Interpolate unknown sample concentrations using the fitted model.
# The two sub-plate measurements (from Plate8A and Plate8B) are later combined.

log_OD <-  fgf_8_all$log_OD

plate8_interpolated <- fgf_8_all %>%
  mutate(log_interpolated_conc = predict(model, newdata = data.frame(x = log_OD))) %>%
  mutate(interpolated_conc = exp(log_interpolated_conc)) %>%
  pivot_wider(names_from = plate, values_from = c(OD, log_OD, log_interpolated_conc, interpolated_conc)) %>%
  mutate(average_conc = rowMeans(select(., starts_with("interpolated_conc")), na.rm = TRUE)) %>%
  mutate(sd = rowSds(as.matrix(select(., starts_with("interpolated_conc")), na.rm = TRUE))) %>%
  mutate(cv = sd/average_conc) %>%
  mutate(real_conc = case_when(
    !sampleName %in% c("STD 1", "STD 3", "STD 5", "STD 7") ~ average_conc * 4,
    TRUE ~ average_conc
  ))

# ==============================================================================
# OUTPUT RESULTS
# ==============================================================================

# Write the interpolated results for Plate 8 to an Excel file.

write_xlsx(plate8_interpolated, "/Users/MangeshK/Desktop/Spectramax/Plate 8/plate8int.xlsx")


###  ----------- END OF INTERPOLATION CODE ----------- ###






##-----Psychosocial scores adjustment-----

# Load necessary libraries
library(readxl)
library(dplyr)
library(car)       # For Type 3 ANOVA
library(visreg)    # For visualization
library(emmeans)   # For emtrends and simple slopes
library(lmtest)    # For model assumptions testing
library(nortest)   # For normality testing
library(effectsize)  # For calculating effect sizes
library(openxlsx)  # For writing results to Excel

### Positive Psychosocial Factors and Fasting FGF21-----

# Import data
df <- readxl::read_excel("Age_fat_positive_All.xlsx")  # Wide format

# Clean up variable issues
df$Condition <- as.factor(df$Condition)

# Assemble all MitoDs in one group
df <- df %>%
  mutate(Group = case_when(
    Condition %in% c("3243A>G", "MELAS", "Deletion") ~ "MitoD",
    TRUE ~ Condition
  ))

# Ensure Group is a factor
df$Group <- as.factor(df$Group)

# Set contrasts for Group to contr.sum
contrasts(df$Group) <- contr.sum

# List of variables to iterate through (including 'Couples' as an example)
variable_list <- c("Couples", 
                   "PSS", 
                   "SSQ_S", "
                   SSQ_N", 
                   "Mdes_Mean", 
                   "soc_compre", 
                   "soc_manage", 
                   "soc_meaning",
                   "pwbs_total", 
                   "pwbs_autonomy", 
                   "pwbs_enviromastery", 
                   "pwbs_personal_gro",
                   "pwbs_pos_rel", 
                   "pwbs_purp_life", 
                   "pwbs_self_accept"
                   )

# Loop through each variable and perform analysis
for (var_name in variable_list) {
  
  # Centering continuous variables (including the one from the loop)
  df[[paste0(var_name, ".c")]] <- scale(df[[var_name]], scale=F)[,1]
  df$Age.c <- scale(df$Age, scale=F)[,1]
  df$Fat.c <- scale(df$Percent_Fat, scale=F)[,1]
  
  # Dynamically create the formula using reformulate
  formula <- reformulate(c(paste0(var_name, ".c * Group"),
                           "Fat.c * Group", "Age.c"), response = "Fasting")
  
  # Run the model with the dynamically created formula
  mdl <- lm(formula, data = df, contrasts = list(Group = "contr.sum"))
  
  # Extract simple slopes (to compute associations for each group separately)
  slopes <- emtrends(mdl, ~Group, var=paste0(var_name, ".c"))
  
  # Test the slopes and get p-values and t-ratios
  slopes_test <- test(slopes, adjust="holm")
  
  # Extract p-values for each group (Control and MitoD)
  p_values <- slopes_test$p.value
  
  # Extract r-values for each group (Control and MitoD) from t-values
  r_values <- effectsize::t_to_r(slopes_test$t.ratio, df_error = df.residual(mdl))
  
  # Combine the group names with their corresponding p-values and r-values
  results_df <- data.frame(
    Group = slopes_test$Group,
    P_value = p_values,
    R_value = r_values
  )
  
  # Save results to Excel
  write.xlsx(results_df, file = paste0("results_", var_name, "_Control_vs_MitoD_Fasting.xlsx"))
}



### Positive Psychosocial Factors and Fed FGF21 (stress time 1 - S.1) ----

# Import data
df <- readxl::read_excel("Age_fat_positive_All.xlsx")  # Wide format

# Clean up variable issues
df$Condition <- as.factor(df$Condition)

# Assemble all MitoDs in one group
df <- df %>%
  mutate(Group = case_when(
    Condition %in% c("3243A>G", "MELAS", "Deletion") ~ "MitoD",
    TRUE ~ Condition
  ))

# Ensure Group is a factor
df$Group <- as.factor(df$Group)

# Set contrasts for Group to contr.sum
contrasts(df$Group) <- contr.sum

# List of variables to iterate through (including 'Couples' as an example)
variable_list <- c("Couples", 
                   "PSS", 
                   "SSQ_S", 
                   "SSQ_N", 
                   "Mdes_Mean", 
                   "soc_compre", 
                   "soc_manage", 
                   "soc_meaning",
                   "pwbs_total", 
                   "pwbs_autonomy", 
                   "pwbs_enviromastery", 
                   "pwbs_personal_gro",
                   "pwbs_pos_rel", 
                   "pwbs_purp_life", 
                   "pwbs_self_accept"
                   )

# Loop through each variable and perform analysis
for (var_name in variable_list) {
  
  # Centering continuous variables (including the one from the loop)
  df[[paste0(var_name, ".c")]] <- scale(df[[var_name]], scale=F)[,1]
  df$Age.c <- scale(df$Age, scale=F)[,1]
  df$Fat.c <- scale(df$Percent_Fat, scale=F)[,1]
  
  # Dynamically create the formula using reformulate
  formula <- reformulate(c(paste0(var_name, ".c * Group"),
                           "Fat.c * Group", "Age.c"), response = "S.1")
  
  # Run the model with the dynamically created formula
  mdl <- lm(formula, data = df, contrasts = list(Group = "contr.sum"))
  
  # Extract simple slopes (to compute associations for each group separately)
  slopes <- emtrends(mdl, ~Group, var=paste0(var_name, ".c"))
  
  # Test the slopes and get p-values and t-ratios
  slopes_test <- test(slopes, adjust="holm")
  
  # Extract p-values for each group (Control and MitoD)
  p_values <- slopes_test$p.value
  
  # Extract r-values for each group (Control and MitoD) from t-values
  r_values <- effectsize::t_to_r(slopes_test$t.ratio, df_error = df.residual(mdl))
  
  # Combine the group names with their corresponding p-values and r-values
  results_df <- data.frame(
    Group = slopes_test$Group,
    P_value = p_values,
    R_value = r_values
  )
  
  # Save results to Excel
  write.xlsx(results_df, file = paste0("results_", var_name, "_Control_vs_MitoD_Fed.xlsx"))
}



### Negative Psychosocial Factors and Fasting FGF21-----

# Import data
df <- readxl::read_excel("Age_fat_Negative_All.xlsx")  # Wide format

# Clean up variable issues
df$Condition <- as.factor(df$Condition)

# Assemble all MitoDs in one group
df <- df %>%
  mutate(Group = case_when(
    Condition %in% c("3243A>G", "MELAS", "Deletion") ~ "MitoD",
    TRUE ~ Condition
  ))

# Ensure Group is a factor
df$Group <- as.factor(df$Group)

# Set contrasts for Group to contr.sum
contrasts(df$Group) <- contr.sum

# List of variables to iterate through (including 'Couples' as an example)
variable_list <- c("Loneliness_Scale", 
                   "PSS", 
                   "TICS", 
                   "Daily_Hassles_Scale", 
                   "CTQ_Emotional_Abuse", 
                   "CTQ_Physical_Abuse", 
                   "CTQ_Sexual_Abuse", 
                   "CTQ_Emotional_Neglect",
                   "CTQ_Physical_Neglect", 
                   "STAY_State_Anxiety", 
                   "STAY_Trait_Anxiety", 
                   "BDI",
                   "MBI_Exhaustion", 
                   "MBI_Cynicism", 
                   "MBI_Inefficacy", 
                   "MBI_Total", 
                   "Mdes_Negative", 
                   "LEQ_Stressful_Events", 
                   "LEQ_Positive_Events",
                   "LEQ_Negative_Events",
                   "PCL_Total", 
                   "STRAIN_Total", 
                   "STRAIN_Severity")

# Loop through each variable and perform analysis
for (var_name in variable_list) {
  
  # Centering continuous variables (including the one from the loop)
  df[[paste0(var_name, ".c")]] <- scale(df[[var_name]], scale=F)[,1]
  df$Age.c <- scale(df$Age, scale=F)[,1]
  df$Fat.c <- scale(df$Percent_Fat, scale=F)[,1]
  
  # Dynamically create the formula using reformulate
  formula <- reformulate(c(paste0(var_name, ".c * Group"),
                           "Fat.c * Group", "Age.c"), response = "Fasting")
  
  # Run the model with the dynamically created formula
  mdl <- lm(formula, data = df, contrasts = list(Group = "contr.sum"))
  
  # Extract simple slopes (to compute associations for each group separately)
  slopes <- emtrends(mdl, ~Group, var=paste0(var_name, ".c"))
  
  # Test the slopes and get p-values and t-ratios
  slopes_test <- test(slopes, adjust="holm")
  
  # Extract p-values for each group (Control and MitoD)
  p_values <- slopes_test$p.value
  
  # Extract r-values for each group (Control and MitoD) from t-values
  r_values <- effectsize::t_to_r(slopes_test$t.ratio, df_error = df.residual(mdl))
  
  # Combine the group names with their corresponding p-values and r-values
  results_df <- data.frame(
    Group = slopes_test$Group,
    P_value = p_values,
    R_value = r_values
  )
  
  # Save results to Excel
  write.xlsx(results_df, file = paste0("results_", var_name, "_Control_vs_MitoD.xlsx"))
}


### Negative Psychosocial Factors and Fed FGF21 (stress time 1 - S.1) ----

# Import data
df <- readxl::read_excel("Age_fat_Negative_All.xlsx")  # Wide format

# Clean up variable issues
df$Condition <- as.factor(df$Condition)

# Assemble all MitoDs in one group
df <- df %>%
  mutate(Group = case_when(
    Condition %in% c("3243A>G", "MELAS", "Deletion") ~ "MitoD",
    TRUE ~ Condition
  ))

# Ensure Group is a factor
df$Group <- as.factor(df$Group)

# Set contrasts for Group to contr.sum
contrasts(df$Group) <- contr.sum

# List of variables to iterate through (including 'Couples' as an example)
variable_list <- c("Loneliness_Scale", 
                   "PSS", 
                   "TICS", 
                   "Daily_Hassles_Scale", 
                   "CTQ_Emotional_Abuse", 
                   "CTQ_Physical_Abuse", 
                   "CTQ_Sexual_Abuse", 
                   "CTQ_Emotional_Neglect",
                   "CTQ_Physical_Neglect", 
                   "STAY_State_Anxiety", 
                   "STAY_Trait_Anxiety", 
                   "BDI",
                   "MBI_Exhaustion", 
                   "MBI_Cynicism", 
                   "MBI_Inefficacy", 
                   "MBI_Total", 
                   "Mdes_Negative", 
                   "LEQ_Stressful_Events", 
                   "LEQ_Positive_Events", 
                   "LEQ_Ngative_vents", 
                   "PCL_Total", 
                   "STRAIN_Total", 
                   "STRAIN_Severity"
                   )

# Loop through each variable and perform analysis
for (var_name in variable_list) {
  
  # Centering continuous variables (including the one from the loop)
  df[[paste0(var_name, ".c")]] <- scale(df[[var_name]], scale=F)[,1]
  df$Age.c <- scale(df$Age, scale=F)[,1]
  df$Fat.c <- scale(df$Percent_Fat, scale=F)[,1]
  
  # Dynamically create the formula using reformulate
  formula <- reformulate(c(paste0(var_name, ".c * Group"),
                           "Fat.c * Group", "Age.c"), response = "S.1")
  
  # Run the model with the dynamically created formula
  mdl <- lm(formula, data = df, contrasts = list(Group = "contr.sum"))
  
  # Extract simple slopes (to compute associations for each group separately)
  slopes <- emtrends(mdl, ~Group, var=paste0(var_name, ".c"))
  
  # Test the slopes and get p-values and t-ratios
  slopes_test <- test(slopes, adjust="holm")
  
  # Extract p-values for each group (Control and MitoD)
  p_values <- slopes_test$p.value
  
  # Extract r-values for each group (Control and MitoD) from t-values
  r_values <- effectsize::t_to_r(slopes_test$t.ratio, df_error = df.residual(mdl))
  
  # Combine the group names with their corresponding p-values and r-values
  results_df <- data.frame(
    Group = slopes_test$Group,
    P_value = p_values,
    R_value = r_values
  )
  
  # Save results to Excel
  write.xlsx(results_df, file = paste0("results_", var_name, "_Control_vs_MitoD_Fed.xlsx"))
}



### ----------- END OF PSYCHOSOCIAL CODE ----------- ###



## Roc Curve for 2 predictors FGF21 + GDF15 combined ------

# for FED FGF21 and GDF15 


library(tidyverse)
library(janitor)
library(dplyr)
library(tidyr)
library(stringr)
install.packages("matrixStats")
library(matrixStats)
install.packages("readxl")
packageVersion("dplyr")
install.packages("writexl")
library(writexl)
library(readxl)
library(purrr)


install.packages("pROC")
install.packages("readr")

# Load libraries
library(readr)
library(pROC)

# Load data
data <- read_xlsx("Biomarker_timepoint_sep.xlsx")

head(data)

# Fit logistic regression model with GF and FG as predictors
model <- glm(Outcome ~ FGF21 + GDF15, data = data, family = binomial)

# Get predicted probabilities
data$predicted_prob <- predict(model, type = "response")

# Generate ROC curve
roc_curve <- roc(data$Outcome, data$predicted_prob)

# Plot ROC curve with percentage axes
plot(roc_curve, main = "ROC Curve", col = "blue", lwd = 2, xlab = "False Positive Rate (%)", ylab = "True Positive Rate (%)", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
axis(1, at = seq(0, 1, by = 0.1), labels = seq(0, 100, by = 10))
axis(2, at = seq(0, 1, by = 0.1), labels = seq(0, 100, by = 10))
box()  # Draw the box around the plot

# Print AUC
auc_value <- auc(roc_curve)
print(paste("AUC:", auc_value))

# Save data with predicted probabilities (optional)
write_xlsx(data, "data_with_predictions_5mintimepoint.xlsx")


### ----------- END OF ROC CODE ----------- ###
