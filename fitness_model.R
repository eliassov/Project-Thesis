
library(dplyr)

fitness_data <- read.csv("Fitness.csv", header = TRUE, stringsAsFactors = FALSE, na.strings = "NA", fileEncoding = "Windows-1252")

fitness_data <- fitness_data %>%
  dplyr::select(Location, Year_ID, Year, Sex, Own_survival, N_Recruits, Least_age) %>%
  filter(grepl("^Hestmann", Location))


if (!require("lme4")) install.packages("lme4")
library(lme4)

# 1. Final Data Prep
fitness_data_clean <- fitness_data %>%
  mutate(
    # Extract the ring number by removing everything up to and including the "_"
    ringnr = sub("^.*_", "", Year_ID),
    
    # Ensure survival is strictly numeric (0 and 1)
    Own_survival = as.numeric(Own_survival),
    
    # Creating age_class: 1 = Juvenile, 2+ = Adult
    age_class = if_else(Least_age == 1, "Juv", "Ad"),
    # Ensuring factors match the morphology data
    Sex = factor(Sex),
    age_class = factor(age_class),
    Year = factor(Year)
      )
  

# 2. Fit the Baseline Survival GLMM
# We use family = binomial for 0/1 survival data
surv_baseline <- glmer(Own_survival ~ 1 + (1 | ringnr) + (1 | Year_factor), 
                       data = fitness_data_clean, 
                       family = binomial(link = "logit"))

# 3. Look at the variance components
summary(surv_baseline)









repro_baseline <- glmer(N_Recruits ~ 1 + (1 | ringnr) + (1 | Year), 
                        data = fitness_data_clean, 
                        family = poisson(link = "log"))

# Look at the variance components
summary(repro_baseline)







# 1. Fit the Reduced Model (Remove (1|ringnr))
# For Survival
surv_reduced <- glmer(Own_survival ~ 1 + (1 | Year_factor), 
                      data = fitness_data_clean, 
                      family = binomial)

# 2. Compare them
anova(surv_reduced, surv_baseline)






# 1. Fit the Reduced Model (Remove (1|ringnr))
# For Survival
repro_reduced <- glmer(N_Recruits ~ 1 + (1 | Year_factor), 
                      data = fitness_data_clean, 
                      family = poisson(link="log"))

# 2. Compare them
anova(repro_reduced, repro_baseline)








library(lme4)

# Fit the Poisson model
rep_mod <- glmer(N_Recruits ~ Sex + age_class + (1 | Year), 
                 family = poisson(link = "log"), 
                 data = fitness_data_clean)

# 1. Summary of effects
summary(rep_mod)

# 2. Overdispersion Check (Crucial for Stan)
# If this number is much larger than 1, we need Negative Binomial in Stan
pearson_resid <- resid(rep_mod, type = "pearson")
overdispersion <- sum(pearson_resid^2) / df.residual(rep_mod)

cat("\n--- OVERDISPERSION RATIO ---\n")
print(overdispersion)







library(ggplot2)

# 1. Get the exact percentage of zeros
zero_pct <- sum(fitness_data_clean$N_Recruits == 0) / nrow(fitness_data_clean)
cat("Percentage of birds with zero recruits:", round(zero_pct * 100, 1), "%\n")

# 2. Make a clean bar plot for the meeting
ggplot(fitness_data_clean, aes(x = N_Recruits)) +
  geom_bar(fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Number of Recruits",
       x = "Number of Recruits (Lifetime)",
       y = "Frequency (Number of Birds)") +
  scale_x_continuous(breaks = 0:max(fitness_data_clean$N_Recruits))





# 1. Extract the Pearson residuals
pearson_resid <- resid(repro_baseline, type = "pearson")

# 2. Calculate the ratio (Sum of squared residuals / degrees of freedom)
overdispersion <- sum(pearson_resid^2) / df.residual(repro_baseline)

cat("\n--- OVERDISPERSION RATIO ---\n")
print(overdispersion)



library(DHARMa)

# 1. Simulate the residuals from your fitted reproduction model
# (Make sure 'repro_baseline' is the glmer model you ran earlier)
sim_res_repro <- simulateResiduals(fittedModel = repro_baseline, plot = FALSE)

# 2. Plot the standardized residuals
plot(sim_res_repro)

# 3. The Ultimate Question for Bob: Test specifically for Zero-Inflation
testZeroInflation(sim_res_repro)








library(dplyr)
library(lme4)

# 1. Create a "Body Size" proxy using PCA on your 5 scaled traits
pca_result <- prcomp(trait_subset_lv[, c("ving_h_std", "tars_h_std", "nebb_l_std", "nebb_h_std", "vekt_std")], 
                     center = TRUE, scale. = TRUE)

# Extract PC1 and average it per bird (in case a bird was caught multiple times)
trait_size <- trait_subset_lv %>%
  mutate(Size_Index = pca_result$x[, 1]) %>%
  dplyr::select(ringnr, Size_Index) %>%
  group_by(ringnr) %>%
  summarise(Size_Index = mean(Size_Index, na.rm = TRUE))

# 2. Merge with your clean Fitness Data
survival_test_data <- fitness_data_clean %>%
  inner_join(trait_size, by = "ringnr")

# 3. Fit the Survival Selection Model
# Size_Index is a FIXED effect. Year is the environmental RANDOM effect. 
# We drop (1 | ringnr) because we already proved it has zero variance.
surv_selection_mod <- glmer(Own_survival ~ Size_Index + (1 | Year), 
                            data = survival_test_data, 
                            family = binomial(link = "logit"))

# 4. Check the results for Bob
summary(surv_selection_mod)
