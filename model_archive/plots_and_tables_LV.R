
# Plots for latent variable animal model: 




library(ggplot2)
library(dplyr)
library(rstan)


# ==============================================================================
# PLOT: Factor Loadings (Zoomed In)
# ==============================================================================

df_loadings <- data.frame(
  Trait = c("Wing Length", "Beak Length"),
  Mean = c(mean(loading_wing), mean(loading_beak)),
  Lower = c(quantile(loading_wing, 0.025), quantile(loading_beak, 0.025)),
  Upper = c(quantile(loading_wing, 0.975), quantile(loading_beak, 0.975))
)

p_loadings <- ggplot(df_loadings, aes(x = Trait, y = Mean)) +
  geom_point(size = 4, color = "darkblue") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.15, size = 1, color = "darkblue") +
  labs(title = "Factor Loadings on Latent 'Body Size'",
       subtitle = "Standardized loadings showing correlation with the latent factor",
       y = "Loading Coefficient", x = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        axis.title.y = element_text(margin = margin(r = 10)))

ggsave("Plot_LV_Loadings_Zoomed.png", p_loadings, width = 5, height = 4)
print(p_loadings)








# Plot for the heritability and 1 - h^2 of LV but called V_A and V_E

library(ggplot2)
library(dplyr)
library(tidyr)
library(rstan)

# 1. Extract Posterior Samples
# We extract the variance components of the Latent Variable (psi)
# defined in your Generated Quantities block
out_lv <- readRDS("Output_LV_Animal_Model_Correct_Sorting.rds")


post_samples_lv <- rstan::extract(out_lv, pars = c("var_psi_a", "var_psi_e"))

# 2. Organize Data
df_post_lv <- data.frame(
  var_psi_a = post_samples_lv$var_psi_a,
  var_psi_e = post_samples_lv$var_psi_e
) %>%
  pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Variance") %>%
  mutate(Parameter = recode(Parameter, 
                            "var_psi_a" = "Latent Genetic Variance (VA)",
                            "var_psi_e" = "Latent Environmental Variance (VE)"))

# 3. Generate the Plot
ggplot(df_post_lv, aes(x = Variance, fill = Parameter)) +
  # Plot Posterior Densities
  geom_density(alpha = 0.6, color = NA) +
  
  # Formatting
  scale_fill_brewer(palette = "Set1") + # Red/Blue contrast usually works well
  labs(x = "Proportion of Latent Variance (0 to 1)", 
       y = "Density",
       title = "Variance Partitioning of the Latent Variable") +
  
  # Constrain X axis to 0-1 (Since total variance is fixed to 1)
  coord_cartesian(xlim = c(0, 1)) +
  
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )



# Plotting the variances scaled back from the LV scale to each original scale
library(latex2exp)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rstan)

# ==============================================================================
# 1. SETUP & DATA PREPARATION
# ==============================================================================

# Load your slimmed RDS (Population Parameters)
out_lv <- readRDS("Output_LV_Final_Relevant.rds")
samples <- rstan::extract(out_lv)

# DEFINE SCALES: 
# Using the variables already in your environment
scale_wing <- sd_wing^2  
scale_beak <- sd_beak^2 

# Calculate Variances on Original Scale
df_wing <- data.frame(
  Trait = "Wing Length",
  VA = samples$lambda[,1]^2 * samples$var_psi_a * scale_wing,
  VE = samples$lambda[,1]^2 * samples$var_psi_e * scale_wing,
  VR = samples$sd_R[,1]^2 * scale_wing
)

df_beak <- data.frame(
  Trait = "Beak Length",
  VA = samples$lambda[,2]^2 * samples$var_psi_a * scale_beak,
  VE = samples$lambda[,2]^2 * samples$var_psi_e * scale_beak,
  VR = samples$sd_R[,2]^2 * scale_beak
)

# Combine into one Master Data Frame
df_all <- bind_rows(df_wing, df_beak) %>%
  pivot_longer(cols = c(VA, VE, VR), names_to = "Component", values_to = "Variance") %>%
  mutate(Component = factor(Component, levels = c("VA", "VE", "VR")))

# Define the LaTeX Labels for consistency across both plots
my_labels <- c(
  "VA" = TeX(r'($\sigma_{A(LV)}^2$)'),
  "VE" = TeX(r'($\sigma_{E(LV)}^2$)'),
  "VR" = TeX(r'($\sigma_{R(Trait)}^2$)')
)

# ==============================================================================
# 2. PLOT A: DENSITY PLOT (Only VA and VE)
# ==============================================================================

df_density <- df_all %>% filter(Component != "VR")

p1 <- ggplot(df_density, aes(x = Variance, fill = Component)) +
  geom_density(alpha = 0.6, color = NA) +
  facet_wrap(~Trait, scales = "free") + 
  scale_fill_brewer(palette = "Set1", labels = my_labels) +
  labs(
       x = "Variance (original scale)", 
       y = "Density") +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold")
  )

print(p1)


# ==============================================================================
# 3. PLOT B: FOREST PLOT (VA, VE, and VR)
# ==============================================================================

df_forest <- df_all %>%
  group_by(Trait, Component) %>%
  summarise(
    Mean = mean(Variance),
    Lower = quantile(Variance, 0.025),
    Upper = quantile(Variance, 0.975),
    .groups = "drop"
  ) %>%
  mutate(Component = factor(Component, levels = c("VR", "VE", "VA")))

p2 <- ggplot(df_forest, aes(x = Component, y = Mean, color = Component)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, size = 1) +
  geom_point(size = 4) +
  facet_wrap(~Trait, scales = "free_x") + 
  coord_flip() + 
  scale_color_brewer(palette = "Set1", labels = my_labels) +
  scale_x_discrete(labels = my_labels) + 
  labs(
       y = "Variance (original scale)", 
       x = NULL) +
  theme_bw() +
  theme(
    legend.position = "none", 
    axis.text.y = element_text(size = 12, color = "black"), 
    strip.text = element_text(size = 14, face = "bold"),
    panel.grid.major.y = element_blank()
  )

print(p2)





# plotting just heritability for latent variable 

library(ggplot2)
library(dplyr)
library(rstan)

# 1. Extract Posterior Samples
# 'h2_psi' is the heritability of the latent variable
post_samples_h2 <- rstan::extract(out_lv, pars = "h2_psi")$h2_psi

# 2. Organize Data
df_h2 <- data.frame(Heritability = post_samples_h2)

# 3. Generate the Plot
ggplot(df_h2, aes(x = Heritability)) +
  # Plot Posterior Density
  geom_density(fill = "#69b3a2", alpha = 0.7, color = NA) +
  
  # Add a vertical line for the mean (optional but helpful)
  geom_vline(aes(xintercept = mean(Heritability)), 
             color = "black", linetype = "dashed", size = 0.8) +
  
  # Formatting
  labs(x = expression(paste("Heritability (", h^2, ") of Latent Variable")), 
       y = "Density",
       title = "Posterior Density of Latent Variable Heritability") +
  
  # Constrain X axis to valid range 0-1
  coord_cartesian(xlim = c(0, 1)) +
  
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# 4. Print Summary Statistics to Console
cat("Mean Heritability:", round(mean(post_samples_h2), 3), "\n")
cat("95% CI:", round(quantile(post_samples_h2, c(0.025, 0.975)), 3), "\n")




# Plot to compare heritabilities from LV model and univariate model

library(ggplot2)
library(dplyr)
library(tidyr)
library(rstan)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

# A. Load Latent Variable Model
out_lv <- readRDS("Output_LV_Final_Relevant.rds")
samps_lv <- rstan::extract(out_lv)

# B. Load Univariate Models
out_uni_wing <- readRDS("Output_Univariate_Animal_Model_Priors.rds") # <--- REPLACE with actual filename
out_uni_beak <- readRDS("Output_Univariate_Animal_Model_Priors_Beak.rds") # 


samps_uni_wing <- rstan::extract(out_uni_wing)
samps_uni_beak <- rstan::extract(out_uni_beak)

# ==============================================================================
# 2. PREPARE LATENT VARIABLE HERITABILITIES (Shared Genetics)
# ==============================================================================

# Wing (Trait 1)
va_wing_lv <- samps_lv$lambda[,1]^2 * samps_lv$var_psi_a
vp_wing_lv <- (samps_lv$lambda[,1]^2 * 1) + samps_lv$sd_R[,1]^2
h2_wing_lv <- va_wing_lv / vp_wing_lv

# Beak (Trait 2)
va_beak_lv <- samps_lv$lambda[,2]^2 * samps_lv$var_psi_a
vp_beak_lv <- (samps_lv$lambda[,2]^2 * 1) + samps_lv$sd_R[,2]^2
h2_beak_lv <- va_beak_lv / vp_beak_lv

df_lv <- data.frame(
  "Wing Length" = h2_wing_lv,
  "Beak Length" = h2_beak_lv,
  check.names = FALSE
) %>%
  pivot_longer(cols = everything(), names_to = "Trait", values_to = "Heritability") %>%
  mutate(Model = "Latent Variable (Shared)")

# ==============================================================================
# 3. PREPARE UNIVARIATE HERITABILITIES (Total Genetics)
# ==============================================================================

df_uni_wing <- data.frame(
  Heritability = samps_uni_wing$heritability,
  Trait = "Wing Length",
  Model = "Univariate (Total)"
)

df_uni_beak <- data.frame(
  Heritability = samps_uni_beak$heritability,
  Trait = "Beak Length",
  Model = "Univariate (Total)"
)

# Combine all data
df_plot <- bind_rows(df_lv, df_uni_wing, df_uni_beak)

# ==============================================================================
# 4. FIX ORDERING AND PLOT
# ==============================================================================

df_plot$Trait <- factor(df_plot$Trait, levels = c("Wing Length", "Beak Length"))

ggplot(df_plot, aes(x = Heritability, fill = Model)) +
  geom_density(alpha = 0.5, color = NA) +
  facet_wrap(~Trait, scales = "free") +
  
  scale_fill_manual(values = c("Latent Variable (Shared)" = "#E41A1C", 
                               "Univariate (Total)" = "#377EB8")) +
  
  # Removed Title and Subtitle
  labs(x = expression(paste("Heritability (", h^2, ")")),
       y = "Density") +
  
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12)
  )





# plot to compare the covariances between lv model and bivariate model

library(ggplot2)
library(dplyr)
library(tidyr)
library(rstan)


scale_cov <- sd_wing * sd_beak 

# ==============================================================================
# 2. EXTRACT FROM LATENT VARIABLE MODEL (Implied Covariances)
# ==============================================================================
out_lv <- readRDS("Output_LV_Final_Relevant.rds")
samps_lv <- rstan::extract(out_lv)

# Calculate Implied Covariances
# Formula: Lambda1 * Lambda2 * Variance_Component * Scale_Factor
cov_a_lv <- samps_lv$lambda[,1] * samps_lv$lambda[,2] * samps_lv$var_psi_a * scale_cov
cov_e_lv <- samps_lv$lambda[,1] * samps_lv$lambda[,2] * samps_lv$var_psi_e * scale_cov

cov_r_lv <- rep(0, length(cov_a_lv))

df_lv <- data.frame(
  Cov_A = cov_a_lv,
  Cov_E = cov_e_lv,
  Cov_R = cov_r_lv
) %>%
  pivot_longer(everything(), names_to = "Component", values_to = "Covariance") %>%
  mutate(Model = "Latent Variable (Implied)")

# ==============================================================================
# 3. EXTRACT FROM BIVARIATE MODEL (Direct Estimates)
# ==============================================================================
out_biv <- readRDS("Output_Bivariate_Final_Relevant.rds")
samps_biv <- rstan::extract(out_biv)


cov_a_biv <- samps_biv$Sigma_A[, 1, 2]
cov_e_biv <- samps_biv$Sigma_E[, 1, 2]
cov_r_biv <- samps_biv$Sigma_R[, 1, 2]

df_biv <- data.frame(
  Cov_A = cov_a_biv,
  Cov_E = cov_e_biv,
  Cov_R = cov_r_biv
) %>%
  pivot_longer(everything(), names_to = "Component", values_to = "Covariance") %>%
  mutate(Model = "Bivariate (Direct)")

# ==============================================================================
# 4. PLOT COMPARISON
# ==============================================================================
df_plot <- bind_rows(df_lv, df_biv) %>%
  mutate(Component = factor(Component, levels = c("Cov_A", "Cov_E", "Cov_R")))

ggplot(df_plot, aes(x = Covariance, fill = Model)) +
  geom_density(alpha = 0.5, color = NA) +
  facet_wrap(~Component, scales = "free") +
  
  scale_fill_manual(values = c("Latent Variable (Implied)" = "#E41A1C", 
                               "Bivariate (Direct)" = "#377EB8")) +
  
  labs(title = "Covariance Comparison: LV vs. Bivariate Model",
       subtitle = "Checking if the single factor captures the correlation structure",
       x = expression(paste("Covariance (mm"^2, ")")), 
       y = "Density") +
  
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12)
  )





# covariances compared but with cov_R = 0 for lv


library(ggplot2)
library(dplyr)
library(tidyr)
library(rstan)

# ==============================================================================
# 1. PREPARE DATA
# ==============================================================================

my_labeller <- as_labeller(c(
  "Cov_A" = "Genetic Covariance",
  "Cov_E" = "Environmental Covariance",
  "Cov_R" = "Residual Covariance"
))

# ==============================================================================
# 2. SPLIT DATA
# ==============================================================================

# Subset 1: All Densities EXCEPT the fixed LV Residual
df_density <- df_plot %>%
  filter(!(Model == "Latent Variable (Implied)" & Component == "Cov_R"))

# Subset 2: The Fixed Line for LV Residual
df_vline <- data.frame(
  Component = "Cov_R",
  xint = 0,
  Model = "Latent Variable (Implied)"
)

# ==============================================================================
# 3. GENERATE PLOT
# ==============================================================================

ggplot() +
  # A. Draw Density Curves
  geom_density(data = df_density, aes(x = Covariance, fill = Model), 
               alpha = 0.5, color = NA) +
  
  # B. Draw Vertical Line (Only for LV Residual)
  geom_vline(data = df_vline, aes(xintercept = xint, color = Model),
             size = 1.0, linetype = "dashed") +
  
  # C. Faceting
  facet_wrap(~Component, scales = "free", labeller = my_labeller) +
  
  # Styling
  scale_fill_manual(values = c("Latent Variable (Implied)" = "#E41A1C", 
                               "Bivariate (Direct)" = "#377EB8")) +
  scale_color_manual(values = c("Latent Variable (Implied)" = "#E41A1C")) +
  
  labs(x = expression(paste("Covariance (mm"^2, ")")), 
       y = "Density") +
  
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12)
  )






# biplot/ordination plot 



library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# ==============================================================================
# 1. LOAD DATA & PREPARE INDIVIDUAL SCORES
# ==============================================================================
if (!exists("A_hestmannoy")) stop("Error: 'A_hestmannoy' matrix not found.")

# A. Load Summary
df_summary <- read.csv("Output_LV_Summary_Correct_Sorting.csv", row.names = 1)

# B. Extract Individual Scores
psi_rows <- df_summary %>%
  filter(grepl("^psi\\[", rownames(.))) %>% 
  filter(!grepl("psi_a", rownames(.))) %>%   
  filter(!grepl("psi_e", rownames(.))) %>%   
  select(mean, X2.5., X97.5.) %>%           
  rename(Mean = mean, Lower = X2.5., Upper = X97.5.)

psi_rows$ringnr <- rownames(A_hestmannoy)

# C. Prepare Sex Data & Merge
# traitData <- read.csv("morphology.txt", sep = ";", header = TRUE, na.strings=c("NA",""))
clean_sex <- traitData %>%
  mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr)))) %>%
  mutate(sex = case_when(
    !is.na(scriptsex) & scriptsex %in% c("m", "pm") ~ "Male",
    !is.na(scriptsex) & scriptsex %in% c("f", "pf") ~ "Female",
    fieldsex %in% c("m", "pm") ~ "Male",
    fieldsex %in% c("f", "pf") ~ "Female",
    TRUE ~ NA_character_
  )) %>%
  dplyr::select(ringnr, sex) %>%
  distinct(ringnr, .keep_all = TRUE) %>%
  filter(!is.na(sex))

df_plot_all <- psi_rows %>% inner_join(clean_sex, by = "ringnr")

# Sample 50 random birds
set.seed(19) 
df_birds <- df_plot_all %>% 
  sample_n(50) %>%
  mutate(y_jitter = runif(n(), min = -0.08, max = 0.08))

# ==============================================================================
# 2. PREPARE ARROW DATA
# ==============================================================================

# --- A. Extract Loadings ---
loadings_data <- df_summary[c("lambda1", "lambda_free[1]"), c("mean", "X2.5.", "X97.5.")]
colnames(loadings_data) <- c("Mean", "Lower", "Upper")
lambda_wing <- loadings_data["lambda1", ]
lambda_beak <- loadings_data["lambda_free[1]", ]

scale_factor <- 1.0 

# --- B. Trait Arrows (Only these remain) ---
# Positions: -0.14 and -0.20 (Gap is 0.06, much tighter)
df_trait_arrows <- data.frame(
  Label = c("Wing Length", "Beak Length"),
  x_start = 0,
  x_end = c(lambda_wing$Mean, lambda_beak$Mean) * scale_factor,
  x_min = c(lambda_wing$Lower, lambda_beak$Lower) * scale_factor, 
  x_max = c(lambda_wing$Upper, lambda_beak$Upper) * scale_factor, 
  y_pos = c(-0.14, -0.20), 
  Color = "black", 
  fontface = "bold",
  hjust_val = -0.1 
)

# ==============================================================================
# 3. GENERATE THE PLOT
# ==============================================================================

ggplot() +
  # --- A. The Birds ---
  geom_errorbarh(data = df_birds, aes(xmin = Lower, xmax = Upper, y = y_jitter, color = sex), 
                 height = 0, alpha = 0.4, size = 0.4) +
  geom_point(data = df_birds, aes(x = Mean, y = y_jitter, color = sex, shape = sex),
             size = 2.5, alpha = 0.8) +
  
  # --- B. Zero Line ---
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", size=0.5) +
  
  # --- C. Arrows ---
  geom_segment(data = df_trait_arrows, aes(x = x_start, xend = x_end, y = y_pos, yend = y_pos),
               arrow = arrow(length = unit(0.2, "cm")), size = 0.8, color = "black") +
  
  # --- D. Trait Loading Error Bars ---
  geom_errorbarh(data = df_trait_arrows, 
                 aes(xmin = x_min, xmax = x_max, y = y_pos), 
                 height = 0.02, size = 0.5, color = "black") + 
  
  # --- E. Labels ---
  geom_text(data = df_trait_arrows, 
            aes(x = x_end, y = y_pos, label = Label, hjust = hjust_val),
            color = "black", fontface = "bold", vjust = 0.5, size = 3.5) +
  
  # --- F. Formatting ---
  scale_y_continuous(limits = c(-0.25, 0.12), breaks = NULL) + 
  
  scale_x_continuous(breaks = seq(-1.5, 1.5, 0.5)) +
  coord_cartesian(xlim = c(-1.5, 1.5)) +
  
  scale_color_manual(values = c("Female" = "#E41A1C", "Male" = "#377EB8"),
                     breaks = c("Female", "Male")) +
  scale_shape_manual(values = c(16, 17),
                     breaks = c("Female", "Male")) +
  
  labs(title = NULL,
       x = "Latent Variable (Standardized Units)",
       y = NULL,
       color = "Sex", shape = "Sex") +
  
  theme_minimal() +
  theme(
    axis.line.x = element_line(size = 0.8), 
    panel.grid.major.y = element_blank(), 
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(), 
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    axis.title.x = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 10)
  )









# also ordination plot/biplot but sorted the random individuals by the lv value


library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# ==============================================================================
# 1. LOAD DATA & PREPARE INDIVIDUAL SCORES
# ==============================================================================
if (!exists("A_hestmannoy")) stop("Error: 'A_hestmannoy' matrix not found.")

# A. Load Summary
df_summary <- read.csv("Output_LV_Summary_Correct_Sorting.csv", row.names = 1)

# B. Extract Individual Scores
psi_rows <- df_summary %>%
  filter(grepl("^psi\\[", rownames(.))) %>% 
  filter(!grepl("psi_a", rownames(.))) %>%   
  filter(!grepl("psi_e", rownames(.))) %>%   
  select(mean, X2.5., X97.5.) %>%           
  rename(Mean = mean, Lower = X2.5., Upper = X97.5.)

psi_rows$ringnr <- rownames(A_hestmannoy)

# C. Prepare Sex Data & Merge
# Ensure traitData is in your environment
if (!exists("traitData")) stop("Error: 'traitData' not found.")

clean_sex <- traitData %>%
  mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr)))) %>%
  mutate(sex = case_when(
    !is.na(scriptsex) & scriptsex %in% c("m", "pm") ~ "Male",
    !is.na(scriptsex) & scriptsex %in% c("f", "pf") ~ "Female",
    fieldsex %in% c("m", "pm") ~ "Male",
    fieldsex %in% c("f", "pf") ~ "Female",
    TRUE ~ NA_character_
  )) %>%
  dplyr::select(ringnr, sex) %>%
  distinct(ringnr, .keep_all = TRUE) %>%
  filter(!is.na(sex))

df_plot_all <- psi_rows %>% inner_join(clean_sex, by = "ringnr")

# --- MODIFIED SECTION START ---
set.seed(19) 
df_birds <- df_plot_all %>% 
  sample_n(50) %>%
  # 1. Sort by the Mean Latent Value (Low to High)
  arrange(Mean) %>% 
  # 2. Assign Y-values sequentially instead of randomly
  # This spreads them evenly from bottom (-0.08) to top (0.08) of the band
  mutate(y_jitter = seq(from = -0.08, to = 0.08, length.out = n()))
# --- MODIFIED SECTION END ---

# ==============================================================================
# 2. PREPARE ARROW DATA
# ==============================================================================

# --- A. Extract Loadings ---
loadings_data <- df_summary[c("lambda1", "lambda_free[1]"), c("mean", "X2.5.", "X97.5.")]
colnames(loadings_data) <- c("Mean", "Lower", "Upper")
lambda_wing <- loadings_data["lambda1", ]
lambda_beak <- loadings_data["lambda_free[1]", ]

scale_factor <- 1.0 

# --- B. Trait Arrows ---
df_trait_arrows <- data.frame(
  Label = c("Wing Length", "Beak Length"),
  x_start = 0,
  x_end = c(lambda_wing$Mean, lambda_beak$Mean) * scale_factor,
  x_min = c(lambda_wing$Lower, lambda_beak$Lower) * scale_factor, 
  x_max = c(lambda_wing$Upper, lambda_beak$Upper) * scale_factor, 
  y_pos = c(-0.14, -0.20), 
  Color = "black", 
  fontface = "bold",
  hjust_val = -0.1 
)

# ==============================================================================
# 3. GENERATE THE PLOT
# ==============================================================================

ggplot() +
  # --- A. The Birds ---
  geom_errorbarh(data = df_birds, aes(xmin = Lower, xmax = Upper, y = y_jitter, color = sex), 
                 height = 0, alpha = 0.4, size = 0.4) +
  geom_point(data = df_birds, aes(x = Mean, y = y_jitter, color = sex, shape = sex),
             size = 2.5, alpha = 0.8) +
  
  # --- B. Zero Line ---
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", size=0.5) +
  
  # --- C. Arrows ---
  geom_segment(data = df_trait_arrows, aes(x = x_start, xend = x_end, y = y_pos, yend = y_pos),
               arrow = arrow(length = unit(0.2, "cm")), size = 0.8, color = "black") +
  
  # --- D. Trait Loading Error Bars ---
  geom_errorbarh(data = df_trait_arrows, 
                 aes(xmin = x_min, xmax = x_max, y = y_pos), 
                 height = 0.02, size = 0.5, color = "black") + 
  
  # --- E. Labels ---
  geom_text(data = df_trait_arrows, 
            aes(x = x_end, y = y_pos, label = Label, hjust = hjust_val),
            color = "black", fontface = "bold", vjust = 0.5, size = 3.5) +
  
  # --- F. Formatting ---
  scale_y_continuous(limits = c(-0.25, 0.12), breaks = NULL) + 
  
  scale_x_continuous(breaks = seq(-1.5, 1.5, 0.5)) +
  coord_cartesian(xlim = c(-1.5, 1.5)) +
  
  scale_color_manual(values = c("Female" = "#E41A1C", "Male" = "#377EB8"),
                     breaks = c("Female", "Male")) +
  scale_shape_manual(values = c(16, 17),
                     breaks = c("Female", "Male")) +
  
  labs(title = NULL,
       x = "Latent Variable (Standardized Units)",
       y = NULL,
       color = "Sex", shape = "Sex") +
  
  theme_minimal() +
  theme(
    axis.line.x = element_line(size = 0.8), 
    panel.grid.major.y = element_blank(), 
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(), 
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    axis.title.x = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 10)
  )















# plot to compare variances from lv model and univariate animal model


library(ggplot2)
library(dplyr)
library(patchwork) 

clean_theme <- theme_minimal() + 
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5), # Compact title
    axis.title.x = element_blank(),     # Remove X label from individual plots (add back later)
    legend.position = "none",           # Hide legend for now
    panel.grid.minor = element_blank()  # Cleaner look
  )

# 2. Create the 4 Individual Plots

# --- Plot A: Beak (Genetic) ---
p_beak_gen <- df_plot %>% 
  filter(Trait == "Beak Length", Component == "Additive Genetic Variance") %>%
  ggplot(aes(x = Variance, fill = Model)) +
  geom_density(alpha = 0.6, color = NA) +
  labs(title = "Beak: Genetic") + 
  clean_theme

# --- Plot B: Wing (Genetic) ---
p_wing_gen <- df_plot %>% 
  filter(Trait == "Wing Length", Component == "Additive Genetic Variance") %>%
  ggplot(aes(x = Variance, fill = Model)) +
  geom_density(alpha = 0.6, color = NA) +
  labs(title = "Wing: Genetic") + 
  clean_theme

# --- Plot C: Beak (Environment) ---
p_beak_env <- df_plot %>% 
  filter(Trait == "Beak Length", Component == "Permanent Env. Variance") %>%
  ggplot(aes(x = Variance, fill = Model)) +
  geom_density(alpha = 0.6, color = NA) +
  labs(title = "Beak: Environment") + 
  clean_theme

# --- Plot D: Wing (Environment) ---
p_wing_env <- df_plot %>% 
  filter(Trait == "Wing Length", Component == "Permanent Env. Variance") %>%
  ggplot(aes(x = Variance, fill = Model)) +
  geom_density(alpha = 0.6, color = NA) +
  labs(title = "Wing: Environment") + 
  clean_theme

# 3. Combine them
final_plot <- (p_beak_gen | p_wing_gen) / (p_beak_env | p_wing_env) +
  
  # Add the tags automatically
  plot_annotation(tag_levels = 'a') + 
  
  # Collect the legend at the bottom
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom", 
        legend.title = element_blank()) # Remove "Model" title from legend if obvious

# 4. Add a single X-axis label at the very bottom (Optional but cleaner)
final_plot_labeled <- wrap_elements(panel = final_plot) + 
  labs(tag = "Variance (mm²)") +
  theme(
    plot.tag = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0),
    plot.tag.position = "bottom"
  )

print(final_plot_labeled)

ggsave("variance_panel_plot.pdf", final_plot_labeled, width = 7, height = 5)
