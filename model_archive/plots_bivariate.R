

# Bivariate animal model plots: 






# Heritability comparison plot 


# 1. Extract Heritabilities
# Bivariate
h2_wing_biv <- samples_biv$heritability1
h2_beak_biv <- samples_biv$heritability2

# Latent Variable
h2_wing_lv <- samples_lv$h2_traits[,1]
h2_beak_lv <- samples_lv$h2_traits[,2]
h2_size_lv <- samples_lv$h2_psi

# 2. Combine into Data Frame
df_comp <- data.frame(
  Value = c(h2_wing_biv, h2_wing_lv, h2_beak_biv, h2_beak_lv, h2_size_lv),
  Model = c(rep("Bivariate", length(h2_wing_biv)), rep("Latent Var", length(h2_wing_lv)),
            rep("Bivariate", length(h2_beak_biv)), rep("Latent Var", length(h2_beak_lv)),
            rep("Latent Factor", length(h2_size_lv))),
  Trait = c(rep("Wing", length(h2_wing_biv) * 2), 
            rep("Beak", length(h2_beak_biv) * 2),
            rep("Size Factor", length(h2_size_lv)))
)

# 3. Plot
p_comp <- ggplot(df_comp, aes(x = Value, fill = Model)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Trait, scales = "free_y") + # Separate panels for Wing, Beak, and Factor
  scale_fill_manual(values = c("Bivariate" = "#999999", 
                               "Latent Var" = "#E69F00", 
                               "Latent Factor" = "#009E73")) +
  labs(title = "Comparison of Heritability Estimates",
       subtitle = "Latent Variable model recovers similar heritabilities to Bivariate model",
       x = "Heritability (h^2)", y = "Density") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("Plot_Heritability_Comparison.png", p_comp, width = 8, height = 4)
print(p_comp)




# Correlation comparison plot with manual values

library(ggplot2)


df_cor <- data.frame(
  Type = c("Genetic Correlation (r_A)", "Perm. Env. Correlation (r_PE)", "Residual Correlation (r_R)"),
  Mean = c(0.56, 0.56, 0.42),
  Lower = c(0.29, 0.47, 0.39),
  Upper = c(0.76, 0.64, 0.45)
)

# 2. Plot
p_cor <- ggplot(df_cor, aes(x = Type, y = Mean, color = Type)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, size = 1.2) +
  labs(title = "Correlations between Wing Length and Beak Length",
       subtitle = "Comparison of Genetic, Permanent Environmental, and Residual levels",
       y = "Correlation Coefficient (r)", x = NULL) +
  ylim(0, 1) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11, face = "bold"))

# Save
ggsave("Plot_Bivariate_Correlations.png", p_cor, width = 7, height = 5)
print(p_cor)





# Correlation Densities 

library(rstan)
library(ggplot2)
library(bayesplot)
library(dplyr)
library(tidyr)

# 1. Load the Bivariate Model
fit_biv <- readRDS("Output_Bivariate_Animal_Model.rds")
samples_biv <- rstan::extract(fit_biv)

df_cors <- data.frame(
  Genetic = samples_biv$cor_A[,1,2],
  PermEnv = samples_biv$cor_E[,1,2],
  Residual = samples_biv$cor_R[,1,2]
) %>% 
  pivot_longer(cols = everything(), names_to = "Type", values_to = "Correlation")

# Plot
p_cor_dens <- ggplot(df_cors, aes(x = Correlation, fill = Type)) +
  geom_density(alpha = 0.6, color = NA) +
  scale_fill_manual(values = c("Genetic" = "#D55E00",    # Red/Orange
                               "PermEnv" = "#56B4E9",    # Blue
                               "Residual" = "#999999")) + # Grey
  labs(title = "Posterior Distributions of Trait Correlations",
       subtitle = "Genetic (rA) and Permanent Env (rPE) correlations are stronger than Residual (rR)",
       x = "Correlation Coefficient (r)", y = "Density") +
  xlim(-0.2, 1.0) + # Focus on the relevant range
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("Plot_Bivariate_Cor_Density.png", p_cor_dens, width = 7, height = 5)
print(p_cor_dens)

# ==============================================================================
# PLOT 2: Variance Partitioning (Wing vs. Beak)
# ==============================================================================


vars_mean <- data.frame(
  Trait = rep(c("Wing Length", "Beak Length"), each = 3),
  Component = rep(c("Additive Genetic", "Permanent Env", "Residual"), 2),
  Value = c(
    mean(samples_biv$Sigma_A[,1,1]), mean(samples_biv$Sigma_E[,1,1]), mean(samples_biv$Sigma_R[,1,1]),
    mean(samples_biv$Sigma_A[,2,2]), mean(samples_biv$Sigma_E[,2,2]), mean(samples_biv$Sigma_R[,2,2])
  )
)

vars_mean <- vars_mean %>%
  group_by(Trait) %>%
  mutate(Proportion = Value / sum(Value)) %>%
  ungroup()

vars_mean$Component <- factor(vars_mean$Component, 
                              levels = c("Residual", "Permanent Env", "Additive Genetic"))

p_var_stack <- ggplot(vars_mean, aes(x = Trait, y = Proportion, fill = Component)) +
  geom_col(width = 0.5, color = "black", alpha = 0.8) +
  scale_fill_manual(values = c("Additive Genetic" = "#D55E00", 
                               "Permanent Env" = "#E69F00", 
                               "Residual" = "#999999")) +
  labs(title = "Variance Architecture: Wing vs. Beak",
       subtitle = "Comparing the sources of phenotypic variation",
       y = "Proportion of Total Variance", x = NULL) +
  theme_minimal() +
  theme(legend.position = "right")

ggsave("Plot_Bivariate_Variance_Stack.png", p_var_stack, width = 6, height = 5)
print(p_var_stack)





