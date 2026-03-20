
# Univariate animal model plots: 




# Plot without priors of the variance components with density 

library(ggplot2)
library(dplyr)
library(tidyr)

out_1 <- readRDS("Output_Univariate_Animal_Model_Priors.rds")

# 1. Extract Posterior Samples
post_samples <- rstan::extract(out_1, pars = c("var_A", "var_E", "var_R"))

# 2. Organize Data
df_post <- data.frame(
  var_A = post_samples$var_A,
  var_E = post_samples$var_E,
  var_R = post_samples$var_R
) %>%
  pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Variance") %>%
  mutate(Parameter = recode(Parameter, 
                            "var_A" = "VA",
                            "var_E" = "VE",
                            "var_R" = "VR"))

# 3. Generate the Plot (No Priors)
ggplot(df_post, aes(x = Variance, fill = Parameter)) +
  # Plot Posterior Densities
  geom_density(alpha = 0.6, color = NA) +
  
  # Formatting
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Variance", y = "Density") +
  theme_minimal() +
  theme(
    legend.position = "top",          # Puts legend at the top
    legend.title = element_blank(),   # Removes "Parameter" title from legend
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )




# plot comparing variance components of univariate model for wing length and beak length

library(ggplot2)
library(dplyr)
library(tidyr)
library(rstan)

# 1. LOAD DATA
out_wing <- readRDS("Output_Univariate_Animal_Model_Priors.rds")
out_beak <- readRDS("Output_Univariate_Animal_Model_Priors_Beak.rds")

# 2. HELPER FUNCTION
get_variance_df <- function(stan_out, trait_name) {
  # Extract samples
  samps <- rstan::extract(stan_out, pars = c("var_A", "var_E", "var_R"))
  
  # Format into a clean dataframe
  data.frame(
    var_A = samps$var_A,
    var_E = samps$var_E,
    var_R = samps$var_R
  ) %>%
    mutate(Trait = trait_name) %>%
    pivot_longer(cols = c(var_A, var_E, var_R), 
                 names_to = "Component", 
                 values_to = "Variance")
}

# 3. PREPARE DATA
df_wing <- get_variance_df(out_wing, "(a) Wing Length")
df_beak <- get_variance_df(out_beak, "(b) Beak Length")

df_combined <- bind_rows(df_wing, df_beak) %>%
  mutate(Component = recode(Component, 
                            "var_A" = "Additive Genetic",
                            "var_E" = "Permanent Environment",
                            "var_R" = "Residual"))

# 4. PLOT
p_variances <- ggplot(df_combined, aes(x = Variance, fill = Component)) +
  geom_density(alpha = 0.6, color = NA) +
  
  # Create two panels with independent X-axis scales
  facet_wrap(~Trait, scales = "free") +
  
  # Colors and Labels
  scale_fill_brewer(palette = "Set2") + 
  labs(x = expression(Variance~(mm^2)), 
       y = "Posterior Density") +
  
  # Clean Theme
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10)
  )

print(p_variances)

ggsave("univariate_variance_components.pdf", p_variances, width = 8, height = 4)
