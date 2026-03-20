# ==============================================================================
# 1. SETUP & LIBRARIES
# ==============================================================================
personal_lib <- "~/R/library"
if (!dir.exists(personal_lib)) dir.create(personal_lib, recursive = TRUE)
.libPaths(c(personal_lib, .libPaths()))

required_packages <- c("rstan", "tidyverse", "nadiv", "MCMCglmm") 
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, lib = personal_lib, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# ==============================================================================
# 2. LOAD & CLEAN DATA (OPTIMIZED ORDER)
# ==============================================================================

# Pedigree
pedData <- read.csv("pedigree_data.txt", sep = " ", header = TRUE, na.strings = "NA")
pedigree_clean <- pedData %>%
  dplyr::select(ringnr, dam, sire) %>%
  mutate(across(c(ringnr, dam, sire), ~ as.character(trimws(gsub("_.*$", "", .)))),
         dam = if_else(is.na(dam) | dam == "" | dam == "NA", NA_character_, dam),
         sire = if_else(is.na(sire) | sire == "" | sire == "NA", NA_character_, sire))

# External Parents
all_parents <- unique(c(pedigree_clean$dam, pedigree_clean$sire))
external_parents <- setdiff(all_parents, c(pedigree_clean$ringnr, NA_character_))
if (length(external_parents) > 0) {
  founder_df <- data.frame(ringnr = external_parents, dam = NA_character_, sire = NA_character_)
  pedigree_full <- bind_rows(founder_df, pedigree_clean) %>% distinct(ringnr, .keep_all = TRUE)
} else {
  pedigree_full <- pedigree_clean
}

pedigree_prepped <- prepPed(pedigree_full, gender = NULL, check = TRUE)

# Trait Data
traitData <- read.csv("morphology.txt", sep = ";", header = TRUE, 
                      na.strings = c("NA", ""), fileEncoding = "Windows-1252", 
                      stringsAsFactors = FALSE, colClasses = c(ringnr = "character"))

# A. LOAD & LIGHT FORMATTING (Crucial for matching)
trait_raw <- traitData %>%
  filter(sted_r == "hestmannøy", 
         !is.na(ving_h), !is.na(nebb_l), !is.na(tars_h), !is.na(nebb_h), !is.na(vekt)) %>%
  mutate(ringnr = as.character(trimws(gsub("_.*$", "", ringnr))))

# B. ALIGN WITH PEDIGREE (The "Gatekeeper" Step)
valid_ids <- intersect(trait_raw$ringnr, pedigree_prepped$ringnr)
trait_subset_temp <- trait_raw %>% filter(ringnr %in% valid_ids)

cat("Birds in Raw File:    ", nrow(trait_raw), "\n")
cat("Birds with Pedigree:  ", nrow(trait_subset_temp), "\n")

# C. DEEP CLEANING (Sex, Age, Bad IDs) on the SUBSET
bad_ids <- c("8N27558", "8N42527", "8N13933", "8M71874", "8N87712", "8L89527")

trait_subset_temp <- trait_subset_temp %>%
  filter(!ringnr %in% bad_ids) %>%
  mutate(raw_sex = if_else(!is.na(scriptsex) & scriptsex != "u", scriptsex, fieldsex)) %>%
  mutate(sex = case_when(
    raw_sex %in% c("m", "pm") ~ "m", 
    raw_sex %in% c("f", "pf") ~ "f", 
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(sex)) %>%
  mutate(
    s_age = trimws(as.character(scriptage)),
    f_age = trimws(as.character(fieldage)),
    age_class = case_when(
      grepl("^1K", s_age) ~ "juvenile",
      grepl("^[2-9]K", s_age) | grepl("^1[0-9]K", s_age) ~ "adult",
      f_age %in% c("j", "pj") ~ "juvenile",
      f_age %in% c("a", "pa") ~ "adult",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(age_class))

# D. STATISTICAL OUTLIER REMOVAL (Residual Check)
outlier_ids_to_kill <- c()
traits_to_check <- c("ving_h", "tars_h", "nebb_l", "nebb_h", "vekt")

for(t in traits_to_check) {
  f <- as.formula(paste(t, "~ sex + age_class"))
  fit <- lm(f, data = trait_subset_temp, na.action = na.exclude)
  
  resids <- rstudent(fit)
  cutoff <- if(t == "vekt") 6.0 else 4.0 
  
  bad_rows <- which(abs(resids) > cutoff)
  
  if(length(bad_rows) > 0) {
    new_bad_ids <- trait_subset_temp$ringnr[bad_rows]
    outlier_ids_to_kill <- c(outlier_ids_to_kill, new_bad_ids)
    cat("Trait:", t, "- Removing", length(new_bad_ids), "IDs\n")
  }
}

# E. FINAL DATASET
trait_subset_lv <- trait_subset_temp %>%
  filter(!ringnr %in% unique(outlier_ids_to_kill))

cat("Final Clean N: ", nrow(trait_subset_lv), "\n")

# # Z-score standardization
# trait_subset_lv <- trait_subset_lv %>%
#   mutate(
#     ving_h_std = scale(ving_h),
#     nebb_l_std = scale(nebb_l),
#     tars_h_std = scale(tars_h),
#     nebb_h_std = scale(nebb_h),
#     vekt_std   = scale(vekt)
#   )



# NEW: Log-transform FIRST, then scale
trait_subset_lv <- trait_subset_lv %>%
  mutate(
    ving_h_std = as.numeric(scale(log(ving_h))),
    nebb_l_std = as.numeric(scale(log(nebb_l))),
    tars_h_std = as.numeric(scale(log(tars_h))),
    nebb_h_std = as.numeric(scale(log(nebb_h))),
    vekt_std   = as.numeric(scale(log(vekt)))
  )


cat("Data Loaded. N Rows:", nrow(trait_subset_lv), "\n")




# 1. Convert to Date (R handles YYYY-MM-DD automatically)
trait_subset_lv$date_obj <- as.Date(trait_subset_lv$dato)

# 2. Extract month
trait_subset_lv$month <- as.numeric(format(trait_subset_lv$date_obj, "%m"))



# 1. Extract the Year as a factor
trait_subset_lv <- trait_subset_lv %>%
  mutate(
    Year = factor(format(date_obj, "%Y")),
    # Create the 3-level Season we discussed
    Season_3 = factor(case_when(
      month %in% c(5, 6, 7) ~ "Breeding",
      month %in% c(8, 9)    ~ "Molt",
      TRUE                  ~ "Winter"
    ))
  )





# # Assuming trait_subset_lv has: sex, age_class, Season_3, Year, and traits
# traits_to_check <- c("ving_h", "tars_h", "nebb_l", "nebb_h")
# 
# for (t in traits_to_check) {
#   cat("\n\n==============================================\n")
#   cat("ANALYZING TRAIT:", t, "\n")
#   cat("==============================================\n")
#   
#   # 1. Fit the model (using log-transformed values as decided)
#   # We use REML = TRUE for variance/random effect estimation
#   form <- as.formula(paste0("log(", t, ") ~ sex + age_class + Season_3 + (1 | Year)"))
#   m_reml <- lmer(form, data = trait_subset_lv)
#   
#   # 2. Fixed Effects Summary (p-values)
#   # This tells you if Sex, Age, or Season are significant
#   print(summary(m_reml))
#   
#   # 3. Random Effect Significance (Year)
#   # Using the manual Likelihood Ratio Test approach
#   m_lm <- lm(as.formula(paste0("log(", t, ") ~ sex + age_class + Season_3")), data = trait_subset_lv)
#   LRT_stat <- as.numeric(2 * (logLik(m_reml) - logLik(m_lm)))
#   p_val_year <- pchisq(LRT_stat, df = 1, lower.tail = FALSE) / 2
#   cat("\n--- Random Effect (Year) P-Value: ", p_val_year, "\n")
#   
#   # 4. Variance Partitioning (R-Squared)
#   var_fixed  <- var(predict(m_reml, re.form = NA))
#   var_random <- as.numeric(VarCorr(m_reml)$Year)
#   var_resid  <- sigma(m_reml)^2
#   var_total  <- var_fixed + var_random + var_resid
#   
#   cat("\n--- Variance Explained ---\n")
#   cat("Fixed Effects (R2m):    ", round((var_fixed / var_total) * 100, 2), "%\n")
#   cat("Year Effect:            ", round((var_random / var_total) * 100, 2), "%\n")
#   cat("Leftover Individual:    ", round((var_resid / var_total) * 100, 2), "%\n")
# }
# 





# n_measurers <- length(unique(trait_subset_lv$init))
# cat("Total unique measurers:", n_measurers, "\n")

# 2. Check observations per measurer
# measurer_counts <- table(trait_subset_lv$init)
# # print(sort(measurer_counts, decreasing = TRUE))
# 
# # 3. Percentage of data covered by top measurers
# prop_table <- prop.table(sort(measurer_counts, decreasing = TRUE)) * 100
# print(head(prop_table, 10))



# Grouping low-volume measurers using case_when
trait_subset_lv <- trait_subset_lv %>%
  group_by(init) %>%
  mutate(init_clean = case_when(
    n() < 10 ~ "OTHER",
    TRUE     ~ as.character(init)
  )) %>%
  ungroup() %>%
  mutate(init_clean = factor(init_clean))


###############################################################

############### DID I NOT ADD YEAR HERE?????? ##################

##############################################################

############# =============================================




# # Verify it worked
# cat("Number of levels in init_clean:", length(levels(trait_subset_lv$init_clean)), "\n")
# 
# 
# table(trait_subset_lv$init_clean)




# library(lme4)
# 
# # Full model with Measurer (init_clean)
# m_with_init <- lmer(log(nebb_h) ~ sex + age_class + Season_3 + (1 | Year) + (1 | init_clean), 
#                     data = trait_subset_lv, REML = FALSE)
# 
# # Reduced model without Measurer
# m_no_init <- lmer(log(nebb_h) ~ sex + age_class + Season_3 + (1 | Year), 
#                   data = trait_subset_lv, REML = FALSE)
# 
# # Likelihood Ratio Test
# anova(m_no_init, m_with_init)



# # 2. Run the comparison models for Body Mass (vekt)
# # This includes the Year random effect and the new Season fixed effect
# library(lme4)
# 
# # Linear version
# lin_mod <- lmer(vekt ~ sex + age_class + Season_3 + (1 | Year), data = trait_subset_lv)
# 
# # Log-transformed version
# log_mod <- lmer(log(vekt) ~ sex + age_class + Season_3 + (1 | Year), data = trait_subset_lv)
# 
# # 3. Final visual check
# par(mfrow=c(1,2))
# qqnorm(resid(lin_mod), main="Linear Vekt QQ")
# qqline(resid(lin_mod), col="red")
# 
# qqnorm(resid(log_mod), main="Log(Vekt) QQ")
# qqline(resid(log_mod), col="red")
# 
# 
# 
# 
# # 1. Compare Residual Distributions (Histograms)
# par(mfrow=c(2,2))
# 
# # Linear Residuals
# hist(resid(lin_mod), breaks=50, main="Linear Resids: Histogram",
#      xlab="Residual Value", col="lightblue")
# plot(fitted(lin_mod), resid(lin_mod), main="Linear: Resid vs Fit",
#      xlab="Fitted Values", ylab="Residuals", pch=20, col=rgb(0,0,0,0.2))
# abline(h=0, col="red")
# 
# # Log Residuals
# hist(resid(log_mod), breaks=50, main="Log Resids: Histogram",
#      xlab="Residual Value", col="lightgreen")
# plot(fitted(log_mod), resid(log_mod), main="Log: Resid vs Fit",
#      xlab="Fitted Values", ylab="Residuals", pch=20, col=rgb(0,0,0,0.2))
# abline(h=0, col="red")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# traits_remaining <- c("ving_h", "tars_h", "nebb_l", "nebb_h")
# 
# # Set up a 4x2 plotting grid
# par(mfrow=c(4,2), mar=c(4,4,2,1))
# 
# for (t in traits_remaining) {
#   # Fit the log-model with our new fixed/random effects
#   form <- as.formula(paste0("log(", t, ") ~ sex + age_class + Season_3 + (1 | Year)"))
#   mod <- lmer(form, data = trait_subset_lv)
# 
#   # Histogram
#   hist(resid(mod), breaks=30, main=paste("Log-Hist:", t),
#        col="lightgreen", xlab="Residuals")
# 
#   # QQ Plot
#   qqnorm(resid(mod), main=paste("Log-QQ:", t))
#   qqline(resid(mod), col="red")
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# traits_remaining <- c("ving_h", "tars_h", "nebb_l", "nebb_h")
# 
# # Set up plotting grid (4 rows, 4 columns: Hist-Lin, Hist-Log, QQ-Lin, QQ-Log)
# par(mfrow=c(4,4), mar=c(3,3,2,1), oma=c(0,0,2,0))
# 
# for (t in traits_remaining) {
#   # Formulas
#   form_lin <- as.formula(paste0(t, " ~ sex + age_class + Season_3 + (1 | Year)"))
#   form_log <- as.formula(paste0("log(", t, ") ~ sex + age_class + Season_3 + (1 | Year)"))
# 
#   # Fit Models
#   m_lin <- lmer(form_lin, data = trait_subset_lv)
#   m_log <- lmer(form_log, data = trait_subset_lv)
# 
#   # 1. Linear Histogram
#   hist(resid(m_lin), main=paste("Lin-Hist:", t), col="lightblue", breaks=30)
#   # 2. Log Histogram
#   hist(resid(m_log), main=paste("Log-Hist:", t), col="lightgreen", breaks=30)
# 
#   # 3. Linear QQ
#   qqnorm(resid(m_lin), main=paste("Lin-QQ:", t), cex=0.5)
#   qqline(resid(m_lin), col="red")
#   # 4. Log QQ
#   qqnorm(resid(m_log), main=paste("Log-QQ:", t), cex=0.5)
#   qqline(resid(m_log), col="red")
# }
# mtext("Comparison of Linear vs Log Residuals", outer=TRUE, cex=1.2)
# 
# 
# 
# 
# 
# 
# traits_remaining <- c("ving_h", "tars_h", "nebb_l", "nebb_h")
# 
# # Set up a 4x2 grid to see Linear vs Log side-by-side for each trait
# par(mfrow=c(4,2), mar=c(4,4,2,1))
# 
# for (t in traits_remaining) {
#   # Formulas
#   form_lin <- as.formula(paste0(t, " ~ sex + age_class + Season_3 + (1 | Year)"))
#   form_log <- as.formula(paste0("log(", t, ") ~ sex + age_class + Season_3 + (1 | Year)"))
# 
#   # Fit Models
#   m_lin <- lmer(form_lin, data = trait_subset_lv)
#   m_log <- lmer(form_log, data = trait_subset_lv)
# 
#   # Linear Cloud
#   plot(fitted(m_lin), resid(m_lin),
#        main=paste("Linear Cloud:", t),
#        xlab="Fitted Values", ylab="Residuals",
#        pch=20, col=rgb(0,0,0,0.2))
#   abline(h=0, col="red", lwd=2)
# 
#   # Log Cloud
#   plot(fitted(m_log), resid(m_log),
#        main=paste("Log Cloud:", t),
#        xlab="Fitted Values", ylab="Residuals",
#        pch=20, col=rgb(0,0,0,0.2))
#   abline(h=0, col="red", lwd=2)
# }
# 
# 
# 
# 
# 
# library(lme4)
# 
# # ---------------------------------------------------------
# # 1. TEST FIXED EFFECTS (Using Likelihood Ratio Tests)
# # Note: We MUST use REML = FALSE when testing fixed effects
# # ---------------------------------------------------------
# # The Full Model
# mod_full <- lmer(log(vekt) ~ sex + age_class + Season_3 + (1 | Year),
#                  data = trait_subset_lv, REML = FALSE)
# 
# # The Reduced Models (Drop one effect at a time)
# mod_no_sex    <- lmer(log(vekt) ~ age_class + Season_3 + (1 | Year), data = trait_subset_lv, REML = FALSE)
# mod_no_age    <- lmer(log(vekt) ~ sex + Season_3 + (1 | Year), data = trait_subset_lv, REML = FALSE)
# mod_no_season <- lmer(log(vekt) ~ sex + age_class + (1 | Year), data = trait_subset_lv, REML = FALSE)
# 
# cat("\n--- P-VALUES FOR FIXED EFFECTS ---\n")
# print(anova(mod_no_sex, mod_full))    # Tests Sex
# print(anova(mod_no_age, mod_full))    # Tests Age
# print(anova(mod_no_season, mod_full)) # Tests Season
# 
# # ---------------------------------------------------------
# # 2. TEST RANDOM EFFECT (Year)
# # ---------------------------------------------------------
# # We refit the full model using standard REML=TRUE for variances
# mod_reml <- lmer(log(vekt) ~ sex + age_class + Season_3 + (1 | Year), data = trait_subset_lv)
# 
# # A model with NO random effect (just a standard linear model)
# mod_lm <- lm(log(vekt) ~ sex + age_class + Season_3, data = trait_subset_lv)
# 
# # Exact Likelihood Ratio Test for the Random Effect
# LRT_stat <- as.numeric(2 * (logLik(mod_reml) - logLik(mod_lm)))
# # We divide the p-value by 2 because variances cannot be negative (bounded at 0)
# p_val_year <- pchisq(LRT_stat, df = 1, lower.tail = FALSE) / 2
# 
# cat("\n--- P-VALUE FOR YEAR RANDOM EFFECT ---\n")
# cat("Likelihood Ratio Test p-value:", p_val_year, "\n")
# 
# # ---------------------------------------------------------
# # 3. MANUAL R-SQUARED (Variance Explained)
# # ---------------------------------------------------------
# # Extract the three sources of variance:
# var_fixed  <- var(predict(mod_reml, re.form = NA)) # Variance from Sex + Age + Season
# var_random <- as.numeric(VarCorr(mod_reml)$Year)   # Variance from Year
# var_resid  <- sigma(mod_reml)^2                    # Leftover (Individual) variance
# 
# var_total <- var_fixed + var_random + var_resid
# 
# # Calculate exact percentages
# R2_marginal    <- var_fixed / var_total
# R2_conditional <- (var_fixed + var_random) / var_total
# 
# cat("\n--- VARIANCE EXPLAINED ---\n")
# cat("Fixed Effects alone (R2 Marginal):    ", round(R2_marginal * 100, 2), "%\n")
# cat("Fixed + Random Year (R2 Conditional): ", round(R2_conditional * 100, 2), "%\n")
# cat("Leftover Individual Variance:         ", round((var_resid / var_total) * 100, 2), "%\n")
# 
# 
# 




# ==============================================================================
# 3. PEDIGREE PROCESSING
# ==============================================================================
pedigree_working <- pedigree_prepped
colnames(pedigree_working) <- c("id", "dam", "sire") 
pedigree_working$id   <- as.character(pedigree_working$id)
pedigree_working$dam  <- as.character(pedigree_working$dam)
pedigree_working$sire <- as.character(pedigree_working$sire)

ped_pruned <- MCMCglmm::prunePed(pedigree_working, keep = trait_subset_lv$ringnr)
A_inv_obj <- MCMCglmm::inverseA(ped_pruned, nodes = "ALL", scale = FALSE)
ordered_ids <- rownames(A_inv_obj$Ainv)

ped_ordered <- data.frame(
  id = ordered_ids,
  dam = as.character(A_inv_obj$pedigree[, 2]), 
  sire = as.character(A_inv_obj$pedigree[, 3]),
  stringsAsFactors = FALSE
)

dii_vector <- A_inv_obj$dii
Na <- nrow(ped_ordered)

dam_idx_mapped  <- match(ped_ordered$dam, ped_ordered$id)
sire_idx_mapped <- match(ped_ordered$sire, ped_ordered$id)

dam_idx_mapped[is.na(dam_idx_mapped)]   <- Na + 1
sire_idx_mapped[is.na(sire_idx_mapped)] <- Na + 1

id_map <- setNames(1:Na, ped_ordered$id)

# # ==============================================================================
# # 4. DATA MAPPING & SAFETY CHECK
# # ==============================================================================
# animal_idx <- id_map[trait_subset_lv$ringnr]
# missing_mask <- is.na(animal_idx)
# 
# if (sum(missing_mask) > 0) {
#   cat("WARNING: Removing", sum(missing_mask), "rows. ID mismatch persists.\n")
#   trait_subset_lv <- trait_subset_lv[!missing_mask, ]
#   animal_idx <- id_map[trait_subset_lv$ringnr]
# }
# 
# Y_mat <- as.matrix(trait_subset_lv[, c("tars_h_std", "ving_h_std", "nebb_l_std", "nebb_h_std", "vekt_std")])
# X_mat <- model.matrix(~ sex + age_class, data = trait_subset_lv)
# 
# cat("==============================================\n")
# cat("        FINAL DATA VERIFICATION CHECK          \n")
# cat("==============================================\n")
# cat("1. Total Data Points (N_obs):     ", nrow(trait_subset_lv), "\n")
# cat("   (This is the number of rows Stan analyzes)\n\n")
# cat("2. Unique Phenotyped Birds:       ", length(unique(trait_subset_lv$ringnr)), "\n")
# cat("   (This is the actual number of birds you caught)\n\n")
# cat("3. Total Pedigree Size (Na):      ", Na, "\n")
# cat("   (This is Phenotyped Birds + Ancestors)\n")
# cat("==============================================\n")
# 
# # ==============================================================================
# # 5. STAN EXECUTION
# # ==============================================================================
# dataset_lv_final <- list(
#   No = nrow(trait_subset_lv),
#   Nt = 5, 
#   Nlv = 2,
#   Y = Y_mat,
#   X = X_mat,
#   K = ncol(X_mat),
#   animal = as.array(as.integer(animal_idx)),
#   Na = Na,
#   dam = as.array(as.integer(dam_idx_mapped)), 
#   sire = as.array(as.integer(sire_idx_mapped)),
#   dii = as.array(dii_vector)
# )
# 
# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())
# 
# tomonitor_lv <- c(
#   "beta",         
#   "Lambda",       
#   "sd_R",         
#   "h2_lv",
#   "rho",
#   "LV",
#   "h2_trait",
#   "pe2_trait",
#   "var_G_trait",
#   "var_PE_trait",
#   "var_R_trait",
#   "var_P_total",
#   "lp__"          
# )
# 
# if (any(dam_idx_mapped[dam_idx_mapped <= Na] >= (1:Na)[dam_idx_mapped <= Na])) {
#   stop("CRITICAL ERROR: A Dam appears in the list AFTER her child. Recursion will fail.")
# }
# if (any(sire_idx_mapped[sire_idx_mapped <= Na] >= (1:Na)[sire_idx_mapped <= Na])) {
#   stop("CRITICAL ERROR: A Sire appears in the list AFTER his child.")
# }
# 
# init_fn_shrinkage <- function() {
#   L_init <- matrix(rnorm(dataset_lv_final$Nt * dataset_lv_final$Nlv, 0, 0.1), 
#                    nrow = dataset_lv_final$Nt, 
#                    ncol = dataset_lv_final$Nlv)
#   
#   L_init[1, 1] <- runif(1, 0.8, 1.1) 
#   L_init[1, 2] <- runif(1, -0.1, 0.1) 
#   
#   list(
#     lambda_raw = L_init, 
#     rho = runif(1, 0.4, 0.6),        
#     h2_lv = runif(2, 0.4, 0.6),      
#     sd_R = runif(dataset_lv_final$Nt, 0.8, 1.2), 
#     beta = matrix(rnorm(dataset_lv_final$K * dataset_lv_final$Nt, 0, 1),
#                   nrow = dataset_lv_final$K),
#     w_a = matrix(rnorm(dataset_lv_final$Nlv * dataset_lv_final$Na, 0, 0.1),
#                  nrow = dataset_lv_final$Nlv),
#     w_pe = matrix(rnorm(dataset_lv_final$Nlv * dataset_lv_final$Na, 0, 0.1),
#                   nrow = dataset_lv_final$Nlv)
#   )
# }
# 
# out_lv <- stan(
#   file = 'good_shrink_prior.stan', 
#   data = dataset_lv_final,
#   init = init_fn_shrinkage,       
#   chains = 4, 
#   pars = tomonitor_lv,
#   control = list(adapt_delta = 0.9),
#   iter = 8000,          
#   warmup = 4000
# )
# 
# saveRDS(out_lv, 'Output_shrink_prior.rds')





# ==============================================================================
# 4. DATA MAPPING & FIXED/RANDOM EFFECTS SETUP
# ==============================================================================
animal_idx <- id_map[trait_subset_lv$ringnr]

# 1. Update the Fixed Effects Matrix (Added Season_3)
X_mat <- model.matrix(~ sex + age_class + Season_3, data = trait_subset_lv)

# 2. Extract Random Effects Arrays
year_factor <- as.integer(factor(trait_subset_lv$Year))
init_factor <- as.integer(factor(trait_subset_lv$init_clean))

# 3. Trait Matrix (Make sure they are the LOG SCALED ones)
Y_mat <- as.matrix(trait_subset_lv[, c("tars_h_std", "ving_h_std", "nebb_l_std", "nebb_h_std", "vekt_std")])

# ==============================================================================
# 5. STAN EXECUTION (UPDATED LIST)
# ==============================================================================
dataset_lv_final <- list(
  No = nrow(trait_subset_lv),
  Nt = 5, 
  Nlv = 2,
  Y = Y_mat,
  X = X_mat,
  K = ncol(X_mat),
  animal = as.array(as.integer(animal_idx)),
  Na = Na,
  dam = as.array(as.integer(dam_idx_mapped)), 
  sire = as.array(as.integer(sire_idx_mapped)),
  dii = as.array(dii_vector),
  
  # NEW: Pass the Random Effects data
  N_year = max(year_factor),
  year_idx = as.array(year_factor),
  N_init = max(init_factor),
  init_idx = as.array(init_factor)
)


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Update Initialization to include the new parameters
init_fn_shrinkage <- function() {
  L_init <- matrix(rnorm(dataset_lv_final$Nt * dataset_lv_final$Nlv, 0, 0.1), 
                   nrow = dataset_lv_final$Nt, 
                   ncol = dataset_lv_final$Nlv)
  L_init[1, 1] <- runif(1, 0.8, 1.1) 
  L_init[1, 2] <- runif(1, -0.1, 0.1) 
  
  list(
    lambda_raw = L_init, 
    rho = runif(1, 0.4, 0.6),        
    h2_lv = runif(2, 0.4, 0.6),      
    sd_R = runif(dataset_lv_final$Nt, 0.5, 1.0), 
    beta = matrix(rnorm(dataset_lv_final$K * dataset_lv_final$Nt, 0, 1),
                  nrow = dataset_lv_final$K),
    w_a = matrix(rnorm(dataset_lv_final$Nlv * dataset_lv_final$Na, 0, 0.1),
                 nrow = dataset_lv_final$Nlv),
    w_pe = matrix(rnorm(dataset_lv_final$Nlv * dataset_lv_final$Na, 0, 0.1),
                  nrow = dataset_lv_final$Nlv),
    
    # NEW: Initialize Random Effects
    sd_year = runif(dataset_lv_final$Nt, 0.05, 0.2),
    sd_init = runif(dataset_lv_final$Nt, 0.05, 0.2),
    z_year = matrix(rnorm(dataset_lv_final$N_year * dataset_lv_final$Nt, 0, 0.1),
                    nrow = dataset_lv_final$N_year),
    z_init = matrix(rnorm(dataset_lv_final$N_init * dataset_lv_final$Nt, 0, 0.1),
                    nrow = dataset_lv_final$N_init)
  )
}

# Run the updated model
out_lv_shaved <- stan(
  file = 'shrink_prior_fixed_random.stan', 
  data = dataset_lv_final,
  init = init_fn_shrinkage,        
  chains = 4, 
  pars = c("beta", "Lambda", "sd_R", "sd_year", "sd_init", "h2_lv", "rho", 
           "h2_trait", "pe2_trait", "var_G_trait", "var_PE_trait", 
           "var_year_trait", "var_init_trait", "var_R_trait", "var_P_total", "lp__"),
  control = list(adapt_delta = 0.95), # Bumped slightly for stability with extra params
  iter = 10000,          
  warmup = 5000,
  sample_file = "output_chain"
)

saveRDS(out_lv_shaved, 'Output_shrink_improved.rds')





# ==============================================================================
# 6. RESULTS & DIAGNOSTICS
# ==============================================================================
out_lv <- readRDS('Output_shrink_improved.rds')

print(out_lv,
      pars = c("beta", "Lambda", "sd_R", "sd_year", "sd_init", "h2_lv", "rho", 
               "h2_trait", "pe2_trait", "var_G_trait", "var_PE_trait", 
               "var_year_trait", "var_init_trait", "var_R_trait", "var_P_total", "lp__"),      
      probs = c(0.025, 0.975),
      digits_summary = 3)

# library(rstan)
# 
# draws <- as.array(out_lv)
# 
# for (chain in 1:4) {
#   chain_mean <- mean(draws[, chain, "Lambda[4,2]"])
#   
#   if (chain_mean < 0) {
#     cat("Flipping signs for Factor 2 in Chain", chain, "\n")
#     draws[, chain, "Lambda[1,2]"] <- draws[, chain, "Lambda[1,2]"] * -1
#     draws[, chain, "Lambda[2,2]"] <- draws[, chain, "Lambda[2,2]"] * -1
#     draws[, chain, "Lambda[3,2]"] <- draws[, chain, "Lambda[3,2]"] * -1
#     draws[, chain, "Lambda[4,2]"] <- draws[, chain, "Lambda[4,2]"] * -1
#     draws[, chain, "Lambda[5,2]"] <- draws[, chain, "Lambda[5,2]"] * -1
#   }
# }
# 
# corrected_summary <- rstan::monitor(draws, warmup = 0, print = FALSE)
# summary_df <- as.data.frame(corrected_summary)
# lambdas_only <- summary_df[grep("Lambda", rownames(summary_df)), ]
# print(lambdas_only[, c("mean", "sd", "2.5%", "97.5%", "n_eff", "Rhat")], digits = 3)
# 
# beta_means   <- summary_df[grep("^beta\\[", rownames(summary_df)), "mean"]
# Lambda_means <- summary_df[grep("^Lambda\\[", rownames(summary_df)), "mean"]
# LV_means     <- summary_df[grep("^LV\\[", rownames(summary_df)), "mean"]
# 
# B_hat  <- matrix(beta_means,   nrow = dataset_lv_final$K,  ncol = dataset_lv_final$Nt)
# L_hat  <- matrix(Lambda_means, nrow = dataset_lv_final$Nt, ncol = dataset_lv_final$Nlv)
# LV_hat <- matrix(LV_means,     nrow = dataset_lv_final$Na, ncol = dataset_lv_final$Nlv)
# 
# LV_obs <- LV_hat[dataset_lv_final$animal, ]
# Mu_hat <- (dataset_lv_final$X %*% B_hat) + (LV_obs %*% t(L_hat))
# Residuals <- dataset_lv_final$Y - Mu_hat
# 
# cat("Residual Correlation Matrix:\n")
# resid_cor_matrix <- cor(Residuals)
# print(round(resid_cor_matrix, 2))
# 
# trait_names <- c("Tarsus", "Wing", "Beak L", "Beak H", "Mass")
# par(mfrow = c(2, 3)) 
# 
# for(i in 1:5) {
#   qqnorm(Residuals[, i], main = paste("QQ Plot:", trait_names[i]))
#   qqline(Residuals[, i], col = "red", lwd = 2)
# }
# 
# par(mfrow = c(2, 3))
# 
# for(i in 1:5) {
#   plot(Mu_hat[, i], Residuals[, i], 
#        main = paste("Resid vs Fit:", trait_names[i]),
#        xlab = "Predicted Value", 
#        ylab = "Residual",
#        pch = 16, col = rgb(0,0,0,0.2)) 
#   abline(h = 0, col = "red", lwd = 2)
# }



# out_old <- readRDS('Output_Joint_Starter.rds')
