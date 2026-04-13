library(tidyverse)
library(lme4)
library(glmmTMB)
library(ggtext)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: Load raw data
# ══════════════════════════════════════════════════════════════════════════════

clean_coral <- read_csv(
  "~/MEDS/capstone/eds-coral-figs-storage/vedika-r-files/data/updated_coral_tidy_2013-2024.csv"
)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2: Filter and clean transect names
# ══════════════════════════════════════════════════════════════════════════════

recruit_raw <- clean_coral %>%
  filter(!grepl("^P", transect)) %>%
  dplyr::select(site, habitat, transect, dyn_recruitment, year, taxa) %>%
  mutate(
    # Standardize duplicate/variant transect names
    # NOTE: TRUE ~ transect is required to preserve non-duplicate names
    transect = case_when(
      transect == "T02_OLD"     ~ "T02",
      transect == "T03_PRE2020" ~ "T03",
      transect == "T04_OLD"     ~ "T04",
      transect == "T02_2021"    ~ "T02",
      transect == "TO1"         ~ "T01",
      TRUE                      ~ transect   # keep all others unchanged
    )
  )

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3: Summarize at transect level, preserving site and habitat
#
# KEY FIX: group_by must include site and habitat so they survive summarize()
# and are available as fixed effects in the model.
#
# KEY FIX: Create transect_id = site_transect (e.g., "LTER1_T01") so the
# random effect recognizes that T01 at LTER1 ≠ T01 at LTER2.
# This implements the nested structure Max described: transects within sites.
# ══════════════════════════════════════════════════════════════════════════════

recruit_transect <- recruit_raw %>%
  group_by(taxa, site, habitat, transect, year) %>%
  summarize(
    # dyn_recruitment is 0/1 per coral row — sum gives total recruits present
    # in this transect-site-year-taxa combination
    n_recruits = sum(dyn_recruitment, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Standardize to recruits per 5 m² by dividing raw count by 5
    recruits_std = n_recruits / 5,
    recruits_int = round(recruits_std),        # integer for negative binomial
    pre_post     = case_when(                  # binary year period (per Max)
      year >= 2013 & year <= 2018 ~ 0L,
      year >= 2019 & year <= 2024 ~ 1L
    ),
    year_s       = (year - min(year)) / (max(year) - min(year)),  # scaled 0–1
    transect_id  = paste(site, transect, sep = "_")               # unique ID
  )

# Confirm unique transect IDs look correct
sort(unique(recruit_transect$transect_id))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4: Sparsity check — identify taxa too sparse to model
# ══════════════════════════════════════════════════════════════════════════════

sparsity_check <- recruit_transect %>%
  group_by(taxa) %>%
  summarize(n_nonzero = sum(recruits_int > 0), .groups = "drop")

print(sparsity_check)
# Acr: very sparse — degenerate model, exclude from GLMM curves
# Mil: zero nonzero — no model possible, exclude entirely
# Poc: sufficient signal
# Por: sufficient signal

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5: Full GLMM — answers BOTH research questions in one model
#
# FIXED EFFECTS:
#   habitat  — 2 levels (BR, FR); fixed because you want a direct comparison
#              between back reef and fore reef. Random effects are unreliable
#              with only 2 levels (insufficient variance estimation).
#   site     — fixed because your Q2 is specifically about named-site
#              differences (LTER1 vs LTER2 etc.) and you want those contrasts.
#   taxa     — fixed; specific taxa comparisons are of interest
#   factor(year) — fixed; fewer assumptions than random year, avoids
#              confounding with intercept (per Max's notes)
#
# RANDOM EFFECTS:
#   (1 | transect_id) — transects are the sampling units nested within sites.
#              Using transect_id (e.g., LTER1_T01) ensures the model treats
#              each transect as unique, implementing the nested structure
#              Transect ⊂ Site that Max described.
#
# FAMILY: nbinom2 — negative binomial for overdispersed count data.
#         Check overdispersion after fitting (see Step 6).
# ══════════════════════════════════════════════════════════════════════════════

mod_recruit <- glmmTMB(
  recruits_int ~ habitat + site + taxa + factor(year) + (1 | transect_id),
  family    = nbinom2,
  data      = recruit_transect,
  na.action = na.exclude
)

summary(mod_recruit)

# ── Marginal effects (ggeffects recommended by Max) ──────────────────────────
# install.packages("ggeffects")  # if not installed
library(ggeffects)

# Habitat effect (averaged across sites, taxa, years)
ggeffect(mod_recruit, terms = "habitat")

# Site effect
ggeffect(mod_recruit, terms = "site")

# Year trend
ggeffect(mod_recruit, terms = "year")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 6: Model diagnostics
# Check these before trusting results (per Max's notes)
# ══════════════════════════════════════════════════════════════════════════════

# ── Overdispersion check ─────────────────────────────────────────────────────
# Ratio of sum(residuals²) / residual df; values >> 1 = overdispersed
pearson_resid <- residuals(mod_recruit, type = "pearson")
overdisp_ratio <- sum(pearson_resid^2) / df.residual(mod_recruit)
cat("Overdispersion ratio:", round(overdisp_ratio, 3), "\n")
# If >> 1: consider ziNBinom2 (zero-inflated) or check for outliers
# If << 1: underdispersed — consider nbinom1 or Poisson

# ── Temporal autocorrelation check (ACF) ─────────────────────────────────────
# Plot residuals vs. year to visually inspect for trend
recruit_transect$resid <- NA
recruit_transect$resid[!is.na(recruit_transect$recruits_int)] <- residuals(mod_recruit)

ggplot(recruit_transect, aes(x = year, y = resid)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ taxa) +
  labs(title = "Residuals vs. Year — check for temporal autocorrelation",
       x = "Year", y = "Pearson Residual") +
  theme_light()

# ACF per transect_id (looking for decay structure)
# If autocorrelation is present, consider adding ar1() correlation structure in glmmTMB

# ══════════════════════════════════════════════════════════════════════════════
# STEP 7: Spaghetti plots — per-taxa GLMMs with transect_id random effect
#
# These are the original single-taxa trend models.
# Now using transect_id (unique) instead of bare transect.
# ══════════════════════════════════════════════════════════════════════════════

dat_poc <- recruit_transect %>% filter(taxa == "Poc")
dat_por <- recruit_transect %>% filter(taxa == "Por")

mod_poc <- glmer.nb(
  recruits_int ~ year_s + (1 | transect_id),
  data = dat_poc
)
summary(mod_poc)

mod_por <- glmer.nb(
  recruits_int ~ year_s + (1 | transect_id),
  data = dat_por
)
summary(mod_por)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 8: Population-level predictions (re.form = NA = marginalizes over
# random effects, giving the average curve across all transects)
# ══════════════════════════════════════════════════════════════════════════════

make_preds <- function(model, dat, taxa_name) {
  years   <- sort(unique(dat$year))
  years_s <- (years - min(dat$year)) / (max(dat$year) - min(dat$year))
  tibble(
    taxa = taxa_name,
    year = years,
    fit  = predict(model,
                   newdata = data.frame(year_s = years_s),
                   re.form = NA,
                   type    = "response")
  )
}

pred_poc <- make_preds(mod_poc, dat_poc, "Poc")
pred_por <- make_preds(mod_por, dat_por, "Por")
pred_all <- bind_rows(pred_poc, pred_por)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 9: Bar chart summary for plotting
# ══════════════════════════════════════════════════════════════════════════════

recruit_bar <- recruit_transect %>%
  group_by(taxa, year) %>%
  summarize(mean_recruits = mean(recruits_std, na.rm = TRUE), .groups = "drop")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 10: Plot
# ══════════════════════════════════════════════════════════════════════════════

coral_colors <- c(
  Acr = "blue",
  Poc = "red",
  Por = "darkgreen",
  Mil = "darkorange"
)

taxa_labels <- c(
  Acr = "*Acropora*",
  Mil = "*Millepora*",
  Poc = "*Pocillopora*",
  Por = "*Porites*"
)

ggplot() +
  geom_col(
    data  = recruit_bar,
    aes(x = factor(year), y = mean_recruits, fill = taxa),
    alpha = 0.5, width = 0.7
  ) +
  geom_point(
    data  = recruit_transect,
    aes(x = factor(year), y = recruits_std, color = taxa),
    shape = 1, size = 1.5, alpha = 0.5
  ) +
  geom_line(
    data      = pred_all,
    aes(x = factor(year), y = fit, color = taxa, group = taxa),
    linewidth = 1.2
  ) +
  geom_point(
    data  = pred_all,
    aes(x = factor(year), y = fit, color = taxa),
    size  = 2.5
  ) +
  facet_wrap(
    ~ taxa,
    scales   = "free_y",
    ncol     = 2,
    labeller = labeller(taxa = taxa_labels)
  ) +
  scale_fill_manual(values  = coral_colors) +
  scale_color_manual(values = coral_colors) +
  scale_x_discrete(breaks = as.character(c(2013, 2015, 2017, 2019, 2021, 2023))) +
  labs(
    title   = "Recruitment Over Time",
    x       = "Year",
    y       = "Recruits per 5 m²",
    caption = "GLMM curves shown for Poc and Por only (glmer.nb, negative binomial).\nAcr and Mil excluded (insufficient recruitment signal)."
  ) +
  theme_light() +
  theme(
    axis.title       = element_text(size = 12),
    axis.text        = element_text(size = 10),
    axis.text.x      = element_text(angle = 45, vjust = 0.5),
    strip.text       = ggtext::element_markdown(size = 12),
    strip.background = element_rect(fill = "#325963"),
    legend.position  = "none",
    plot.caption     = element_text(size = 9, color = "grey50")
  )







# library(tidyverse)
# library(lme4)
# library(glmmTMB)
# # ══════════════════════════════════════════════════════════════════════════════
# # Step 1: Load and filter
# clean_coral <- read_csv(
#   "~/MEDS/capstone/eds-coral-figs-storage/vedika-r-files/data/updated_coral_tidy_2013-2024.csv")
# 
# #######_________________ recruitment raw
# 
# recruit_raw <- clean_coral %>%
#   filter(
#   # year >= 2013 & year <= 2018,
#     # site == "LTER1",
#     # habitat == "BR",
#     !grepl("^P", transect)
#   ) %>%  dplyr::select(site, habitat, transect, dyn_recruitment, year, taxa)
# 
# ### First question: Differences between sites (site is a fixed effect?)
# #recruits_pre_lter1 <- recruit_transect_pre %>% filter(site == "LTER1")
# 
# # Clean up duplicate transects
# unique(recruit_transect$transect)
# #[1] "T01"         "T02"         "T03"         "T04"         "T02_2021"    "T04_OLD"     "T02_OLD"     "T03_PRE2020" "TO1"        
# 
# # Add on duplicates
# transect_duplicates <-  c("T02_OLD","T03_PRE2020", "T04_OLD" ,"T02_2021" ,"TO1") 
# 
# recruit_raw  %>% # casewhen() ---> str_detect("1", then "transect 1")
# mutate(transect = case_when(
#   transect == "T02_OLD" ~ "T02", 
#   transect == "T03_PRE2020" ~ "T03", 
#   transect == "T04_OLD" ~"T04", 
#   transect == "T02_2021" ~ "T02", 
#   transect == "TO1" ~ "T01"
# ))
# 
# ?str_detect()
# 
# # Create unique transect IDs to truly account for nested spatial structure (transects should be spatially unique to each site)
# recruit_transect<- recruit_raw %>% mutate(
#   pre_post = case_when(
#     year >= 2013 & year <= 2018 ~ "0", 
#     year >= 2019 & year <= 2024 ~ "1"
#   )
# ) %>%
#   group_by(taxa, transect, year) %>%
#            #, site, habitat) %>%
#   summarize(n_recruits = sum(dyn_recruitment, na.rm = TRUE),
#     recruits_std = n_recruits * 5 / 25, 
#     recruits_int = round(recruits_std), 
#     .groups = "drop"
#   ) 
# 
# 
# # Step 3: mean per year for bar chart
# recruit_bar_pre <- recruit_transect %>%
#   group_by(taxa, year, transect) %>%
#   summarize(mean_recruits = mean(recruits_std, na.rm = TRUE), .groups = "drop")
# 
# # Step 4: Sparsity check to identify potential zero inflation
# sparsity_check <- recruit_transect %>%
#   group_by(taxa) %>%
#   summarize(n_nonzero = sum(recruits_int > 0), .groups = "drop") 
# 
# print(sparsity_check)
# 
# # Step 5: Fit model
# 
# # mod_recruit_lter1 <- glmer.nb(
# #   recruits_int ~ taxa + factor(year) + (1 | transect),
# #   data = recruit_transect_pre
# #   # na.action = na.exclude()
# # )
# 
# mod_recruit <- glmmTMB(
#         recruits_int ~ factor(year) + taxa + (1|transect) + (1|habitat) + site, # Update to include transect ID (no habitat, showing variance between sites)
#       # dispformula = ~1,
#       family      = nbinom2, # Both habitat and site in same model... average over another... site differences, then either ex. pick one habitat, average across sites to compare difference habitats
#       # family = negative.binomial(link = "log"),
#         data        = recruit_transect,
#         na.action   = na.exclude
#       )
# 
# summary(mod_recruit)
# 
# 
# ### Second question: Differences in habitat
# 
# # ══════════════════════════════════════════════════════════════════════════════
# # Step 2: Summarize at transect level 
# recruit_transect <- recruit_raw %>%
#   group_by(taxa, transect, year) %>%
#   summarize(n_recruits = sum(dyn_recruitment, na.rm = TRUE), 
#     recruits_std = n_recruits * 5 / 25,
#     recruits_int = round(recruits_std), 
#     .groups = "drop"
#   )
# 
# 
# # ══════════════════════════════════════════════════════════════════════════════
# # Step 3: mean per year for bar chart
# recruit_bar <- recruit_transect %>%
#   group_by(taxa, year) %>%
#   mutate(mean_recruits = mean(recruits_std, na.rm = TRUE), .groups = "drop")
# 
# # ══════════════════════════════════════════════════════════════════════════════
# # Step 4: Sparsity check to identify potential zero inflation
# sparsity_check <- recruit_transect %>%
#   group_by(taxa) %>%
#   summarize(n_nonzero = sum(recruits_int > 0), .groups = "drop") 
# 
# print(sparsity_check)
# # Acr: 2 nonzero years = model degenerate, exclude from curves
# # Mil: 0 nonzero years = no model possible, exclude from curves
# # Poc: enough signal to fit glmm
# # Por: enough signal to fit glmm
# 
# 
# 
# 
# # ══════════════════════════════════════════════════════════════════════════════
# # Step 5: Prep data with scaled year (0-1) for numerical stability
# dat_poc <- recruit_transect %>%
#   filter(taxa == "Poc") %>%
#   mutate(year_s = (year - min(year)) / (max(year) - min(year)))
# 
# dat_por <- recruit_transect %>%
#   filter(taxa == "Por") %>%
#   mutate(year_s = (year - min(year)) / (max(year) - min(year)))
# 
# dat_acr <- recruit_transect %>%
#   filter(taxa == "Acr") %>%
#   mutate(year_s = (year - min(year)) / (max(year) - min(year)))
# 
# # ══════════════════════════════════════════════════════════════════════════════
# # Step 6: Fit GLMMs (negative binomial, scaled year, random transect effect)
# recruit_transect<- recruit_transect %>% mutate(year_s = (year - min(year)) / (max(year) - min(year)))
# glmer.nb(recruits_int ~ factor(year) + (1 | transect) + (1 | site) + (1 | transect) )
# # Pocillopora 
# mod_poc <- glmer.nb(
#   recruits_int ~ year_s + (1 | transect),
#   data = dat_poc
# )
# summary(mod_poc)
# 
# # Porites 
# mod_por <- glmer.nb(
#   recruits_int ~ year_s + (1 | transect),
#   data = dat_por
#  # na.action = na.exclude()
# )
# summary(mod_por)
# 
# # Acr 
# mod_acr <- glmer.nb(
#   recruits_int ~ year_s + (1 | transect),
#   data = dat_acr
# )
# summary(mod_acr)
# # # Models with factor(year
# # 
# # # --- Pocillopora ---
# # mod_poc <- glmer.nb(
# #   recruits_int ~ factor(year) + (1 | transect),
# #   data = dat_poc
# # )
# # summary(mod_poc)
# # 
# # # --- Porites ---
# # mod_por <- glmer.nb(
# #   recruits_int ~ factor(year) + (1 | transect),
# #   data = dat_por
# # )
# # summary(mod_por)
# # 
# # # Predictions with factor(year
# # 
# # years_poc <- sort(unique(dat_poc$year))
# # 
# # pred_poc <- tibble(
# #   taxa = "Poc",
# #   year = years_poc,
# #   fit  = predict(mod_poc,
# #                  newdata = data.frame(year = years_poc),
# #                  re.form = NA,
# #                  type    = "response")
# # )
# # 
# # years_por <- sort(unique(dat_por$year))
# # 
# # pred_por <- tibble(
# #   taxa = "Por",
# #   year = years_por,
# #   fit  = predict(mod_por,
# #                  newdata = data.frame(year = years_por),
# #                  re.form = NA,
# #                  type    = "response")
# # )
# # 
# # pred_all <- bind_rows(pred_poc, pred_por)
# # ══════════════════════════════════════════════════════════════════════════════
# # STEP 7: Predictions to then back-transform year scaling for plotting
# 
# # Pocillopora 
# years_poc   <- sort(unique(dat_poc$year))
# years_poc_s <- (years_poc - min(dat_poc$year)) / (max(dat_poc$year) - min(dat_poc$year))
# 
# pred_poc <- tibble(
#   taxa = "Poc",
#   year = years_poc,
#   fit  = predict(mod_poc,
#                  newdata = data.frame(year_s = years_poc_s),
#                  re.form = NA,
#                  type    = "response")
# )
# 
# # --- Porites ---
# years_por   <- sort(unique(dat_por$year))
# years_por_s <- (years_por - min(dat_por$year)) / (max(dat_por$year) - min(dat_por$year))
# 
# pred_por <- tibble(
#   taxa = "Por",
#   year = years_por,
#   fit  = predict(mod_por,
#                  newdata = data.frame(year_s = years_por_s),
#                  re.form = NA,
#                  type    = "response")
# )
# 
# ##___ Do for Acr
# years_acr   <- sort(unique(dat_acr$year))
# years_acr_s <- (years_acr - min(dat_acr$year)) / (max(dat_acr$year) - min(dat_acr$year))
# 
# pred_acr <- tibble(
#   taxa = "Acr",
#   year = years_acr,
#   fit  = predict(mod_acr,
#                  newdata = data.frame(year_s = years_acr_s),
#                  re.form = NA,
#                  type    = "response")
# )
# 
# 
# 
# # combine — Acr and Mil excluded (too sparse, degenerate models)
# pred_all <- bind_rows(pred_poc, pred_por, pred_acr)
# 
# # ══════════════════════════════════════════════════════════════════════════════
# # Step 8: Plot
# 
# coral_colors <- c(
#   Acr = "blue",
#   Poc = "red",
#   Por = "darkgreen",
#   Mil = "darkorange"
# )
# 
# taxa_labels <- c(
#   Acr = "*Acropora*",
#   Mil = "*Millepora*",
#   Poc = "*Pocillopora*",
#   Por = "*Porites*"
# )
# 
# ggplot() +
#   # mean bars — all taxa
#   geom_col(
#     data  = recruit_bar,
#     aes(x = factor(year), y = mean_recruits, fill = taxa),
#     alpha = 0.5, width = 0.7
#   ) +
#   # raw transect points  all taxa
#   geom_point(
#     data  = recruit_transect,
#     aes(x = factor(year), y = recruits_std, color = taxa),
#     shape = 1, size = 1.5, alpha = 0.5
#   ) +
#   # GLMM fitted curves — Poc and Por only
#   geom_line(
#     data      = pred_all,
#     aes(x = factor(year), y = fit, color = taxa, group = taxa),
#     linewidth = 1.2
#   ) +
#   geom_point(
#     data  = pred_all,
#     aes(x = factor(year), y = fit, color = taxa),
#     size  = 2.5
#   ) +
#   facet_wrap(
#     ~ taxa,
#     scales   = "free_y",
#     ncol     = 2,
#     labeller = labeller(taxa = taxa_labels)
#   ) +
#   scale_fill_manual(values  = coral_colors) +
#   scale_color_manual(values = coral_colors) +
#   scale_x_discrete(breaks = c(2013, 2015, 2017, 2019, 2021, 2023)) +
#   labs(
#     title   = "Recruitment Over Time at LTER1 Back Reef",
#     x       = "Year",
#     y       = "Recruits per 5 m²",
#     caption = "GLMM curves shown for Poc and Por only (lme4, negative binomial).\nAcr and Mil excluded (insufficient recruitment signal)."
#   ) +
#   theme_light() +
#   theme(
#     axis.title       = element_text(size = 12),
#     axis.text        = element_text(size = 10),
#     axis.text.x      = element_text(angle = 45, vjust = 0.5),
#     strip.text       = ggtext::element_markdown(size = 12),
#     strip.background = element_rect(fill = "#325963"),
#     legend.position  = "none",
#     plot.caption     = element_text(size = 9, color = "grey50")
#   )
# 
# 
# ############_________________________
# # ══════════════════════════════════════════════════════════════════════════════
# # Step 2: Summarize at transect level 
# recruit_transect <- recruit_raw %>%
#   group_by(taxa, transect, year) %>%
#   summarize(n_recruits = sum(dyn_recruitment, na.rm = TRUE), .groups = "drop") %>%
#   mutate(
#     recruits_std = n_recruits * 5 / 25,
#     recruits_int = round(recruits_std)
#   )
# 
# # ══════════════════════════════════════════════════════════════════════════════
# # Step 3: mean per year for bar chart
# recruit_bar <- recruit_transect %>%
#   group_by(taxa, year) %>%
#   summarize(mean_recruits = mean(recruits_std, na.rm = TRUE), .groups = "drop")
# 
# # ══════════════════════════════════════════════════════════════════════════════
# # Step 4: Sparsity check to identify potential zero inflation
# sparsity_check <- recruit_transect %>%
#   group_by(taxa) %>%
#   summarize(n_nonzero = sum(recruits_int > 0), .groups = "drop")
# 
# print(sparsity_check)
# # Acr: 2 nonzero years = model degenerate, exclude from curves
# # Mil: 0 nonzero years = no model possible, exclude from curves
# # Poc: enough signal to fit glmm
# # Por: enough signal to fit glmm
# 
# # ══════════════════════════════════════════════════════════════════════════════
# # Step 5: Prep data with scaled year (0-1) for numerical stability
# dat_poc <- recruit_transect %>%
#   filter(taxa == "Poc") %>%
#   mutate(year_s = (year - min(year)) / (max(year) - min(year)))
# 
# dat_por <- recruit_transect %>%
#   filter(taxa == "Por") %>%
#   mutate(year_s = (year - min(year)) / (max(year) - min(year)))
# 
# dat_acr <- recruit_transect %>%
#   filter(taxa == "Acr") %>%
#   mutate(year_s = (year - min(year)) / (max(year) - min(year)))
# 
# # ══════════════════════════════════════════════════════════════════════════════
# # Step 6: Fit GLMMs (negative binomial, scaled year, random transect effect)
# 
# # Pocillopora 
# mod_poc <- glmer.nb(
#   recruits_int ~ year_s + (1 | transect),
#   data = dat_poc
# )
# summary(mod_poc)
# 
# # Porites 
# mod_por <- glmer.nb(
#   recruits_int ~ year_s + (1 | transect),
#   data = dat_por
# )
# summary(mod_por)
# 
# # Acr 
# mod_acr <- glmer.nb(
#   recruits_int ~ year_s + (1 | transect),
#   data = dat_acr
# )
# 
# summary(mod_acr)
# 
# # # Models with factor(year
# # 
# # # --- Pocillopora ---
# # mod_poc <- glmer.nb(
# #   recruits_int ~ factor(year) + (1 | transect),
# #   data = dat_poc
# # )
# # summary(mod_poc)
# # 
# # # --- Porites ---
# # mod_por <- glmer.nb(
# #   recruits_int ~ factor(year) + (1 | transect),
# #   data = dat_por
# # )
# # summary(mod_por)
# # 
# # # Predictions with factor(year
# # 
# # years_poc <- sort(unique(dat_poc$year))
# # 
# # pred_poc <- tibble(
# #   taxa = "Poc",
# #   year = years_poc,
# #   fit  = predict(mod_poc,
# #                  newdata = data.frame(year = years_poc),
# #                  re.form = NA,
# #                  type    = "response")
# # )
# # 
# # years_por <- sort(unique(dat_por$year))
# # 
# # pred_por <- tibble(
# #   taxa = "Por",
# #   year = years_por,
# #   fit  = predict(mod_por,
# #                  newdata = data.frame(year = years_por),
# #                  re.form = NA,
# #                  type    = "response")
# # )
# # 
# # pred_all <- bind_rows(pred_poc, pred_por)
# # ══════════════════════════════════════════════════════════════════════════════
# # STEP 7: Predictions to then back-transform year scaling for plotting
# 
# # Pocillopora 
# years_poc   <- sort(unique(dat_poc$year))
# years_poc_s <- (years_poc - min(dat_poc$year)) / (max(dat_poc$year) - min(dat_poc$year))
# 
# pred_poc <- tibble(
#   taxa = "Poc",
#   year = years_poc,
#   fit  = predict(mod_poc,
#                  newdata = data.frame(year_s = years_poc_s),
#                  re.form = NA,
#                  type    = "response")
# )
# 
# # --- Porites ---
# years_por   <- sort(unique(dat_por$year))
# years_por_s <- (years_por - min(dat_por$year)) / (max(dat_por$year) - min(dat_por$year))
# 
# pred_por <- tibble(
#   taxa = "Por",
#   year = years_por,
#   fit  = predict(mod_por,
#                  newdata = data.frame(year_s = years_por_s),
#                  re.form = NA,
#                  type    = "response")
# )
# 
# # combine — Acr and Mil excluded (too sparse, degenerate models)
# pred_all <- bind_rows(pred_poc, pred_por)
# 
# # ══════════════════════════════════════════════════════════════════════════════
# # STEP 8: Plot
# 
# coral_colors <- c(
#   Acr = "blue",
#   Poc = "red",
#   Por = "darkgreen",
#   Mil = "darkorange"
# )
# 
# taxa_labels <- c(
#   Acr = "*Acropora*",
#   Mil = "*Millepora*",
#   Poc = "*Pocillopora*",
#   Por = "*Porites*"
# )
# 
# ggplot() +
#   # mean bars — all taxa
#   geom_col(
#     data  = recruit_bar,
#     aes(x = factor(year), y = mean_recruits, fill = taxa),
#     alpha = 0.5, width = 0.7
#   ) +
#   # raw transect points — all taxa
#   geom_point(
#     data  = recruit_transect,
#     aes(x = factor(year), y = recruits_std, color = taxa),
#     shape = 1, size = 1.5, alpha = 0.5
#   ) +
#   # GLMM fitted curves — Poc and Por only
#   geom_line(
#     data      = pred_all,
#     aes(x = factor(year), y = fit, color = taxa, group = taxa),
#     linewidth = 1.2
#   ) +
#   geom_point(
#     data  = pred_all,
#     aes(x = factor(year), y = fit, color = taxa),
#     size  = 2.5
#   ) +
#   facet_wrap(
#     ~ taxa,
#     scales   = "free_y",
#     ncol     = 2,
#     labeller = labeller(taxa = taxa_labels)
#   ) +
#   scale_fill_manual(values  = coral_colors) +
#   scale_color_manual(values = coral_colors) +
#   scale_x_discrete(breaks = c(2013, 2015, 2017, 2019, 2021, 2023)) +
#   labs(
#     title   = "Recruitment Over Time at LTER1 Back Reef",
#     x       = "Year",
#     y       = "Recruits per 5 m²",
#     caption = "GLMM curves shown for Poc and Por only (lme4, negative binomial).\nAcr and Mil excluded (insufficient recruitment signal)."
#   ) +
#   theme_light() +
#   theme(
#     axis.title       = element_text(size = 12),
#     axis.text        = element_text(size = 10),
#     axis.text.x      = element_text(angle = 45, vjust = 0.5),
#     strip.text       = ggtext::element_markdown(size = 12),
#     strip.background = element_rect(fill = "#325963"),
#     legend.position  = "none",
#     plot.caption     = element_text(size = 9, color = "grey50")
#   )
# 
# ######_____________ GLMM for OR LTER1
# library(tidyverse)
# library(lme4)
# library(MASS)
# 
# # ══════════════════════════════════════════════════════════════════════════════
# # Step 1: Load and filter —
# clean_coral <- read_csv(
#   "~/MEDS/capstone/eds-coral-figs-storage/vedika-r-files/data/updated_coral_tidy_2013-2024.csv")
# 
# recruit_raw_or <- clean_coral %>%
#   filter(
#     #year >= 2013 & year <= 2018,
#   #  year >= 2019,
#     site    == "LTER1",
#     habitat == "OR",
#     !grepl("^P", transect)
#   )
# 
# ## pre or post bleaching status
# recruit_raw_or<- recruit_raw_or %>% mutate(bleaching_status = 
#                                                                case_when(
#                                                                  year >= 2013 & year <= 2018 ~ "pre", 
#                                                                  year >= 2019 & year <= 2024 ~ "post"
#                                                                ))
# recruit_raw_or <- recruit_raw_or %>% mutate(bleaching_status = factor(bleaching_status))
# 
# 
# # ══════════════════════════════════════════════════════════════════════════════
# # Step 2: Summarize at transect level
# # recruit_transect_or <- recruit_raw_or %>%
# #   group_by(taxa, transect, year) %>%
# #   summarize(n_recruits = sum(dyn_recruitment, na.rm = TRUE), .groups = "drop") %>%
# #   mutate(
# #     recruits_std = n_recruits * 5 / 25,
# #     recruits_int = round(recruits_std)
# #   ) 
# 
# recruit_transect_or <- recruit_raw_or %>%
#   group_by(taxa, transect, bleaching_status, year) %>%
#   summarize(n_recruits = sum(dyn_recruitment, na.rm = TRUE), .groups = "drop") %>%
#   mutate(
#     recruits_std = n_recruits * 5 / 25,
#     recruits_int = round(recruits_std)
#   )
# 
# # ══════════════════════════════════════════════════════════════════════════════
# # Step 3: mean per year for bar chart
# # recruit_bar_or <- recruit_transect_or %>%
# #   group_by(taxa, year) %>%
# #   summarize(mean_recruits = mean(recruits_std, na.rm = TRUE), .groups = "drop")
# recruit_bar_or <- recruit_transect_or %>%
#   group_by(taxa, bleaching_status, year) %>%
#   summarize(mean_recruits = mean(recruits_std, na.rm = TRUE), .groups = "drop")
# 
# # ══════════════════════════════════════════════════════════════════════════════
# # Step 4: Sparsity check
# sparsity_check_or <- recruit_transect_or %>%
#   group_by(taxa) %>%
#   summarize(n_nonzero = sum(recruits_int > 0), .groups = "drop")
# 
# print(sparsity_check_or)
# 
# # ══════════════════════════════════════════════════════════════════════════════
# # Step 5: Prep data with scaled year
# dat_poc_or <- recruit_transect_or %>%
#   filter(taxa == "Poc") %>%
#   mutate(year_s = (year - min(year)) / (max(year) - min(year)))
# 
# dat_por_or <- recruit_transect_or %>%
#   filter(taxa == "Por")%>%
#   mutate(year_s = (year - min(year)) / (max(year) - min(year)))
# 
# dat_acr_or <- recruit_transect_or %>%
#   filter(taxa == "Acr")%>%
# mutate(year_s = (year - min(year)) / (max(year) - min(year)))
# 
# # ══════════════════════════════════════════════════════════════════════════════
# # STEP 6: Fit models
# # --- Pocillopora: full GLMM ---
# mod_poc_or <- glmer.nb(
#   recruits_int ~ bleaching_status + year_s + (1 | transect),
#   data = dat_poc_or
# )
# summary(mod_poc_or)
# 
# # --- Porites: full GLMM ---
# mod_por_or <- glmer.nb(
#   recruits_int ~ bleaching_status + year_s + (1 | transect),
#   data = dat_por_or
# )
# summary(mod_por_or)
# 
# # Acropora: simple GLM, Poisson (counts 0/1, no overdispersion) ---
# mod_acr_or <- glmer.nb(
#   recruits_int ~ bleaching_status + year_s + (1 | transect),
#  # family = poisson,
#   data   = dat_acr_or
# )
# 
# summary(mod_acr_or)
# 
# # mod_acr_or <- glm(
# #   recruits_int ~ year_s,
# #   family = poisson,
# #   data   = dat_acr_or
# # )
# # summary(mod_acr_or)
# 
# # --- Millepora: check sparsity first ---
# recruit_transect_or %>% filter(taxa == "Mil") %>% count(year, recruits_int)
# # skip if all/mostly zeros
# 
# # ══════════════════════════════════════════════════════════════════════════════
# # Step 7: Predictions
# # --- Pocillopora ---
# years_poc_or   <- sort(unique(dat_poc_or$year))
# years_poc_or_s <- (years_poc_or - min(dat_poc_or$year)) / (max(dat_poc_or$year) - min(dat_poc_or$year))
#  
# pred_poc_or <- tibble(
#   taxa = "Poc",
#   year = years_poc_or,
#   bleaching_status = dat_poc_or$bleaching_status, 
#   fit  = predict(mod_poc_or,
#                  newdata = data.frame(year_s = years_poc_or_s),
#                  re.form = NA,
#                  type    = "response")
# )
# 
# # --- Porites ---
# years_por_or   <- sort(unique(dat_por_or$year))
# years_por_or_s <- (years_por_or - min(dat_por_or$year)) / (max(dat_por_or$year) - min(dat_por_or$year))
# 
# pred_por_or <- tibble(
#   taxa = "Por",
#   year = years_por_or,
#   fit  = predict(mod_por_or,
#                  newdata = data.frame(year_s = years_por_or_s),
#                  re.form = NA,
#                  type    = "response")
# )
# 
# # --- Acropora ---
# years_acr_or   <- sort(unique(dat_acr_or$year))
# years_acr_or_s <- (years_acr_or - min(dat_acr_or$year)) / (max(dat_acr_or$year) - min(dat_acr_or$year))
# 
# # pred_acr_or <- tibble(
# #   taxa = "Acr",
# #   year = years_acr_or,
# #   fit  = predict(mod_acr_or,
# #                  newdata = data.frame(year_s = years_acr_or_s),
# #                  type    = "response")   # no re.form for glm
# # )
# 
# # if poisson
# pred_acr_or <- tibble(
#   taxa = "Acr",
#   year = years_acr_or,
#   fit  = predict(mod_acr_or,
#                  newdata = data.frame(year_s = years_acr_or_s),
#                  re.form = NA,
#                  type    = "response")
# )
# 
# 
# # combine all predictions
# pred_all_or <- bind_rows(pred_poc_or, pred_por_or, pred_acr_or)
# 
# # ══════════════════════════════════════════════════════════════════════════════
# # Step 8: Plot
# 
# coral_colors <- c(
#   Acr = "blue",
#   Poc = "red",
#   Por = "darkgreen",
#   Mil = "darkorange"
# )
# 
# taxa_labels <- c(
#   Acr = "*Acropora*",
#   Mil = "*Millepora*",
#   Poc = "*Pocillopora*",
#   Por = "*Porites*"
# )
# 
# ggplot() +
#   # mean bars for all taxa
#   geom_col(
#     data  = recruit_bar_or,
#     aes(x = factor(year), y = mean_recruits, fill = taxa),
#     alpha = 0.8, width = 0.7
#   ) +
#   # raw transect points for all taxa
#   geom_point(
#     data  = recruit_transect_or,
#     aes(x = factor(year), y = recruits_std, color = taxa),
#     shape = 1, size = 1.5, alpha = 0.5
#   ) +
#   # GLMM/GLM fitted curves — Poc, Por, Acr
#   geom_line(
#     data      = pred_all_or,
#     aes(x = factor(year), y = fit, color = taxa, group = taxa),
#     linewidth = 1.2
#   ) +
#   geom_point(
#     data  = pred_all_or,
#     aes(x = factor(year), y = fit, color = taxa),
#     size  = 2.5
#   ) +
#   facet_wrap(
#     ~ taxa,
#     scales   = "free_y",
#     ncol     = 2,
#     labeller = labeller(taxa = taxa_labels)
#   ) +
#   scale_fill_manual(values  = coral_colors) +
#   scale_color_manual(values = coral_colors) +
#   scale_x_discrete(breaks = c(2013, 2015, 2017, 2019, 2021, 2023)) +
#   labs(
#     title   = "Recruitment Over Time at LTER1 Fore Reef",
#     x       = "Year",
#     y       = "Recruits per 5 m²",
#     caption = "GLMM curves: Poc and Por (lme4, negative binomial), Acr (GLM, Poisson).\nMil excluded (insufficient recruitment signal)."
#   ) +
#   theme_light() +
#   theme(
#     axis.title       = element_text(size = 12),
#     axis.text        = element_text(size = 10),
#     axis.text.x      = element_text(angle = 45, vjust = 0.5),
#     strip.text       = ggtext::element_markdown(size = 12),
#     strip.background = element_rect(fill = "#325963"),
#     legend.position  = "none",
#     plot.caption     = element_text(size = 9, color = "grey50")
#   )

############_____________ function to plot glmm with averaged 
# plot_recruitment_1 <- function(coral_data,
#                                site_name = NULL,
#                                habitat_name = NULL,
#                                facet_level = c("transect", "site", "habitat"),
#                                add_glmm = TRUE) {
#   
#   facet_level <- match.arg(facet_level)
#   
#   coral_colors <- c(
#     Acr = "blue",
#     Poc = "red",
#     Por = "darkgreen",
#     Mil = "darkorange"
#   )
#   
#   habitat_labels <- c(BR = "Back Reef", OR = "Fore Reef")
#   
#   taxa_labels <- c(
#     Acr = "*Acropora*",
#     Mil = "*Millepora*",
#     Poc = "*Pocillopora*",
#     Por = "*Porites*"
#   )
#   
#   transect_labels <- c(T01 = "T1", T02 = "T2", T03 = "T3", T04 = "T4")
#   
#   facet_labeller <- labeller(
#     habitat  = habitat_labels,
#     taxa     = taxa_labels,
#     transect = transect_labels,
#     .default = label_value
#   )
#   
#   # ── Data prep ───────────────────────────────────────────────────────────────
#   plot_data <- coral_data %>%
#     {if (!is.null(site_name))    dplyr::filter(., site %in% site_name) else .} %>%
#     {if (!is.null(habitat_name)) dplyr::filter(., habitat == habitat_name) else .} %>%
#     dplyr::filter(!grepl("^P", transect)) %>%
#     dplyr::group_by(transect, site, habitat, taxa, year) %>%
#     dplyr::summarize(
#       n_recruits = sum(dyn_recruitment, na.rm = TRUE),
#       .groups = "drop"
#     ) %>%
#     dplyr::mutate(recruits_std = n_recruits * 5 / 25)  # per 5 m²
#   
#   # ── Mean per year for bars (replacing sum) ───────────────────────────────── 
#   bar_data <- plot_data %>%
#     dplyr::group_by(site, habitat, taxa, year) %>%
#     dplyr::summarize(
#       mean_recruits = mean(recruits_std, na.rm = TRUE),
#       .groups = "drop"
#     )
#   
#   # ── GLMM fits ────────────────────────────────────────────────────────────────
#   pred_data <- NULL
#   
#   if (add_glmm) {
#     # determine which taxa have enough signal to model
#     sparse_check <- plot_data %>%
#       dplyr::group_by(site, habitat, taxa) %>%
#       dplyr::summarize(
#         nonzero = sum(recruits_std > 0),
#         .groups = "drop"
#       ) %>%
#       dplyr::filter(nonzero >= 3)  # need at least 3 nonzero years to fit
#     
#     pred_data <- sparse_check %>%
#       purrr::pmap_dfr(function(site, habitat, taxa, nonzero) {
#         dat <- plot_data %>%
#           dplyr::filter(
#             site    == !!site,
#             habitat == !!habitat,
#             taxa    == !!taxa
#           ) %>%
#           dplyr::mutate(recruits_int = round(recruits_std))
#         
#         mod <- tryCatch(
#           glmmTMB::glmmTMB(
#             recruits_int ~ factor(year) + (1 | transect),
#             dispformula = ~1,
#             family      = glmmTMB::nbinom2,
#             data        = dat,
#             na.action   = na.exclude
#           ),
#           error = function(e) NULL
#         )
#         
#         if (is.null(mod)) return(NULL)
#         
#         years <- sort(unique(dat$year))
#         tibble::tibble(
#           site    = site,
#           habitat = habitat,
#           taxa    = taxa,
#           year    = years,
#           fit     = predict(mod,
#                             newdata = data.frame(year = years),
#                             re.form = NA,
#                             type    = "response")
#         )
#       })
#   }
#   
#   facet_formula <- switch(
#     facet_level,
#     "transect" = transect ~ taxa,
#     "site"     = site ~ taxa,
#     "habitat"  = habitat ~ taxa
#   )
#   
#   theme_fn <- if (facet_level == "site") ggthemes::theme_igray() else theme_light()
#   
#   # ── Plot ─────────────────────────────────────────────────────────────────────
#   p <- ggplot() +
#     # mean bars
#     geom_col(
#       data = bar_data,
#       aes(x = factor(year), y = mean_recruits, fill = taxa, group = taxa),
#       alpha = 0.8
#     ) +
#     # raw transect points
#     geom_point(
#       data  = plot_data,
#       aes(x = factor(year), y = recruits_std, color = taxa),
#       shape = 1, size = 1.5, alpha = 0.5
#     ) +
#     facet_grid(
#       facet_formula,
#       scales   = "free",
#       labeller = facet_labeller
#     ) +
#     scale_fill_manual(values  = coral_colors) +
#     scale_color_manual(values = coral_colors) +
#     scale_x_discrete(breaks = c(2014, 2016, 2018, 2020, 2022, 2024)) +
#     labs(
#       x    = "Year",
#       y    = "Recruits per 5 m²",
#       fill = "Taxa"
#     ) +
#     theme_fn +
#     theme(
#       axis.title       = element_text(size = 12),
#       axis.text        = element_text(size = 10),
#       legend.text      = element_text(size = 10),
#       legend.title     = element_text(size = 10),
#       strip.text       = ggtext::element_markdown(size = 14),
#       strip.background = element_rect(fill = "#325963"),
#       panel.background = element_rect(fill = "transparent", color = NA),
#       plot.background  = element_rect(fill = "transparent", color = NA),
#       axis.text.x      = element_text(angle = 45, vjust = 0.5),
#       legend.position  = "none"
#     )
#   
#   # add GLMM curve if available
#   if (!is.null(pred_data) && nrow(pred_data) > 0) {
#     p <- p +
#       geom_line(
#         data      = pred_data,
#         aes(x = factor(year), y = fit, color = taxa, group = taxa),
#         linewidth = 1.2
#       ) +
#       geom_point(
#         data  = pred_data,
#         aes(x = factor(year), y = fit, color = taxa),
#         size  = 2.5
#       )
#   }
#   
#   return(p)
# }
# 
# 
# plot_recruitment_1(clean_coral, site_name = "LTER1", facet_level = "habitat")
# 
# ##########_________________ loess fit
# plot_recruitment_loess <- function(coral_data,
#                                site_name = NULL,
#                                habitat_name = NULL,
#                                facet_level = c("transect", "site", "habitat"),
#                                add_smoother = TRUE) {
#   
#   facet_level <- match.arg(facet_level)
#   
#   coral_colors <- c(Acr = "blue", Poc = "red", Por = "darkgreen", Mil = "darkorange")
#   habitat_labels <- c(BR = "Back Reef", OR = "Fore Reef")
#   taxa_labels <- c(
#     Acr = "*Acropora*", Mil = "*Millepora*",
#     Poc = "*Pocillopora*", Por = "*Porites*"
#   )
#   transect_labels <- c(T01 = "T1", T02 = "T2", T03 = "T3", T04 = "T4")
#   
#   facet_labeller <- labeller(
#     habitat = habitat_labels, taxa = taxa_labels,
#     transect = transect_labels, .default = label_value
#   )
#   
#   # ── Data prep ────────────────────────────────────────────────────────────────
#   plot_data <- coral_data %>%
#     {if (!is.null(site_name))    dplyr::filter(., site %in% site_name) else .} %>%
#     {if (!is.null(habitat_name)) dplyr::filter(., habitat == habitat_name) else .} %>%
#     dplyr::filter(!grepl("^P", transect)) %>%
#     dplyr::group_by(transect, site, habitat, taxa, year) %>%
#     dplyr::summarize(n_recruits = sum(dyn_recruitment, na.rm = TRUE), .groups = "drop") %>%
#     dplyr::mutate(recruits_std = n_recruits * 5 / 25)
#   
#   # ── Mean per year for bars ───────────────────────────────────────────────────
#   bar_data <- plot_data %>%
#     dplyr::group_by(site, habitat, taxa, year) %>%
#     dplyr::summarize(mean_recruits = mean(recruits_std, na.rm = TRUE), .groups = "drop")
#   
#   facet_formula <- switch(
#     facet_level,
#     "transect" = transect ~ taxa,
#     "site"     = site ~ taxa,
#     "habitat"  = habitat ~ taxa
#   )
#   
#   theme_fn <- if (facet_level == "site") ggthemes::theme_igray() else theme_light()
#   
#   # ── Plot ─────────────────────────────────────────────────────────────────────
#   p <- ggplot() +
#     geom_col(
#       data  = bar_data,
#       aes(x = factor(year), y = mean_recruits, fill = taxa, group = taxa),
#       alpha = 0.8
#     ) +
#     geom_point(
#       data  = plot_data,
#       aes(x = factor(year), y = recruits_std, color = taxa),
#       shape = 1, size = 1.5, alpha = 0.5
#     ) +
#     facet_grid(facet_formula, scales = "free", labeller = facet_labeller) +
#     scale_fill_manual(values  = coral_colors) +
#     scale_color_manual(values = coral_colors) +
#     scale_x_discrete(breaks = c(2014, 2016, 2018, 2020, 2022, 2024)) +
#     labs(x = "Year", y = "Recruits per 5 m²", fill = "Taxa") +
#     theme_fn +
#     theme(
#       axis.title       = element_text(size = 12),
#       axis.text        = element_text(size = 10),
#       strip.text       = ggtext::element_markdown(size = 14),
#       strip.background = element_rect(fill = "#325963"),
#       panel.background = element_rect(fill = "transparent", color = NA),
#       plot.background  = element_rect(fill = "transparent", color = NA),
#       axis.text.x      = element_text(angle = 45, vjust = 0.5),
#       legend.position  = "none"
#     )
#   
#   # ── Loess smoother ───────────────────────────────────────────────────────────
#   if (add_smoother) {
#     p <- p +
#       geom_smooth(
#         data    = bar_data,                          # smooth over yearly means
#         aes(x   = factor(year), y = mean_recruits,
#             color = taxa, group = taxa),
#         method  = "loess",                           # matches reference figure
#         se      = FALSE,                             # no confidence ribbon
#         linewidth = 1.2,
#         span    = 0.75                               # adjust 0.5-1 for more/less smoothing
#       )
#   }
#   
#   return(p)
# }


# call it
#plot_recruitment_loess(clean_coral, site_name = "LTER1", facet_level = "habitat")
# # ── 1. Load & prep ALL taxa ───────────────────────────────────────────────────
# recruit_data <- clean_coral %>%
#   filter(
#     site     == "LTER1",
#     habitat  == "BR",
#     transect %in% c("T01", "T02", "T03", "T04")
#   ) %>%
#   group_by(taxa, transect, year) %>%
#   summarize(n_recruits = sum(dyn_recruitment, na.rm = TRUE), .groups = "drop") %>%
#   mutate(recruits_std = round(n_recruits / 5))
# 
# # ── 2. Fit GLMM only for taxa with enough signal ──────────────────────────────
# model_taxa <- c("Poc", "Por")  # skip Acr and Mil — too sparse
# 
# models <- model_taxa %>%
#   set_names() %>%
#   map(~ {
#     dat <- recruit_data %>% filter(taxa == .x)
#     glmmTMB(
#       recruits_std ~ factor(year) + (1 | transect),
#       dispformula = ~1,
#       family      = nbinom2,
#       data        = dat,
#       na.action   = na.exclude
#     )
#   })
# 
# # ── 3. Predictions only for modeled taxa ──────────────────────────────────────
# pred_df <- model_taxa %>%
#   set_names() %>%
#   map_dfr(~ {
#     dat   <- recruit_data %>% filter(taxa == .x)
#     years <- sort(unique(dat$year))
#     tibble(
#       taxa = .x,
#       year = years,
#       fit  = predict(
#         models[[.x]],
#         newdata = data.frame(year = years),
#         re.form = NA,
#         type    = "response"
#       )
#     )
#   })
# 
# # ── 4. Bar chart summary for ALL taxa ─────────────────────────────────────────
# bar_df <- recruit_data %>%
#   group_by(taxa, year) %>%
#   summarize(mean_recruits = mean(recruits_std, na.rm = TRUE), .groups = "drop")
# 
# # ── 5. Color palette ──────────────────────────────────────────────────────────
# taxa_colors <- c(
#   "Acr" = "#4472C4",
#   "Mil" = "#FFC000",
#   "Poc" = "#FF0000",
#   "Por" = "#2E7D32"
# )
# 
# # ── 6. Plot — bars for all 4, curve only for Poc & Por ────────────────────────
# ggplot() +
#   # raw transect points — all taxa
#   geom_point(
#     data  = recruit_data,
#     aes(x = year, y = recruits_std, color = taxa),
#     alpha = 0.4, size = 2, shape = 1
#   ) +
#   # mean bars — all taxa
#   geom_col(
#     data  = bar_df,
#     aes(x = year, y = mean_recruits, fill = taxa),
#     alpha = 0.5, width = 0.7
#   ) +
#   # GLMM fitted curve — Poc & Por only
#   geom_line(
#     data      = pred_df,
#     aes(x = year, y = fit, color = taxa),
#     linewidth = 1.2
#   ) +
#   geom_point(
#     data  = pred_df,
#     aes(x = year, y = fit, color = taxa),
#     size  = 2.5
#   ) +
#   scale_color_manual(values = taxa_colors) +
#   scale_fill_manual(values  = taxa_colors) +
#   scale_x_continuous(breaks = unique(recruit_data$year)) +
#   facet_wrap(~ taxa, scales = "free_y", ncol = 2) +
#   labs(
#     title   = "Recruitment Over Time — LTER1 Back Reef",
#     x       = "Year",
#     y       = "Recruits per m²",
#     caption = "Fitted curves shown for Poc and Por only (Acr and Mil too sparse to model)"
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.x      = element_text(angle = 45, hjust = 1),
#     legend.position  = "none",
#     strip.text       = element_text(face = "italic", size = 11),
#     strip.background = element_rect(fill = "grey90")
#   )

#######___________________________
# library(tidyverse)
# library(dplyr)
# library(gganimate)
# 
# # Load in data
# coral <- read_csv(
#   "~/MEDS/capstone/eds-coral-figs-storage/vedika-r-files/data/updated_coral_tidy_2013-2024.csv")
# 
# library(MASS)
# library(glmmTMB)
# 
# # # Load in data 
# # clean_coral <- read_csv(here::here("vedika-r-files", "data", "volume_corals.csv"))
# clean_coral <- coral %>% filter(transect %in% c("T01", "T02", "T03", "T04"))
# 
# 
# # Standardize recruitment
# clean_coral_standard <- clean_coral %>% 
#   filter(taxa == "Por", site == "LTER1", habitat == "BR") %>% 
#   group_by(transect, year) %>% 
#   summarize(n_recruits = sum(dyn_recruitment, na.rm = TRUE),
#             .groups = "drop") %>%
#   mutate(recruits_std = n_recruits / 5) # Standardize per 5 m^2
# 
# 
# clean_coral_standard <- clean_coral_standard %>%
#   mutate(log_area = log(5))
# 
# recruit_time <- glmmTMB(
#   n_recruits ~ factor(year) + (1 | transect) + offset(log_area),
#   dispformula = ~1,
#   family = nbinom2,
#   data = clean_coral_standard,
#   na.action = na.exclude
# )
# 
# 
# 
# ###########_____________________
# # recruit_time <-glmmTMB(
# #   recruits_std ~ factor(year), dispformula=~1 | transect, family= negative.binomial, data = clean_coral_standard, na.action=na.exclude)
# # 
# # #recruit_summary <- summary(recruit_time)
# # 
# # clean_coral_standard$fit.recruit<- predict(recruit_time, type = "response")
# # 
# # view(clean_coral_standard)
# # clean_coral_standard %>% group_by(year) %>% mutate(recruits_year = sum(fit.recruit, na.rm = T)) %>%
# #   ggplot(aes(factor(year), recruits_std)) + 
# #   geom_col() + 
# #   geom_point(aes(x = factor(year), y = recruits_year), color = "blue", size = 1)
# # 
# 
# 
# 
# 
# # coral_death <- coral %>% 
# #   group_by(coral_number) %>%
# #   arrange(year) %>% 
# #   mutate(survival = case_when(
# #     dyn_death == 1 ~ 0,
# #     is.na(length) ~ NA,
# #     TRUE ~ 1
# #   ), volume = case_when(
# #     is.na(size_t0) ~ size_t_minus1,
# #     TRUE~ size_t0
# #   )
# #   ) %>% 
# #   ungroup()
# # 
# # library(MASS)
# # library(glmmTMB)
# # acr_site_1 <- coral_death %>% filter(site == "LTER1", taxa == "Acr", habitat == "OR")
# # surv_acr_site_1 <- glmmTMB(survival ~ volume, dispformula=~1|coral_number, family = binomial, data = acr_site_1, na.action = na.exclude)
# # sum <- summary(surv_acr_site_1)
# # 
# # #acr_site_1$fit.Surv <- 
# # exp(sum$coefficients[,1])/(1+exp(sum$fitted[,1])) #back transform the model fit and add to data table
# # 
# # 
# # 
# # view(acr_site_1)
# 
# 
# 
# 
