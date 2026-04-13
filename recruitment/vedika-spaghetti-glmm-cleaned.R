library(tidyverse)
library(lme4)
library(glmmTMB)
library(ggtext)
library(ggeffects)

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
    # Standardize duplicate/variant transect names.
    # TRUE ~ transect preserves non-duplicate names (without it, unmatched
    # transects become NA and drop silently from the model).
    transect = case_when(
      transect == "T02_OLD"     ~ "T02",
      transect == "T03_PRE2020" ~ "T03",
      transect == "T04_OLD"     ~ "T04",
      transect == "T02_2021"    ~ "T02",
      transect == "TO1"         ~ "T01",
      TRUE                      ~ transect
    )
  )

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3: Summarize at transect level
#
# dyn_recruitment is 0/1 per coral row. Summing gives total recruits present
# in a given transect x site x habitat x taxa x year combination.
# Dividing by 5 standardizes to recruits per 5 m².
#
# group_by includes site and habitat so they survive summarize() and remain
# available as fixed effects in the model.
#
# transect_id = paste(site, transect) creates a globally unique transect
# identifier (e.g., "LTER1_T01") so the random effect correctly treats
# T01 at LTER1 as distinct from T01 at LTER2 — the nested Transect ⊂ Site
# structure described in meeting notes.
# ══════════════════════════════════════════════════════════════════════════════

recruit_transect <- recruit_raw %>%
  group_by(taxa, site, habitat, transect, year) %>%
  summarize(
    n_recruits = sum(dyn_recruitment, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(
    recruits_std = n_recruits / 5,
    recruits_int = round(recruits_std),
    pre_post     = case_when(
      year >= 2013 & year <= 2018 ~ 0L,
      year >= 2019 & year <= 2024 ~ 1L
    ),
    period      = ifelse(pre_post == 0, "Pre (2013-2018)", "Post (2019-2024)"),
    year_s      = (year - min(year)) / (max(year) - min(year)),
    transect_id = paste(site, transect, sep = "_")
  )

# Confirm unique transect IDs look correct
sort(unique(recruit_transect$transect_id))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4: Sparsity check — identify taxa too sparse to model
# ══════════════════════════════════════════════════════════════════════════════

sparsity_check <- recruit_transect %>%
  group_by(taxa, period) %>%
  summarize(
    n_rows    = n(),
    n_nonzero = sum(recruits_int > 0),
    .groups   = "drop"
  )
print(sparsity_check)
# Acr and Mil are too sparse in both periods — complete separation confirmed
# by diagnose() (non-finite SEs, flat Hessian directions). Excluded from GLMMs.
# Poc and Por have sufficient signal in both periods.

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5: Split into pre/post periods, filter to modelable taxa (Poc, Por)
#
# Acr and Mil excluded: complete separation in multiple site x habitat combos
# (all-zero counts) makes their fixed effect coefficients inestimable and
# destabilizes the Hessian. Confirmed via diagnose() in both periods.
# ══════════════════════════════════════════════════════════════════════════════

dat_pre  <- recruit_transect %>%
  filter(pre_post == 0, taxa %in% c("Poc", "Por"))

dat_post <- recruit_transect %>%
  filter(pre_post == 1, taxa %in% c("Poc", "Por"))

range(dat_pre$year)   # should be 2013-2018
range(dat_post$year)  # should be 2019-2024

# ══════════════════════════════════════════════════════════════════════════════
# STEP 6: Complete separation check before fitting
#
# Any site x habitat x taxa combo with all-zero counts pushes a fixed effect
# coefficient toward -Inf, destabilizing the Hessian.
# ══════════════════════════════════════════════════════════════════════════════

sep_pre <- dat_pre %>%
  group_by(site, habitat, taxa) %>%
  summarize(total = sum(recruits_int), .groups = "drop") %>%
  filter(total == 0)

sep_post <- dat_post %>%
  group_by(site, habitat, taxa) %>%
  summarize(total = sum(recruits_int), .groups = "drop") %>%
  filter(total == 0)

cat("Pre-bleaching zero combos:\n");  print(sep_pre)
cat("Post-bleaching zero combos:\n"); print(sep_post)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 7: Full GLMMs — one per period, answers both research questions
#
# FIXED EFFECTS (both models):
#   habitat      — 2 levels (BR vs OR); fixed for direct comparison.
#                  Random effects unreliable with only 2 levels.
#   site         — fixed; Q2 asks about named-site differences (LTER1 vs
#                  LTER2 etc.) requiring explicit contrasts.
#   taxa         — fixed (Poc vs Por).
#   factor(year) — fixed within each period; fewer assumptions than random
#                  year, avoids confounding with intercept (per Max).
#
# PRE-BLEACHING: nbinom2 + random transect effect
#   - Data is overdispersed pre-bleaching (confirmed by overdispersion check).
#   - Transect-level clustering is estimable (sufficient between-transect
#     variance in recruitment counts).
#   - (1 | transect_id) accounts for repeated measures on the same transect
#     across years within the pre period.
#   - Large Z-statistics on some coefficients = well-estimated parameters,
#     not a model failure. Use profile CIs for inference (see Step 9).
#
# POST-BLEACHING: Poisson, no random effect
#   Two findings from diagnose() drove these choices:
#   (1) Overdispersion ratio = 1.009 — data is NOT overdispersed post-
#       bleaching. nbinom2 dispersion blew up to ~21, confirming Poisson
#       is the correct family.
#   (2) Random effect variance collapsed (theta = -17.5) — transect-level
#       clustering disappeared post-bleaching. Recruitment became so low
#       and spatially homogeneous that there is no between-transect variance
#       left to estimate.
#   Ecological interpretation: post-bleaching recruitment is rare enough to
#   be Poisson-distributed, and spatial patchiness has been lost. This
#   asymmetry between models is a finding, not a modeling inconsistency.
# ══════════════════════════════════════════════════════════════════════════════

mod_pre <- glmmTMB(
  recruits_int ~ habitat + site + taxa + factor(year) + (1 | transect_id),
  family    = nbinom2,
  data      = dat_pre,
  na.action = na.exclude
)

mod_post <- glmmTMB(
  recruits_int ~ habitat + site + taxa + factor(year),
  family    = poisson,
  data      = dat_post,
  na.action = na.exclude
)

summary(mod_pre)
summary(mod_post)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 8: Model diagnostics
# ══════════════════════════════════════════════════════════════════════════════

# ── Overdispersion ───────────────────────────────────────────────────────────
check_overdisp <- function(model, label) {
  pr    <- residuals(model, type = "pearson")
  ratio <- sum(pr^2) / df.residual(model)
  cat(label, "overdispersion ratio:", round(ratio, 3), "\n")
  # >> 1: overdispersed (nbinom2 appropriate)
  # ~= 1: not overdispersed (Poisson appropriate)
  # << 1: underdispersed
}

check_overdisp(mod_pre,  "Pre-bleaching  (nbinom2)")
check_overdisp(mod_post, "Post-bleaching (Poisson)")

# ── Convergence ──────────────────────────────────────────────────────────────
diagnose(mod_pre)
diagnose(mod_post)

# ── Temporal autocorrelation — residuals vs year ─────────────────────────────
dat_pre$resid  <- NA
dat_post$resid <- NA
dat_pre$resid[!is.na(dat_pre$recruits_int)]   <- residuals(mod_pre,  type = "pearson")
dat_post$resid[!is.na(dat_post$recruits_int)] <- residuals(mod_post, type = "pearson")

resid_all <- bind_rows(dat_pre, dat_post)

ggplot(resid_all, aes(x = year, y = resid)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(taxa ~ period, scales = "free_x") +
  labs(
    title = "Pearson Residuals vs. Year — temporal autocorrelation check",
    x     = "Year",
    y     = "Pearson Residual"
  ) +
  theme_light()
# Systematic trend in residuals (linear drift, all negative near end of period)
# would indicate autocorrelation. If present, add ar1() structure in glmmTMB.

# ══════════════════════════════════════════════════════════════════════════════
# STEP 9: Profile confidence intervals and marginal effects
#
# Profile CIs used instead of Wald because large Z-statistics on some
# pre-bleaching coefficients indicate the Wald approximation may be
# unreliable. Profile CIs are valid regardless.
# ══════════════════════════════════════════════════════════════════════════════

# Profile CIs — comment out during iteration (slow), run for final results
# confint(mod_pre,  method = "profile")
# confint(mod_post, method = "profile")

# Q1: Habitat effect (BR vs OR), averaged over sites, taxa, years
ggeffect(mod_pre,  terms = "habitat")
ggeffect(mod_post, terms = "habitat")

# Q2: Site effect, averaged over habitats, taxa, years
ggeffect(mod_pre,  terms = "site")
ggeffect(mod_post, terms = "site")

# Year trend within each period
ggeffect(mod_pre,  terms = "year")
ggeffect(mod_post, terms = "year")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 10: Spaghetti GLMMs — fit one model per site × habitat × taxa × period
#
# Each panel in the final plot = one site × habitat × taxa combination,
# split into pre and post periods. We fit a separate glmer.nb per group
# so the GLMM curve reflects the recruitment trajectory at that specific
# location and taxon.
#
# Model: recruits_int ~ year_s + (1 | transect_id)
#   - year_s scaled 0-1 WITHIN each period independently (0 = first year
#     of that period, 1 = last). This prevents the pre-period range from
#     being compressed when fitting period-specific models.
#   - (1 | transect_id) accounts for repeated measures on the same transect
#     across years within the period.
#
# Groups with insufficient data (< 2 non-zero observations or only 1
# transect) are skipped — a fitted curve requires some variance to estimate.
# ══════════════════════════════════════════════════════════════════════════════

# All site x habitat x taxa x period combinations to model
plot_data <- recruit_transect %>%
  filter(taxa %in% c("Poc", "Por"))

combinations <- plot_data %>%
  distinct(site, habitat, taxa, period, pre_post)

# Fit a glmer.nb for each combination, skip if insufficient data
fit_panel_model <- function(site_name, habitat_name, taxa_name,
                            period_label, pre_post_val, data) {
  
  dat <- data %>%
    filter(
      site    == site_name,
      habitat == habitat_name,
      taxa    == taxa_name,
      pre_post == pre_post_val
    ) %>%
    mutate(year_s = (year - min(year)) / (max(year) - min(year)))
  
  # Skip if fewer than 2 non-zero observations or only 1 transect
  # (model cannot estimate a trend or random effect)
  if (sum(dat$recruits_int > 0) < 2 || n_distinct(dat$transect_id) < 2) {
    return(NULL)
  }
  
  tryCatch(
    glmer.nb(recruits_int ~ year_s + (1 | transect_id), data = dat),
    error   = function(e) NULL,
    warning = function(w) {
      # Still return model on warnings (e.g. large Z-stats) but print notice
      cat("Warning for", site_name, habitat_name, taxa_name, period_label, ":\n",
          conditionMessage(w), "\n")
      tryCatch(
        glmer.nb(recruits_int ~ year_s + (1 | transect_id), data = dat),
        error = function(e) NULL
      )
    }
  )
}

# Run all combinations
panel_models <- combinations %>%
  mutate(
    model = pmap(
      list(site, habitat, taxa, period, pre_post),
      fit_panel_model,
      data = plot_data
    )
  )

cat("Models fitted:", sum(!map_lgl(panel_models$model, is.null)), "/",
    nrow(panel_models), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 11: Generate predictions from each panel model
#
# re.form = NA marginalizes over random transect effect → population-level
# curve representing the average trajectory at that site x habitat x taxa.
# ══════════════════════════════════════════════════════════════════════════════

make_panel_preds <- function(model, site_name, habitat_name,
                             taxa_name, period_label, pre_post_val, data) {
  if (is.null(model)) return(NULL)
  
  dat     <- data %>%
    filter(site == site_name, habitat == habitat_name,
           taxa == taxa_name, pre_post == pre_post_val) %>%
    mutate(year_s = (year - min(year)) / (max(year) - min(year)))
  
  years   <- sort(unique(dat$year))
  years_s <- (years - min(dat$year)) / (max(dat$year) - min(dat$year))
  
  tibble(
    site    = site_name,
    habitat = habitat_name,
    taxa    = taxa_name,
    period  = period_label,
    year    = years,
    fit     = predict(
      model,
      newdata = data.frame(year_s = years_s),
      re.form = NA,
      type    = "response"
    )
  )
}

pred_panels <- panel_models %>%
  pmap_dfr(function(site, habitat, taxa, period, pre_post, model) {
    make_panel_preds(model, site, habitat, taxa, period, pre_post, plot_data)
  })

# ══════════════════════════════════════════════════════════════════════════════
# STEP 12: Summarize raw data for bars and points
# ══════════════════════════════════════════════════════════════════════════════

# Mean recruits per site x habitat x taxa x year x period (for bars)
recruit_bar <- plot_data %>%
  group_by(site, habitat, taxa, year, period) %>%
  summarize(mean_recruits = mean(recruits_std, na.rm = TRUE), .groups = "drop")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 13: Plot — one panel per site x habitat x taxa, split by period
#
# Layout: rows = site x habitat x taxa, columns = period (pre | post)
# Each panel shows:
#   - Grey bars:   mean recruits per year across transects
#   - Open points: raw transect-level recruits_std (individual observations)
#   - Colored line + points: GLMM fitted population-level trend curve
#
# Panels without a curve had insufficient data to fit a model (< 2 non-zero
# obs or only 1 transect).
# ══════════════════════════════════════════════════════════════════════════════

coral_colors <- c(Poc = "red", Por = "darkgreen")

taxa_labels <- c(
  Poc = "*Pocillopora*",
  Por = "*Porites*"
)

# Create a combined facet label: "LTER1 — BR — Pocillopora"
plot_data <- plot_data %>%
  mutate(panel_label = paste(site, habitat, taxa, sep = " \u2014 "))

recruit_bar <- recruit_bar %>%
  mutate(panel_label = paste(site, habitat, taxa, sep = " \u2014 "))

pred_panels <- pred_panels %>%
  mutate(panel_label = paste(site, habitat, taxa, sep = " \u2014 "))

ggplot() +
  # Mean bars
  geom_col(
    data  = recruit_bar,
    aes(x = factor(year), y = mean_recruits),
    fill  = "grey70", alpha = 0.7, width = 0.7
  ) +
  # Raw transect points
  geom_point(
    data  = plot_data,
    aes(x = factor(year), y = recruits_std, color = taxa),
    shape = 1, size = 1.5, alpha = 0.6
  ) +
  # GLMM fitted curve
  geom_line(
    data      = pred_panels,
    aes(x = factor(year), y = fit, color = taxa, group = 1),
    linewidth = 1.1
  ) +
  geom_point(
    data  = pred_panels,
    aes(x = factor(year), y = fit, color = taxa),
    size  = 2
  ) +
  facet_grid(
    panel_label ~ period,
    scales = "free"
  ) +
  scale_color_manual(values = coral_colors) +
  labs(
    title   = "Coral Recruitment by Site, Habitat, and Taxa",
    x       = "Year",
    y       = "Recruits per 5 m²",
    caption = paste0(
      "Bars = mean across transects. Open circles = individual transect values.\n",
      "Curves = glmer.nb GLMM fit (1 | transect_id), population-level prediction.\n",
      "Panels without a curve had insufficient data to fit a model."
    )
  ) +
  theme_light() +
  theme(
    axis.title        = element_text(size = 11),
    axis.text         = element_text(size = 8),
    axis.text.x       = element_text(angle = 45, vjust = 0.5),
    strip.text.y      = element_text(size = 8),
    strip.text.x      = element_text(size = 10, face = "bold"),
    strip.background  = element_rect(fill = "#325963"),
    legend.position   = "none",
    plot.caption      = element_text(size = 8, color = "grey50")
  )

# ── Save individual site plots if needed ─────────────────────────────────────
# To produce one figure per site (less crowded), loop over sites:
#
for (s in unique(plot_data$site)) {
  p <- ggplot() + ... + filter(site == s) ...
  ggsave(paste0("recruitment_", s, ".png"), p, width = 10, height = 8)
}
