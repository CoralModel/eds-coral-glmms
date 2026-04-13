library(tidyverse)
library(lme4)
library(glmmTMB)
library(ggtext)
library(ggeffects)

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Load raw data
# ══════════════════════════════════════════════════════════════════════════════

clean_coral <- read_csv(
  "~/MEDS/capstone/eds-coral-figs-storage/vedika-r-files/data/updated_coral_tidy_2013-2024.csv"
)


# ══════════════════════════════════════════════════════════════════════════════
# Step 2: Filter and clean transect names
# ══════════════════════════════════════════════════════════════════════════════

recruit_raw <- clean_coral %>%
  filter(!grepl("^P", transect)) %>%
  dplyr::select(site, habitat, transect, dyn_recruitment, year, taxa) %>%
  mutate(
    # Standardize duplicate and variant transect names
    transect = case_when(
      transect == "T02_OLD"     ~ "T02",
      transect == "T03_PRE2020" ~ "T03",
      transect == "T04_OLD"     ~ "T04",
      transect == "T02_2021"    ~ "T02",
      transect == "TO1"         ~ "T01",
      TRUE                      ~ transect
    ),
    # Recode variant taxa labels
    taxa = case_when(
      taxa == "A" ~ "Acr",
      TRUE        ~ taxa
    )
  )


# ══════════════════════════════════════════════════════════════════════════════
# Step 3: Summarize at transect level
#
# dyn_recruitment is binary per coral row; sum gives total recruits per
# transect x site x habitat x taxa x year combo.
# Dividing by 5 standardizes to recruits per 5 m².
#
# transect_id = paste(site, transect) creates a unique transect identifier
# so T01 at LTER1 is treated as distinct from T01 at LTER2.
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
    period       = ifelse(pre_post == 0, "Pre (2013-2018)", "Post (2019-2024)"),
    # year_s scaled 0-1 GLOBALLY across full 2013-2024 range so slopes are
    # comparable across periods
    year_s       = (year - min(year)) / (max(year) - min(year)),
    transect_id  = paste(site, transect, sep = "_"),
    pre_post_f   = factor(pre_post, levels = c(0, 1),
                          labels = c("Pre", "Post")),
    taxa         = factor(taxa),
    site         = factor(site),
    habitat      = factor(habitat)
  ) %>%
  filter(!is.na(year), !is.na(pre_post))

# Sanity checks
cat("Years in data:", sort(unique(recruit_transect$year)), "\n")
cat("Any NA years:", anyNA(recruit_transect$year), "\n")
cat("Unique transect IDs:\n")
print(sort(unique(recruit_transect$transect_id)))


# ══════════════════════════════════════════════════════════════════════════════
# Step 4: Sparsity check — identify taxa too sparse to model reliably
#
# Taxa with very few non-zero observations in a given period will cause
# complete separation and inestimable coefficients. Flag these before modeling.
# ══════════════════════════════════════════════════════════════════════════════

sparsity_check <- recruit_transect %>%
  group_by(taxa, period) %>%
  summarize(
    n_rows       = n(),
    n_nonzero    = sum(recruits_int > 0),
    pct_nonzero  = round(100 * mean(recruits_int > 0), 1),
    .groups      = "drop"
  )
print(sparsity_check)
# Acr and Mil are typically too sparse (all-zero in multiple combos).
# Poc and Por have sufficient signal. Sparse taxa are skipped automatically
# in the modeling loop below (Step 5) when they fail convergence.


# ══════════════════════════════════════════════════════════════════════════════
# Step 5: Overdispersion helper
# ══════════════════════════════════════════════════════════════════════════════

check_overdisp <- function(model, label) {
  pr    <- residuals(model, type = "pearson")
  ratio <- sum(pr^2) / df.residual(model)
  cat(label, "overdispersion ratio:", round(ratio, 3), "\n")
  # >> 1: overdispersed (nbinom2 appropriate)
  # ~= 1: not overdispersed (Poisson appropriate)
  # << 1: underdispersed
}


# ══════════════════════════════════════════════════════════════════════════════
# Step 6: Per-taxa GLMM loop
#
# Two parallel questions, each answered with a separate model per taxa:
#
#   Q1 — Habitat differences (BR vs OR) over time:
#        Fixed:  habitat * year_s
#        Random: (1 | site/transect_id)
#        → site accounts for between-site clustering so habitat estimate is clean.
#          Transects nest within sites as the physical structure dictates.
#
#   Q2 — Site differences over time:
#        Fixed:  site * year_s
#        Random: (1 | habitat/transect_id)
#        → habitat accounts for BR vs OR clustering so site estimates are clean.
#
# year_s * fixed_effect interaction asks: does the TRAJECTORY of recruitment
# over time differ between habitats (Q1) or between sites (Q2)?
#
# Family chosen per taxa per period based on overdispersion ratio:
#   nbinom2 if ratio >> 1, poisson if ratio ~= 1.
#   Default here is nbinom2; switch to poisson if diagnose() flags it.
#
# Models that fail to converge or have insufficient data are stored as NULL
# and reported at the end.
# ══════════════════════════════════════════════════════════════════════════════

taxa_list <- levels(recruit_transect$taxa)  # Acr, Mil, Poc, Por

# Storage lists
q1_models <- list()   # habitat differences
q2_models <- list()   # site differences

for (tx in taxa_list) {

  cat("\n══════════════════════════════════════════════\n")
  cat("Taxa:", tx, "\n")
  cat("══════════════════════════════════════════════\n")

  dat <- recruit_transect %>% filter(taxa == tx)

  # Quick sparsity gate: skip if fewer than 10 non-zero observations total
  if (sum(dat$recruits_int > 0) < 10) {
    cat("  Skipping", tx, "— too sparse (<10 non-zero obs)\n")
    q1_models[[tx]] <- NULL
    q2_models[[tx]] <- NULL
    next
  }

  # ── Q1: Habitat differences over time ──────────────────────────────────────
  cat("  Fitting Q1 (habitat * year_s | site/transect_id)...\n")

  q1_models[[tx]] <- tryCatch(
    glmmTMB(
      recruits_int ~ habitat * year_s + (1 | site/transect_id),
      family    = nbinom2,
      data      = dat,
      na.action = na.exclude
    ),
    error   = function(e) { cat("  Q1 ERROR:", conditionMessage(e), "\n"); NULL },
    warning = function(w) { cat("  Q1 WARNING:", conditionMessage(w), "\n"); NULL }
  )

  if (!is.null(q1_models[[tx]])) {
    cat("  Q1 summary:\n")
    print(summary(q1_models[[tx]]))
    check_overdisp(q1_models[[tx]], paste("  Q1", tx))
    diagnose(q1_models[[tx]])
  }

  # ── Q2: Site differences over time ─────────────────────────────────────────
  cat("  Fitting Q2 (site * year_s | habitat/transect_id)...\n")

  q2_models[[tx]] <- tryCatch(
    glmmTMB(
      recruits_int ~ site * year_s + (1 | habitat/transect_id),
      family    = nbinom2,
      data      = dat,
      na.action = na.exclude
    ),
    error   = function(e) { cat("  Q2 ERROR:", conditionMessage(e), "\n"); NULL },
    warning = function(w) { cat("  Q2 WARNING:", conditionMessage(w), "\n"); NULL }
  )

  if (!is.null(q2_models[[tx]])) {
    cat("  Q2 summary:\n")
    print(summary(q2_models[[tx]]))
    check_overdisp(q2_models[[tx]], paste("  Q2", tx))
    diagnose(q2_models[[tx]])
  }
}

# Convergence report
cat("\n══════════════════════════════════════════════\n")
cat("Convergence summary:\n")
cat("Q1 (habitat) converged:", paste(names(Filter(Negate(is.null), q1_models)), collapse = ", "), "\n")
cat("Q1 (habitat) failed:   ", paste(names(Filter(is.null,         q1_models)), collapse = ", "), "\n")
cat("Q2 (site)    converged:", paste(names(Filter(Negate(is.null), q2_models)), collapse = ", "), "\n")
cat("Q2 (site)    failed:   ", paste(names(Filter(is.null,         q2_models)), collapse = ", "), "\n")


# ══════════════════════════════════════════════════════════════════════════════
# Step 7: Marginal effects — predicted recruitment over time
#
# ggpredict() marginalizes over the random effects, giving the population-
# average trajectory for each habitat (Q1) or each site (Q2).
# ══════════════════════════════════════════════════════════════════════════════

for (tx in taxa_list) {

  # Q1 — habitat trajectories
  if (!is.null(q1_models[[tx]])) {
    cat("\nQ1 marginal effects —", tx, "(habitat x year_s):\n")
    print(ggpredict(q1_models[[tx]], terms = c("year_s [all]", "habitat")))
  }

  # Q2 — site trajectories
  if (!is.null(q2_models[[tx]])) {
    cat("\nQ2 marginal effects —", tx, "(site x year_s):\n")
    print(ggpredict(q2_models[[tx]], terms = c("year_s [all]", "site")))
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# Step 8: Spaghetti plots — one per taxa
#
# Each plot shows recruitment over time per transect (thin lines = raw data
# per transect) with a population-level GLMM fitted curve overlaid, faceted
# by habitat (Q1 framing).
#
# Layout: facet_wrap(~ habitat) so BR and OR are side by side.
# Color:  site, so between-site spread is visible within each habitat panel.
# ══════════════════════════════════════════════════════════════════════════════

coral_colors <- c(
  LTER1 = "#1f77b4",
  LTER2 = "#ff7f0e",
  LTER3 = "#2ca02c",
  LTER4 = "#d62728",
  LTER5 = "#9467bd",
  LTER6 = "#8c564b"
)

output_dir <- "~/MEDS/capstone/eds-coral-figs-storage/vedika-r-files/figs/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

for (tx in taxa_list) {

  dat <- recruit_transect %>% filter(taxa == tx)

  # Per-transect raw means for spaghetti lines
  transect_means <- dat %>%
    group_by(site, habitat, transect, transect_id, year, year_s) %>%
    summarize(mean_recruits = mean(recruits_std, na.rm = TRUE), .groups = "drop")

  # Population-level predictions from Q1 model (habitat framing)
  pred_q1 <- if (!is.null(q1_models[[tx]])) {
    preds <- ggpredict(q1_models[[tx]], terms = c("year_s [all]", "habitat"))
    # Back-convert year_s to calendar year for x-axis readability
    year_range <- range(dat$year)
    as_tibble(preds) %>%
      rename(year_s = x, habitat = group, fit = predicted) %>%
      mutate(year = year_s * (year_range[2] - year_range[1]) + year_range[1])
  } else {
    NULL
  }

  p <- ggplot() +
    # Spaghetti: one line per transect
    geom_line(
      data  = transect_means,
      aes(x = year, y = mean_recruits, group = transect_id, color = site),
      alpha = 0.45, linewidth = 0.6
    ) +
    geom_point(
      data  = transect_means,
      aes(x = year, y = mean_recruits, color = site),
      alpha = 0.5, size = 1.5
    ) +
    # Population-level GLMM curve (Q1 model)
    { if (!is.null(pred_q1))
        geom_line(
          data      = pred_q1,
          aes(x = year, y = fit),
          color     = "black",
          linewidth = 1.2,
          linetype  = "dashed"
        )
    } +
    facet_wrap(~ habitat, labeller = label_both) +
    scale_color_manual(values = coral_colors, name = "Site") +
    scale_x_continuous(breaks = seq(2013, 2024, by = 2)) +
    labs(
      title   = paste0(tx, " recruitment over time by habitat"),
      x       = "Year",
      y       = "Recruits per 5 m²",
      caption = paste0(
        "Thin lines = individual transect trajectories. ",
        "Dashed black = population-level GLMM fit\n",
        "(habitat * year_s, random: 1 | site/transect_id).",
        if (is.null(q1_models[[tx]])) " Model did not converge." else ""
      )
    ) +
    theme_light() +
    theme(
      axis.title       = element_text(size = 11),
      axis.text        = element_text(size = 9),
      axis.text.x      = element_text(angle = 45, vjust = 0.5),
      strip.text       = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "#325963"),
      plot.title       = element_text(size = 13, face = "bold"),
      plot.caption     = element_text(size = 8, color = "grey50"),
      legend.position  = "right"
    )

  filename <- paste0(output_dir, "recruitment_spaghetti_", tx, ".png")
  ggsave(filename, p, width = 10, height = 5, dpi = 300)
  cat("Saved:", filename, "\n")
  print(p)
}
