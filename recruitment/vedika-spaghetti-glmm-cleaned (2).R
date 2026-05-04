library(tidyverse)
library(glmmTMB)
library(ggeffects)

# Q1: Island-wide habitat differences for pre vs post bleaching
#     Formula = habitat * pre_post_f + (1|site/transect_id) + (1|year)
#

# Families: Acr = Poisson, Poc/Por = nbinom2
# Reference levels: habitat = BR, period = Pre, site = LTER1


# ── Load data ──────────────────────────────────────────────────────────────

# clean_coral <- read_csv(
#   "~/MEDS/capstone/eds-coral-figs-storage/vedika-r-files/data/updated_coral_tidy_2013-2024.csv"
# )


clean_coral <- read_csv(
  "~/MEDS/capstone/eds-coral-data-storage-management/analysisdata/volume_corals.csv"
)


# raw_data <- read_csv(
#   "~/MEDS/capstone/eds-coral-data-storage-management/data/coral_tidy_dyn_2013-2024.csv")
# 
# unique(raw_data)
# raw_data %>% filter(habitat == "LTER4")
# 
# view(raw_data)

# ── Step 2: Clean ─────────────────────────────────────────────────────────────

recruit_raw <- clean_coral %>%
  filter(!grepl("^P", transect)) %>%
  dplyr::select(site, habitat, transect, dyn_recruitment, year, taxa) %>%
  mutate(
    transect = case_when(
      transect == "T02_OLD"     ~ "T02",
      transect == "T03_PRE2020" ~ "T03",
      transect == "T04_OLD"     ~ "T04",
      transect == "T02_2021"    ~ "T02",
      transect == "TO1"         ~ "T01",
      TRUE                      ~ transect
    ),
    taxa = if_else(taxa == "A", "Acr", taxa)
  )


# ── Step 3: Summarize to transect and year level ──────────────────────────────────────

recruit_transect <- recruit_raw %>%
  group_by(taxa, site, habitat, transect, year) %>%
  summarize(n_recruits = sum(dyn_recruitment, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    recruits_std = n_recruits / 5,
    recruits_int = as.integer(round(recruits_std)),
    pre_post_f   = factor(
      case_when(
        year >= 2013 & year <= 2018 ~ "Pre",
        year >= 2019 & year <= 2024 ~ "Post"
      ),
      levels = c("Pre", "Post")
    ),
    transect_id = paste(site, transect, sep = "_"),
    taxa        = factor(taxa),
    site        = relevel(factor(site),    ref = "LTER1"),
    habitat     = relevel(factor(habitat), ref = "BR")
  ) %>%
  filter(!is.na(pre_post_f))

# Per-taxa and per period subsets (a for loop or function for cleaner production?)
dat_acr      <- recruit_transect %>% filter(taxa == "Acr")
dat_poc      <- recruit_transect %>% filter(taxa == "Poc")
dat_por      <- recruit_transect %>% filter(taxa == "Por")
# dat_acr_pre  <- recruit_transect %>% filter(taxa == "Acr", pre_post_f == "Pre")
# dat_acr_post <- recruit_transect %>% filter(taxa == "Acr", pre_post_f == "Post")
# dat_poc_pre  <- recruit_transect %>% filter(taxa == "Poc", pre_post_f == "Pre")
# dat_poc_post <- recruit_transect %>% filter(taxa == "Poc", pre_post_f == "Post")
# dat_por_pre  <- recruit_transect %>% filter(taxa == "Por", pre_post_f == "Pre")
# dat_por_post <- recruit_transect %>% filter(taxa == "Por", pre_post_f == "Post")
# 

# ── Step 4: Data structure checks ─────────────────────────────────────────────

cat("Sparsity: nonzero observations per taxa × period combo:")
recruit_transect %>%
  group_by(taxa, pre_post_f) %>%
  summarize(
    n           = n(),
    n_nonzero   = sum(recruits_int > 0),
    pct_nonzero = round(100 * mean(recruits_int > 0), 1),
    .groups     = "drop"
  ) %>% print()


cat("Nonzero obs per taxa × site × habitat × period:")
cat("zero = structural zero, so breaks site*habitat interaction\n\n")
recruit_transect %>%
  filter(taxa != "Mil") %>%
  group_by(taxa, site, habitat, pre_post_f) %>%
  summarize(
    n_nonzero = sum(recruits_int > 0),
    mean_rec  = round(mean(recruits_std), 3),
    .groups   = "drop"
  ) %>%
  mutate(flag = if_else(n_nonzero == 0, "zero", "")) %>%
  print(n = Inf)


# ── Step 5: Exploratory bar charts ───────────────────────────────────────────

recruit_transect %>%
  filter(taxa != "Mil") %>%
  ggplot(aes(x = factor(year), y = recruits_std, fill = pre_post_f)) +
  geom_col() +
  facet_grid(habitat ~ taxa) +
  scale_fill_manual(values = c(Pre = "#4e9af1", Post = "#f17a4e"),
                    name = "Period") +
  theme_bw() +
  labs(title = "Recruitment by habitat and taxa",
       x = "Year", y = "Recruits per 5 m²")

recruit_transect %>%
  filter(taxa != "Mil") %>%
  mutate(site_habitat = paste(site, habitat, sep = "-")) %>%
  ggplot(aes(x = factor(year), y = recruits_std, fill = pre_post_f)) +
  geom_col() +
  facet_grid(site_habitat ~ taxa) +
  scale_fill_manual(values = c(Pre = "#4e9af1", Post = "#f17a4e"),
                    name = "Period") +
  theme_bw() +
  labs(title = "Recruitment by site-habitat combination and taxa",
       x = "Year", y = "Recruits per 5 m²")



## Model validation and assumptions 
# 1. for negative binomial, we will need to test for the overdispersion (qqplot), heteroscedacity, and zero-inflation
# 2. Temporal autocorrelation??


#_______________________________________________________________________
# Q1: Island wide habitat x period 
# Formula: habitat * pre_post_f + (1|site/transect_id) + (1|year)
# site/transect_id: sites pooled as random... island-wide habitat question
# year crossed random effect: absorbs annual variation

# ── Q1 Acr ───────────────────────────────────────────────────────────────────

cat("Q1 Acr ~ habitat * pre_post_f | site/transect_id + year\n")
cat("Family: Poisson\n")

q1_acr <- glmmTMB(
  recruits_int ~ habitat * pre_post_f + (1 | site/transect_id) + (1 | year),
  family    = nbinom2,
  data      = dat_acr,
  na.action = na.exclude
)

summary(q1_acr)
## Pre-bleaching, OR recruited much more than BR (habitat OR is 3.08 on log scale and significant). 
# Gap shrank post-bleaching. 

cat("Overdispersion ratio:\n")
pr_q1_acr <- residuals(q1_acr, type = "pearson")
cat(round(sum(pr_q1_acr^2) / df.residual(q1_acr), 3), "\n")

cat("\nQ1 Acr marginal effects: habitat × period:\n")
print(ggpredict(q1_acr, terms = c("habitat", "pre_post_f"), bias_correction = T)) # Adjust for Jensen's inequality based on console feedback


# ── Q1 Poc ───────────────────────────────────────────────────────────────────

cat("Q1 Poc ~ habitat * pre_post_f | site/transect_id + year\n")
cat("Family: nbinom2\n")

q1_poc <- glmmTMB(
  recruits_int ~ habitat * pre_post_f + (1 | site/transect_id) + (1 | year),
  family    = nbinom2,
  data      = dat_poc,
  na.action = na.exclude
)

summary(q1_poc)
testDispersion(q1_poc)
# Pre-bleaching OR recruited more than BR. Post-bleaching had flipped pattern: or dropped sharply while BR held (interaction is p < .0001). 
# or recruitment seemed to drop significantly post bleaching 

cat("\nOverdispersion ratio:\n")
pr_q1_poc <- residuals(q1_poc, type = "pearson")
cat(round(sum(pr_q1_poc^2) / df.residual(q1_poc), 3), "\n")

cat("\nQ1 Poc marginal effects — habitat × period:\n")
print(ggpredict(q1_poc, terms = c("habitat", "pre_post_f")))




# ── Q1 Por ───────────────────────────────────────────────────────────────────

cat("Q1 Por ~ habitat * pre_post_f | site/transect_id + year\n")
cat("Family: nbinom2\n")


q1_por <- glmmTMB(
  recruits_int ~ habitat * pre_post_f + (1 | site/transect_id) + (1 | year),
  family    = nbinom2,
  data      = dat_por,
  na.action = na.exclude
)

summary(q1_por)
# Opposit pattern seen: BR recruited much more than OR pre-bleaching (habitatOR = -1.91 and is significant). 
# Both habitats droppped post bleaching (pre_post is significant); however, the gap is not significant since interaction is not significant. 
# Observing marginal effects show that BR droppped from 2.43 to .25, Or from .36 to .06

cat("\nOverdispersion ratio:\n")
pr_q1_por <- residuals(q1_por, type = "pearson")
cat(round(sum(pr_q1_por^2) / df.residual(q1_por), 3), "\n")

cat("\nQ1 Por marginal effects: habitat × period:\n")
print(ggpredict(q1_por, terms = c("habitat", "pre_post_f")))



## Then general trend is that all three taxa show strong bleaching effects; however, the habitat story differs by taxa. 
# Acr and Poc are both favored OR pre-bleaching and saw OR collapse post-bleaching. 
#Por saw that BR is favored through and both habitats have crashed recruitment following the bleaching event

#_________________________________
# Q2A: Site x habitat x period for a three-way interaction
# Formula = site * habitat * pre_post_f + (1|transect_id) + (1|year)

# site and habitat are fixed... transect_id only in random structure
# Poc: converges cleanly, use this as primary Q2 model for Poc
# Por: two three-way terms have blown SEs (structural zeros LTER2/4-OR post)
#      so use Q2B for Por instead
# Acr: not attempted here since structural zeros make this not estimable

# ── Q2A Poc ──────────────────────────────────────────────────────────────────
cat("Q2A Poc ~ site * habitat * pre_post_f | transect_id + year\n")
cat("Family: nbinom2\n")
cat("Primary Q2 model for Poc — converges cleanly\n")

q2a_poc <- glmmTMB(
  recruits_int ~ site * habitat * pre_post_f + (1 | transect_id) + (1 | year),
  family    = nbinom2,
  data      = dat_poc,
  na.action = na.exclude
)

summary(q2a_poc)

cat("\nOverdispersion ratio:\n")
pr_q2a_poc <- residuals(q2a_poc, type = "pearson")
cat(round(sum(pr_q2a_poc^2) / df.residual(q2a_poc), 3), "\n")

cat("\nQ2A Poc marginal effects — site × habitat × period:\n")
print(ggpredict(q2a_poc, terms = c("site", "habitat", "pre_post_f")))


## This is the cleanest model for the three way...significant three-way interactions show that site-habitat ranking shifted post-bleaching. 
# Key coefficients are: habitatOR:pre_post_fPost = -2.98 (p<0.001) (OR recruitment crashed post-bleaching at reference site LTER1)...
# but siteLTER2:habitatOR:pre_post_fPost = 2.24 (p=0.049) and siteLTER5:habitatOR:pre_post_fPost = 3.37 (p=0.006) meaning LTER2-OR and LTER5-OR did NOT crash as badly as LTER1-OR. 
#  marginal means show that: 
# LTER1-OR went from 1.46 to 0.07, but LTER5-OR only went from 0.90 to 0.41. site-habitat ranking seems to have shifted







# ── Q2A Por ──────────────────────────────────────────────────────────────────

cat("Q2A Por ~ site * habitat * pre_post_f | transect_id + year\n")
cat("Family: nbinom2\n")
cat("Note: LTER2:habitatOR:pre_post_fPost and LTER4:habitatOR:pre_post_fPost\n")
cat("have blown SEs (structural zeros post-bleaching). Use Q2B for Por.\n")
cat("════════════════════════════════════════════════════\n")

q2a_por <- glmmTMB(
  recruits_int ~ site * habitat * pre_post_f + (1 | transect_id) + (1 | year),
  family    = nbinom2,
  data      = dat_por,
  na.action = na.exclude
)

summary(q2a_por)

cat("\nOverdispersion ratio:\n")
pr_q2a_por <- residuals(q2a_por, type = "pearson")
cat(round(sum(pr_q2a_por^2) / df.residual(q2a_por), 3), "\n")

cat("\nQ2A Por marginal effects — site × habitat × period:\n")
print(ggpredict(q2a_por, terms = c("site", "habitat", "pre_post_f")))


#_____________________
# Q2B: sitex habitat (separate pre and post) S
# Formula: site * habitat + (1|transect_id) + (1|year)
# Used for Por (avoids three-way terms from structural zeros)
# and Acr (separate pre/post shown for completeness- coefficients with
# huge SEs on BR cells reported as not estimable, this could be a biological finding)
# ══════════════════════════════════════════════════════════════════════════════

# ── Q2B Por PRE ──────────────────────────────────────────────────────────────
cat("Q2B Por PRE ~ site * habitat | transect_id + year\n")
cat("Family: nbinom2 | Period: Pre-bleaching (2013-2018)\n")

q2b_por_pre <- glmmTMB(
  recruits_int ~ site * habitat + (1 | transect_id) + (1 | year),
  family    = nbinom2,
  data      = dat_por_pre,
  na.action = na.exclude
)

summary(q2b_por_pre)

cat("\nOverdispersion ratio:\n")
pr_q2b_por_pre <- residuals(q2b_por_pre, type = "pearson")
cat(round(sum(pr_q2b_por_pre^2) / df.residual(q2b_por_pre), 3), "\n")

cat("\nQ2B Por PRE marginal effects — site × habitat:\n")
print(ggpredict(q2b_por_pre, terms = c("site", "habitat")))


# ── Q2B Por POST ─────────────────────────────────────────────────────────────

cat("Q2B Por POST ~ site * habitat | transect_id + year\n")
cat("Family: nbinom2 | Period: Post-bleaching (2019-2024)\n")
cat("Note: LTER2:habitatOR and LTER4:habitatOR have large SEs —\n")
cat("structural zeros (no Por recruitment in OR at those sites post-bleaching)\n")

q2b_por_post <- glmmTMB(
  recruits_int ~ site * habitat + (1 | transect_id) + (1 | year),
  family    = nbinom2,
  data      = dat_por_post,
  na.action = na.exclude
)

summary(q2b_por_post)

cat("\nOverdispersion ratio:\n")
pr_q2b_por_post <- residuals(q2b_por_post, type = "pearson")
cat(round(sum(pr_q2b_por_post^2) / df.residual(q2b_por_post), 3), "\n")

cat("\nQ2B Por POST marginal effects — site × habitat:\n")
ggpredict(q2b_por_post, terms = c("site", "habitat")) %>%
  as_tibble() %>%
  mutate(conf.high = pmin(conf.high, 10)) %>%
  print()

## General summary: siteLTER2:habitatOR:pre_post_fPost and siteLTER4:habitatOR:pre_post_fPost have SEs in the thousands...those are structural zeros (LTER2-OR and LTER4-OR had zero Por recruits post-bleaching). 
# The siteLTER5:habitatOR:pre_post_fPost term is fine (p=0.005)...the  Q2B pre/post split handles this better, pre-bleaching shows LTER4-OR significantly lower than LTER1-OR (siteLTER4:habitatOR = -2.25, p<0.001), LTER5-OR also lower (p=0.012). 
# Post-bleaching, LTER5-OR recovered relative to others (p=0.047) while LTER2-OR and LTER4-OR had near zero recruitment



#___________________________________________________________
# ── Q2B Acr pre shown here for justification of completeness but BR  not estimable
cat("Q2B Acr PRE ~ site * habitat | transect_id + year\n")
cat("Family: Poisson | Period: Pre-bleaching\n")
cat("Note: BR cells have structural zeros at LTER1 and LTER5.\n")
cat("Coefficients involving BR will have NaN SEs — not estimable.\n")
cat("OR cells are interpretable.\n")
cat("════════════════════════════════════════════════════\n")

q2b_acr_pre <- glmmTMB(
  recruits_int ~ site * habitat + (1 | transect_id) + (1 | year),
  family    = poisson(),
  data      = dat_acr_pre,
  na.action = na.exclude
)

summary(q2b_acr_pre)

cat("\nOverdispersion ratio:\n")
pr_q2b_acr_pre <- residuals(q2b_acr_pre, type = "pearson")
cat(round(sum(pr_q2b_acr_pre^2) / df.residual(q2b_acr_pre), 3), "\n")

cat("\nQ2B Acr PRE marginal effects — OR sites only (BR not estimable):\n")
print(ggpredict(q2b_acr_pre, terms = c("site", "habitat")))


# ── Q2B Acr Post shown for completeness but some BR not estimable
cat("Q2B Acr POST ~ site * habitat | transect_id + year\n")
cat("Family: Poisson | Period: Post-bleaching\n")
cat("Note: LTER2-BR and LTER5-BR are structural zeros post-bleaching.\n")
cat("LTER1-BR and LTER4-BR are estimable. All OR cells estimable.\n")
q2b_acr_post <- glmmTMB(
  recruits_int ~ site * habitat + (1 | transect_id) + (1 | year),
  family    = poisson(),
  data      = dat_acr_post,
  na.action = na.exclude
)

summary(q2b_acr_post)

cat("\nOverdispersion ratio:\n")
pr_q2b_acr_post <- residuals(q2b_acr_post, type = "pearson")
cat(round(sum(pr_q2b_acr_post^2) / df.residual(q2b_acr_post), 3), "\n")

cat("\nQ2B Acr post marginal effects — estimable cells only:\n")
print(ggpredict(q2b_acr_post, terms = c("site", "habitat")))

## Generally, all models failed for Acr... Acr may not be estimable from the data... Acr sitexhabitat interactions are not estimable from the data... not enough signal in BR... 
# This may mean that recruitment disappeared from BR afterpost-bleaching at these sites

# Model convergence fails here without splitting.... also SE is huge??
summary(glmmTMB(
  recruits_int ~ site * habitat * pre_post_f + (1 | transect_id) + (1 | year),
  family    = nbinom2,
  data      = dat_acr,
  na.action = na.exclude
))











# Plots
# Raw data as time series backdrop with GLMM marginal means per period overlaid

period_colors <- c(Pre = "#4e9af1", Post = "#f17a4e")
site_colors   <- c(LTER1 = "#1f77b4", LTER2 = "#ff7f0e",
                   LTER4 = "#2ca02c", LTER5 = "#9467bd")
taxa_labels   <- c(Acr = "Acropora", Poc = "Pocillopora", Por = "Porites")

# ── Plot Q1 ───────────────────────────────────────────────────────────────────

raw_q1 <- recruit_transect %>%
  filter(taxa %in% c("Acr", "Poc", "Por")) %>%
  group_by(taxa, habitat, year, pre_post_f) %>%
  summarize(mean_std = mean(recruits_std, na.rm = TRUE), .groups = "drop") %>%
  mutate(taxa = factor(taxa, levels = c("Acr", "Poc", "Por")))

preds_q1_acr <- ggpredict(q1_acr, terms = c("habitat", "pre_post_f")) %>%
  as_tibble() %>%
  rename(habitat = x, pre_post_f = group) %>%
  mutate(taxa = "Acr",
         year_mid = if_else(pre_post_f == "Pre", 2015.5, 2021.5))

preds_q1_poc <- ggpredict(q1_poc, terms = c("habitat", "pre_post_f")) %>%
  as_tibble() %>%
  rename(habitat = x, pre_post_f = group) %>%
  mutate(taxa = "Poc",
         year_mid = if_else(pre_post_f == "Pre", 2015.5, 2021.5))

preds_q1_por <- ggpredict(q1_por, terms = c("habitat", "pre_post_f")) %>%
  as_tibble() %>%
  rename(habitat = x, pre_post_f = group) %>%
  mutate(taxa = "Por",
         year_mid = if_else(pre_post_f == "Pre", 2015.5, 2021.5))

preds_q1 <- bind_rows(preds_q1_acr, preds_q1_poc, preds_q1_por) %>%
  mutate(taxa = factor(taxa, levels = c("Acr", "Poc", "Por")))

p_q1 <- ggplot() +
  geom_vline(xintercept = 2018.5, linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  annotate("text", x = 2019, y = Inf, label = "← Bleaching",
           hjust = 0, vjust = 2, size = 2.8, color = "grey40") +
  geom_line(
    data = raw_q1,
    aes(x = year, y = mean_std, color = pre_post_f,
        group = interaction(habitat, pre_post_f)),
    alpha = 0.5, linewidth = 0.6
  ) +
  geom_point(
    data = raw_q1,
    aes(x = year, y = mean_std, color = pre_post_f),
    alpha = 0.6, size = 1.8
  ) +
  geom_pointrange(
    data = preds_q1,
    aes(x = year_mid, y = predicted,
        ymin = conf.low, ymax = conf.high,
        color = pre_post_f),
    size = 0.9, linewidth = 1.3
  ) +
  facet_grid(
    habitat ~ taxa,
    labeller = labeller(taxa = taxa_labels),
    scales   = "free_y"
  ) +
  scale_color_manual(values = period_colors, name = "Period") +
  scale_x_continuous(breaks = seq(2013, 2024, by = 3)) +
  labs(
    title    = "Q1: Island-wide habitat recruitment pre vs post bleaching",
    subtitle = "Small points/lines = raw annual means  |  Large points + CI = GLMM marginal mean per period",
    x = "Year", y = "Recruits per 5 m²",
    caption  = "Model: habitat × period + (1|site/transect_id) + (1|year). Acr = Poisson; Poc, Por = nbinom2."
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#2c5f6e"),
    strip.text       = element_text(color = "white", face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "bottom"
  )

print(p_q1)
ggsave(
  "~/MEDS/capstone/eds-coral-figs-storage/vedika-r-files/figs/Q1_habitat_period.png",
  p_q1, width = 10, height = 6, dpi = 300
)


# ── Plot Q2 ───────────────────────────────────────────────────────────────────
# Poc: Q2A combined predictions split by period
# Por: Q2B separate pre/post predictions
# Acr: Q2B predictions shown with caveat that BR cells are not estimable

raw_q2 <- recruit_transect %>%
  filter(taxa %in% c("Acr", "Poc", "Por")) %>%
  group_by(taxa, site, habitat, pre_post_f) %>%
  summarize(mean_std = mean(recruits_std, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    taxa         = factor(taxa, levels = c("Acr", "Poc", "Por")),
    pre_post_f   = factor(pre_post_f, levels = c("Pre", "Post")),
    site_habitat = paste(site, habitat, sep = "-")
  )

# Poc predictions from Q2A, split into pre and post rows
preds_q2_poc <- ggpredict(q2a_poc,
                          terms = c("site", "habitat", "pre_post_f")) %>%
  as_tibble() %>%
  rename(site = x, habitat = group, pre_post_f = facet) %>%
  mutate(taxa = "Poc",
         site_habitat = paste(site, habitat, sep = "-"),
         pre_post_f   = factor(pre_post_f, levels = c("Pre", "Post")),
         conf.high    = pmin(conf.high, 10))

# Por predictions from Q2B pre and post separately
preds_q2_por_pre <- ggpredict(q2b_por_pre,
                              terms = c("site", "habitat")) %>%
  as_tibble() %>%
  rename(site = x, habitat = group) %>%
  mutate(taxa = "Por", pre_post_f = "Pre",
         site_habitat = paste(site, habitat, sep = "-"),
         pre_post_f   = factor(pre_post_f, levels = c("Pre", "Post")),
         conf.high    = pmin(conf.high, 10))

preds_q2_por_post <- ggpredict(q2b_por_post,
                               terms = c("site", "habitat")) %>%
  as_tibble() %>%
  rename(site = x, habitat = group) %>%
  mutate(taxa = "Por", pre_post_f = "Post",
         site_habitat = paste(site, habitat, sep = "-"),
         pre_post_f   = factor(pre_post_f, levels = c("Pre", "Post")),
         conf.high    = pmin(conf.high, 10))

# Acr predictions from Q2B pre and post — BR cells will have missing CIs
preds_q2_acr_pre <- ggpredict(q2b_acr_pre,
                              terms = c("site", "habitat")) %>%
  as_tibble() %>%
  rename(site = x, habitat = group) %>%
  mutate(taxa = "Acr", pre_post_f = "Pre",
         site_habitat = paste(site, habitat, sep = "-"),
         pre_post_f   = factor(pre_post_f, levels = c("Pre", "Post")),
         conf.high    = pmin(conf.high, 10))

preds_q2_acr_post <- ggpredict(q2b_acr_post,
                               terms = c("site", "habitat")) %>%
  as_tibble() %>%
  rename(site = x, habitat = group) %>%
  mutate(taxa = "Acr", pre_post_f = "Post",
         site_habitat = paste(site, habitat, sep = "-"),
         pre_post_f   = factor(pre_post_f, levels = c("Pre", "Post")),
         conf.high    = pmin(conf.high, 10))

preds_q2 <- bind_rows(
  preds_q2_poc,
  preds_q2_por_pre,
  preds_q2_por_post,
  preds_q2_acr_pre,
  preds_q2_acr_post
) %>%
  mutate(taxa = factor(taxa, levels = c("Acr", "Poc", "Por")))

p_q2 <- ggplot() +
  geom_point(
    data = raw_q2,
    aes(x = site_habitat, y = mean_std, color = site),
    size = 2.5, alpha = 0.5,
    position = position_jitter(width = 0.15, seed = 42)
  ) +
  geom_pointrange(
    data = preds_q2,
    aes(x = site_habitat, y = predicted,
        ymin = conf.low, ymax = conf.high),
    color = "black", size = 0.6, linewidth = 1
  ) +
  facet_grid(
    taxa ~ pre_post_f,
    labeller = labeller(taxa = taxa_labels),
    scales   = "free_y"
  ) +
  scale_color_manual(values = site_colors, name = "Site") +
  labs(
    title    = "Q2: Site × habitat recruitment pre vs post bleaching",
    subtitle = "Colored points = raw period means  |  Black points + CI = GLMM marginal mean",
    x = "Site-habitat", y = "Recruits per 5 m²",
    caption  = paste(
      "Poc: Q2A site × habitat × period three-way (nbinom2).",
      "Por: Q2B separate pre/post site × habitat (nbinom2); LTER2/4-OR post = structural zeros.",
      "Acr: Q2B separate pre/post; BR cells at LTER1/2/5 not estimable (structural zeros).",
      sep = "\n"
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#2c5f6e"),
    strip.text       = element_text(color = "white", face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "bottom"
  )

print(p_q2)
ggsave(
  "~/MEDS/capstone/eds-coral-figs-storage/vedika-r-files/figs/Q2_site_habitat_period.png",
  p_q2, width = 11, height = 8, dpi = 300
)