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

######____________ predict
library(tidyverse)
library(glmmTMB)

# ── Prediction Grids 
# Formula: recruits_int ~ habitat * pre_post_f + (1|site/transect_id) + (1|year)
# re.form = NA to marginalise over all random effects island-wide
# type = "response"  to back transform from log scale to count scale


# ── Acropora 

acr_pred_grid_q1 <- dat_acr %>%
  distinct(habitat, pre_post_f)

acr_pred_grid_q1$predicted <- predict(q1_acr, newdata = acr_pred_grid_q1,
                                      type = "response", re.form = NA)

acr_pred_grid_q1


### manual math
# BR Pre  = intercept only (reference level)
exp(-3.8370)                                    # = 0.0216 (match)

# BR Post = intercept + post effect
exp(-3.8370 + 0.8367)                           # = 0.0498 (match)

# OR Pre  = intercept + OR effect
exp(-3.8370 + 3.0602)                           # = 0.460  (match)

# OR Post = intercept + OR + post + interaction
exp(-3.8370 + 3.0602 + 0.8367 + (-1.5736))     # = 0.220  (match)


# ── Pocillopora 
poc_pred_grid_q1 <- dat_poc %>%
  distinct(habitat, pre_post_f)

poc_pred_grid_q1$predicted <- predict(q1_poc, newdata = poc_pred_grid_q1,
                                      type = "response", re.form = NA)

poc_pred_grid_q1


# ── Porites 

por_pred_grid_q1 <- dat_por %>%
  distinct(habitat, pre_post_f)

por_pred_grid_q1$predicted <- predict(q1_por, newdata = por_pred_grid_q1,
                                      type = "response", re.form = NA)

por_pred_grid_q1


# ── Combine for plotting 

q1_pred_all <- bind_rows(
  acr_pred_grid_q1 %>% mutate(taxa = "Acropora"),
  poc_pred_grid_q1 %>% mutate(taxa = "Pocillopora"),
  por_pred_grid_q1 %>% mutate(taxa = "Porites")
) %>%
  mutate(
    taxa       = factor(taxa, levels = c("Acropora", "Pocillopora", "Porites")),
    pre_post_f = factor(pre_post_f, levels = c("Pre", "Post"))
  )


# ── Bar plot 
# ******predict() gives point estimates only... so no se (how to get??), use emmeans??
# For se bars need emmeans() or bootstrap 

period_colors <- c(Pre = "#4e9af1", Post = "#f17a4e")

ggplot(q1_pred_all,
       aes(x = habitat, y = predicted, fill = pre_post_f)) +
  geom_col(position = position_dodge(width = 0.75),
           width = 0.65, colour = "white", linewidth = 0.3) +
  facet_wrap(~ taxa, scales = "free_y", ncol = 3) +
  scale_fill_manual(
    values = period_colors,
    name   = "Period",
    labels = c(Pre = "Pre-bleaching (2013–2018)",
               Post = "Post-bleaching (2019–2024)")
  ) +
  labs(
    title    = "Island-wide coral recruitment by habitat and bleaching period",
    subtitle = "Marginal predicted counts (re.form = NA), no SE available from predict()",
    x        = "Habitat",
    y        = "Predicted recruits per 5 m²"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background   = element_rect(fill = "grey92", colour = NA),
    strip.text         = element_text(face = "bold"),
    legend.position    = "bottom",
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(face = "bold"),
    plot.subtitle      = element_text(colour = "grey40")
  )


# Add se to grid
acr_se <- predict(q1_acr, newdata = acr_pred_grid_q1,
                  type = "response", re.form = NA, se.fit = TRUE)
acr_pred_grid_q1$se <- acr_se$se.fit

# Same for poc and por
poc_se <- predict(q1_poc, newdata = poc_pred_grid_q1,
                  type = "response", re.form = NA, se.fit = TRUE)
poc_pred_grid_q1$se <- poc_se$se.fit

por_se <- predict(q1_por, newdata = por_pred_grid_q1,
                  type = "response", re.form = NA, se.fit = TRUE)
por_pred_grid_q1$se <- por_se$se.fit

# Combine
q1_pred_all <- bind_rows(
  acr_pred_grid_q1 %>% mutate(taxa = "Acropora"),
  poc_pred_grid_q1 %>% mutate(taxa = "Pocillopora"),
  por_pred_grid_q1 %>% mutate(taxa = "Porites")
) %>%
  mutate(
    taxa       = factor(taxa, levels = c("Acropora", "Pocillopora", "Porites")),
    pre_post_f = factor(pre_post_f, levels = c("Pre", "Post"))
  )

# Plot
ggplot(q1_pred_all,
       aes(x = habitat, y = predicted, fill = pre_post_f)) +
  geom_col(position = position_dodge(width = 0.75),
           width = 0.65, colour = "white", linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = predicted - 1.96 * se,
        ymax = predicted + 1.96 * se),
    position = position_dodge(width = 0.75),
    width = 0.25, linewidth = 0.6
  ) +
  facet_wrap(~ taxa, scales = "free_y", ncol = 3) +
  scale_fill_manual(
    values = c(Pre = "#4e9af1", Post = "#f17a4e"),
    name   = "Period",
    labels = c(Pre = "Pre-bleaching (2013–2018)",
               Post = "Post-bleaching (2019–2024)")
  ) +
  labs(
    title = "Island-wide coral recruitment by habitat and bleaching period",
    x     = "Habitat",
    y     = "Predicted recruits per 5 m²"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background   = element_rect(fill = "grey92", colour = NA),
    strip.text         = element_text(face = "bold"),
    legend.position    = "bottom",
    panel.grid.major.x = element_blank(),
    plot.title         = element_text(face = "bold")
  )

