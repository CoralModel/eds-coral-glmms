# MCR LTER Coral Demography — GLMM Analysis

<img width="1880" height="714" alt="acr" src="https://github.com/user-attachments/assets/44d13c08-0049-44a8-a915-5b3769fb4c35" />

<br>

This repository contains taxa-specific generalized linear mixed models (GLMMs) estimating rates of **recruitment**, **growth**, and **survival** for three coral taxa (*Acropora*, *Pocillopora*, *Porites*) across backreef and forereef habitats at the Moorea Coral Reef (MCR) LTER site. Models are fit before and after the 2019 bleaching event to assess how this disturbance affected coral population dynamics island-wide.

## Research questions

Each vital rate is examined under two model structures:

- **Model 1 (Habitat):** Do rates differ between forereef and backreef habitats before and after the 2019 bleaching event?
- **Model 2 (Site × Habitat):** Do rates differ across all site-habitat combinations (LTER1 backreef, LTER1 forereef, LTER2 backreef, LTER2 forereef, etc.) before and after the 2019 bleaching event?

## Model details

**Taxa:** *Acropora*, *Pocillopora*, *Porites* (separate models per taxa)

**Growth** models the size-dependence of proportional growth rate using `log10(size_t1)` as the response and `log10(size_t0)` as a covariate. Random effects account for colonies nested in transects nested in sites, plus inter-annual variation.

**Recruitment** and **survival** scripts follow the same two-model structure across the same taxa, sites, and bleaching periods.

**Reference groups:** Backreef / before bleaching (Model 1); LTER1 backreef / before bleaching (Model 2).

## Repository structure

```
├── growth_glmm.qmd          # Size-dependent growth rate models
├── recruitment_glmm.qmd     # Recruitment rate models
├── survival_glmm.qmd        # Survival rate models
├── eds-coral-glmms.Rproj
├── README.md
└── LICENSE
```

> **Note:** Input data (`volume_corals.csv`) is read from the companion data cleaning repository [`eds-coral-data-storage-management`](https://github.com/) and is not stored here.

## Data access

Raw data are sourced from the **Moorea Coral Reef (MCR) LTER** and are not yet publicly available. Processed input data are produced by the companion data cleaning repository. 

## Authors

- [Kylie Newcomer](https://github.com/kylienewcomer)
- [Joaquin Sandoval](https://github.com/sandovaljoaquin)
- [Vedika Shirtekar](https://github.com/vedikaS-byte)

## References & acknowledgements

Data are from the **Moorea Coral Reef Long Term Ecological Research** (MCR LTER) site, funded by the National Science Foundation.
