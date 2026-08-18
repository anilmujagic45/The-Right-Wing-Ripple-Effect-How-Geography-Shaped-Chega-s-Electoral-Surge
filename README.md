# The Right-Wing Ripple Effect: How Geography Shaped Chega's Electoral Surge

**Team Basta** — Manuel Maria Fróis, Flemming Kastrop, Anil Mujagic, Hugo Miguel Simões, Michele Rinaldi
*ISEG, March 2026*

## Overview

This project investigates the spatial socio-economic drivers behind the rise of the Portuguese right-wing populist party **Chega**, using a dataset of 278 mainland Portuguese municipalities (2016–2024). We combine PCA, clustering, and spatial econometrics to test whether Chega's electoral surge (2019–2024) is driven by demographic decline, economic frustration, digital mobilisation, and/or spatial diffusion across neighbouring municipalities.

## Methodology

PCA (4 components: *Ghost Town*, *Utopia*, *Chega Campaign*, *TikTok Gen*) → K-Means clustering → Moran's I / LISA spatial diagnostics → spatial regression models (OLS, SEM, SLM, **SDM** preferred, SARAR robustness) → Geographically Weighted Regression (GWR) for spatial non-stationarity.

## Key Findings

- Strong positive spatial autocorrelation in Chega's vote change (Global Moran's I = 0.460).
- The SDM removes residual spatial autocorrelation (ρ = 0.700, p < 0.001) and is the preferred model.
- **TikTok Gen** (digital/youth mobilisation) shows the largest total effect.
- **Ghost Town** (demographic decline) has a *negative* effect on Chega growth.
- **Chega Campaign** (economic precarity) is insignificant globally but a strong localized driver in the Alentejo/south (via GWR).
- Strong spatial dependence confirms spatial diffusion in Chega's electoral growth.

## References

Anselin (1988, 2019, 2020), Kaiser (1960), Kelejian & Prucha (2010), Moran (1950). Full bibliography in the report.
