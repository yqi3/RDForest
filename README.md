# Replication files for Liu & Qi (2025)
This repo contains replication files for the figures and tables in [Liu & Qi (2025)](https://arxiv.org/abs/2303.11721).

For replication of the simulations, please download and unzip the file "Data.zip" from [here](https://drive.google.com/file/d/1_enCyBQjVbAcFvbLCOin4J2D42NkepOk/view?usp=sharing), and place the unzipped folder under the `Replications` directory.

For the simulations, implementations for local linear regressions using [`rdrobust`](https://rdpackages.github.io/rdrobust/) and for random/local linear forests using [`grf`](https://grf-labs.github.io/grf/reference/index.html) are done via
```
# R version 4.3.2 (2023-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux trixie/sid
```

Implementations for the minimax-optimal estimator using [`optrdd`](https://github.com/swager/optrdd) are done via
 ```
R version 4.1.3 (2022-03-10)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Monterey 12.5
```

For the empirical application using the Covid19 Funding data ([data source 1](https://drive.google.com/file/d/1_enCyBQjVbAcFvbLCOin4J2D42NkepOk/view?usp=sharing), [data source 2](https://healthdata.gov/Hospital/COVID-19-Reported-Patient-Impact-and-Hospital-Capa/uqq2-txqb/about_data)), please download `rand_hcris_ffy_hosp_a_2024_11_01.dta` and `COVID-19_Reported_Patient_Impact_and_Hospital_Capacity_by_Facility_--_RAW_20240927.csv` and place them under the `Replications/Figures/Figure 4 & 5 (Covid19 Funding)` directory. Alternatively, users could directly use `Covid19_funding_data_cleaned.csv`, which is our combined and cleaned data for visualization and estimation. Implementations for random/local linear forests using [`grf`](https://grf-labs.github.io/grf/reference/index.html) are done via
```
R version 4.4.1 (2024-06-14)
Platform: x86_64-apple-darwin20
Running under: macOS Monterey 12.7.6
```

Note: Reproducibility of [`grf`](https://grf-labs.github.io/grf/reference/index.html) depends on both seed and computing platform. Exact reproducibility is therefore not guaranteed. See [here](https://grf-labs.github.io/grf/REFERENCE.html#forests-predict-different-values-depending-on-the-platform-even-though-the-seed-is-the-same) for more details.
