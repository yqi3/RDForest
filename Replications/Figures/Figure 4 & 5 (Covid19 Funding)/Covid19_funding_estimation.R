######################################
# Figure 5: Estimated Treatment Effect
# with CIs Using RFs & LLFs
######################################
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(grf) # NOTE: Requires packageVersion("grf")>=2.4.0. Results may not be exactly reproducible across platforms; see: https://grf-labs.github.io/grf/REFERENCE.html#forests-predict-different-values-depending-on-the-platform-even-though-the-seed-is-the-same
library(gridExtra)
source("../../Helper Fns/KT_Helper_Fns.R")
set.seed(1234)

#### ---- Estimation ---- ####
data <- read.csv("Covid19_funding_data_cleaned.csv")

## Separate into treated and control samples. Running variables: DPP, UCC per bed, and profit margin.
trt <- data[data$safety_net==1,] # n_trt = 656
ctrl <- data[data$safety_net==0,] # n_ctrl = 3199
feature_trt <- as.data.frame(trt[,c(2,6,5)])
feature_ctrl <- as.data.frame(ctrl[,c(2,6,5)])

## Build RFs
B <- 1000 # number of trees
forestT <- tuned_RDForest(feature_trt, trt$total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_sum, 
                          tune.parameters = c(), 
                          type = "trt",
                          num.trees=B)
forestC <- tuned_RDForest(feature_ctrl, ctrl$total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_sum, 
                          tune.parameters = c(), 
                          type = "ctrl",
                          num.trees=B)

## Function for estimating and obtaining CI for treatment effect with given focal point
rf_rdd <- function(focal) {
  predictionT <- predict(forestT, focal, estimate.variance = TRUE)
  predictionC <- predict(forestC, focal, estimate.variance = TRUE)
  prediction_rf <- predictionT$predictions - predictionC$predictions
  se_BLB <- sqrt(predictionT$variance.estimates+predictionC$variance.estimates)
  CI_rf <- c(prediction_rf-1.96*se_BLB, prediction_rf+1.96*se_BLB)
  return(list(prediction_rf, CI_rf))
}


## Generate treatment effect plots: 1) Fix UCC per bed and profit margin, vary DPP
n_grid_points <- 100
DPP_grid <- seq(0.202, quantile(data$sum_pctg_ssi_mdcd_days, 0.9), length.out = 100)
point_estimates <- numeric(n_grid_points)
lower_CI <- numeric(n_grid_points)
upper_CI <- numeric(n_grid_points)

# Loop over the DPP grid and apply the rf_rdd function
for (i in 1:n_grid_points) {
  focal_input <- matrix(c(DPP_grid[i], 25000, 0.03), nrow = 1)
  result <- rf_rdd(focal_input)
  point_estimates[i] <- result[[1]]
  lower_CI[i] <- result[[2]][1]
  upper_CI[i] <- result[[2]][2]
}

# Combine results into a data frame for plotting
DPP_results_rf <- data.frame(
  DPP = DPP_grid,
  point_estimate = point_estimates,
  lower_CI = lower_CI,
  upper_CI = upper_CI
)

## Generate treatment effect plots: 2) Fix DPP and profit margin, vary UCC per bed
UCC_per_bed_grid <- seq(25000, quantile(data$ucc_per_bed, 0.9), length.out = 100)
point_estimates <- numeric(n_grid_points)
lower_CI <- numeric(n_grid_points)
upper_CI <- numeric(n_grid_points)

# Loop over the UCC per bed grid and apply the rf_rdd function
for (i in 1:n_grid_points) {
  focal_input <- matrix(c(0.202, UCC_per_bed_grid[i], 0.03), nrow = 1)
  result <- rf_rdd(focal_input)
  point_estimates[i] <- result[[1]]
  lower_CI[i] <- result[[2]][1]
  upper_CI[i] <- result[[2]][2]
}

# Combine results into a data frame for plotting
UCC_per_bed_results_rf <- data.frame(
  UCC_per_bed = UCC_per_bed_grid,
  point_estimate = point_estimates,
  lower_CI = lower_CI,
  upper_CI = upper_CI
)

## Generate treatment effect plots: 3) Fix DPP and UCC per bed, vary profit margin
profit_margin_grid <- seq(quantile(data$total_margin, 0.1), 0.03, length.out = 100)
point_estimates <- numeric(n_grid_points)
lower_CI <- numeric(n_grid_points)
upper_CI <- numeric(n_grid_points)

# Loop over the profit margin grid and apply the rf_rdd function
for (i in 1:n_grid_points) {
  focal_input <- matrix(c(0.202, 25000, profit_margin_grid[i]), nrow = 1)
  result <- rf_rdd(focal_input)
  point_estimates[i] <- result[[1]]
  lower_CI[i] <- result[[2]][1]
  upper_CI[i] <- result[[2]][2]
}

# Combine results into a data frame for plotting
profit_margin_results_rf <- data.frame(
  profit_margin = profit_margin_grid,
  point_estimate = point_estimates,
  lower_CI = lower_CI,
  upper_CI = upper_CI
)


### Estimation using LLF
## Set up parameters for local linear forests
d <- 3
pi_ <- 1
lower_bound_default <- 1-(1+ (d/(1.3*pi_))*log(0.05)/log(1-0.05))^(-1)  # subsample rate for d=2, alpha=0.05; modified for local linear forests with improved rates
buffer <- 1e-30 # buffer to adjust the focal point to be inside the interior of the support

## Build LLFs
forestT <- ll_regression_forest(feature_trt, trt$total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_sum, 
                                num.trees=B, mtry=1,
                                honesty = TRUE, 
                                enable.ll.split = TRUE,
                                ll.split.weight.penalty = TRUE,
                                sample.fraction = 0.4*ceiling(nrow(trt)^lower_bound_default)/nrow(trt))

forestC <- ll_regression_forest(feature_ctrl, ctrl$total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_sum, 
                                num.trees=B, mtry=1,
                                honesty = TRUE, 
                                enable.ll.split = TRUE,
                                ll.split.weight.penalty = TRUE,
                                sample.fraction = 0.4*ceiling(nrow(ctrl)^lower_bound_default)/nrow(ctrl))

## Function for estimating and obtaining CI for treatment effect with given focal point
llf_rdd <- function(focal) {
  predT <- predict(forestT, focal+c(buffer, buffer, -buffer), 
                   estimate.variance = TRUE,
                   ll.weight.penalty = TRUE
  )
  predictionT <- as.numeric(predT$predictions)
  predC <- predict(forestC, focal+c(-buffer, -buffer, buffer), 
                   estimate.variance = TRUE,
                   ll.weight.penalty = TRUE
  )
  predictionC <- as.numeric(predC$predictions)
  prediction_llf <- predictionT - predictionC
  se_BLB <- sqrt(predT$variance.estimates+predC$variance.estimates)
  CI_llf <- c(prediction_llf-1.96*se_BLB, prediction_llf+1.96*se_BLB)
  return(list(prediction_llf, CI_llf))
}

## Generate treatment effect plots: 1) Fix UCC per bed and profit margin, vary DPP
point_estimates <- numeric(n_grid_points)
lower_CI <- numeric(n_grid_points)
upper_CI <- numeric(n_grid_points)

# Loop over the DPP grid and apply the llf_rdd function
for (i in 1:n_grid_points) {
  focal_input <- matrix(c(DPP_grid[i], 25000, 0.03), nrow = 1)
  result <- llf_rdd(focal_input)
  point_estimates[i] <- result[[1]]
  lower_CI[i] <- result[[2]][1]
  upper_CI[i] <- result[[2]][2]
}

# Combine results into a data frame for plotting
DPP_results_llf <- data.frame(
  DPP = DPP_grid,
  point_estimate = point_estimates,
  lower_CI = lower_CI,
  upper_CI = upper_CI
)

## Generate treatment effect plots: 2) Fix DPP and profit margin, vary UCC per bed
point_estimates <- numeric(n_grid_points)
lower_CI <- numeric(n_grid_points)
upper_CI <- numeric(n_grid_points)

# Loop over the UCC per bed grid and apply the llf_rdd function
for (i in 1:n_grid_points) {
  focal_input <- matrix(c(0.202, UCC_per_bed_grid[i], 0.03), nrow = 1)
  result <- llf_rdd(focal_input)
  point_estimates[i] <- result[[1]]
  lower_CI[i] <- result[[2]][1]
  upper_CI[i] <- result[[2]][2]
}

# Combine results into a data frame for plotting
UCC_per_bed_results_llf <- data.frame(
  UCC_per_bed = UCC_per_bed_grid,
  point_estimate = point_estimates,
  lower_CI = lower_CI,
  upper_CI = upper_CI
)

## Generate treatment effect plots: 3) Fix DPP and UCC per bed, vary profit margin
point_estimates <- numeric(n_grid_points)
lower_CI <- numeric(n_grid_points)
upper_CI <- numeric(n_grid_points)

# Loop over the profit margin grid and apply the llf_rdd function
for (i in 1:n_grid_points) {
  focal_input <- matrix(c(0.202, 25000, profit_margin_grid[i]), nrow = 1)
  result <- llf_rdd(focal_input)
  point_estimates[i] <- result[[1]]
  lower_CI[i] <- result[[2]][1]
  upper_CI[i] <- result[[2]][2]
}

# Combine results into a data frame for plotting
profit_margin_results_llf <- data.frame(
  profit_margin = profit_margin_grid,
  point_estimate = point_estimates,
  lower_CI = lower_CI,
  upper_CI = upper_CI
)

## RF plots
plot_DPP_range <- c(min(DPP_results_rf$lower_CI, DPP_results_llf$lower_CI),
                    max(DPP_results_rf$upper_CI, DPP_results_llf$upper_CI))
plot_UCC_range <- c(min(UCC_per_bed_results_rf$lower_CI, UCC_per_bed_results_llf$lower_CI),
                    max(UCC_per_bed_results_rf$upper_CI, UCC_per_bed_results_llf$upper_CI))

plot_PM_range <- c(min(profit_margin_results_rf$lower_CI, profit_margin_results_llf$lower_CI),
                   max(profit_margin_results_rf$upper_CI, profit_margin_results_llf$upper_CI))

DPP_plot_rf <- ggplot(DPP_results_rf, aes(x = DPP)) +
  geom_line(aes(y = point_estimate), color = "blue") +  # Point estimates line
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "blue", alpha = 0.2) +  # CI bands
  labs(x = "DPP", y = "Treatment Effect Estimate",
       title = "Varying DPP") +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(plot_DPP_range)

UCC_per_bed_plot_rf <- ggplot(UCC_per_bed_results_rf, aes(x = UCC_per_bed)) +
  geom_line(aes(y = point_estimate), color = "blue") +  # Point estimates line
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "blue", alpha = 0.2) +  # CI bands
  labs(x = "UCC per Bed", y = "Treatment Effect Estimate",
       title = "Varying UCC per Bed") +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(plot_UCC_range)

profit_margin_plot_rf <- ggplot(profit_margin_results_rf, aes(x = profit_margin)) +
  geom_line(aes(y = point_estimate), color = "blue") +  # Point estimates line
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "blue", alpha = 0.2) +  # CI bands
  labs(x = "Profit Margin", y = "Treatment Effect Estimate",
       title = "Varying Profit Margin") +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(plot_PM_range)

ggsave(paste0("emp_rf_temp.png"),
       grid.arrange(DPP_plot_rf, UCC_per_bed_plot_rf, profit_margin_plot_rf, ncol = 3),
       width = 16, height = 6)

## LLF plots
DPP_plot_llf <- ggplot(DPP_results_llf, aes(x = DPP)) +
  geom_line(aes(y = point_estimate), color = "blue") +  # Point estimates line
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "blue", alpha = 0.2) +  # CI bands
  labs(x = "DPP", y = "Treatment Effect Estimate",
       title = "Varying DPP") +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(plot_DPP_range)

UCC_per_bed_plot_llf <- ggplot(UCC_per_bed_results_llf, aes(x = UCC_per_bed)) +
  geom_line(aes(y = point_estimate), color = "blue") +  # Point estimates line
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "blue", alpha = 0.2) +  # CI bands
  labs(x = "UCC per Bed", y = "Treatment Effect Estimate",
       title = "Varying UCC per Bed") +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(plot_UCC_range)

profit_margin_plot_llf <- ggplot(profit_margin_results_llf, aes(x = profit_margin)) +
  geom_line(aes(y = point_estimate), color = "blue") +  # Point estimates line
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "blue", alpha = 0.2) +  # CI bands
  labs(x = "Profit Margin", y = "Treatment Effect Estimate",
       title = "Varying Profit Margin") +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5)) + ylim(plot_PM_range)

ggsave(paste0("emp_llf_temp.png"),
       grid.arrange(DPP_plot_llf, UCC_per_bed_plot_llf, profit_margin_plot_llf, ncol = 3),
       width = 16, height = 6)
