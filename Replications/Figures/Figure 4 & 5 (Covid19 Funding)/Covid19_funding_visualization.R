#############################
# Figure 4: Visualization of 
# the CARES Act data.
#############################
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(plotly)

#### ---- Visualization ---- ####
data <- read.csv("Covid19_funding_data_cleaned.csv")

# For better display, adjust the variable ranges
dpp_range <- c(quantile(data$sum_pctg_ssi_mdcd_days, .025),
               quantile(data$sum_pctg_ssi_mdcd_days, .975))
ucc_range <- c(quantile(data$ucc_per_bed, .025),
               quantile(data$ucc_per_bed, .975))
margin_range <- c(quantile(data$total_margin, .025),
                  quantile(data$total_margin, .975))

# Create the plot
p <- plot_ly(data = data) %>%
  add_trace(
    type = "scatter3d",
    mode = "markers",
    x = ~sum_pctg_ssi_mdcd_days,
    y = ~ucc_per_bed,
    z = ~total_margin,
    color = ~safety_net,  # color depends on this variable
    colors = c("black", "red"),  # choose any two colors
    marker = list(
      size = 2,
      opacity = 0.5
    )
  ) %>%
  layout(
    scene = list(
      xaxis = list(title = "DPP", range = dpp_range),
      yaxis = list(title = "UCC per bed", range = ucc_range),
      zaxis = list(title = "Profit Margin", range = margin_range)
    ),
    showlegend = FALSE,
    title = 'Illustration of the Treatment Boundary'
  )

# Plane 1: DPP fixed at 0.202
p <- p %>%
  add_trace(
    type = "mesh3d",
    x = c(0.202, 0.202, 0.202, 0.202),
    y = c(25000, max(ucc_range), max(ucc_range), 25000),
    z = c(min(margin_range), min(margin_range), 0.03, 0.03),
    i = c(0, 1, 2, 3),
    j = c(1, 2, 3, 0),
    k = c(2, 3, 0, 1),
    facecolor = rep("pink", 4),
    opacity = 0.25
  )

# Plane 2: UCC per bed fixed at 25000
p <- p %>%
  add_trace(
    type = "mesh3d",
    x = c(0.202, max(dpp_range), max(dpp_range), 0.202),
    y = c(25000, 25000, 25000, 25000),
    z = c(min(margin_range), min(margin_range), 0.03, 0.03),
    i = c(0, 1, 2, 3),
    j = c(1, 2, 3, 0),
    k = c(2, 3, 0, 1),
    facecolor = rep("pink", 4),
    opacity = 0.25
  )

# Plane 3: Profit Margin fixed at 0.03
p <- p %>%
  add_trace(
    type = "mesh3d",
    x = c(0.202, max(dpp_range), max(dpp_range), 0.202),
    y = c(25000, 25000, max(ucc_range), max(ucc_range)),
    z = c(0.03, 0.03, 0.03, 0.03),
    i = c(0, 1, 2, 3),
    j = c(1, 2, 3, 0),
    k = c(2, 3, 0, 1),
    facecolor = rep("pink", 4),
    opacity = 0.25
  )

# Line 1: DPP fixed at 0.202, UCC per bed fixed at 25000
p <- p %>%
  add_trace(
    type = "scatter3d",
    mode = "lines",
    x = c(0.202, 0.202),
    y = c(25000, 25000),
    z = c(min(margin_range), 0.03),
    line = list(color = "deepskyblue", width = 5)
  )

# Line 2: DPP fixed at 0.202, Profit Margin fixed at 0.03
p <- p %>%
  add_trace(
    type = "scatter3d",
    mode = "lines",
    x = c(0.202, 0.202),
    y = c(25000, max(ucc_range)),
    z = c(0.03, 0.03),
    line = list(color = "deepskyblue", width = 5)
  )

# Line 3: UCC per bed fixed at 25000, Profit Margin fixed at 0.03
p <- p %>%
  add_trace(
    type = "scatter3d",
    mode = "lines",
    x = c(0.202, max(dpp_range)),
    y = c(25000, 25000),
    z = c(0.03, 0.03),
    line = list(color = "deepskyblue", width = 5)
  )
p
