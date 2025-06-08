#########################################################
# Script for cleaning th HCRIS and COVID Reported Patient 
# Impact and Hospital Capacity by Facility datasets.
# Produces Covid19_funding_data_cleaned.csv, which can 
# be used for visualization and estimation directly.
#########################################################
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(haven)
library(dplyr)
library(lubridate)

#### ---- Data Prep ---- ####
### Treatment eligibility
rand_full <- read_dta("rand_hcris_ffy_hosp_a_2024_11_01.dta") # data source: https://www.hospitaldatasets.org/

# Subset to federal fiscal year 2019
rand <- subset(rand_full, ffy==2019)

# We focus on acute care and critical access hospitals that are either located in one of the 50 states or Washington DC
rand <- rand[!(rand$state_name %in% c("Puerto Rico", "US Virgin Islands")), ]

# Relevant variables: hospital IDs, DPP, UCC, beds, and profit margin
rand <- rand[c("prvdr_num", "sum_pctg_ssi_mdcd_days", "nonmdcr_uncomp_costs_only10", "beds", "total_margin")]

# Impute 0 for missing values, as per Appendix H of Narita & Yata (2023)
rand$sum_pctg_ssi_mdcd_days[is.na(rand$sum_pctg_ssi_mdcd_days)] <- 0
rand$nonmdcr_uncomp_costs_only10[is.na(rand$nonmdcr_uncomp_costs_only10)] <- 0
rand$total_margin[is.na(rand$total_margin)] <- 0

# Removed rows with NA beds, as per Appendix H of Narita & Yata (2023)
rand <- rand[!is.na(rand$beds), ]

# Delete rows with negative UCC, then compute UCC per bed
rand <- rand[rand$nonmdcr_uncomp_costs_only10>=0,]
rand$ucc_per_bed <- rand$nonmdcr_uncomp_costs_only10/rand$beds

### Outcome variable: Total reports of patients currently hospitalized in an adult inpatient bed who have laboratory-confirmed or suspected COVID-19, including those in observation beds reported during the week spanning August 2nd to August 8th 2020
Y <- read.csv("COVID-19_Reported_Patient_Impact_and_Hospital_Capacity_by_Facility_--_RAW_20240927.csv")
Y <- Y %>% 
  mutate(collection_week = ymd(collection_week))
week <- Y %>% 
  filter(collection_week == "2020-08-02") %>% 
  select(ccn, total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_sum, total_adult_patients_hospitalized_confirmed_covid_7_day_sum) %>% 
  rename(prvdr_num = ccn)

# As per Appendix H of Narita & Yata (2023), we impute the value # Confirmed / Suspected COVID Patients and # Confirmed COVID Patients to missing if the latter is greater than the former. Negative values are also imputed tp missing.
week$total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_sum[week$total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_sum<0] <- NA

week$total_adult_patients_hospitalized_confirmed_covid_7_day_sum[week$total_adult_patients_hospitalized_confirmed_covid_7_day_sum<0] <- NA

week$total_adult_patients_hospitalized_confirmed_covid_7_day_sum[is.na(week$total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_sum)] <- NA

week$total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_sum[!is.na(week$total_adult_patients_hospitalized_confirmed_covid_7_day_sum) & week$total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_sum<week$total_adult_patients_hospitalized_confirmed_covid_7_day_sum] <- NA

# The following step deletes 1 hospital (prvdr_num = 050515), which has 2 occurrences in the "week" data frame. The 2 occurrences have quite different values for total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_sum and total_adult_patients_hospitalized_confirmed_covid_7_day_sum. This step also deletes 14 rows with blank prvdr_num, as these rows cannot be matched with the RAND data anyway.
week <- week %>%
  group_by(prvdr_num) %>%
  filter(n() == 1) %>%
  ungroup()

### Combine treatment eligibility and outcome by provider number
df <- rand %>%
  left_join(week[,c(1,2)], by = "prvdr_num")

### Add treatment indicator
df$safety_net <- df$sum_pctg_ssi_mdcd_days>=0.202 & df$ucc_per_bed>=25000 & df$total_margin<=0.03

### Keep complete entries only
data <- df[complete.cases(df$sum_pctg_ssi_mdcd_days, df$ucc_per_bed, df$total_margin, df$safety_net, df$total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_sum),] # n = 3855

### Write combined & cleaned data
write.csv(data, "Covid19_funding_data_cleaned.csv", row.names = F)
