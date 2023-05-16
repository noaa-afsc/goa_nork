# 2023 goa northern rockfish assessment code
# ben.williams@noaa.gov

library(afscdata)

# globals ----
year = 2023
species = "NORK"
region = "goa"

# setup 
setup_folders(year)
accepted_model(2022, "m22.1", 2023)

# query data ----
goa_dusk(year, off_yr = TRUE)
