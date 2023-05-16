# 2023 goa northern rockfish assessment code
# ben.williams@noaa.gov

library(afscdata)

# globals ----
year = 2023
species = "NORK"
region = "goa"

# setup 
setup_folders(year)

# query data ----
akfin <- connect()
goa_dusk(year, off_yr = TRUE)