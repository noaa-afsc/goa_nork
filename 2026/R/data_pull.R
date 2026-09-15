# 2026 northern rockfish data pulls
# update assessment
# ben.williams@noaa.gov
# 2026-09

# notes: 
# catch and survey age composition data are updated in the model
# weight and size at age are additionally updated
# the model is transitioned to use `RTMButils`` functions

# load ----
library(afscdata)
library(afscassess)

# globals ----
year = 2026
species = "NORK"
area = "goa"
rec_age = 2
plus_age = 45
# norpac_species = 303
ages = 2:45
length_bins = lenbins = 15:45

# setup 
# setup_folders(year)
# accepted_model(2024, "m24", 2026)
# dir.create(here::here(year, "m24_2026"))
# file.copy(here::here(year, "base", "model.RDS"),
#           here::here(year, "m24_2026")) # copy the RTMB model 

# data ----
goa_nork(year) 

# tables 
db = afscdata::connect()
afscdata::q_psc(year=2026, target='k', area='goa', db=db, save = TRUE) 
afscdata::q_nontarget(year=2026, target='k', area='goa', db=db, save = TRUE) 
afscdata::q_specs(year=year, species = 'NORK', area ='goa', db=db)
afscdata::disconnect(db)

# data processing
tac = afscassess::get_tac(year, area = area)
afscassess::clean_catch(year, species, TAC = tac)
afscassess::fish_age_comp(year, rec_age = rec_age, plus_age = plus_age,
  lenbins = lengths)
afscassess::fish_length_comp(year = year, rec_age = rec_age, lenbins = length_bins)
afscassess::bts_gap_age_comp(year, rec_age = rec_age, plus_age = plus_age)
afscassess::bts_gap_length_comp(year = year, lenbins = length_bins)

age_data <- read.csv(here::here(year, 'data', 'raw', "goa_bts_gap_specimen_data.csv"))
length_data <- read.csv(here::here(year, 'data', 'raw', "goa_bts_gap_length_data.csv"))
afscassess::saa_waa(year = year, age_data = age_data,
  length_data = length_data,  rec_age = rec_age, 
  len_bins = length_bins,  save = T)
