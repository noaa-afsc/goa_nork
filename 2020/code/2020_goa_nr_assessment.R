# run 2020 Gulf of Alaska northern rockfish assessment

# ben.williams@noaa.gov
# 2020-10

# see README.md in the 2020 folder for complete model descriptions

# load ----
library(rockfishr)
# detach("package:rockfishr", unload=TRUE)

# globals ----

year = 2020
species <- "NORK"
rec_age <- 2
plus_age <- 45
TAC <- c(3786, 3681, 4528)
admb_home <- "C:/Program Files (x86)/ADMB-12.1"
model_name <- "updated_nr"
dat_name <- "goa_nr_2020"

afsc_user = "WILLIAMSB"
afsc_pwd = "H@ines2010_2021"
akfin_user = "bwilliams"
akfin_pwd = "jnu$6350"

mcmc = 1e+07
mcsave = 2000

# data ----
# setup folders 
modeldir(year)

# query databases 
raw_data(species, year, afsc_user, afsc_pwd, akfin_user, akfin_pwd)

# clean and process datatdata ----
# catch
clean_catch(year, TAC)

# biological data 
# old age-error matrix ----
# using the old reader_tester file
ageage(reader_tester = "reader_tester_old.csv", species, year, admb_home) # file provide in "user_input" folder
fish_age_comp(year, rec_age, plus_age)
survey_age_comp(year, rec_age, plus_age)
fish_size_comp(year, rec_age)
survey_size_comp(year)
size_at_age(year, admb_home, rec_age)
weight_at_age(year, admb_home, rec_age)

# note: must provide file for VAST or DB estimates are output
# design-based model
survey_biomass(year) 
concat_dat(year, "db", rec_age, plus_age)

# VAST "bridge model"
survey_biomass(year, "VAST_estimate_mesa.csv") 
concat_dat(year, "m18.2", rec_age, plus_age)

# VAST GAP assessment
survey_biomass(year, "VAST_estimates.csv") 
concat_dat(year, "m18.2a", rec_age, plus_age)

# new age-error matrix ----
# rerun functions as some are dependent on age error matrix output
# updated (as of 2020-10) reader_tester file
ageage(reader_tester = NULL, species, year, admb_home) # read_tester file is provide in the "user_input" folder
fish_age_comp(year, rec_age, plus_age)
survey_age_comp(year, rec_age, plus_age)
fish_size_comp(year, rec_age)
survey_size_comp(year)
size_at_age(year, admb_home, rec_age)
weight_at_age(year, admb_home, rec_age)
concat_dat(year, "m18.2b", rec_age, plus_age)

# run models ----

# currently running through the command line - otherwise get odd output file names that do not work with functions

# run_admb <- function(year, model, model_name, mcmc, mcsave){
  setwd(here::here(year, model))
  
  # compile
  R2admb::compile_admb(model_name, verbose = TRUE)
  
  #Run ADMB Model with MCMC
  system(paste(model_name, "-mcmc", mcmc, "-mcsave", mcsave))
  system(paste(model_name, "-mceval"))  # to get "evalout.proj"
  
  # reset environment 
  setwd(here::here())


model = "db"
# cleanup output ----
process_results(year, model, model_name, data_name, rec_age, plus_age, mcmc = mcmc, mcsave = mcsave)

base_plots(year, model, model_name, rec_age)


model = "m18.2"
# cleanup output ----
process_results(year, model, model_name, data_name, rec_age, plus_age, mcmc = mcmc, mcsave = mcsave, survey = "VAST_estimate_mesa.csv")

base_plots(year, model, model_name, rec_age)

model = "m18.2a"
# cleanup output ----
process_results(year, model, model_name, data_name, rec_age, plus_age, mcmc = mcmc, mcsave = mcsave, survey = "VAST_estimates.csv")

base_plots(year, model, model_name, rec_age)


model = "m18.2b"
# cleanup output ----
library(rockfishr)
process_results(year, model, model_name, data_name, rec_age, plus_age, mcmc = mcmc, mcsave = mcsave, survey = "VAST_estimates2.csv")

base_plots(year, model, model_name, rec_age)


# run retro ----
run_retro(year, model = "db", tpl_name = "updated_nr", n_retro = 10, admb_home = admb_home, mcmc = mcmc, mcsave = mcsave)


# best f ----
data <- proj_data(year, model, model_name)  

best_f(data$best, data$m, last_ofl = data$last_ofl, type = 1)


# plot retros ----
plot_retro(year, model)
plot_retro_survey(year, model)
plot_re("Y:/ABL_MESA/SAFES2020/Apportionment/Plot data - NR")




