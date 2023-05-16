library(tidyverse)
library(afscassess)
year = 2022
model_name = "nr"
dat_name = "goa_nr_2022"
rec_age = 2
plus_age = 45
mcmc = 150000


model = "m18.2b"
process_results(year, model, model_name, dat_name, rec_age, plus_age, mcmc=5000000, mcsave=1000, len_bins = "lbins.csv")
run_retro(year, model, model_name, dat_name, mcmc)
plot_params(year, model, model_name)
plot_retro(year, model)
afscassess::plot_biomass(year, model)
afscassess::plot_catch(year, model)
afscassess::plot_phase(year, model, model_name)
afscassess::plot_rec(year, model, rec_age)
afscassess::plot_rec_ssb(year, model, rec_age)
afscassess::plot_selex(year, model)
afscassess::plot_survey(year, model)
afscassess::plot_swath(year, model)

model = "m22"
process_results(year, model, model_name, dat_name, rec_age, plus_age, mcmc=5000000, mcsave=1000, len_bins = "lbins.csv")
run_retro(year, model, model_name, dat_name, mcmc)
plot_params(year, model, model_name)
plot_retro(year, model)
afscassess::plot_biomass(year, model)
afscassess::plot_catch(year, model)
afscassess::plot_phase(year, model, model_name)
afscassess::plot_rec(year, model, rec_age)
afscassess::plot_rec_ssb(year, model, rec_age)
afscassess::plot_selex(year, model)
afscassess::plot_survey(year, model)
afscassess::plot_swath(year, model)


model = "m22.1"
process_results(year, model, model_name, dat_name, rec_age, plus_age, mcmc=5000000, mcsave=1000, len_bins = "lbins2.csv")
run_retro(year, model, model_name, dat_name, mcmc)
plot_params(year, model, model_name)
plot_retro(year, model)
afscassess::plot_comps(year, model)
afscassess::plot_biomass(year, model)
afscassess::plot_catch(year, model)
afscassess::plot_phase(year, model, model_name)
afscassess::plot_rec(year, model, rec_age)
afscassess::plot_rec_ssb(year, model, rec_age)
afscassess::plot_selex(year, model)
afscassess::plot_survey(year, model)
afscassess::plot_swath(year, model)
plot_catch_bio(year, model, model_name)
afscassess::fac_table(year, model)
afscassess::sac_table(year, model)
afscassess::fsc_table(year, model)
afscassess::ssc_table(year, model)



model = "m22.1a"
process_results(year, model, model_name, dat_name, rec_age, plus_age, mcmc=200000, mcsave=200, len_bins = "lbins2.csv")
run_retro(year, model, model_name, dat_name, mcmc, n_retro=10)
plot_params(year, model, model_name)
plot_retro(year, model)
afscassess::plot_comps(year, model)
afscassess::plot_biomass(year, model)
afscassess::plot_catch(year, model)
afscassess::plot_phase(year, model, model_name)
afscassess::plot_rec(year, model, rec_age)
afscassess::plot_rec_ssb(year, model, rec_age)
afscassess::plot_selex(year, model)
afscassess::plot_survey(year, model)
afscassess::plot_swath(year, model)
afscassess::fac_table(year, model)
afscassess::sac_table(year, model)
afscassess::fsc_table(year, model)
afscassess::ssc_table(year, model)


model = "m22.1b"
process_results(year, model, model_name, dat_name, rec_age, plus_age, mcmc=200000, mcsave=200, len_bins = "lbins2.csv")
run_retro(year, model, model_name, dat_name, mcmc)
plot_params(year, model, model_name)
plot_retro(year, model)
afscassess::plot_comps(year, model)
afscassess::plot_biomass(year, model)
afscassess::plot_catch(year, model)
afscassess::plot_phase(year, model, model_name)
afscassess::plot_rec(year, model, rec_age)
afscassess::plot_rec_ssb(year, model, rec_age)
afscassess::plot_selex(year, model)
afscassess::plot_survey(year, model)
afscassess::plot_swath(year, model)
params_table(year, model, model_name)

model = "db"
process_results(year, model, model_name, dat_name, rec_age, plus_age, mcmc=100000, mcsave=100, len_bins = "lbins.csv")

afscassess::plot_comps(year, model)
afscassess::plot_biomass(year, model)
afscassess::plot_catch(year, model)
afscassess::plot_phase(year, model, model_name)
afscassess::plot_rec(year, model, rec_age)
afscassess::plot_rec_ssb(year, model, rec_age)
afscassess::plot_selex(year, model)
afscassess::plot_survey(year, model)
afscassess::plot_swath(year, model)
params_table(year, model, model_name)


plot_compare_survey(year, 
                    models = c("2022, base", "2022, m18.2b", "2022, m22"))

plot_compare_survey(year, 
                    models = c("2022, m22", "2022, m22.1", "2022, m22.1a", "2022, m22.1b"))

plot_compare_biomass(year, 
                    models = c("2022, m22", "2022, m22.1", "2022, m22.1a", "2022, m22.1b"))




png(filename=here::here(year, "figs", "fcomparesrv.png"), width = 6.5, height = 6.5, 
    units = "in", type ="cairo", res = 200)
plot_compare_survey(year, 
                    models = c("2022, base", "2022, m18.2b", "2022, m22"))
dev.off()

png(filename=here::here(year, "figs", "fcomparebio.png"), width = 6.5, height = 6.5, 
    units = "in", type ="cairo", res = 200)

plot_compare_biomass(year, 
                     models = c("2022, base", "2022, m18.2b", "2022, m22"))
dev.off()


png(filename=here::here(year, "figs", "fcomparesrv2.png"), width = 6.5, height = 6.5, 
    units = "in", type ="cairo", res = 200)
plot_compare_survey(year, 
                    models = c("2022, m22", "2022, m22.1", "2022, m22.1a", "2022, m22.1b"))
dev.off()

png(filename=here::here(year, "figs", "fcomparebio2.png"), width = 6.5, height = 6.5, 
    units = "in", type ="cairo", res = 200)

plot_compare_biomass(year, 
                     models = c("2022, m22", "2022, m22.1", "2022, m22.1a", "2022, m22.1b"))
dev.off()


