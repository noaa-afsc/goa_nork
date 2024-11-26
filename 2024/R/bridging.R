# load ----
library(RTMB)
library(tidyverse)
library(Matrix)
library(tmbstan)
library(shinystan)
library(here)
library(scico)
theme_set(afscassess::theme_report())
# devtools::install_github('Cole-Monnahan-NOAA/adnuts', ref='sparse_M')
library(adnuts)
# install.packages('StanEstimators', repos = c('https://andrjohns.r-universe.dev', 'https://cloud.r-project.org'))
library(StanEstimators)

year = 2024


source(here::here(year, 'rtmb', 'bridge_model.r'))
source(here::here(year, 'r', 'utils.r'))


# inputs 
ages = 2:45
years = 1961:2022
length_bins = 15:45

# bridge data from 2022 base model m22.1
REP = readLines(here::here(year, "base", "nr.rep"))
PAR = readLines(here::here(year, "base", "nr.par"))
DAT = readLines(here::here(year, "base", "goa_nr_2022.dat"))
CTL = readLines(here::here(year, "base", "goa_nr_2022.ctl"))

maa = as.numeric(stringr::str_split(REP[grep('Maturity', REP)], " ")[[1]][3:52])
waa = as.numeric(stringr::str_split(REP[grep('Weight', REP)], " ")[[1]][3:52])
wt_mature = waa * maa / 2
catch_obs = c(as.numeric(stringr::str_split(REP[grep('Obs_Catch', REP)], " ")[[1]][3:19]),
              as.numeric(stringr::str_split(REP[grep('Obs_Catch_Later', REP)], " ")[[1]][22:66]))

srv_yrs = as.numeric(stringr::str_split(DAT[grep('Trawl survey years', DAT)+1], " ")[[1]]) 
srv_ind = ifelse(years %in% srv_yrs, 1, 0)
srv_obs = as.numeric(stringr::str_split(DAT[grep('Trawl survey years', DAT)+3], " ")[[1]]) 
srv_sd = as.numeric(stringr::str_split(DAT[grep('Trawl survey years', DAT)+5], " ")[[1]]) 

fish_age_yrs = as.numeric(stringr::str_split(DAT[grep('Fishery Age Composition', DAT)+5], " ")[[1]]) 
as.data.frame(t(mapply(rbind, stringr::str_split(REP[grep('Obs_P_fish_age', REP)+1:length(fish_age_yrs)], " ")))) -> fish_ac
names(fish_ac) <- c('year', 'b', paste0('age-', ages), 'b1', 'eff_N', 'b2', 'N', 'b3', 'sdnr')
fish_ac %>% 
  mutate(across(c(year, starts_with('age'), N), as.numeric)) -> fish_ac

fish_age_ind = ifelse(years %in% fish_age_yrs, 1, 0)
fish_age_iss = fish_ac$N
fish_age_obs = fish_ac %>% 
  select(starts_with('age')) %>% 
  as.matrix() %>% 
  t() %>% 
  unname()

srv_age_yrs = as.numeric(stringr::str_split(DAT[grep('Trawl Survey Age Composition', DAT)+5], " ")[[1]]) 
as.data.frame(t(mapply(rbind, stringr::str_split(REP[grep('Obs_P_srv1_age', REP)+1:length(srv_age_yrs)], " ")))) -> srv_ac
names(srv_ac) <- c('year', 'b', paste0('age-', ages), 'b1', 'eff_N', 'b2', 'N', 'b3', 'sdnr')
srv_ac %>% 
  mutate(across(c(year, starts_with('age'), N), as.numeric)) -> srv_ac

srv_age_ind = ifelse(years %in% srv_age_yrs, 1, 0)
srv_age_iss = srv_ac$N
srv_age_obs = srv_ac %>% 
  select(starts_with('age')) %>% 
  as.matrix() %>% 
  t() %>% 
  unname()

fish_size_yrs = as.numeric(stringr::str_split(DAT[grep('Fishery Size Composition', DAT)+5], " ")[[1]]) 
as.data.frame(t(mapply(rbind, stringr::str_split(REP[grep('Obs_P_fish_size', REP)+1:length(fish_size_yrs)], " ")))) -> fish_sc
names(fish_sc) <- c('year', 'b', paste0('l-', length_bins), 'b1', 'eff_N', 'b2', 'N', 'b3', 'sdnr')
fish_sc %>% 
  mutate(across(c(year, starts_with('l'), N), as.numeric)) -> fish_sc

fish_size_ind = ifelse(years %in% fish_size_yrs, 1, 0)
fish_size_iss = fish_sc$N
fish_size_obs = fish_sc %>% 
  select(starts_with('l')) %>% 
  as.matrix() %>% 
  t() %>% 
  unname()

as.data.frame(t(mapply(rbind, stringr::str_split(DAT[grep('Size-age transition matrix:', DAT)+2:(1+length(maa))], " ")))) %>% 
  mutate(across(everything(), as.numeric)) %>% 
  as.matrix() -> saa

as.data.frame(t(mapply(rbind, stringr::str_split(DAT[grep('age error transition matrix:', DAT)+2:(1+length(maa))], " ")))) %>% 
  mutate(across(everything(), as.numeric)) %>% 
  as.matrix() -> ae

yield_ratio = as.numeric(stringr::str_split(CTL[grep("yield", CTL)], "\t")[[1]][1])

data = list(ages = 2:45,
            years = years,
            length_bins = 15:45,
            waa = waa,
            maa = maa,
            wt_mature = wt_mature,
            spawn_mo = 5,
            catch_obs = catch_obs,
            catch_ind = rep(1, length(catch_obs)),
            catch_wt = c(rep(5, 17), rep(50, 45)),
            srv_yrs = srv_yrs,
            srv_ind = srv_ind,
            srv_obs = srv_obs,
            srv_sd = srv_sd,
            srv_wt = 0.25,
            fish_age_yrs = fish_age_yrs,
            fish_age_ind = fish_age_ind,
            fish_age_iss = fish_age_iss,
            fish_age_obs = fish_age_obs,
            fish_age_wt = 0.5,
            srv_age_yrs = srv_age_yrs,
            srv_age_ind = srv_age_ind,
            srv_age_iss = srv_age_iss,
            srv_age_obs = srv_age_obs,
            srv_age_wt = 0.5,
            fish_size_yrs = fish_size_yrs,
            fish_size_ind = fish_size_ind,
            fish_size_iss = fish_size_iss,
            fish_size_obs = fish_size_obs,
            fish_size_wt = 0.5,
            age_error = ae,
            size_age = saa,
            wt_fmort_reg = 0.1,
            wt_rec_var = 1,
            mean_M = 0.06,
            cv_M = 0.05,
            mean_q = 1,
            cv_q = 0.45,
            mean_sigmaR = 1.5,
            cv_sigmaR = 0.01,
            yield_ratio = yield_ratio
)

pars = list(log_M = as.numeric(PAR[grep('logm', PAR)+1]),
            log_a50C = log(as.numeric(PAR[grep('\\ba50:', PAR)+1])),
            deltaC = as.numeric(PAR[grep('\\bdelta:', PAR)+1]),
            log_a50S = log(as.numeric(PAR[grep('\\ba50_srv1:', PAR)+1])),
            deltaS = as.numeric(PAR[grep('\\bdelta_srv1:', PAR)+1]),
            log_q = as.numeric(PAR[grep('\\blog_q_srv1:', PAR)+1]),
            log_mean_R = as.numeric(PAR[grep('\\blog_mean_rec:', PAR)+1]),
            init_log_Rt = rev(as.numeric(stringr::str_split(PAR[grep('\\log_rec_dev:', PAR)+1], " ")[[1]])[2:49]),
            log_Rt = as.numeric(stringr::str_split(PAR[grep('\\log_rec_dev:', PAR)+1], " ")[[1]])[50:111],
            log_mean_F = as.numeric(PAR[grep('\\blog_avg_F:', PAR)+1]),
            log_Ft = as.numeric(stringr::str_split(PAR[grep('\\blog_F_devs:', PAR)+1], " ")[[1]])[2:63],
            log_F35 = log(as.numeric(PAR[grep('\\bmF35:', PAR)+1])),
            log_F40 = log(as.numeric(PAR[grep('\\bmF40:', PAR)+1])),
            log_F50 = log(as.numeric(PAR[grep('\\bmF50:', PAR)+1])),
            sigmaR = as.numeric(PAR[grep('\\bsigr:', PAR)+1])
)

# match the base model
obj <- RTMB::MakeADFun(bridge, 
                       pars,
                       map = list(sigmaR = factor(NA)))  

base <- obj$report(obj$env$last.par.best)
proj_bio(base)

# get same values after fitting
fit <- nlminb(start = obj$par,
              objective = obj$fn,
              gradient = obj$gr,
              control = list(iter.max=100000,
                             eval.max=20000))
base <- obj$report(obj$env$last.par.best)
proj_bio(base)

saveRDS(base, here::here(year, 'base', 'base.rds'))

sdbase = sdreport(obj, getJointPrecision = TRUE)
summary(sdbase, "report") %>% 
  as.data.frame() %>% 
  mutate(lci = Estimate - 1.96*`Std. Error`,
         uci = Estimate + 1.96*`Std. Error`) %>% 
  tibble::rownames_to_column("item") %>% 
  mutate(item = gsub("\\..*", "", item),
         year = c(data$srv_yrs, rep(data$years, 3))) %>% 
  dplyr::select(year, item, value = Estimate, se = `Std. Error`, lci, uci) -> base_se

data.frame(Item = c('Catch', 
                    'Survey', 
                    "Fish age", 
                    'Survey age', 
                    'Fish size', 
                    'Recruitment', 
                    'F regularity', 
                    'SPR penalty', 
                    'M prior', 
                    'q prior', 
                    'Sub total',
                    'L maturity',
                    'C maturity',
                    'Sum maturity'), 
           RTMB = c(base$ssqcatch, 
                    base$like_srv,
                    base$like_fish_age,
                    base$like_srv_age, 
                    base$like_fish_size, 
                    base$like_rec, 
                    base$f_regularity, 
                    base$sprpen, 
                    base$nll_M, 
                    base$nll_q, 
                    round(sum(base$ssqcatch, 
                              base$like_srv, 
                              base$like_fish_age, 
                              base$like_srv_age, 
                              base$like_fish_size, 
                              base$like_rec, 
                              base$f_regularity, 
                              base$sprpen, 
                              base$nll_M, 
                              base$nll_q), 4),
                    NA, NA, NA),
           
           ADMB = c(
             stringr::str_split(REP[grep("SSQ_Catch_Likelihood" , REP) ], " ")[[1]][2],
             stringr::str_split(REP[grep("Survey_Abundance_Index_Likelihood" , REP) ], " ")[[1]][2],
             stringr::str_split(REP[grep("Fishery_Age_Composition_Likelihood" , REP) ], " ")[[1]][2],
             stringr::str_split(REP[grep("Survey_Age_Composition_Likelihood" , REP) ], " ")[[1]][2],
             stringr::str_split(REP[grep("Fishery_Size_Composition_Likelihood" , REP) ], " ")[[1]][2],
             stringr::str_split(REP[grep("Recruitment_Deviations_Likelihood" , REP) ], " ")[[1]][2],
             stringr::str_split(REP[grep("Fishing_Mortality_Regularity_Penalty" , REP) ], " ")[[1]][2],
             stringr::str_split(REP[grep("SPR penalty" , REP) ], " ")[[1]][2],
             stringr::str_split(REP[grep("priors M" , REP) ], " ")[[1]][2],
             stringr::str_split(REP[grep("priors q" , REP)], " ")[[1]][2],
             sum(as.numeric(stringr::str_split(REP[grep("SSQ_Catch_Likelihood" , REP) ], " ")[[1]][2]),
                 as.numeric(stringr::str_split(REP[grep("Survey_Abundance_Index_Likelihood" , REP) ], " ")[[1]][2]),
                 as.numeric(stringr::str_split(REP[grep("Fishery_Age_Composition_Likelihood" , REP) ], " ")[[1]][2]),
                 as.numeric(stringr::str_split(REP[grep("Survey_Age_Composition_Likelihood" , REP) ], " ")[[1]][2]),
                 as.numeric(stringr::str_split(REP[grep("Fishery_Size_Composition_Likelihood" , REP) ], " ")[[1]][2]),
                 as.numeric(stringr::str_split(REP[grep("Recruitment_Deviations_Likelihood" , REP) ], " ")[[1]][2]),
                 as.numeric(stringr::str_split(REP[grep("Fishing_Mortality_Regularity_Penalty" , REP) ], " ")[[1]][2]),
                 as.numeric(stringr::str_split(REP[grep("SPR penalty" , REP) ], " ")[[1]][2]),
                 as.numeric(stringr::str_split(REP[grep("priors M" , REP) ], " ")[[1]][2]),
                 as.numeric(stringr::str_split(REP[grep("priors q" , REP)], " ")[[1]][2])),
             
             as.numeric(stringr::str_split(REP[grep("L mat like" , REP)], " ")[[1]][2]),
             as.numeric(stringr::str_split(REP[grep("C mat like" , REP)], " ")[[1]][2]),
             sum(as.numeric(stringr::str_split(REP[grep("L mat like" , REP)], " ")[[1]][2]),
                 as.numeric(stringr::str_split(REP[grep("C mat like" , REP)], " ")[[1]][2])))#,
           #    stringr::str_split(REP[grep("obj_fun" , REP) ], " ")[[1]][2])
           
) %>% 
  mutate(ADMB = as.numeric(ADMB),
         difference = as.numeric(ADMB) - as.numeric(RTMB)) %>% 
  mutate(across(2:4, round, 4))  %>% 
  saveRDS(., file = here::here(year, 'data', 'user_input', "rtmb-like.RDS"))


pr = proj_bio(base)

data.frame(Item = c('M', 
                    'q', 
                    "Log mean recruitment", 
                    'Log mean F', 
                    'A50 fishery', 
                    'Delta fishery', 
                    'A50 survey', 
                    'Delta survey', 
                    '2023 Total biomass', 
                    '2023 Spawning biomass', 
                    '2023 OFL',
                    '2023 FOFL',
                    '2023 ABC',
                    '2023 FABC'), 
           RTMB = c(base$M, 
                    base$q,
                    base$log_mean_R,
                    base$log_mean_F, 
                    base$a50C, 
                    base$deltaC, 
                    base$a50S, 
                    base$deltaS, 
                    pr[1,3], 
                    pr[1,2], 
                    pr[1,5],
                    pr[1,7],
                    pr[1,4],
                    pr[1,6]),
           
           ADMB = c(
             stringr::str_split(REP[grep('nat_mort' , REP)+1], " ")[[1]][1],
             stringr::str_split(REP[grep('q_trawl' , REP)+1], " ")[[1]][1],
             stringr::str_split(REP[grep('log_mean_rec' , REP)+1], " ")[[1]][1],
             stringr::str_split(PAR[grep('log_avg_F' , PAR)+1], " ")[[1]][1],
             stringr::str_split(PAR[grep('a50:' , PAR)+1], " ")[[1]][1],
             stringr::str_split(PAR[grep('delta:' , PAR)+1], " ")[[1]][1],
             stringr::str_split(PAR[grep('a50_srv1:' , PAR)+1], " ")[[1]][1],
             stringr::str_split(PAR[grep('delta_srv1:' , PAR)+1], " ")[[1]][1],
             stringr::str_split(REP[grep('tot_biom' , REP)+1], " ")[[1]][1],
             stringr::str_split(REP[grep('spawn_biom' , REP)+1], " ")[[1]][1],
             stringr::str_split(REP[grep('\\bOFL for 2023' , REP)+1], " ")[[1]][1],
             stringr::str_split(REP[grep('F_OFL for 2023' , REP)+1], " ")[[1]][1],
             stringr::str_split(REP[grep('\\bABC for 2023' , REP)+1], " ")[[1]][1],
             stringr::str_split(REP[grep('F_ABC for 2023' , REP)+1], " ")[[1]][1])) %>% 
  mutate(ADMB = as.numeric(ADMB),
         difference = as.numeric(ADMB) - as.numeric(RTMB)) %>% 
  mutate(across(2:4, round, 4)) %>% 
  saveRDS(., file = here::here(year, 'data', 'user_input', "rtmb-compare.RDS"))


# comparison plots ----
png(filename=here::here(year, "compare_models", "compare-slx.png"), width = 6.5, height = 5.5,
    units = "in", type ="cairo", res = 200)
read.csv(here::here(year, 'base', 'processed', 'selex.csv')) %>% 
  mutate(model = 'ADMB') %>% 
  bind_rows(data.frame(age = 2:51,
            fish = base$slx[,1],
            srv1 = base$slx[,2],
            maturity = base$maa,
            model = 'RTMB')) %>% 
  tidyr::pivot_longer(-c(age, model)) %>% 
  tidyr::pivot_wider(names_from=model, values_from = value) %>% 
  dplyr::mutate(pd = (RTMB - ADMB) / RTMB) %>% 
  ggplot(aes(age, pd, color = name)) + 
  geom_point() +
  scale_y_continuous(labels=scales::percent) +
  scico::scale_color_scico_d("item", palette='roma') +
  xlab('Age') +
  ylab('Percent difference')
dev.off()


png(filename=here::here(year, "compare_models", "compare-bio.png"), width = 6.5, height = 6.5,
    units = "in", type ="cairo", res = 200)
read.csv(here::here(year, 'base', 'processed', 'bio_rec_f.csv')) %>% 
  mutate(model = 'ADMB') %>% 
  bind_rows(data.frame(year = base$years,
                       tot_biom = base$tot_bio,
                       sp_biom = base$spawn_bio,
                       F = base$Ft,
                       recruits = base$recruits,
                       model = 'RTMB')) %>% 
  tidyr::pivot_longer(-c(year, model)) %>% 
  tidyr::pivot_wider(names_from=model, values_from = value) %>% 
  dplyr::mutate(pd = (RTMB - ADMB) / RTMB) %>% 
  ggplot(aes(year, pd, color = name)) + 
  geom_point() +
  scale_y_continuous(labels=scales::percent) +
  facet_wrap(~name) +
  scico::scale_color_scico_d("item", palette='roma') +
  xlab('Year') +
  ylab('Percent difference')
dev.off()

