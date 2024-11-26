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


source(here::here(year, 'r', 'models.r'))
source(here::here(year, 'r', 'utils.r'))


# inputs 
ages = 2:45
years = 1961:2022
length_bins = 15:45


# base model with 2024 data 

fishery = grep("fish", list.files("C:/Users/Ben.Williams/Work/assessments/goa_nork/2024/data/output"), value=TRUE)
survey = grep("bts_", list.files("C:/Users/Ben.Williams/Work/assessments/goa_nork/2024/data/output"), value=TRUE)
yld = read.csv("C:/Users/Ben.Williams/Work/assessments/goa_nork/2024/data/output/yld_rat.csv")
catch = read.csv(paste0("C:/Users/Ben.Williams/Work/assessments/goa_nork/2024/data/output/", grep("catch", fishery, value=TRUE)))
wta = read.csv("C:/Users/Ben.Williams/Work/assessments/goa_nork/2024/data/output/waa.csv")
names(wta) <- c('age', 'wt')
saa = as.matrix(read.csv("C:/Users/Ben.Williams/Work/assessments/goa_nork/2024/data/output/saa.csv"))
ae = as.matrix(read.csv("C:/Users/Ben.Williams/Work/assessments/goa_nork/2024/data/output/ae_model.csv"))
fishac = read.csv(paste0("C:/Users/Ben.Williams/Work/assessments/goa_nork/2024/data/output/", grep("age", fishery, value=TRUE))) %>% 
  filter(year>=1990)
fishlc = read.csv(paste0("C:/Users/Ben.Williams/Work/assessments/goa_nork/2024/data/output/", grep("length", fishery, value=TRUE))) %>% 
  filter(year>=1991, year!=2024)
tsac = read.csv(paste0("C:/Users/Ben.Williams/Work/assessments/goa_nork/2024/data/output/", grep("age", survey, value=TRUE))) %>% 
  filter(year>=1990)
# tsb = read.csv("C:/Users/Ben.Williams/Work/assessments/goa_nork/2024/data/user_input/vast_lognormal.csv"))
tsb = read.csv("C:/Users/Ben.Williams/Work/assessments/goa_nork/2024/data/user_input/vast_2_1.csv")
REP = readLines("C:/Users/Ben.Williams/Work/assessments/goa_nork/2024/m22.1/nr.rep")
maa = as.numeric(stringr::str_split(REP[grep('Maturity', REP)], " ")[[1]][3:52])
wt_mature = wta$wt * maa / 2

data = list(ages = 2:45,
            years = catch$year,
            length_bins = 15:45,
            waa = wta$wt,
            maa = maa,
            wt_mature = wt_mature,
            spawn_mo = 5,
            catch_obs = catch$catch,
            catch_ind = rep(1, length(catch$year)),
            catch_wt = c(rep(5, 17), rep(50, 47)),
            srv_yrs = tsb$year,
            srv_ind = ifelse(catch$year %in% tsb$year, 1, 0),
            srv_obs = tsb$biomass,
            srv_sd = tsb$se,
            srv_wt = 0.25,
            fish_age_yrs = fishac$year,
            fish_age_ind = ifelse(catch$year %in% fishac$year, 1, 0),
            fish_age_iss = c(25.2467,	39.2224,	54.5595,	44.3967,	55.2952,	74.6048,	48.9655,	53.8815,	68.6773,	69.8816,	87.4294,	100,	86.6187,	84.3934,	62.816,	76.1955),
            fish_age_obs = unname(t(as.matrix(fishac[,-(1:4)]))),
            fish_age_wt = 0.5,
            srv_age_yrs = tsac$year,
            srv_age_ind = ifelse(catch$year %in% tsac$year, 1, 0),
            # srv_age_iss = c(49, 33, 91, 43, 86, 25, 61, 127, 171, 94, 101, 117, 123, 99, 155, 50),
            srv_age_iss = c(28.2956716,28.79701882,42.0642524,38.89731882,89.35478366,30.94949752,77.79463763,100,95.15473983,80.08763593,82.37064037,72.44958071,86.31402473,68.90159616,83.77314565, 74.3933),
            srv_age_obs = unname(t(as.matrix(tsac[,-(1:4)]))),
            srv_age_wt = 0.5,
            fish_size_yrs = fishlc$year ,
            fish_size_ind = ifelse(catch$year %in% fishlc$year, 1, 0),
            fish_size_iss = c(75.2673,69.1482,50.3742,42.9224,59.2575,46.4495,31.0813,75.2867,93.4583,83.2692,72.05,89.8439,100,57.9263,77.1133,64.5234,54.1795),
            fish_size_obs = unname(t(as.matrix(fishlc[,-c(1:4)]))),
            fish_size_wt = 0.5,
            age_error = unname(ae),
            size_age = unname(as.matrix(saa[,-1])),
            wt_fmort_reg = 0.1,
            wt_rec_var = 1,
            mean_M = 0.06,
            cv_M = 0.05,
            mean_q = 1,
            cv_q = 0.45,
            mean_sigmaR = 1.5,
            cv_sigmaR = 0.01,
            yield_ratio = yld$yld
)

# bridge model with estimation
pars = list(log_M = log(0.06),
            log_a50C = log(7.5),
            deltaC = 3.0,
            log_a50S = log(7.3),
            deltaS = 3.8,
            log_q = log(1),
            log_mean_R = 4.3,
            init_log_Rt =rep(0, nrow(data$age_error)-2),
            log_Rt = rep(0, length(data$years)),
            log_mean_F = 0,
            log_Ft =  rep(0, length(data$years)),
            log_F35 = 0,
            log_F40 = 0,
            log_F50 = 0,
            sigmaR = 1.5)

# m22.1 - base model in 2024 ----
obj22.1 <- RTMB::MakeADFun(bridge, 
                           pars,
                           map = list(sigmaR = factor(NA)))  

fit22.1 <- nlminb(start = obj22.1$par,
                  objective = obj22.1$fn,
                  gradient = obj22.1$gr,
                  control = list(iter.max=100000,
                                 eval.max=20000))
m22.1 <- obj22.1$report(obj22.1$env$last.par.best)
proj_bio(m22.1)
sd22.1 = sdreport(obj22.1)
summary(sd22.1, "report") %>% 
  as.data.frame() %>% 
  mutate(lci = Estimate - 1.96*`Std. Error`,
         uci = Estimate + 1.96*`Std. Error`) %>% 
  tibble::rownames_to_column("item") %>% 
  mutate(item = gsub("\\..*", "", item),
         year = c(data$srv_yrs, rep(data$years, 3))) %>% 
  dplyr::select(year, item, value = Estimate, se = `Std. Error`, lci, uci) -> m22.1_se
saveRDS(m22.1, here::here(year, 'm22.1', 'm22.1.rds'))


# m22.1 - rewt - base model iwith survey biomass weighted at 1----
data$srv_wt = 1
obj22.1_rewt <- RTMB::MakeADFun(bridge, 
                           pars,
                           map = list(sigmaR = factor(NA)))  

fit22.1_rewt <- nlminb(start = obj22.1_rewt$par,
                  objective = obj22.1_rewt$fn,
                  gradient = obj22.1_rewt$gr,
                  control = list(iter.max=100000,
                                 eval.max=20000))
m22.1_rewt <- obj22.1_rewt$report(obj22.1_rewt$env$last.par.best)
proj_bio(m22.1_rewt)
saveRDS(m22.1_rewt, here::here(year, 'm22.1', 'm22.1_rewt.rds'))

# 22.1a - lognormal----
REP = readLines("C:/Users/Ben.Williams/Work/assessments/goa_nork/2024/m22.1b/nr.rep")
maa = as.numeric(stringr::str_split(REP[grep('Maturity', REP)], " ")[[1]][3:52])
data$wt_mature = wta$wt * maa / 2
data$maa = maa
data$waa = wta$wt
data$srv_wt = 0.25
obj22.1a <- RTMB::MakeADFun(f, 
                            pars,
                            map = list(sigmaR = factor(NA)))  

fit22.1a <- nlminb(start = obj22.1a$par,
                   objective = obj22.1a$fn,
                   gradient = obj22.1a$gr,
                   control = list(iter.max=100000,
                                  eval.max=20000))
m22.1a <- obj22.1a$report(obj22.1a$env$last.par.best)
proj_bio(m22.1a)
saveRDS(m22.1a, here::here(year, 'm22.1a', 'm22.1a.rds'))
sd22.1a = sdreport(obj22.1a)
summary(sd22.1a, "report") %>% 
  as.data.frame() %>% 
  mutate(lci = Estimate - 1.96*`Std. Error`,
         uci = Estimate + 1.96*`Std. Error`) %>% 
  tibble::rownames_to_column("item") %>% 
  mutate(item = gsub("\\..*", "", item),
         year = c(data$years, data$srv_yrs, rep(data$years, 3)),
         lci = ifelse(item=='recruits' & lci<0, 0, lci)) %>% 
  dplyr::select(year, item, value = Estimate, se = `Std. Error`, lci, uci) -> m22.1a_se

m22.1a_se %>% 
  mutate(item = case_when(item=='log_Rt' ~ 'log recruitment devs',
                          item=='srv_pred' ~ 'predicted survey biomass',
                          item=='recruits' ~ 'recruitment',
                          item=='spawn_bio' ~ 'spawning biomass',
                          item=='tot_bio' ~ 'total_biomass')) %>% 
  ggplot(aes(year, value)) + 
  geom_errorbar(aes(ymin=lci, ymax=uci), width=0.2) +
  geom_line() +
  geom_point() +
  facet_wrap(~item, scales = 'free') +
  expand_limits(y=0) +
  scale_y_continuous(labels = scales::comma) +
  xlab('Year') +
  ylab("")

# m22.1a -rewt
data$srv_wt = 1
obj22.1a_rewt <- RTMB::MakeADFun(f, 
                            pars,
                            map = list(sigmaR = factor(NA)))  

fit22.1a_rewt <- nlminb(start = obj22.1a_rewt$par,
                   objective = obj22.1a_rewt$fn,
                   gradient = obj22.1a_rewt$gr,
                   control = list(iter.max=100000,
                                  eval.max=20000))
m22.1a_rewt <- obj22.1a_rewt$report(obj22.1a_rewt$env$last.par.best)
proj_bio(m22.1a_rewt)
saveRDS(m22.1a_rewt, here::here(year, 'm22.1a', 'm22.1a_rewt.rds'))

# m22.1b ISS---- 
REP = readLines("C:/Users/Ben.Williams/Work/assessments/goa_nork/2024/m22.1b/nr.rep")
maa = as.numeric(stringr::str_split(REP[grep('Maturity', REP)], " ")[[1]][3:52])
data$wt_mature = wta$wt * maa / 2
data$maa = maa
data$srv_age_iss = c(49, 33, 91, 43, 86, 25, 61, 127, 171, 94, 101, 117, 123, 99, 155, 50)
data$srv_wt = 0.25
obj22.1b <- RTMB::MakeADFun(f, 
                            pars, 
                            map = list(sigmaR = factor(NA)))  
fit22.1b <- nlminb(start = obj22.1b$par,
                   objective = obj22.1b$fn,
                   gradient = obj22.1b$gr,
                   control = list(iter.max=100000,
                                  eval.max=20000))
m22.1b <- obj22.1b$report(obj22.1b$env$last.par.best)
proj_bio(m22.1b)
saveRDS(m22.1b, here::here(year, 'm22.1b', 'm22.1b.rds'))
sd22.1b = sdreport(obj22.1b)
summary(sd22.1b, "report") %>% 
  as.data.frame() %>% 
  mutate(lci = Estimate - 1.96*`Std. Error`,
         uci = Estimate + 1.96*`Std. Error`) %>% 
  tibble::rownames_to_column("item") %>% 
  mutate(item = gsub("\\..*", "", item),
         year = c(data$years, data$srv_yrs, rep(data$years, 3))) %>% 
  dplyr::select(year, item, value = Estimate, se = `Std. Error`, lci, uci) -> m22.1b_se


# m22.1b -rewtISS---- 
data$srv_wt = 1
obj22.1b_rewt <- RTMB::MakeADFun(f, 
                            pars, 
                            map = list(sigmaR = factor(NA)))  
fit22.1b_rewt <- nlminb(start = obj22.1b_rewt$par,
                   objective = obj22.1b_rewt$fn,
                   gradient = obj22.1b_rewt$gr,
                   control = list(iter.max=100000,
                                  eval.max=20000))
m22.1b_rewt <- obj22.1b_rewt$report(obj22.1b_rewt$env$last.par.best)
proj_bio(m22.1b_rewt)
saveRDS(m22.1b_rewt, here::here(year, 'm22.1b', 'm22.1b_rewt.rds'))

# m22.1c  VAST  lognormal ----
REP = readLines("C:/Users/Ben.Williams/Work/assessments/goa_nork/2024/m22.1c/nr.rep")
maa = as.numeric(stringr::str_split(REP[grep('Maturity', REP)], " ")[[1]][3:52])
data$wt_mature = wta$wt * maa / 2
data$maa = maa
tsb = read.csv("C:/Users/Ben.Williams/Work/assessments/goa_nork/2024/data/user_input/vast_lognormal.csv")
data$srv_obs = tsb$biomass
data$srv_sd = tsb$se
data$srv_wt = 0.25
obj22.1c <- RTMB::MakeADFun(f, 
                            pars, 
                            map = list(sigmaR = factor(NA)))  
fit22.1c <- nlminb(start = obj22.1c$par,
                   objective = obj22.1c$fn,
                   gradient = obj22.1c$gr,
                   control = list(iter.max=100000,
                                  eval.max=20000))
m22.1c <- obj22.1c$report(obj22.1c$env$last.par.best)
proj_bio(m22.1c)
sdm22.1c = sdreport(obj22.1c)
summary(sdm22.1c, "report") %>% 
  as.data.frame() %>% 
  mutate(lci = Estimate - 1.96*`Std. Error`,
         uci = Estimate + 1.96*`Std. Error`) %>% 
  tibble::rownames_to_column("item") %>% 
  mutate(item = gsub("\\..*", "", item),
         year = c(data$years, data$srv_yrs, rep(data$years, 3))) %>% 
  dplyr::select(year, item, value = Estimate, se = `Std. Error`, lci, uci) -> m22.1c_se

saveRDS(m22.1c, here::here(year, 'm22.1c', 'm22.1c.rds'))

# m22.1c-rewt  VAST  lognormal ----
data$srv_wt = 1
obj22.1c_rewt <- RTMB::MakeADFun(f, 
                            pars, 
                            map = list(sigmaR = factor(NA)))  
fit22.1c_rewt <- nlminb(start = obj22.1c_rewt$par,
                   objective = obj22.1c_rewt$fn,
                   gradient = obj22.1c_rewt$gr,
                   control = list(iter.max=100000,
                                  eval.max=20000))
m22.1c_rewt <- obj22.1c_rewt$report(obj22.1c_rewt$env$last.par.best)
proj_bio(m22.1c_rewt)
saveRDS(m22.1c_rewt, here::here(year, 'm22.1c', 'm22.1c_rewt.rds'))

# m24 ----
# update maturity
REP = readLines("C:/Users/Ben.Williams/Work/assessments/goa_nork/2024/m22.1d/nr.rep")
maa = as.numeric(stringr::str_split(REP[grep('Maturity', REP)], " ")[[1]][3:52])
data$wt_mature = wta$wt * maa / 2
data$maa = maa
data$srv_wt = 0.25
obj24 <- RTMB::MakeADFun(f, 
                         pars, 
                         map = list(sigmaR = factor(NA)))  
fit24 <- nlminb(start = obj24$par,
                objective = obj24$fn,
                gradient = obj24$gr,
                control = list(iter.max=100000,
                               eval.max=20000))

m24 <- obj24$report(obj24$env$last.par.best)
proj_bio(m24)
sd24 = sdreport(obj24, getJointPrecision = TRUE)
summary(sd24, "report") %>% 
  as.data.frame() %>% 
  mutate(lci = Estimate - 1.96*`Std. Error`,
         uci = Estimate + 1.96*`Std. Error`) %>% 
  tibble::rownames_to_column("item") %>% 
  filter(!(item %in% c('q', "M"))) %>% 
  mutate(item = gsub("\\..*", "", item),
         year = c(data$years, data$srv_yrs, rep(data$years, 3))) %>% 
  dplyr::select(year, item, value = Estimate, se = `Std. Error`, lci, uci) -> m24_se


# m24-rewt ----
# update maturity
data$srv_wt = 1.0
obj24_rwt <- RTMB::MakeADFun(f, 
                         pars, 
                         map = list(sigmaR = factor(NA)))  
fit24_rwt <- nlminb(start = obj24_rwt$par,
                objective = obj24_rwt$fn,
                gradient = obj24_rwt$gr,
                control = list(iter.max=100000,
                               eval.max=20000))

m24_rewt <- obj24_rwt$report(obj24_rwt$env$last.par.best)
proj_bio(m24_rewt)
saveRDS(m24_rewt, here::here(year, 'm24', 'm24_rewt.rds'))


plot(m22.1$years, m22.1$spawn_bio) 
lines(m22.1_rewt$years, m22.1_rewt$spawn_bio)
lines(m24_rewt$years, m24_rewt$spawn_bio)
# afscassess::proj_ak(year = 2024, last_full_assess=2022, alt=NULL, folder='m24', species='nork', region='goa', rec_age = 2,
#         off_yr = FALSE, run_name = "Standard",
#         tac_abc = 1, srr = 1, rec_proj = 1, srr_cond = 0, srr_prior = 0,
#         nyrs_proj = 14, nsims = 1000, n_species = 1, abc_mult = 1, pop_scalar = 1000,
#         tac_categories = 1, tac_models = 1, admb_home = NULL)

Q <- sd24$jointPrecision
## if no random effectrs this will be NULL
if(!is.null(Q)){
  M <- solve(Q)
} else {
  M <- sd24$cov.fixed
  Q <- solve(M)
}   
globals <- list(data=data) # NULL
chains <- 5
# Cole recommends 5 chains, 1000 iters, 250 warmup
mcmc <- sample_sparse_tmb(obj24, iter=2000, warmup=400,
                          init='random', chains=chains,
                          cores=chains, metric='dense',
                          Qinv=M, Q=Q,
                          globals=globals, skip_optimization=TRUE)


saveRDS(obj24, here::here(2024, 'm24', 'obj.rds') )
saveRDS(data, here::here(2024, 'm24', 'dat.rds'))
saveRDS(pars, here::here(2024, 'm24', 'pars.rds'))
saveRDS(f, here::here(2024, 'm24', 'model.rds'))
saveRDS(m24_se, file = here::here(2024, 'm24','m24_se.rds'))
saveRDS(mcmc, file = here::here(2024, 'm24','mcmc.rds'))



data = readRDS(here::here(2024, 'm24', 'dat.rds'))
pars = readRDS(here::here(2024, 'm24', 'pars.rds'))
f = readRDS(here::here(2024, 'm24', 'model.rds'))
m24_se = readRDS(file = here::here(2024, 'm24','m24_se.rds'))
mcmc = readRDS(file = here::here(2024, 'm24','mcmcmc.rds'))


summary(mcmc) ## gives you minESS and maxRhat
summary(mcmc$monitor$n_eff)
summary(mcmc$monitor$Rhat)
post <- (adnuts::extract_samples(mcmc))
mpost = as.matrix(post)
postlist <- coda::mcmc.list(lapply(post, coda::mcmc))
coda::traceplot(postlist)
glimpse(post)
sp <- extract_sampler_params(mcmc)
glimpse(sp)
plot_marginals(mcmc, pars=1:6)
## diagnose issues
pairs_admb(mcmc, pars=1:5, order='slow')
pairs_admb(mcmc, pars=1:5, order='mismatch')

# process mcmc - takes a while
mcmcout = list()
for(i in 1:nrow(mpost)){
  mcmcout[[i]] <- obj24$report(mpost[i,])
}

saveRDS(mcmcout, here::here(year, folder, 'mcmcout.RDS'))
mcmcout = readRDS(file = here::here(2024, 'm24','mcmcout.rds'))
purrr::map(mcmcout, 'spawn_bio') %>% 
  purrr::map(., ~ as.data.frame(t(.))) %>% 
  bind_rows(.id = 'sim') -> sp 

purrr::map(mcmcout, 'tot_bio') %>% 
  purrr::map(., ~ as.data.frame(t(.))) %>% 
  bind_rows(.id = 'sim') -> tp 

purrr::map(mcmcout, 'recruits') %>% 
  purrr::map(., ~ as.data.frame(t(.))) %>% 
  bind_rows(.id = 'sim') -> recs 

sp %>% 
  purrr::discard(~is.character(.x)) %>% 
  {
    list(mean = purrr::map(., mean),
         q1 = purrr::map(., ~ quantile(.x, probs = 0.025)),
         q2 = purrr::map(., ~ quantile(.x, probs = 0.975)))
  } %>% 
  bind_rows(.id = 'item') %>% 
  pivot_longer(-item) %>% 
  mutate(year = as.numeric(gsub("V", "", name)) + 1960) %>% 
  filter(year>=1977) %>% 
  pivot_wider(names_from=item, values_from = value) %>% 
  dplyr::select(year, mean_sp=mean, q1_sp=q1, q2_sp=q2) -> sp1


tp %>% 
  purrr::discard(~is.character(.x)) %>% 
  {
    list(mean = purrr::map(., mean),
         q1 = purrr::map(., ~ quantile(.x, probs = 0.025)),
         q2 = purrr::map(., ~ quantile(.x, probs = 0.975)))
  } %>% 
  bind_rows(.id = 'item') %>% 
  pivot_longer(-item) %>% 
  mutate(year = as.numeric(gsub("V", "", name)) + 1960) %>% 
  filter(year>=1977) %>% 
  pivot_wider(names_from=item, values_from = value) %>% 
  dplyr::select(year, mean_tp=mean, q1_tp=q1, q2_tp=q2) -> tp1
  
recs %>% 
  purrr::discard(~is.character(.x)) %>% 
  {
    list(mean = purrr::map(., mean),
         q1 = purrr::map(., ~ quantile(.x, probs = 0.025)),
         q2 = purrr::map(., ~ quantile(.x, probs = 0.975)))
  } %>% 
  bind_rows(.id = 'item') %>% 
  pivot_longer(-item) %>% 
  mutate(year = as.numeric(gsub("V", "", name)) + 1960,
         value = value * 1000) %>% 
  filter(year>=1977) %>% 
  pivot_wider(names_from=item, values_from = value) %>% 
  dplyr::select(-name) -> r1


left_join(r1, tp1) %>% 
  left_join(sp1) %>% 
  saveRDS(., here::here(year, 'data', 'output', 'mcmc_tbl.rds'))






head(post) %>% 
  








# retro rtmb
# setup
if (!dir.exists(here::here(year, folder, "retro"))){
  dir.create(here::here(year, folder, "retro"), recursive=TRUE)
}

data = d0 = readRDS(here::here(year, folder, 'dat.RDS'))
pars = p0 = readRDS(here::here(year, folder, 'pars.RDS'))
lwr = readRDS(here::here(year, folder, 'lower.RDS'))
upr = readRDS(here::here(year, folder, 'upper.RDS'))
model = readRDS(here::here(year, folder, 'model.RDS'))

sdreport(obj)-> sdm24
summary(sdm24) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("item") %>% 
  filter(item == 'log_q') -> qsd

summary(sdm24) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("item") %>% 
  filter(item == 'log_M') -> msd

summary(sdm24) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("item") %>% 
  filter(item == 'F40') -> f40


purrr::map(mcmcout, "q") %>% 
  unlist() -> qs

purrr::map(mcmcout, "M") %>% 
  unlist() -> ms

  data.frame(name = 'q',
             value = m24$q,
             mean = mean(qs),
             median = median(qs),
             std.dev = qsd$`Std. Error`,
             sd = sd(qs),
             lci = quantile(qs, 0.025),
             lci = quantile(qs, 0.975)
             )
  read.csv(here::here(2022, "m22.1", "tables", "mcmc_pars.csv"))

reps = list()
N = sum(catch_ind)

for(i in 1:10) {
  
  data$years = head(d0$years, -i)
  data$catch_ind = head(d0$catch_ind, -i)
  n = sum(data$catch_ind)
  
  data$catch_obs = d0$catch_obs[1:n]
  data$catch_wt = d0$catch_wt[1:n]
  
  data$srv_ind[n:N] <- 0
  data$fish_age_ind[(n-1):N] <- 0
  data$srv_age_ind[(n-1):N] <- 0 
  data$fish_size_ind[(n-1):N] <- 0 
  data = lapply(data, unname)
  
  pars$log_Ft = p0$log_Ft[1:n]
  pars$log_Rt = p0$log_Rt[1:n]
  pars = lapply(pars, unname)
  
  obj = RTMB::MakeADFun(model, 
                        pars, 
                        map = list(sigmaR = factor(NA))) 
  fit = nlminb(start = obj$par,
               objective = obj$fn,
               gradient = obj$gr,
               control = list(iter.max=100000,
                              eval.max=20000))
  out = obj$report(obj$env$last.par.best)
  sds = sdreport(obj)
  reps[[paste0('out',i)]] <- out
  reps[[paste0('sd',i)]] <- sds
}


saveRDS(reps, here::here(year, folder, 'retro', 'reps.rds'))








purrr::map(mcmcout, 'spawn_bio') %>% 
  purrr::map(., ~ as.data.frame(t(.))) %>% 
  bind_rows(.id = 'sim') -> sp 

sp[,ncol(sp)] -> sp_hist

purrr::map(r, 'tot_bio') %>% 
  purrr::map(., ~ as.data.frame(t(.))) %>% 
  bind_rows(.id = 'sim') -> tp 

tp[,ncol(tp)] -> tp_hist


plot_par <- function(item, post=NULL, report, rep_item) {
  
  if(!is.null(post)){
  post[[item]] %>% 
    as.data.frame() %>% 
    ggplot(aes(.)) +
    geom_histogram(ggplot2::aes(y=(after_stat(count))/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                   fill = "lightgray", color = "darkgray", bins = 50) -> p1
  
  max_y <- ggplot_build(p1)$data[[1]] %>%
    dplyr::summarize(max_y = max(y)) %>%
    pull(max_y)
  
  p1 + ylab("") +
    annotate("segment", x = log(report[[rep_item]]), xend = log(report[[rep_item]]), y = 0, yend = max_y,
             linetype = 1, linewidth = 1.25)
  } else {
    item %>% 
      as.data.frame() %>% 
      ggplot(aes(.)) +
      geom_histogram(ggplot2::aes(y=(after_stat(count))/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                     fill = "lightgray", color = "darkgray", bins = 50) -> p1
    
    max_y <- ggplot_build(p1)$data[[1]] %>%
      dplyr::summarize(max_y = max(y)) %>%
      pull(max_y)
    
    report[[rep_item]][max(length(report[[rep_item]]))] -> ri
    
    
    p1 + ylab("") +
      annotate("segment", x = ri, xend = ri, y = 0, yend = max_y,
               linetype = 1, linewidth = 1.25)
  }
  
  }

plot_par(item = "log_q", post=post, report=m24, rep_item='q') +
  xlab('Trawl survey catchability')
plot_par(item = "log_M", post=post, report=m24, rep_item='M') +
  xlab('Trawl survey catchability')
plot_par(item = "log_F40", post=post, report=m24, rep_item='F40') +
  xlab(expression(italic(F)["40%"]))
plot_par(item=sp_hist, report=m24, rep_item='spawn_bio') 
plot_par(item=tp_hist, report=m24, rep_item='tot_bio') 
plot_par(item=tp_hist, report=m24, rep_item='tot_bio') +
  scale_x_continuous(labels=scales::comma) +
  xlab('Total biomass (t)')


hist_dat = pp %>% 
ggplot(aes(.)) +
  geom_histogram(ggplot2::aes(y=(after_stat(count))/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                              fill = "lightgray", color = "darkgray", bins = 50) 


max_y <- ggplot_build(hist_dat)$data[[1]] %>%
  dplyr::summarize(max_y = max(y)) %>%
  pull(max_y)

hist_dat +
  annotate("segment", x = log(m24$M), xend = log(m24$M), y = 0, yend = max_y,
           linewidth = 1.25)


?proj_bio
  
tt = list()
for(i in 1:800){
  tt[[i]] = proj_bio(obj$report(mpost[i,]))
}

tt %>% 
  bind_rows(.id = 'sim') %>% 
  filter(year==2025) %>% 
  summarise(sb_mean = mean(spawn_bio),
            sb_median = median(spawn_bio),
            sb_sd = sd(spawn_bio),
            sb_lci = quantile(spawn_bio, .025),
            sb_uci = quantile(spawn_bio, .975),
            f40_mean = mean(F40),
            f40_median = median(F40),
            f40_sd = sd(F40),
            f40_lci = quantile(F40, .025),
            f40_uci = quantile(F40, .975),
            abc_mean = mean(catch_abc),
            abc_median = median(catch_abc),
            abc_sd = sd(catch_abc),
            abc_lci = quantile(catch_abc, .025),
            abc_uci = quantile(catch_abc, .975)) -> pjs


  data.frame(name = c('F40', "spawn_bio", "ABC"),
             value = c(proj_bio(m24)$F40[1], proj_bio(m24)$spawn_bio[1], proj_bio(m24)$catch_abc[1]),
              mean = c(pjs$f40_mean, pjs$sb_mean, pjs$abc_mean),
             median = c(pjs$f40_median, pjs$sb_median, pjs$abc_median),
             sd = c(pjs$f40_sd, pjs$sb_sd, pjs$abc_sd),
             lci = c(pjs$f40_lci, pjs$sb_lci, pjs$abc_lci),
             uci = c(pjs$f40_uci, pjs$sb_uci, pjs$abc_uci)  ) %>% 
    bind_rows(  summary(sd24, "report") %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column('name') %>% 
    filter(name %in% c('M', 'q')) %>% 
    select(name, value = Estimate) %>% 
    mutate(mean = c(mean(exp(post$log_M)), mean(exp(post$log_q))),
           median = c(median(exp(post$log_M)), median(exp(post$log_q))),
           sd = c(sd(exp(post$log_M)), sd(exp(post$log_q))),
           lci = c(quantile(exp(post$log_M), .025), quantile(exp(post$log_q), .025)),
           uci = c(quantile(exp(post$log_M), .975), quantile(exp(post$log_q), .975)))
    ) %>% 
    saveRDS(here::here(year, 'm24', 'mcmcpars.rds'))
 
  
  pivot_longer()


# establish quantiles
q_name <- tidytable::map_chr(seq(.025,.975,.05), ~ paste0("q", .x*100))
q_fun <- tidytable::map(seq(.025,.975,.05), ~ purrr::partial(quantile, probs = .x, na.rm = TRUE)) %>%
  purrr::set_names(nm = q_name)
yr = year



  
sp %>% 
  dplyr::rename_with(~ as.character(data$years), -sim) %>% 
  tidyr::pivot_longer(-sim) %>% 
  mutate(year = as.numeric(name),
         sim = as.numeric(sim) ) %>%
  bind_rows(filter(tt, id=='spawn_bio') %>% 
              mutate(id = 'projected')) %>% 
  group_by(year) %>% 
  dplyr::summarise_at(dplyr::vars(value), tibble::lst(!!!q_fun, median)) %>% 
  tidytable::pivot_longer(-c(year, median)) %>%
  tidytable::mutate(grouping = tidytable::case_when(name == q_name[1] | name == q_name[20] ~ 1,
                                                    name == q_name[2] | name == q_name[19] ~ 2,
                                                    name == q_name[3] | name == q_name[18] ~ 3,
                                                    name == q_name[4] | name == q_name[17] ~ 4,
                                                    name == q_name[5] | name == q_name[16] ~ 5,
                                                    name == q_name[6] | name == q_name[15] ~ 6,
                                                    name == q_name[7] | name == q_name[14] ~ 7,
                                                    name == q_name[8] | name == q_name[13] ~ 8,
                                                    name == q_name[9] | name == q_name[12] ~ 9,
                                                    name == q_name[10] | name == q_name[11] ~ 10)) %>%
  tidytable::mutate(min = min(value),
                    max = max(value),
                    fill = as.character(ifelse(year>yr, 0, 1)),
                    .by = c(year, grouping)) -> dat


dat %>% 
  ggplot2::ggplot(ggplot2::aes(year, group = grouping)) +
  ggplot2::geom_ribbon(data = dplyr::filter(dat, year<yr),
                       ggplot2::aes(ymin = min, ymax = max), alpha = 0.07) +
  ggplot2::geom_line(ggplot2::aes(y = median)) +
  ggplot2::guides(linetype = "none", fill = "none") +
  ggplot2::geom_ribbon(data = dat, ggplot2::aes(ymin = min, ymax = max), alpha = 0.07) +
  afscassess::scale_x_tickr(data=dat, var=year, to=10, start = 1960) +
  ggplot2::ylab("Spawning biomass (kt)") +
  ggplot2::xlab("Year") +
  ggplot2::expand_limits(y = 0) +
  scale_y_continuous(labels = scales::comma) +
  ggplot2::geom_hline(yintercept = c(m24$B40, m24$B35),
                      lty = c(3, 1))


  ggplot2::ggsave(here::here(year, 'm24', "figs", "swath.png"),
                  width = 6.5, height = 5.5, units = "in", dpi = 200)


  purrr::map(r, 'M') %>% 
    purrr::map(., ~ as.data.frame(t(.))) %>% 
    bind_rows(.id = 'sim') %>% 
    ggplot(aes(V1)) +
    geom_histogram(ggplot2::aes(y=(after_stat(count))/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                   fill = "lightgray", color = "darkgray", bins = 50) +
    geom_vline(xintercept = m24$M)
  
  purrr::map(., bind_rows(.id = 'rep'))
str(post[1,])
obj24$report(post[1,])
mcmc$samples

post[['log_M']]
names(post)

obj

base_se %>% 
  mutate(model = 'base') %>% 
  bind_rows(m22.1_se %>% 
  mutate(model = 'm22.1')) %>% 
  bind_rows(m24_se %>% 
              mutate(model = "m24")) %>% 
  filter(item == 'spawn_bio') %>% 
  ggplot(aes(year, value, color = model, fill = model)) + 
  geom_ribbon(aes(ymin = lci, ymax =uci), alpha = 0.2, color=NA) +
  geom_line() +
  scico::scale_color_scico_d(palette = 'roma', end = 0.8) +
  scico::scale_fill_scico_d(palette = 'roma', end = 0.8) +
  expand_limits(y=0)




# mcmc output ----
obj = readRDS(here::here(year, "m24", "obj.rds"))
mcmc = readRDS(here::here(year, "m24", "mcmc.rds"))


png(filename=here::here(year, "compare_models", "srv_bio_wts.png"), width = 6.5, height = 6.5,
    units = "in", type ="cairo", res = 200)
data.frame(year = data$srv_yrs,
           obs = data$srv_obs,
           m22.1 = m22.1$srv_pred,
           m22.1a = m22.1a$srv_pred,
           m22.1b = m22.1b$srv_pred,
           m22.1c = m22.1c$srv_pred,
           m24 = m24$srv_pred,
           wt = 0.25) %>% 
  bind_rows(
  data.frame(year = data$srv_yrs,
             obs = data$srv_obs,
             m22.1 = m22.1_rewt$srv_pred,
             m22.1a = m22.1a_rewt$srv_pred,
             m22.1b = m22.1b_rewt$srv_pred,
             m22.1c = m22.1c_rewt$srv_pred,
             m24 = m24_rewt$srv_pred,
             wt = 1)) %>%   
  tidyr::pivot_longer(-c(year, obs, wt)) %>% 
  ggplot(aes(year, value, color = name)) + 
geom_line() + 
  geom_point(aes(y=obs), color = 'gray') +
  expand_limits(y=1) +
  facet_wrap(~wt) +
  scale_y_continuous(labels = scales::comma) +
  ylab('Survey biomass (t)') +
  scico::scale_color_scico_d(palette = 'roma')
dev.off()

png(filename=here::here(year, "compare_models", "spawn_bio_wts.png"), width = 6.5, height = 6.5,
    units = "in", type ="cairo", res = 200)
data.frame(year = data$years,
           m22.1 = m22.1$spawn_bio,
           m22.1a = m22.1a$spawn_bio,
           m22.1b = m22.1b$spawn_bio,
           m22.1c = m22.1c$spawn_bio,
           m24 = m24$spawn_bio,
           wt = 0.25) %>% 
  bind_rows(
    data.frame(year = data$years,
               m22.1 = m22.1_rewt$spawn_bio,
               m22.1a = m22.1a_rewt$spawn_bio,
               m22.1b = m22.1b_rewt$spawn_bio,
               m22.1c = m22.1c_rewt$spawn_bio,
               m24 = m24_rewt$spawn_bio,
               wt = 1)) %>%   
  tidyr::pivot_longer(-c(year, wt)) %>% 
  ggplot(aes(year, value, color = name)) + 
  geom_line() + 
  expand_limits(y=1) +
  facet_wrap(~wt) +
  scale_y_continuous(labels = scales::comma) +
  ylab('Spawning biomass (t)') +
  scico::scale_color_scico_d(palette = 'roma')
dev.off()

png(filename=here::here(year, "compare_models", "tot_bio_wts.png"), width = 6.5, height = 6.5,
    units = "in", type ="cairo", res = 200)
data.frame(year = data$years,
           m22.1 = m22.1$tot_bio,
           m22.1a = m22.1a$tot_bio,
           m22.1b = m22.1b$tot_bio,
           m22.1c = m22.1c$tot_bio,
           m24 = m24$tot_bio,
           wt = 0.25) %>% 
  bind_rows(
    data.frame(year = data$years,
               m22.1 = m22.1_rewt$tot_bio,
               m22.1a = m22.1a_rewt$tot_bio,
               m22.1b = m22.1b_rewt$tot_bio,
               m22.1c = m22.1c_rewt$tot_bio,
               m24 = m24_rewt$tot_bio,
               wt = 1)) %>%   
  tidyr::pivot_longer(-c(year, wt)) %>% 
  ggplot(aes(year, value, color = name)) + 
  geom_line() + 
  expand_limits(y=1) +
  facet_wrap(~wt) +
  scale_y_continuous(labels = scales::comma) +
  ylab('Total biomass (t)') +
  scico::scale_color_scico_d(palette = 'roma')
dev.off()
