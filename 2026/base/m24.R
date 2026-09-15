# load ----
library(RTMB)
library(tidyverse)
library(Matrix)
library(tmbstan)
library(shinystan)
library(here)
library(scico)
theme_set(theme_bw())
# devtools::install_github('Cole-Monnahan-NOAA/adnuts', ref='sparse_M')
library(adnuts)
# install.packages('StanEstimators', repos = c('https://andrjohns.r-universe.dev', 'https://cloud.r-project.org'))
library(StanEstimators)
source(here::here(2024, 'rtmb', 'bridge_model.r'))
source(here::here(2024, 'rtmb', 'utils.r'))

# globals
year = 2024
# grab maturity from updated ogive in ADMB output  

REP = readLines(here::here(2024, 'm22.1c', 'nr.rep'))
maa = as.numeric(stringr::str_split(REP[grep('Maturity', REP)], " ")[[1]][3:52])

fishery = grep("fish", list.files(here::here(year, 'data', "output")), value=TRUE)
survey = grep("bts_", list.files(here::here(year, 'data', "output")), value=TRUE)


yld = read.csv(here::here(year, "data", "output", 'yld_rat.csv'))
catch = read.csv(here::here(year, "data", "output", grep("catch", fishery, value=TRUE)))
wta = read.csv(here::here(year, "data", "output", "waa.csv"))
names(wta) <- c('age', 'wt')
saa = read.csv(here::here(year, "data", "output", "saa.csv"))
ae = read.csv(here::here(year, "data", "output", "ae_model.csv"))
fishac = read.csv(here::here(year, "data", "output", grep("age", fishery, value=TRUE))) %>% 
  filter(year>=1990)
fishlc = read.csv(here::here(year, "data", "output", grep("length", fishery, value=TRUE))) %>% 
  filter(year>=1991, year!=2024)
tsac = read.csv(here::here(year, "data", "output", grep("age", survey, value=TRUE))) %>% 
  filter(year>=1990)
tsb = read.csv(here::here(year, "data", "user_input", 'vast_lognormal.csv'))

wt_mature = wta$wt * maa / 2

data = list(ages = 2:45,
            years = catch$year,
            length_bins = 15:45,
            waa = wta$wt,
            maa = maa,
            wt_mature = wta$wt * maa / 2,
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
            srv_age_iss = c(49, 33, 91, 43, 86, 25, 61, 127, 171, 94, 101, 117, 123, 99, 155, 50),
            srv_age_obs = unname(t(as.matrix(tsac[,-(1:4)]))),
            srv_age_wt = 0.5,
            fish_size_yrs = fishlc$year ,
            fish_size_ind = ifelse(catch$year %in% fishlc$year, 1, 0),
            fish_size_iss = c(75.2673,69.1482,50.3742,42.9224,59.2575,46.4495,31.0813,75.2867,93.4583,83.2692,72.05,89.8439,100,57.9263,77.1133,64.5234,54.1795),
            fish_size_obs = unname(t(as.matrix(fishlc[,-c(1:4)]))),
            fish_size_wt = 0.5,
            age_error = as.matrix(ae),
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

saveRDS(data, here::here(2024, 'm24', 'dat.RDS'))

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
saveRDS(pars, here::here(2024, 'm24', 'pars.RDS'))
# parameter bounds - same as ADMB
lower = c(-Inf, # log M
          -Inf, #log a50C
          -Inf, # delta C
          -Inf, # log_a50S
          -Inf, # delta S
          -Inf, # logq
          -15, # log mean R
          rep(-10, length(pars$init_log_Rt)), # init rec devs
          rep(-10, length(data$years)), # rec devs
          -15, # log mean F
          rep(-15, length(data$years)), # Fdevs
          rep(-4.605,3)) # Fspr
saveRDS(lower, here::here(2024, 'm24', 'lower.RDS'))

upper = c(Inf, # log M
          Inf, #log a50C
          Inf, # delta C
          Inf, # log_a50S
          Inf, # delta S
          Inf, # logq
          10, # log mean R
          rep(10,  length(pars$init_log_Rt)), # init rec devs
          rep(10, length(data$years)), # rec devs
          15, # log mean F
          rep(15, length(data$years)), # Fdevs
          rep(0,3)) # Fspr
saveRDS(upper, here::here(2024, 'm24', 'upper.RDS'))
saveRDS(f, here::here(2024, 'm24', 'model.RDS'))
obj <- RTMB::MakeADFun(f, 
                       pars, 
                       map = list(sigmaR = factor(NA)))  
saveRDS(obj, here::here('2024', 'm24', "obj.RDS"))

fit <- nlminb(start = obj$par,
              objective = obj$fn,
              gradient = obj$gr,
              control = list(iter.max=100000,
                             eval.max=20000),
              lower=lower,
              upper=upper)

m24 <- obj$report(obj$env$last.par.best)
saveRDS(m24, file = here::here(2024, 'm24','m24.rds'))

proj_bio(m24)
rep = sdreport(obj)
summary(rep, "report") %>% 
  as.data.frame() %>% 
  mutate(lci = Estimate - 1.96*`Std. Error`,
         uci = Estimate + 1.96*`Std. Error`) %>% 
  tibble::rownames_to_column("item") %>% 
  mutate(item = gsub("\\..*", "", item),
         year = c(data$srv_yrs, rep(data$years, 3))) %>% 
  dplyr::select(year, item, value = Estimate, se = `Std. Error`, lci, uci) -> m24_se
  vroom::vroom_write(m24_se, here::here(year, 'm24', 'm24_se.csv'), delim = ",")
  

  m24_se %>% 
  filter(item == 'srv_pred') %>% 
  ggplot() +
     geom_point(data = tsb, aes(x=year, y = biomass), color = 'lightgray') +
    geom_errorbar(data = tsb, aes(x = year, ymin=lci, ymax = uci), color = 'lightgray', width = 0.2) +
  geom_ribbon(aes(x=year, ymin=lci, ymax=uci), alpha = 0.15) +
  geom_line(aes(year, value)) +
  expand_limits(y=0) +
    scale_y_continuous(labels = scales::comma) +
    ylab('Survey biomass')
  
  m24_se %>% 
    filter(item == 'spawn_bio') %>% 
    ggplot() +
    geom_ribbon(aes(x=year, ymin=lci, ymax=uci), alpha = 0.15) +
    geom_line(aes(year, value)) +
    expand_limits(y=0) +
    scale_y_continuous(labels = scales::comma) +
    ylab('Spawning biomass') +
    geom_hline(yintercept = c(m24$B35, m24$B40), lty = c(1,3))
  
  m24_se %>% 
    filter(item == 'tot_bio') %>% 
    ggplot() +
    geom_ribbon(aes(x=year, ymin=lci, ymax=uci), alpha = 0.15) +
    geom_line(aes(year, value)) +
    expand_limits(y=0) +
    scale_y_continuous(labels = scales::comma) +
    ylab('Total biomass')
  
  
  m24_se %>% 
    filter(item == 'recruits') %>% 
    mutate(lci = ifelse(lci<0, 0, lci)) %>% 
    ggplot(aes(year, value)) +
    geom_point() +
    geom_errorbar(aes(ymin=lci, ymax=uci), alpha = 0.15) +
    expand_limits(y=0) +
    ylab('Recruits')

  
  obj = readRDS(here::here(2024, 'm24', 'obj.RDS'))
  dat = readRDS(here::here(2024, 'm24', 'dat.RDS'))
  
  
  head(dat$catch_ind, -1)
  
  mcmc = readRDS(here::here(year, 'm24', 'mcmc.rds'))
  
