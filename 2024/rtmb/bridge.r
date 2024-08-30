# notes ----
# bridging GOA northern rockfish from ADMB to RTMB
# ben.williams@noaa.gov
# last update
# 2024-08

# load ----
library(RTMB)
library(tidyverse)
library(Matrix)
library(tmbstan)
library(shinystan)
theme_set(afscassess::theme_report())
# devtools::install_github('Cole-Monnahan-NOAA/adnuts', ref='sparse_M')
library(adnuts)
# install.packages('StanEstimators', repos = c('https://andrjohns.r-universe.dev', 'https://cloud.r-project.org'))
library(StanEstimators)

source(here::here(2024, 'rtmb', 'bridge_data.r'))
source(here::here(2024, 'rtmb', 'bridge_pars.r'))
source(here::here(2024, 'rtmb', 'bridge_model.r'))

# original model output 
catch <- read.csv(here::here(2024, 'base_srv_like_iss', 'processed', 'catch.csv'))
srv <- read.csv(here::here(2024, 'base_srv_like_iss', 'processed', 'survey.csv'))
fac <- read.csv(here::here(2024, 'base_srv_like_iss', 'processed', 'fac.csv'))
sac <- read.csv(here::here(2024, 'base_srv_like_iss', 'processed', 'sac.csv'))
fsc <- read.csv(here::here(2024, 'base_srv_like_iss', 'processed', 'fsc.csv'))
slx <- read.csv(here::here(2024, 'base_srv_like_iss', 'processed', 'selex.csv'))
bio <- read.csv(here::here(2024, 'base_srv_like_iss', 'processed', 'bio_rec_f.csv'))
b40 <- read.csv(here::here(2024, 'base_srv_like_iss', 'processed', 'b35_b40_yld.csv'))

# data ----
data <- list(
             ages = ages,
             years = years,
             length_bins = 15:45,
             waa = waa,
             wt_mature = maa * waa / 2,
             spawn_mo = 5,
             catch_ind = rep(1, length(years)),
             catch_obs = catch_obs,
             catch_wt = c(rep(5, 17), rep(50, 45)),
             srv_obs = srv_obs,
             srv_ind = srv_ind,
             srv_sd = srv_sd,
             srv_wt = 0.25,
             fish_age_obs = fish_age_obs,
             fish_age_ind = fish_age_ind,
             fish_age_iss = fish_age_iss,
             fish_age_wt = 0.5,
             srv_age_obs = srv_age_obs,
             srv_age_ind = srv_age_ind,
             srv_age_iss = srv_age_iss,
             srv_age_wt = 0.5,
             fish_size_obs = fish_size_obs,
             fish_size_ind = fish_size_ind,
             fish_size_iss = fish_size_iss,
             fish_size_wt = 0.5,
             age_error = age_error,
             size_age = size_age,
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
pars = list(log_M = log_M,
        log_a50C = log_a50C,
        deltaC = deltaC,
        log_a50S = log_a50S,
        deltaS = deltaS,
        log_q = log_q,
        log_mean_R = log_mean_R,
        init_log_Rt = init_log_Rt,
        log_Rt = log_Rt,
        log_mean_F = log_mean_F,
        log_Ft =  log_Ft,
        log_F35 = log_F35,
        log_F40 = log_F40,
        log_F50 = log_F50,
        sigmaR = sigmaR)

map = list(log_M = factor(NA),
           log_a50C = factor(NA),
           deltaC = factor(NA),
           log_a50S = factor(NA),
           deltaS = factor(NA),
           log_q = factor(NA),
           log_mean_R = factor(NA),
           init_log_Rt = factor(rep(NA, length(init_log_Rt))),
           log_Rt = factor(rep(NA, length(log_Rt))),
           log_mean_F = factor(NA),
           log_Ft = factor(rep(NA, length(log_Ft))),
           log_F35 = factor(NA),
           log_F40 = factor(NA),
           log_F50 = factor(NA),
           sigmaR = factor(NA))

# build the model 
# comparison run to admb model
obj <- RTMB::MakeADFun(f, 
                        pars, 
                        map = map)  
report <- obj$report(obj$env$last.par.best)

# examine model run with same starting point as ADMB
pars2 = list(log_M = log(0.05954247),
             log_a50C = log(7.5),
             deltaC = 3.0,
             log_a50S = log(7.3),
             deltaS = 3.8,
             log_q = log(1),
             log_mean_R = 4.3,
             init_log_Rt = rep(0, length(pars$init_log_Rt)),
             log_Rt = rep(0, length(years)),
             log_mean_F = 0,
             log_Ft =  rep(0, length(years)),
             log_F35 = 0,
             log_F40 = 0,
             log_F50 = 0,
             sigmaR = 1.5)

lower = c(#log(0.001), # log M
          log(3), #log a50C
          .5, # delta C
          log(3), # log_a50S
          0.5, # delta S
          log(0.5), # logq
          -15, # log mean R
          rep(-15, length(pars$init_log_Rt)), # init rec devs
          rep(-15, length(years)), # rec devs
          -15, # log mean F
          rep(-15, length(years)), # Fdevs
          rep(-4.605,3))#,  # Fspr
          0.3) # sigmaR

upper = c(#log(0.15), # log M
          log(12), #log a50C
          8.5, # delta C
          log(12), # log_a50S
          8.5, # delta S
          log(1.5), # logq
          10, # log mean R
          rep(15,  length(pars$init_log_Rt)), # init rec devs
          rep(15, length(years)), # rec devs
          15, # log mean F
          rep(15, length(years)), # Fdevs
          rep(0,3))#,  # Fspr
          10) # sigmaR

obj1 <- RTMB::MakeADFun(f, 
                        pars2,
                        map = list(log_M = factor(NA),
                                    sigmaR = factor(NA)))  
fit1 <- nlminb(start = obj1$par,
               objective = obj1$fn,
               gradient = obj1$gr,
               control = list(iter.max=100000,
                              eval.max=20000),
                lower = lower,
                upper = upper)

## for a stock assessment, inputs depend on whether RE exist or not
sdrep <- sdreport(obj1, getJointPrecision = TRUE)
Q <- sdrep$jointPrecision
## if no random effectrs this will be NULL
if(!is.null(Q)){
    M <- solve(Q)
} else {
    M <- sdrep$cov.fixed
    Q <- solve(M)
}               
report1 <- obj1$report(obj1$env$last.par.best)
report$q
report1$q
report$a50C
report1$a50C
report$deltaC
report1$deltaC
report1$a50S
report1$deltaS
plot(years, report$spawn_bio)
lines(years, report1$spawn_bio)
report$B40
report1$B40

proj_bio(report)
proj_bio(report1)

## if using RTMB, need to add data input to globals
globals <- list(data=data) # NULL
chains <- 3
mcmc <- sample_sparse_tmb(obj1, iter=4000, warmup=500,
                          init='random', chains=chains,
                          cores=chains, metric='dense',
                          Qinv=M, Q=Q,
                          globals=globals, skip_optimization=TRUE)
summary(mcmc) ## gives you minESS and maxRhat

## diagnose issues
pairs_admb(mcmc, pars=1:5, order='slow')
pairs_admb(mcmc, pars=1:5, order='mismatch')


launch_shinyadmb(mcmc)


# update model with priors on selectivity
data$mean_a50C = 7.5
data$cv_a50C = 1
data$mean_deltaC = 3.8
data$cv_deltaC = 1
data$mean_q
data$cv_q = 0.25

obj2 <- RTMB::MakeADFun(f1, 
                        pars2,
                        map = list(log_M = factor(NA),
                                   sigmaR = factor(NA))) 

fit2 <- nlminb(obj2$par,
               obj2$fn,
               obj2$gr,
               control = list(iter.max=100000,
                              eval.max=20000),
                              lower = lower,
                              upper = upper)
## for a stock assessment, inputs depend on whether RE exist or not
sdrep <- sdreport(obj2, getJointPrecision = TRUE)
Q <- sdrep$jointPrecision
## if no random effectrs this will be NULL
if(!is.null(Q)){
    M <- solve(Q)
} else {
    M <- sdrep$cov.fixed
    Q <- solve(M)
}               
report2 <- obj2$report(obj2$env$last.par.best)
report2$log_mean_R
report2$spawn_bio
report2$B40
report2$q

## if using RTMB, need to add data input to globals
globals <- list(data=data) # NULL
chains <- 3
mcmc <- sample_sparse_tmb(obj2, iter=5000, warmup=400,
                          init='random', chains=chains,
                          cores=chains, metric='dense',
                          Qinv=M, Q=Q,
                          globals=globals, skip_optimization=TRUE)
summary(mcmc) ## gives you minESS and maxRhat

## diagnose issues
pairs_admb(mcmc, pars=1:5)
pairs_admb(mcmc, pars=1:5, order='slow')
pairs_admb(mcmc, pars=1:5, order='mismatch')


launch_shinyadmb(mcmc)
report2$a50C
report2$a50S
report2$deltaS
report2$deltaC
report2$slx
plot(1:50, report2$slx[,1], type= 'l')
lines(1:50, report2$slx[,2], col=4)
lines(1:50, c(0.003730717, 0.007114453, 0.013525548, 0.025565162, 0.047802384, 0.087642403, 0.155271892, 0.260204279, 0.402279271, 0.56290318, 0.711336288, 0.825031111, 0.900226422, 0.945250043, 0.970619488, 0.984427181, 0.991800615, 0.995698124, 0.997747193, 0.998821405, 0.999383714, 0.999677831, 0.999831607, 0.99991199, 0.999954004, 0.999975962, 0.999987437, 0.999993435, 0.999996569, 0.999998207, 0.999999063, 0.99999951, 0.999999744, 0.999999866, 0.99999993, 0.999999963, 0.999999981, 0.99999999, 0.999999995, 0.999999997, 0.999999999, 0.999999999, 1, 1, 1, 1, 1, 1, 1, 1)       
, col=2)
lines(1:50, slx[,2])
lines(1:50, slx[,3], col=4)

a50C
[1] 8.093552
a50S
[1] 8.565348
deltaS
[1] 3.905736
deltaC
[1] 0.9329942








### Fishery size compositions

data.frame(length = length_bins, 
           report$fish_size_pred) %>% 
  pivot_longer(-length) %>% 
  mutate(year = rep(fish_size_yrs, length(length_bins)),
         groups = 'rtmb') %>% 
  bind_rows(data.frame(length = length_bins, 
             report$fish_size_pred) %>% 
  pivot_longer(-length) %>% 
  mutate(year = rep(fish_size_yrs, length(length_bins)),
         groups = 'rtmb.1')) %>% 
  bind_rows(fsc) -> df4

df4 %>% 
  mutate(groups = ifelse(groups=='pred', 'base.1a', groups)) %>% 
  ggplot(aes(length, value, color = groups, shape = groups, linetype=groups)) + 
  geom_point(show.legend=FALSE) +
  geom_line(show.legend=T) +
  scale_shape_manual(values = c(NA,19,NA,NA),guide = guide_none()) +
  scale_linetype_manual(values = c(1,0,1,1),guide = guide_none()) +
  scico::scale_color_scico_d(name="",palette = 'roma', begin=0.2) +
  facet_wrap(~year)

ggsave('fsc.png', width=6.5, height=6.5, units='in')

df4 %>% 
  select(length, value, year, groups) %>% 
  pivot_wider(names_from=groups, values_from = value) %>% 
  mutate(diff = (rtmb - pred) / pred,
         Length = factor(length)) %>% 
  pivot_longer(diff) %>% 
  ggplot(aes(year, value, color = Length)) + 
  geom_point() +
  # facet_wrap(~year) +
  scale_y_continuous( labels = scales::percent) +
  scico::scale_color_scico_d(palette = 'roma') +
  ylab('percent difference') +
  geom_hline(yintercept=0, lty=3)

ggsave(here::here(2024, 'base_srv_like_iss', 'figs', "fsc-diff.png"), width=6.5, height=6.5, units='in')


# comparisons ----
## selectivity
as.data.frame(report$slx) %>% 
  rename(fishery = V1, survey = V2) %>% 
  bind_cols(as.data.frame(report1$slx) %>% 
              rename(fishery.1 = V1, survey.1 = V2)) %>% 
  bind_cols(slx) -> slxs 
  
slxs %>% 
  rename(`fishery-rtmb`=fishery, `srv-rtmb`=survey, `fishery-rtmb.1`=fishery.1, `srv-rtmb.1`=survey.1,`fishery-base.1a`=fish,`srv-base.1a`=srv1) %>% 
  pivot_longer(-c(age, maturity)) %>% 
  ggplot(aes(age, value, color = name)) + 
  geom_line() +
  scico::scale_color_scico_d(name="", palette = 'roma') +
  theme(legend.position = c(0.8, 0.2))

ggsave(here::here(2024, 'base_srv_like_iss', 'figs', "slx.png"), width=6.5, height=4.5, units='in')

slxs %>% 
  mutate(fishery = (fishery - fish)/fish,
         survey = (survey - srv1)/srv1) %>% 
  dplyr::select(age, fishery, survey) %>% 
  tidyr::pivot_longer(-age) %>% 
  ggplot(aes(age, value, color = name)) + 
  geom_point() +
  scale_y_continuous(labels = scales::percent) +
  scico::scale_color_scico_d(name="",palette = 'roma') +
  ylab('percent difference') +
  geom_hline(yintercept=0, lty=3) +
  theme(legend.position = c(0.8, 0.2))
ggsave(here::here(2024, 'base_srv_like_iss', 'figs', "slx-diff.png"), width=6.5, height=4.5, units='in')

## biomass
data.frame('ssb-rtmb' = report$spawn_bio,
           'tot-rtmb' = report$tot_bio) %>% 
  bind_cols(data.frame('ssb-rtmb.1' = report1$spawn_bio,
                       'tot-rtmb.1' = report1$tot_bio)) %>% 
  bind_cols(bio) %>% 
  rename(`tot-base.1a`=tot_biom, `ssb-base.1a`=sp_biom, `tot-rtmb`=tot.rtmb, `ssb-rtmb`=ssb.rtmb, `tot-rtmb.1`=tot.rtmb.1, `ssb-rtmb.1`=ssb.rtmb.1) %>% 
  pivot_longer(-c(year, F, recruits)) -> Bs
  
Bs %>% 
  ggplot(aes(year, value, color = name)) + 
  geom_line() +
  scico::scale_color_scico_d(name="",palette = 'roma') +
  expand_limits(y=0) +
  theme(legend.position = c(0.8, 0.2)) +
  ylab('Biomass') +
  xlab('Year')
ggsave(here::here(2024, 'base_srv_like_iss', 'figs', "bio.png"), width=6.5, height=4.5, units='in')


Bs %>% 
  pivot_wider(values_from=value, names_from=name) %>% 
  mutate(tot = (`tot-rtmb` - `tot-base.1a`) / `tot-base.1a`,
         ssb = (`ssb-rtmb` - `ssb-base.1a`) / `ssb-base.1a`) %>% 
  pivot_longer(c(tot, ssb)) %>% 
  ggplot(aes(year, value, color = name)) + 
  geom_point() +
  # facet_wrap(~year) +
  scale_y_continuous(labels = scales::percent) +
  scico::scale_color_scico_d(name="", palette = 'roma') +
  ylab('percent difference') +
  geom_hline(yintercept=0, lty=3)

ggsave(here::here(2024, 'base_srv_like_iss', 'figs', "bio-diff.png"), width=6.5, height=4.5, units='in')

#F

 Bs %>% 
   mutate(id = case_when(grepl('ssb', name) ~ sub('ssb-', "", name),
                         grepl('tot', name) ~ sub('tot-', "", name))) -> Fs
   
 Fs %>% 
   ggplot(aes(year, F, color = id)) +
   geom_line() +
   scico::scale_color_scico_d(name="",palette = 'roma') +
   expand_limits(y=0) +
   theme(legend.position = c(0.8, 0.4)) +
   ylab('F') +
   xlab('Year')

  Fs %>% 
    distinct(year, F, id) %>% 
    pivot_wider(values_from=F, names_from=id) %>% 
    mutate(diff = (rtmb - `base.1a`) / `base.1a`) %>% 
  ggplot(aes(year, diff)) + 
  geom_point() +
  # facet_wrap(~year) +
  scale_y_continuous( labels = scales::percent) +
  scico::scale_color_scico_d(palette = 'roma') +
  ylab('percent difference') +
  geom_hline(yintercept=0, lty=3)
  
  
# survey biomass
data.frame(rtmb = report$srv_pred) %>% 
  bind_cols(rtmb.1 = report1$srv_pred) %>% 
  bind_cols(srv) %>% 
  rename(base.1a = pred) -> ds

  ds %>% 
    pivot_longer(c(rtmb, base.1a, rtmb.1)) %>% 
  ggplot(aes(year, value, color = name)) + 
  geom_point(aes(y=biomass), color = 'black') +
  geom_line() +
  scico::scale_color_scico_d(name="",palette='roma', end = 0.6) +
  expand_limits(y=0) +
  ylab('Biomass') +
  xlab('Years') +
  theme(legend.position = c(0.8, 0.2))

  ggsave(here::here(2024, 'base_srv_like_iss', 'figs', "srvbio.png"), width=6.5, height=4.5, units='in')
  
## fish age comp
  data.frame(length = length_bins, 
             report$fish_size_pred) %>% 
    pivot_longer(-length) %>% 
    mutate(year = rep(fish_size_yrs, length(length_bins)),
           groups = 'rtmb') %>% 
    bind_rows(data.frame(length = length_bins, 
                         report$fish_size_pred) %>% 
                pivot_longer(-length) %>% 
                mutate(year = rep(fish_size_yrs, length(length_bins)),
                       groups = 'rtmb.1')) %>% 
    bind_rows(fsc) -> df4
  
data.frame(age = data$ages, 
           report$fish_age_pred) %>% 
  pivot_longer(-age) %>% 
  mutate(year = rep(fish_age_yrs, length(ages)),
         groups = 'rtmb') %>%
  bind_rows(data.frame(age = data$ages, 
                       report1$fish_age_pred) %>% 
              pivot_longer(-age) %>% 
              mutate(year = rep(fish_age_yrs, length(ages)),
                     groups = 'rtmb.1')) %>% 
  bind_rows(fac)-> df5
  
df5 %>%  
  mutate(groups = ifelse(groups=='pred', 'admb.1a', groups)) %>% 
  ggplot(aes(age, value, color = groups, shape = groups, linetype=groups)) + 
  geom_point(show.legend=FALSE) +
  geom_line(show.legend=T) +
  scale_shape_manual(values = c(NA,19,NA),guide = guide_none()) +
  scale_linetype_manual(values = c(1,0,1),guide = guide_none()) +
  scico::scale_color_scico_d(name="",palette = 'roma', begin=0.2) +
  facet_wrap(~year)
ggsave(here::here(2024, 'base_srv_like_iss', 'figs','fac.png'), width=6.5, height=6.5, units='in')

df5 %>% 
  select(age, value, year, groups) %>% 
  filter(groups!='obs') %>% 
  pivot_wider(names_from=groups, values_from = value) %>% 
  mutate(diff = (rtmb - pred) / pred,
         Age = factor(age)) %>% 
  pivot_longer(diff) %>% 
  ggplot(aes(year, value, color = Age)) + 
  geom_point() +
  # facet_wrap(~year) +
  scale_y_continuous( labels = scales::percent) +
  scico::scale_color_scico_d(palette = 'roma') +
  ylab('percent difference') +
  geom_hline(yintercept=0, lty=3)

ggsave(here::here(2024, 'base_srv_like_iss', 'figs', "fac-diff.png"), width=6.5, height=6.5, units='in')


## survey age comp
data.frame(age = data$ages, 
           report$srv_age_pred) %>% 
  pivot_longer(-age) %>% 
  mutate(year = rep(srv_age_yrs, length(ages)),
         groups = 'rtmb') %>% 
  bind_rows(sac)-> df6

df6 %>%  
  mutate(groups = ifelse(groups=='pred', 'admb', groups)) %>% 
  ggplot(aes(age, value, color = groups, shape = groups, linetype=groups)) + 
  geom_point(show.legend=FALSE) +
  geom_line(show.legend=T) +
  scale_shape_manual(values = c(NA,19,NA),guide = guide_none()) +
  scale_linetype_manual(values = c(1,0,1),guide = guide_none()) +
  scico::scale_color_scico_d(name="",palette = 'roma', begin=0.2) +
  facet_wrap(~year)
ggsave(here::here(2024, 'base_srv_like_iss', 'figs','sac.png'), width=6.5, height=6.5, units='in')

df6 %>% 
  select(age, value, year, groups) %>% 
  filter(groups!='obs') %>% 
  pivot_wider(names_from=groups, values_from = value) %>% 
  mutate(diff = (rtmb - pred) / pred,
         Age = factor(age)) %>% 
  pivot_longer(diff) %>% 
  ggplot(aes(year, value, color = Age)) + 
  geom_point() +
  # facet_wrap(~year) +
  scale_y_continuous( labels = scales::percent) +
  scico::scale_color_scico_d(palette = 'roma') +
  ylab('percent difference') +
  geom_hline(yintercept=0, lty=3)

ggsave(here::here(2024, 'base_srv_like_iss', 'figs', "sac-diff.png"), width=6.5, height=6.5, units='in')


# compare 
report1

## survey age comp
data.frame(age = data$ages, 
           report$srv_age_pred) %>% 
  pivot_longer(-age) %>% 
  mutate(year = rep(srv_age_yrs, length(ages)),
         groups = 'rtmb') %>% 
  bind_rows(data.frame(age = data$ages, 
                       report1$srv_age_pred) %>% 
              pivot_longer(-age) %>% 
              mutate(year = rep(srv_age_yrs, length(ages)),
                     groups = 'rtmb1')) %>% 
  bind_rows(sac) -> df6

df6 %>%  
  mutate(groups = ifelse(groups=='pred', 'admb', groups)) %>% 
  ggplot(aes(age, value, color = groups, shape = groups, linetype=groups)) + 
  geom_point(show.legend=FALSE) +
  geom_line(show.legend=T) +
  scale_shape_manual(values = c(NA,19,NA, NA),guide = guide_none()) +
  scale_linetype_manual(values = c(1,0,1,1),guide = guide_none()) +
  scico::scale_color_scico_d(name="",palette = 'roma', begin=0.2) +
  facet_wrap(~year)


# survey biomass
data.frame(rtmb = report$srv_pred,
           rtmb1 = report1$srv_pred) %>% 
  bind_cols(srv) %>% 
  rename(admb = pred) -> ds

ds %>% 
  pivot_longer(c(rtmb, rtmb1)) %>% 
  ggplot(aes(year, value, color = name)) + 
  geom_point(aes(y=biomass), color = 'black') +
  geom_line() +
  scico::scale_color_scico_d(name="",palette='roma') +
  expand_limits(y=0) +
  ylab('Biomass') +
  xlab('Years') +
  theme(legend.position = c(0.8, 0.2))


# ggsave(here::here(2024, 'base_srv_like_iss', 'figs', "srvbio.png"), width=6.5, height=4.5, units='in')

## fish age comp
data.frame(age = data$ages, 
           report$fish_age_pred) %>% 
  pivot_longer(-age) %>% 
  mutate(year = rep(fish_age_yrs, length(ages)),
         groups = 'rtmb') %>% 
  bind_rows(data.frame(age = data$ages, 
                       report1$fish_age_pred) %>% 
              pivot_longer(-age) %>% 
              mutate(year = rep(fish_age_yrs, length(ages)),
                     groups = 'rtmb.1')) %>% 
  bind_rows(fac %>% 
              mutate(groups = 'admb.1a'))-> df5

df5 %>%  
  # mutate(groups = ifelse(groups=='pred', 'admb.1a', groups)) %>% distinct(groups)
  ggplot(aes(age, value, color = groups, shape = groups, linetype=groups)) + 
  geom_point(show.legend=FALSE) +
  geom_line(show.legend=T) +
  scale_shape_manual(values = c(19,19, 19),guide = guide_none()) +
  # scale_linetype_manual(values = c(1,0,1,1),guide = guide_none()) +
  scico::scale_color_scico_d(name="",palette = 'roma', begin=0.2) +
  facet_wrap(~year)

### Fishery size compositions

data.frame(age = data$ages, 
           report$fish_age_pred) %>% 
  pivot_longer(-age) %>% 
  mutate(year = rep(fish_age_yrs, length(ages)),
         groups = 'rtmb') %>% 
  bind_rows(data.frame(age = data$ages, 
                       report1$fish_age_pred) %>% 
              pivot_longer(-age) %>% 
              mutate(year = rep(fish_age_yrs, length(ages)),
                     groups = 'rtmb.1')) %>% 
  bind_rows(fac) -> df5

df5 %>% 
  mutate(groups = ifelse(groups=='pred', 'base.1a', groups)) %>% 
  ggplot(aes(age, value, color = groups, shape = groups, linetype=groups)) + 
  geom_point(show.legend=FALSE) +
  geom_line(show.legend=T) +
  scale_shape_manual(values = c(NA,19,NA,NA),guide = guide_none()) +
  scale_linetype_manual(values = c(1,0,1,1),guide = guide_none()) +
  scico::scale_color_scico_d(name="",palette = 'roma', begin=0.2) +
  facet_wrap(~year)

ggsave('fsc.png', width=6.5, height=6.5, units='in')

df5 %>% 
  select(age, value, year, groups) %>% 
  pivot_wider(names_from=groups, values_from = value) %>% 
  mutate(diff = (rtmb - pred) / pred,
         Age = factor(age)) %>% 
  pivot_longer(diff) %>% 
  ggplot(aes(year, value, color = Age)) + 
  geom_point() +
  # facet_wrap(~year) +
  scale_y_continuous( labels = scales::percent) +
  scico::scale_color_scico_d(palette = 'roma') +
  ylab('percent difference') +
  geom_hline(yintercept=0, lty=3)

ggsave(here::here(2024, 'base_srv_like_iss', 'figs', "fac-diff.png"), width=6.5, height=6.5, units='in')


par_name = 'log_q'
par_values = seq(-5, 5, length.out = 50)
par_like = single_likelihood(par_name, par_values, obj1, fit1)
plot(par_values, par_like, type="l", xlab=par_name, ylab="Log-Likelihood")

data.frame(pars = par_values, 
           log_q = par_like) -> likes

par_name = 'log_a50S'
par_values = seq(-5, 5, length.out = 50)
par_like = single_likelihood(par_name, par_values, obj1, fit1)
lines(par_values, par_like, col = 4)

likes$log_a50S = par_like
likes %>% 
  # mutate(log_delta_catch = ifelse(log_delta_catch>265, NA, log_delta_catch)) %>% 
  pivot_longer(-c(pars)) %>% 
  filter(value<500) %>%
  ggplot(aes(pars, value, color = name)) + 
  geom_line()
glimpse(likes)
