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
             maa = maa,
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
proj_bio(report)
report$B40

# examine model run with same starting point as ADMB
pars2 = list(log_M = log(0.06),
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
obj1$report(obj1$env$last.par.best)$spawn_bio
fit1 <- nlminb(start = obj1$par,
               objective = obj1$fn,
               gradient = obj1$gr,
               control = list(iter.max=100000,
                              eval.max=20000),
                lower = lower,
                upper = upper)
names(obj1$par)
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
proj_bio(report1)
report1$M
report1$q
vals = seq(0, 1, length.out = 50)
plot((vals), single_likelihood(param_name='log_q', 
                  param_values = vals,
                  fit = fit1,
                  obj = obj1),
     type = 'l')

data.frame(vals = (seq(log(3), log(12), length.out = 50)),
           a50C = single_likelihood(param_name='log_a50C', 
                             param_values = vals,
                             fit = fit2,
                             obj = obj2),
           a50S = single_likelihood(param_name='log_a50S', 
                                    param_values = vals,
                                    fit = fit2,
                                    obj = obj2)) %>% 
  pivot_longer(-vals, names_to = 'parameter', values_to = 'likelihood') %>% 
  ggplot(aes(vals, likelihood, color = parameter)) + 
  geom_line()

vals = seq(-1.2, .5, length.out = 50)
plot(exp(vals), single_likelihood(param_name='log_q', 
                             param_values = vals,
                             fit = fit2,
                             obj = obj2),
     type = 'l')


report1$q
report1$M
report$a50C
report1$a50C
report$deltaC
report1$deltaC
report1$a50S
report1$deltaS
report$deltaS
plot(years, report$spawn_bio)
lines(years, report1$spawn_bio)
report$B40
report1$B40

proj_bio(report)
proj_bio(report1)

## if using RTMB, need to add data input to globals
globals <- list(data=data) # NULL
chains <- 3
# Cole recommends 5 chains, 1000 iters, 250 warmup
mcmc <- sample_sparse_tmb(obj1, iter=4000, warmup=500,
                          init='random', chains=chains,
                          cores=chains, metric='dense',
                          Qinv=M, Q=Q,
                          globals=globals, skip_optimization=TRUE)
summary(mcmc) ## gives you minESS and maxRhat
summary(mcmc$monitor$n_eff)
summary(mcmc$monitor$Rhat)
post <- extract_samples(mcmc, as.list=TRUE)
postlist <- coda::mcmc.list(lapply(post, coda::mcmc))
coda::traceplot(postlist)
glimpse(post)
sp <- extract_sampler_params(mcmc)
glimpse(sp)
plot_marginals(mcmc, pars=1:5)
## What if you want a posterior for derived quantities in the report? Just
## loop through each posterior sample (row) and call the report function
## which returns a list. The last column is the log-posterior density (lp__)
## and needs to be dropped

obj1$report(post[1,-ncol(post)])         # sd0 is only element
sd0 <- rep(NA, len=nrow(post))
for(i in 1:nrow(post)){
  r <- obj2$report(post[1,-ncol(post)])
  sd0[i] <- r$B40
}
hist(sd0)

get_item <- function(item='spawn_bio', post, obj, reps = NULL) {
  mapit = function(i, post, obj) {
    surv_i = paste0(obj$report(post[i,-ncol(post)]),"$", item)
    list(surv_i = surv_i)
  }
  results = purrr::map(1:reps, mapit, post=post, obj=obj)
  do.call(cbind, map(results, "surv_i")) 
}


# also try plot_uncertainties()

## diagnose issues
pairs_admb(mcmc2, pars=1:5, order='slow')
pairs_admb(mcmc2, pars=1:5, order='mismatch')


launch_shinyadmb(mcmc)


# update model with priors on selectivity
# 2nd try ----
data$mean_a50C = 7.5
data$cv_a50C = 1
data$mean_deltaC = 3.8
data$cv_deltaC = 1
data$mean_q
data$cv_q = 0.45
data$cv_M = 0.15
obj2 <- RTMB::MakeADFun(f1, 
                        pars2,
                        map = list(#log_M = factor(NA),
                                   sigmaR = factor(NA))) 
obj2$report(obj2$env$last.par.best)$spawn_bio
lower = c(log(0.001), # log M
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

upper = c(log(0.15), # log M
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
report2$M
report2$q
report2$B40

report2$a50C
report2$deltaC
report2$a50S
report2$deltaS
report2$log_mean_R
report2$F40
report2$B0
report2$B40
proj_bio(report2)
report2$like_srv
report2$like_fish_age
report2$like_srv_age
report2$like_fish_size
report2$like_rec
report2$sprpen
report2$f_regularity
report2$nll_M
report2$nll_q
report2$nll_a50C
report2$nll_deltaC
report2$nll_a50S
report2$nll_deltaS

valsq = seq(-1.2, 1.2, length.out = 50)
valsM = seq(log(0.01), log(0.15), length.out = 50)
data.frame(valsq,
           par = single_likelihood(param_name='log_q', 
                             param_values = valsq,
                             fit = fit1,
                             obj = obj1)) %>% 
  # dplyr::mutate(vals = exp(vals)) %>% 
  ggplot(aes(valsq, par)) + 
  geom_line()

valsM = seq(log(0.00001), log(0.15), length.out = 50)
data.frame(vals=valsM,
           M = single_likelihood(param_name='log_M', 
                                 param_values = valsM,
                                 fit = fit1,
                                 obj = obj1)) %>% 
  # dplyr::mutate(vals = exp(vals)) %>% 
  ggplot(aes(vals, M)) + 
  geom_line() +
  xlab('Natural mortality (M)') +
  ylab('Likelihood')


vals = seq(log(5), log(12), length.out = 50)
data.frame(vals,
           a50C = single_likelihood(param_name='log_a50C', 
                                    param_values = vals,
                                    fit = fit1,
                                    obj = obj1),
           a50S = single_likelihood(param_name='log_a50S', 
                                    param_values = vals,
                                    fit = fit1,
                                    obj = obj1)) %>% 
  # mutate(vals = exp(vals)) %>% 
  pivot_longer(-vals, names_to = 'parameter', values_to = 'likelihood') %>% 
  ggplot(aes(vals, likelihood, color = parameter)) + 
  geom_line()

pairs_admb(mcmc, order='mismatch', pars=1:5) 
## if using RTMB, need to add data input to globals
globals <- list(data=data) # NULL
chains <- 5
# Cole recommends 5 chains, 1000 iters, 250 warmup
mcmc2 <- sample_sparse_tmb(obj2, iter=500, warmup=250,
                          init='random', chains=chains,
                          cores=chains, metric='dense',
                          Qinv=M, Q=Q,control=list(adapt_delta=0.95),
                          globals=globals, skip_optimization=TRUE)

saveRDS(mcmc2, file = here::here(2024, 'rtmb', 'mcmc.RData'))
readRDS(here::here(2024, 'rtmb', 'mcmc.RData'))

summary(mcmc2) ## gives you minESS and maxRhat
summary(mcmc2$monitor$n_eff)
summary(mcmc2$monitor$Rhat)
post2 <- extract_samples(mcmc2, as.list=TRUE)
postlist2 <- coda::mcmc.list(lapply(post2, coda::mcmc))
coda::traceplot(postlist2)
glimpse(post2)
sp2 <- extract_sampler_params(mcmc2)

np = extract_sampler_params(mcmc2) %>%
  pivot_longer(-c(chain, iteration), names_to='Parameter', values_to='Value') %>% 
  select(Iteration=iteration, Parameter, Value, Chain=chain) %>%
  mutate(Parameter=factor(Parameter),
         Iteration=as.integer(Iteration),
         Chain=as.integer(Chain)) 
  

adnuts::plot_uncertainties(mcmc2)
aaa = adnuts::extract_samples(mcmc2)

mcmc_nuts_acceptance(np) + ggtitle("NUTS Energy Diagnostic") + theme_minimal()

glimpse(sp2)
plot_marginals(mcmc2, pars=1:6)
## diagnose issues
pairs_admb(mcmc2, pars=1:5, order='slow')
pairs_admb(mcmc2, pars=1:5, order='mismatch')

post = as.matrix(mcmc2)
rstan::traceplot(mcmc2, pars = c("log_M", "log_q"), inc_warmup = FALSE, nrow = 2)
# # 3nd try ----

mcmc2$mle$est

data$srv_age_wt <- 1
data$srv_wt <- 1
data$fish_age_wt <- 1
data$fish_size_wt <- 1

obj3 <- RTMB::MakeADFun(f1, 
                        pars2,
                        map = list(log_M = factor(NA),
                                   sigmaR = factor(NA))) 
fit3 <- nlminb(obj3$par,
               obj3$fn,
               obj3$gr,
               control = list(iter.max=100000,
                              eval.max=20000),
               lower = lower,
               upper = upper)
## for a stock assessment, inputs depend on whether RE exist or not
sdrep3 <- sdreport(obj3, getJointPrecision = TRUE)

fit3$par["log_q"]
# Q <- sdrep$jointPrecision
# ## if no random effectrs this will be NULL
# if(!is.null(Q)){
#   M <- solve(Q)
# } else {
#   M <- sdrep$cov.fixed
#   Q <- solve(M)
# }               
report3 <- obj3$report(obj3$env$last.par.best)
report3$M
report3$q
report3$sigmaR
proj_bio(report3)
report1$log_mean_R
report2$log_mean_R
report2$spawn_bio
report3$spawn_bio
report2$B40
report3$B40

report2$q
report3$q
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

proj_bio(report2)

adnuts::launch_shinyadmb(mcmc)
report2$a50C
report2$a50S
report2$deltaS
report2$deltaC
report2$q
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
             report2$fish_size_pred) %>% 
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
  # bind_cols(as.data.frame(report2$slx) %>% 
  #             rename(fishery.1 = V1, survey.1 = V2)) %>% 
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
data.frame('ssb-RTMB' = report$spawn_bio,
           'tot-RTMB' = report$tot_bio) %>% 
  # bind_cols(data.frame('ssb-RTMB.1' = report2$spawn_bio,
  #                      'tot-RTMB.1' = report2$tot_bio)) %>% 
  bind_cols(bio) %>% 
  rename(`tot-Base.1a`=tot_biom, `ssb-Base.1a`=sp_biom) %>% 
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
  mutate(tot = (`tot.RTMB` - `tot-Base.1a`) / `tot-Base.1a`,
         ssb = (`ssb.RTMB` - `ssb-Base.1a`) / `ssb-Base.1a`) %>% 
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
                         grepl('tot', name) ~ sub('tot-', "", name))) %>% 
   distinct(year, F, id) %>% 
   filter(!(id %in% c('ssb.RTMB', 'ssb.RTMB.1'))) -> Fs
   
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
    mutate(diff = (tot.RTMB - `Base.1a`) / `Base.1a`) %>% 
  ggplot(aes(year, diff)) + 
  geom_point() +
  # facet_wrap(~year) +
  scale_y_continuous( labels = scales::percent) +
  scico::scale_color_scico_d(palette = 'roma') +
  ylab('percent difference') +
  geom_hline(yintercept=0, lty=3)
  
  
# survey biomass
data.frame(RTMB = report$srv_pred) %>% 
  # bind_cols(RTMB.1 = report2$srv_pred) %>% 
  bind_cols(srv) %>% 
  rename(Base.1a = pred) -> ds

  ds %>% 
    pivot_longer(c(RTMB, Base.1a, RTMB.1)) %>% 
  ggplot(aes(year, value, color = name)) + 
  geom_point(aes(y=biomass), color = 'darkgray') +
    geom_errorbar(aes(ymin=lci, ymax=uci), color = 'darkgray', width = 0.2) +
  geom_line() +
  scico::scale_color_scico_d(name="",palette='roma', end = 0.6) +
  expand_limits(y=0) +
  ylab('Biomass') +
  xlab('Years') +
  theme(legend.position = c(0.7, 0.2))

  ggsave(here::here(2024, 'base_srv_like_iss', 'figs', "srvbio.png"), width=6.5, height=4.5, units='in')
  
  
  ds %>% 
    select(year, RTMB, Base.1a) %>% 
    mutate(diff = (RTMB - Base.1a) / Base.1a) %>% 
    ggplot(aes(year, diff)) + 
    geom_point() +
    xlab('Year') +
    scale_y_continuous(labels = scales::percent) +
    # scico::scale_color_scico_d(name="", palette = 'roma') +
    ylab('percent difference') +
    geom_hline(yintercept=0, lty=3)
  
  ggsave(here::here(2024, 'base_srv_like_iss', 'figs', "srvbio-diff.png"), width=6.5, height=4.5, units='in')
  
## fish age comp
  data.frame(length = length_bins, 
             report$fish_size_pred) %>% 
    pivot_longer(-length) %>% 
    mutate(year = rep(fish_size_yrs, length(length_bins)),
           groups = 'rtmb') %>% 
    # bind_rows(data.frame(length = length_bins, 
    #                      report2$fish_size_pred) %>% 
    #             pivot_longer(-length) %>% 
    #             mutate(year = rep(fish_size_yrs, length(length_bins)),
    #                    groups = 'rtmb.1')) %>% 
    bind_rows(fsc) -> df4
  
data.frame(age = data$ages, 
           report$fish_age_pred) %>% 
  pivot_longer(-age) %>% 
  mutate(year = rep(fish_age_yrs, length(ages)),
         groups = 'RTMB') %>%
  # bind_rows(data.frame(age = data$ages, 
  #                      report2$fish_age_pred) %>% 
  #             pivot_longer(-age) %>% 
  #             mutate(year = rep(fish_age_yrs, length(ages)),
  #                    groups = 'RTMB.1')) %>% 
  bind_rows(fac) -> df5
  
df5 %>%  
  mutate(groups = ifelse(groups=='pred', 'Base.1a', groups)) %>% 
  ggplot(aes(age, value, color = groups, shape = groups, linetype=groups)) + 
  geom_point(show.legend=FALSE) +
  geom_line(show.legend=T) +
  # scale_shape_manual(values = c(NA,19,NA, NA),guide = guide_none()) +
  # scale_linetype_manual(values = c(1,0,1, 1),guide = guide_none()) +
  scico::scale_color_scico_d(name="",palette = 'roma', begin=0.2) +
  facet_wrap(~year)
ggsave(here::here(2024, 'base_srv_like_iss', 'figs','fac.png'), width=6.5, height=6.5, units='in')

df5 %>% 
  select(age, value, year, groups) %>% 
  filter(groups!='obs') %>% 
  pivot_wider(names_from=groups, values_from = value) %>% 
  mutate(diff = (RTMB - pred) / pred,
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

data.frame(age = data$ages, 
           report$srv_age_pred) %>% 
  pivot_longer(-age) %>% 
  mutate(year = rep(srv_age_yrs, length(ages)),
         groups = 'RTMB') %>%
  bind_rows(data.frame(age = data$ages, 
                       report2$srv_age_pred) %>% 
              pivot_longer(-age) %>% 
              mutate(year = rep(srv_age_yrs, length(ages)),
                     groups = 'RTMB.1')) %>% 
  bind_rows(sac) -> df6


df6 %>%  
  mutate(groups = ifelse(groups=='pred', 'Base.1a', groups)) %>% 
  ggplot(aes(age, value, color = groups, shape = groups, linetype=groups)) + 
  geom_point(show.legend=FALSE) +
  geom_line(show.legend=T) +
  scale_shape_manual(values = c(NA,19,NA,NA),guide = guide_none()) +
  scale_linetype_manual(values = c(1,0,1,1),guide = guide_none()) +
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
data.frame(RTMB = report$srv_pred,
           RTMB.1 = report2$srv_pred) %>% 
  bind_cols(srv) %>% 
  rename(ADMB.1a = pred) -> ds

ds %>% 
  pivot_longer(c(RTMB, RTMB.1, ADMB.1a)) %>% 
  ggplot(aes(year, value, color = name)) + 
  geom_point(aes(y=biomass), color = 'black') +
  geom_line() +
  scico::scale_color_scico_d(name="",palette='roma') +
  expand_limits(y=0) +
  ylab('Biomass') +
  xlab('Years') +
  theme(legend.position = c(0.8, 0.2))


ggsave(here::here(2024, 'base_srv_like_iss', 'figs', "srvbio.png"), width=6.5, height=4.5, units='in')

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
         groups = 'RTMB') %>% 
  bind_rows(data.frame(age = data$ages, 
                       report2$fish_age_pred) %>% 
              pivot_longer(-age) %>% 
              mutate(year = rep(fish_age_yrs, length(ages)),
                     groups = 'RTMB.1')) %>% 
  bind_rows(fac) -> df5

df5 %>% 
  mutate(groups = ifelse(groups=='pred', 'Base.1a', groups)) %>% 
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


par_name = 'log_a50S'
par_values = seq(0, 2.5, length.out = 50)
par_like = single_likelihood(par_name, par_values, obj2, fit2)


par_name = 'log_a50C'
par_values = seq(0, 2.5, length.out = 50)
par_like2 = single_likelihood(par_name, par_values, obj2, fit2)

par_name = 'deltaC'
par_values = seq(3, 12, length.out = 50)
dc = single_likelihood(par_name, par_values, obj2, fit2)


par_name = 'deltaS'
par_values= seq(3, 12, length.out = 50)
ds = single_likelihood(par_name, par_values, obj2, fit2)



data.frame(pars = (par_values), 
           log_a50S = par_like,
           log_a50C = par_like2) %>% 
  pivot_longer(-pars) %>% 
ggplot(aes(pars, value, color = name)) + 
  geom_line()

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

report$slx[,2] - report2$slx[,2]



load(here::here('2024','sep_pt','figs','a50C.RData'), a50C <- new.env())
load(here::here('2024','sep_pt','figs','deltaC.RData'), deltaC <- new.env())
load(here::here('2024','sep_pt','figs','a50S.RData'), a50S <- new.env())
load(here::here('2024','sep_pt','figs','deltaS.RData'), deltaS <- new.env())
load(here::here('2024','sep_pt','figs','q.RData'), q <- new.env())

a50C$shinystan_density_gg$data %>% 
  ggplot(aes(x,y)) +
  geom_area(alpha = 0.2) +
  stat_function(fun = dnorm, args = list(mean = (data$mean_a50C), sd = (data$cv_a50C))) +
  expand_limits(x=c(5,10)) +
  ylab('') +
  xlab(expression('Fishery selectivity A'[50]))

deltaC$shinystan_density_gg$data %>% 
  ggplot(aes(x,y)) +
  geom_area(alpha = 0.2) +
  stat_function(fun = dnorm, args = list(mean = (data$mean_deltaC), sd = (data$cv_deltaC))) +
  expand_limits(x=c(1,6)) +
  ylab('') +
  xlab(expression('Fishery delta'))

a50S$shinystan_density_gg$data %>% 
  ggplot(aes(x,y)) +
  geom_area(alpha = 0.2) +
  stat_function(fun = dnorm, args = list(mean = (data$mean_a50C), sd = (data$cv_a50C))) +
  expand_limits(x=c(5,10)) +
  ylab('') +
  xlab(expression('Survey selectivity A'[50]))

deltaS$shinystan_density_gg$data %>% 
  ggplot(aes(x,y)) +
  geom_area(alpha = 0.2) +
  stat_function(fun = dnorm, args = list(mean = (data$mean_deltaC), sd = (data$cv_deltaC))) +
  # expand_limits(x=c(5,10)) +
  ylab('') +
  xlab(expression('Survey delta'))

q$shinystan_density_gg$data %>% 
  ggplot(aes(x,y)) +
  geom_area(alpha = 0.2) +
  stat_function(fun = dnorm, args = list(mean = (data$mean_q), sd = (data$cv_q))) +
  # expand_limits(x=c(5,2.5)) +
  ylab('') +
  xlab(expression('Survey delta'))

bayesplot::mcmc_areas(mcmc, pars = c('log_q', 'log_a50C'))
bayesplot::plot_density(shinystan_density_gg$data$y)
data$mean_q

par_name = 'log_a50C'
par_values = seq(0, 3, length.out = 50)
par_like = single_likelihood(par_name, par_values, obj, fit)
# Plot the profile
plot(par_values, par_like, type="l", xlab=par_name, ylab="Log-Likelihood")
