# 2024 GOA northern SAFE figures
# switching from ADMB to RTMB so figure functions will not work without some tweaks
# ben.williams@noaa.gov

# load ----
library(afscassess)
library(tidyverse)
theme_set(theme_report())
library(vroom)
library(here)
library(patchwork)
# globals ----

year = 2024
folder = 'm24'
source(here::here(year, 'r', 'utils.r'))

data = readRDS(here::here(year, folder, 'dat.rds'))
base = readRDS(here::here(year, 'base', 'base.rds'))
m22.1 = readRDS(here::here(year, 'm22.1', 'm22.1.rds'))
m22.1a = readRDS(here::here(year, 'm22.1a', 'm22.1a.rds'))
m22.1b = readRDS(here::here(year, 'm22.1b', 'm22.1b.rds'))
m22.1c = readRDS(here::here(year, 'm22.1c', 'm22.1c.rds'))
m24 = readRDS(here::here(year, folder, 'm24.rds'))
m24_se = readRDS(here::here(year, folder, 'm24_se.rds'))
mcmc = readRDS(file = here::here(2024, 'm24','mcmc.rds'))
post = adnuts::extract_samples(mcmc)
mpost = as.matrix(post)
prior_ssb = vroom(here(year, 'data', 'user_input', 'prior_sp_bio.csv')) 
r = list()
# process posterior
for(i in 1:nrow(mpost)){
  r[[i]] <- obj24$report(mpost[i,])
}

pb = purrr::map(r, proj_bio)
purrr::map(pb, item = 'catch_abc') %>% 
  purrr::map(., ~ as.data.frame(t(.))) %>% 
  bind_rows(.id = 'sim') %>% 
  pull(V1) -> abc_hist


purrr::map(r, 'Ft') %>% 
  purrr::map(., ~ as.data.frame(t(.))) %>% 
  bind_rows(.id = 'sim') -> fp 



purrr::map(r, 'spawn_bio') %>% 
  purrr::map(., ~ as.data.frame(t(.))) %>% 
  bind_rows(.id = 'sim') -> sp 

sp[,ncol(sp)] -> sp_hist

purrr::map(r, 'tot_bio') %>% 
  purrr::map(., ~ as.data.frame(t(.))) %>% 
  bind_rows(.id = 'sim') -> tp 

tp[,ncol(tp)] -> tp_hist

purrr::map(r, 'recruits') %>% 
  purrr::map(., ~ as.data.frame(t(.))) %>% 
  bind_rows(.id = 'sim') -> recs 

# catch ---- 
png(filename=here::here("2024", "m24", "figs", "fcatch.png"), width = 6.5, height = 4.5,
    units = "in", type ="cairo", res = 200)

data.frame(year = dat$years,
           observed = dat$catch_obs,
           m24 = m24$catch_pred) %>% 
  tidyr::pivot_longer(-year) %>%
  ggplot(aes(year, value, color = name, linetype = name, shape = name)) + 
  geom_point() +
  geom_line() +
  scale_y_continuous(labels = scales::comma) +
  scale_linetype_manual("", values = c(1,0)) +
  scale_shape_manual("", values = c(NA,19)) +
  afscassess::scale_x_tickr('Year', data=m24_se, var=year) +
  scico::scale_color_scico_d("", palette = 'grayC', end=0.3) +
  ylab('Catch (t)') +
  theme(legend.position = c(x=0.8, y=0.8)) 

dev.off()


vroom::vroom(here::here(2024, 'data', 'raw', 'fish_catch_data.csv')) %>% 
  mutate(day = lubridate::yday(week_end_date)) %>% 
  group_by(year, day) %>% 
  summarise(catch = sum(weight_posted)) %>% 
  group_by(year) %>% 
  mutate(cumsum = cumsum(catch)) %>% 
  ungroup() -> df
df %>% 
  filter(year==max(year)) %>% 
  filter(day==max(day)) -> pt


png(filename=here::here("2024", "m24", "figs", "fcatch1.png"), width = 6.5, height = 6.5,
    units = "in", type ="cairo", res = 200)  
df %>% 
  ggplot(aes(day, cumsum, color = year, group = year)) + 
  geom_line() +
  geom_point(data=pt, size = 2) +
  scico::scale_color_scico("Year", palette = 'roma') +
  scale_y_continuous(labels=scales::comma) +
  xlab('Julian Day') +
  ylab('Cumulative catch (t)') +
  theme(legend.position = c(x=0.2, y=0.8)) 
dev.off()  


# past model ssb ----
png(filename=here::here("2024", "m24", "figs", "past_ssb.png"), width = 6.5, height = 5.5,
    units = "in", type ="cairo", res = 200)
prior_ssb %>% 
  tidyr::pivot_longer(-year) %>% 
  ggplot(aes(year, value, color = name)) +
  geom_line() +
  scico::scale_color_scico_d("", palette='roma') +
  expand_limits(y=0) +
  theme(legend.position = c(x=0.8, y=0.2)) +
  ylab('Female spawning biomass (t)') +
  scale_y_continuous(labels=scales::comma) +
  afscassess::scale_x_tickr(data=prior_ssb, var=year, to=10, start=1970)
dev.off()  

# histograms ----
p1 = plot_par(item = "log_q", mcmc=mcmc) +
  xlab(expression("Log trawl survey catchability ("*italic(q)*")"))
p2 = plot_par(item = "log_M", mcmc=mcmc) +
  xlab(expression("Log natural mortality ("*italic(M)*")")) +
  ylab("Probability density") 
p3 = plot_par(item = "log_F40", mcmc=mcmc) +
  xlab(expression("Log "*italic(F)["40%"]*""))
p4 = plot_par(item=abc_hist, report=proj_bio(m24), rep_item='catch_abc', proj=TRUE)  +
  xlab("ABC (t)") + scale_x_continuous(labels = scales::comma)
p5 = plot_par(item=sp_hist, report=m24, rep_item='spawn_bio')  +
  scale_x_continuous(labels=scales::comma) +
  xlab('Current spawning biomass (t)')
p6 = plot_par(item=tp_hist, report=m24, rep_item='tot_bio') +
  scale_x_continuous(labels=scales::comma) +
  xlab('Current total biomass (t)')

png(filename=here::here(year, folder, "figs", "hists.png"), width = 6.5, height = 6.5,
    units = "in", type ="cairo", res = 200)
print((cowplot::plot_grid(p1, p4, p2,  p5, p3, p6, align = "v", ncol = 2, rel_heights = c(0.5, 0.5))))
dev.off()


# biomass estimates ----
base_se %>% 
  filter(item %in% c('spawn_bio', 'tot_bio')) %>% 
  mutate(name = case_when(item == 'spawn_bio' ~ 'Spawning biomass',
                          TRUE ~ 'Total biomass'),
         model = "base") %>% 
  dplyr::select(year, name, biomass=value, model) -> base_df

m22.1_se %>% 
  filter(item %in% c('spawn_bio', 'tot_bio')) %>% 
  mutate(name = case_when(item == 'spawn_bio' ~ 'Spawning biomass',
                          TRUE ~ 'Total biomass'),
         model = "m22.1") %>% 
  dplyr::select(year, name, biomass=value, model) -> m22.1_df
           
sp %>% 
  dplyr::rename_with(~ as.character(data$years), -sim) %>% 
  tidyr::pivot_longer(-sim) %>% 
  mutate(year = as.numeric(name),
         sim = as.numeric(sim),
         name = 'Spawning biomass',
         model = 'm24') %>% 
  bind_rows(tp %>% 
  dplyr::rename_with(~ as.character(data$years), -sim) %>% 
  tidyr::pivot_longer(-sim) %>% 
  mutate(year = as.numeric(name),
         sim = as.numeric(sim),
         name = 'Total biomass',
         model = 'm24')) %>% 
  tidytable::summarise(median = median(value),
                       lci = quantile(value, 0.025),
                       uci = quantile(value, 0.975),
                       .by = c(year, name, model)) %>% 
  left_join(data.frame(year = dat$years,
                       tot = m24$tot_bio,
                       bio = m24$spawn_bio)) %>%
  tidytable::mutate(biomass = ifelse(name == "Total biomass", tot, bio)) %>%
  tidytable::select(-tot, -bio) %>% 
  bind_rows(m22.1_df, base_df) -> df

png(filename=here::here(year, folder, "figs", "est_biomass.png"), width = 6.5, height = 6.5,
    units = "in", type ="cairo", res = 200)
df %>% 
  ggplot2::ggplot(ggplot2::aes(year, biomass, color = model, fill = model)) +
  ggplot2::geom_line() +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lci, ymax = uci), alpha = 0.1, color = NA) +
  ggplot2::facet_wrap(~name, dir = "v", scales = "free_y") +
  ggplot2::scale_y_continuous(name = "Biomass (t)", labels = scales::comma) +
  ggplot2::expand_limits(y = 0) +
  scico::scale_color_scico_d(palette = 'roma') +
  afscassess::scale_x_tickr(name = "Year", data=df, var=year, to=10, start = 1960) +
  theme(legend.position = c(x=0.2, y=0.8)) 
  
dev.off()

# fishing mortality ----

png(filename=here::here(year, folder, "figs", "Ft.png"), width = 6.5, height = 5.5,
    units = "in", type ="cairo", res = 200)
fp %>% 
  dplyr::rename_with(~ as.character(data$years), -sim) %>% 
  tidyr::pivot_longer(-sim) %>% 
  mutate(year = as.numeric(name),
         sim = as.numeric(sim)) %>% 
  select(-name) %>% 
  tidytable::summarise(lci = quantile(value, 0.025),
                       uci = quantile(value, 0.975),
                       .by = c(year)) %>% 
  left_join(data.frame(year = dat$years,
                       M24 = m24$Ft,
                       m22.1 = m22.1$Ft)) %>% 
  left_join(data.frame(year = head(dat$years, -2),
                       base = base$Ft)) %>% 
  tidyr::pivot_longer(-c(year, lci, uci)) %>% 
  ggplot(aes(year, value, color = name, fill = name)) + 
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.1, color = NA) +
  geom_line() + 
  ggplot2::expand_limits(y = 0) +
  scico::scale_color_scico_d('model', palette = 'roma') +
  scico::scale_fill_scico_d('model', palette = 'roma') +
  afscassess::scale_x_tickr(name = "Year", data=df, var=year, to=10, start = 1960) +
  theme(legend.position = c(x=0.8, y=0.8)) +
  ggplot2::ylab("Fishing mortality rate (F)\n")
dev.off()





# plot catch/biomass ----
catch <- vroom::vroom(here::here(year, "data", "output", "fish_catch.csv"))
pc = vroom::vroom(here::here(year, 'data', 'output', 'yld_rat.csv')) %>% 
  select(year, catch = proj_catch)
catch %>% 
  mutate(catch = ifelse(year==pc$year, pc$catch, catch)) -> catch 

std = readRDS(here::here(year, "m24", "m24_se.rds")) 


std %>% 
filter(item=='tot_bio') %>% 
  left_join(catch) %>% 
  filter(year >= 1991) %>%
  mutate(perc = catch / value,
         lci = catch / lci,
         uci = catch / uci,
         mean = mean(perc)) %>%
  dplyr::select(year, value, mean, perc, lci, uci) -> df

png(filename=here::here(year, "m24", "figs", "catch_bio.png"), width = 6.5, height = 6.5,
    units = "in", type ="cairo", res = 200)
df %>%
  ggplot2::ggplot(ggplot2::aes(year, perc)) +
  ggplot2::geom_line() +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lci, ymax = uci), alpha = 0.2) +
  ggplot2::geom_hline(yintercept = df$mean, lty = 3) +
  ggplot2::expand_limits(y = c(0, 0.08)) +
  afscassess::scale_x_tickr(data=df, var=year, start = 1990) +
  afscassess::theme_report() +
  ggplot2::xlab("\nYear") +
  ggplot2::ylab("Catch/Biomass\n")


dev.off()


# recruits ----
png(filename=here::here(year, folder, "figs", "recruits.png"), width = 6.5, height = 5.5,
    units = "in", type ="cairo", res = 200)

recs %>% 
  dplyr::rename_with(~ as.character(data$years), -sim) %>% 
  tidyr::pivot_longer(-sim) %>% 
  mutate(year = as.numeric(name),
         sim = as.numeric(sim)) %>% 
  select(-name) %>% 
  tidytable::summarise(lci = quantile(value, 0.025),
                       uci = quantile(value, 0.975),
                       .by = c(year)) %>% 
  left_join(data.frame(year = dat$years,
                       M24 = m24$recruits,
                       m22.1 = m22.1$recruits)) %>% 
  left_join(data.frame(year = head(dat$years, -2),
                       base = base$recruits)) %>% 
  tidyr::pivot_longer(-c(year, lci, uci)) %>% 
  ggplot(aes(year, value, color = name, fill = name)) + 
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.1, color = NA) +
  geom_line() + 
  ggplot2::expand_limits(y = 0) +
  scico::scale_color_scico_d('model', palette = 'roma') +
  scico::scale_fill_scico_d('model', palette = 'roma') +
  afscassess::scale_x_tickr(name = "Year", data=df, var=year, to=10, start = 1960) +
  theme(legend.position = c(x=0.8, y=0.8)) +
  ggplot2::ylab("Age-2+ Recruitment (millions)")

dev.off()

# spawn-recr ----
png(filename=here::here(year, folder, "figs", "recr-ssb.png"), width = 6.5, height = 5.5,
    units = "in", type ="cairo", res = 200)

data.frame(year = head(m24$years,-2),
            rec = tail(m24$recruits,-2),
            bio = head(m24$spawn_bio, -2)) %>% 
  mutate( label = stringr::str_sub(year, 3),
          decade = (floor(year / 10) * 10)) %>% 
  ggplot2::ggplot(ggplot2::aes(bio, rec, color = decade)) +
  ggplot2::geom_label(ggplot2::aes(label=label, color = decade),
                      label.size = 0, show.legend = FALSE, size = 3, family="Times", alpha = 0.85) +
  ggplot2::expand_limits(x = 0, y = 0) +
  scico::scale_color_scico(palette = "roma") +
  geom_path(show.legend = FALSE) +
  scale_x_continuous(labels = scales::comma) +
  ggplot2::xlab("SSB (kt)") +
  ggplot2::ylab("Recruitment (millions)")

dev.off()

           
# selectivity ----

data.frame(age = 2:51,
           Fishery = m24$slx[,1],
           Survey = m24$slx[,2],
           `Original maturity` = m22.1$maa,
           `Updated maturity` = m24$maa) %>% 
  tidytable::pivot_longer(-age) %>% 
  mutate(name = factor(name, levels = c('Fishery', 'Survey', 'Original.maturity', 'Updated.maturity')))  -> slx
  
png(filename=here::here(year, folder, "figs", "selex.png"), width = 6.5, height = 5.5,
    units = "in", type ="cairo", res = 200)

slx %>% 
  ggplot2::ggplot(ggplot2::aes(age, value, color = name)) +
  ggplot2::geom_line() +
  scico::scale_color_scico_d(name = "", palette = 'roma', begin = 0.0)+
  ggplot2::ylab("Selectivity/Maturity\n") +
  afscassess::scale_x_tickr(name = "Age", data=slx, var=age, to=10, start = 0) +
  ggplot2::theme(legend.justification=c(1,0),
                 legend.position=c(0.9,0.2)) +
  coord_cartesian(xlim = c(0, 45))

dev.off()


pairs_admb(mcmc, pars=1:5, order='slow')
pairs_admb(mcmc, pars=1:5, order='mismatch')
plot_marginals(mcmc, pars=1:6)
obj24$mle

posterior <- extract_samples(mcmc, inc_lp = FALSE)
par.names <- names(posterior)
mcmc$mle$est[1]


hist(posterior[, 1], breaks = 30, yaxs = "i", 
     freq = FALSE, col = gray(0.8))
lines(x1, y1, col = "red", lwd = 2)


mle <- (filter(m24_se, item=='spawn_bio') %>% tail(1))$value
se <- (filter(m24_se, item=='spawn_bio') %>% tail(1))$se
x1 <- seq(qnorm(0.001, mle, se), qnorm(0.999, mle, 
                                       se), len = 100)
y1 <- dnorm(x1, mle, se)


# phase plane ----
Fabc <- 0.35/0.4

segs = data.frame(x1 = rep(c(0.05, 0.4/0.35), 2),
                  x2 = rep(c(0.4/0.35, 2.8), 2),
                  y1 = c(0, 1, 0, Fabc),
                  y2 = c(1, 1, Fabc, Fabc),
                  group = factor(c("ofl", "ofl", "abc", "abc"),
                                 levels = c("ofl", "abc")))


png(filename=here::here(year, folder, "figs", "phase_plane.png"), width = 6.5, height = 5.5,
    units = "in", type ="cairo", res = 200)

data.frame(year = min(m24$years):(max(m24$years)+2),
           x = c(m24$spawn_bio, proj_bio(m24)$spawn_bio[1],  
                 proj_bio(m24)$spawn_bio[2]) / m24$B35,
           y = c(m24$Ft, proj_bio(m24)$F40 * yld$yld) / proj_bio(m24)$F35[1]) %>% 
  tidytable::mutate(label = stringr::str_sub(year, 3),
                    decade = (floor(year / 10) * 10)) %>% 
  ggplot2::ggplot(ggplot2::aes(x, y)) +
  geom_path(aes(color = decade), show.legend = FALSE) +
  ggplot2::geom_label(ggplot2::aes(label=label, color = decade), label.size = 0,
                      show.legend = FALSE, size = 3, family="Times", alpha = 0.5) +
  ggplot2::geom_segment(data = segs, ggplot2::aes(x1, y1, xend = x2, yend = y2, linetype = group)) +
  ggplot2::scale_linetype_manual(values = c(1, 3),
                                 labels = c(expression(italic(F[OFL])),
                                            expression(italic(F[ABC]))),
                                 name = "") +
  scico::scale_color_scico(palette = "roma") +
  ggplot2::ylab(expression(italic(F/F["35%"]))) +
  ggplot2::xlab(expression(italic(SSB/B["35%"]))) +
  ggplot2::theme(legend.justification=c(1,0),
                 legend.position=c(0.9,0.85)) 


dev.off()


# comps ----
fsc = resids(obs = data$fish_size_obs, pred = m24$fish_size_pred, yrs = data$fish_size_yrs, iss=data$fish_size_iss, ind = length_bins, label = 'Length (cm)')
fac = resids(obs = data$fish_age_obs, pred = m24$fish_age_pred, yrs = data$fish_age_yrs, iss=data$fish_age_iss, ind = ages, label = 'Age')
sac = resids(obs = data$srv_age_obs, pred = m24$srv_age_pred, yrs = data$srv_age_yrs, iss=data$srv_age_iss, ind = ages, label = 'Age')

# fishery size comp
png(filename=here::here(year, folder, "figs", "fsc_resid.png"), width = 6.5, height = 5.5,
    units = "in", type ="cairo", res = 200)

(fsc$pearson +
    fsc$osa) /
  (fsc$agg + fsc$ss) + plot_layout(guides = "collect") 
dev.off()

png(filename=here::here(year, folder, "figs", "fsc.png"), width = 6.5, height = 5.5,
    units = "in", type ="cairo", res = 200)
fsc$annual
dev.off()

# fishery age comp
png(filename=here::here(year, folder, "figs", "fac_resid.png"), width = 6.5, height = 5.5,
    units = "in", type ="cairo", res = 200)

(fac$pearson +
    fac$osa) /
  (fac$agg  + fac$ss) + plot_layout(guides = "collect") 
dev.off()

png(filename=here::here(year, folder, "figs", "fac.png"), width = 6.5, height = 5.5,
    units = "in", type ="cairo", res = 200)
fac$annual
dev.off()

# survey age comp
png(filename=here::here(year, folder, "figs", "sac_resid.png"), width = 6.5, height = 5.5,
    units = "in", type ="cairo", res = 200)

(sac$pearson +
    sac$osa) /
  (sac$agg  + sac$ss) + plot_layout(guides = "collect") 
dev.off()

png(filename=here::here(year, folder, "figs", "sac.png"), width = 6.5, height = 5.5,
    units = "in", type ="cairo", res = 200)
sac$annual
dev.off()


# comp bubbles ----

obs = as.data.frame(data$fish_size_obs)
names(obs) <- data$fish_size_yrs
obs %>% 
  mutate(length = 15:45) %>% 
  tidyr::pivot_longer(-length) %>% 
  mutate(year = as.numeric(name)) %>% 
  ggplot(aes(length, year, size = value)) + 
  geom_point() +
  scale_size_area()

obs = as.data.frame(data$fish_age_obs)
names(obs) <- data$fish_age_yrs
obs %>% 
  mutate(age = 2:45, 
         id = 'fishery') %>% 
  tidyr::pivot_longer(-c(age, id)) %>% 
  mutate(year = as.numeric(name)) -> obs

obs %>% 
ggplot(aes(age, year, size = value)) + 
  geom_point() +
  scale_size_area()

obs1 = as.data.frame(data$srv_age_obs)
names(obs1) <- data$srv_age_yrs
obs1 %>% 
  mutate(age = 2:45, 
         id = 'survey') %>% 
  tidyr::pivot_longer(-c(age, id)) %>% 
  mutate(year = as.numeric(name)) -> obs1 

png(filename=here::here(year, "m24", "figs", "age_comp_in.png"), width = 6.5, height = 5.5,
    units = "in", type ="cairo", res = 200)
obs %>% 
  bind_rows(obs1) %>% 
  ggplot(aes(year, age, size = value, color = id)) + 
  geom_point() +
  scale_size_area() +
  scico::scale_color_scico_d(palette = 'roma', begin = 0.1) +
  ylab('Age') +
  xlab('Year')
dev.off()


srv_size <- read_csv("2024/data/output/goa_bts_sizecomp.csv")


obs = as.data.frame(data$fish_size_obs)
names(obs) <- data$fish_size_yrs
obs %>% 
  mutate(length = 15:45, 
         id = 'fishery') %>% 
  tidyr::pivot_longer(-c(length, id)) %>% 
  mutate(year = as.numeric(name)) -> obs


obs1 = as.data.frame(t(srv_size[-c(1:4)]))
names(obs1) <- data$srv_yrs
obs1 %>% 
  mutate(length = 15:45, 
         id = 'survey') %>% 
  tidyr::pivot_longer(-c(length, id)) %>% 
  mutate(year = as.numeric(name)) -> obs1 

png(filename=here::here(year, "m24", "figs", "size_comp_in.png"), width = 6.5, height = 5.5,
    units = "in", type ="cairo", res = 200)
obs %>% 
  bind_rows(obs1) %>% 
  ggplot(aes(year, length, size = value, color = id)) + 
  geom_point(alpha = 0.6) +
  scale_size_area() +
  scico::scale_color_scico_d(palette = 'roma', begin = 0.1) +
  ylab('Length') +
  xlab('Year')
dev.off()

# survey inputs ----
db = vroom::vroom(here::here('2024', 'data', 'output', 'goa_total_bts_biomass.csv')) %>% 
  mutate(id = 'design-based',
         srv = '2023')
vg = vroom::vroom(here::here('2022', 'data', 'output', 'goa_ts_biomass_m22.1.csv')) %>% 
  mutate(id = 'VAST - Gamma',
         srv = '2021')
vg24 = vroom::vroom(here::here('2024', 'data', 'user_input', 'vast_2_1.csv')) %>% 
  mutate(lci = biomass - se * 1.96,
         uci = biomass + se * 1.96,
         id = 'VAST - Gamma',
         srv = '2023')
vl = vroom::vroom(here::here('2022', 'data', 'user_input', 'vast_lognormal.csv')) %>% 
  mutate(id = 'VAST - lognormal',
         srv = '2021')
vl24 = vroom::vroom(here::here('2024', 'data', 'output', 'goa_total_bts_biomass_vast.csv')) %>% 
  mutate(id = 'VAST - lognormal',
         srv = '2023')


png(filename=here::here(year, "compare_models", "surveys.png"), width = 6.5, height = 5.5,
    units = "in", type ="cairo", res = 200)
bind_rows(db, vg, vg24, vl, vl24) %>% 
  ggplot(aes(year, biomass, color = srv, group = srv)) + 
  geom_point(position=position_dodge(width=0.5)) + 
  geom_line(alpha = 0.2) +
  geom_errorbar(aes(ymin=lci, ymax=uci), alpha=0.3, width = 0.2, 
                position=position_dodge(width=0.5)) +
  scico::scale_color_scico_d("survey year", palette = 'roma') +
  expand_limits(y=0) +
  scale_y_continuous(labels=scales::comma) +
  ylab('Survey biomass (t)') +
  xlab('Year') +
  facet_wrap(~id) +
  ggplot2::theme(legend.justification=c(1,0),
                 legend.position=c(0.9,0.85)) 
dev.off()



data.frame(
  year = data$srv_yrs,
  observed = data$srv_obs,
  sd = data$srv_sd) %>% 
  mutate(olci = observed - 1.96 * sd,
         ouci = observed + 1.96 * sd) %>% 
  left_join(m24_se %>% 
              filter(item=='srv_pred')) %>% 
  rename(predicted = value) %>% 
  pivot_longer(c(observed, predicted)) -> df
  

png(filename=here::here(year, "m24", "figs", "srv_biomass.png"), width = 6.5, height = 5.5,
    units = "in", type ="cairo", res = 200)

df %>% 
  ggplot(aes(year, value, color = name, linetype = name)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.1, color = NA,
              show.legend = FALSE) +
  geom_errorbar(aes(ymin = olci, ymax = ouci), width = 0.2, color = 'gray') +
  geom_point() +
  geom_line() +
  scale_color_manual("", values = c(NA, 1)) +
  scale_linetype_manual("", values = c(0,1)) +
  expand_limits(y = 0) + 
  afscassess::scale_x_tickr(data = df, var = year)  +
  scale_y_continuous(labels = scales::comma) +
  ylab('Survey biomass (t)') +
  xlab('Year') +
  theme(legend.position = c(0.8, 0.85))

dev.off()

png(filename=here::here(year, "compare_models", "surveys_out.png"), width = 6.5, height = 5.5,
    units = "in", type ="cairo", res = 200)
data.frame(year = dat$srv_yrs,
           m24.lognormal = m24$srv_pred,
           m22.1c.lognormal = m22.1c$srv_pred,
           m22.1b.Gamma = m22.1b$srv_pred,
           m22.1a.Gamma = m22.1a$srv_pred,
           m22.1.Gamma = m22.1$srv_pred,
           base.Gamma = c(base$srv_pred, NA)) %>% 
  tidyr::pivot_longer(-year) %>% 
  ggplot(aes(year, value, color = name)) + 
  geom_line() +
  scico::scale_color_scico_d("model", palette = 'roma') +
  expand_limits(y = 0) +
  scale_y_continuous(labels = scales::comma) +
  ylab('Survey biomass (t)') +
  xlab('Year') +
  theme(legend.position = c(0.2, 0.25))
dev.off()
           
# ISS ----
png(filename=here::here(year, "m24", "figs", "iss.png"), width = 3.5, height = 3.5,
    units = "in", type ="cairo", res = 200)
data.frame(year = c(1990, 1993, 1996, 1999, 2001, 2003, 2005, 2007, 2009, 2011, 2013, 2015, 2017, 2019, 2021),
           hybrid = c(28.2956716,28.79701882,42.0642524,38.89731882,89.35478366,30.94949752,77.79463763,100,95.15473983,80.08763593,82.37064037,72.44958071,86.31402473,68.90159616,83.77314565),
           iss = c(52, 32, 90, 42, 88, 26, 66, 126, 168, 91, 104, 117, 118, 97, 153)) %>% 
  tidyr::pivot_longer(-year) %>% 
  ggplot(aes(year, value, color = name)) + 
  geom_point() + 
  afscassess::theme_report() + 
  ylab('Sample size') +
  xlab('Year') +
  scico::scale_color_scico_d("", palette = 'roma', end = 0.8) +
  theme(legend.position = c(0.2, 0.75)) +
  expand_limits(y=0)
dev.off()           
