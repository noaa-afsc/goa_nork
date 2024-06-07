### Cecilia O'Leary, cecilia.oleary@noaa.gov; Madison Hall, madison.hall@noaa.gov
### code to produce model-based indices using VAST software implemented in R
### created by CO March 2020, updated by MBH April 2023
###2023 settings R v4.0.3: VAST v3.10.0, FishStatsUtils v2.12.0, cpp VAST_v14_0_1, TMB v1.9.2, Matrix v1.5-3, DHARMa 0.4.6

# edited ben.williams@noaa.gov 2024-04

# load ----
library(TMB)
library(VAST)
library(DHARMa)
library(FishStatsUtils)
library(Matrix)
library(here)
library(vroom)
library(dplyr)
library(ggplot2)
# data ----
Ndata <- vroom(here::here('2024', 'R', 'vast', "N_production_data_2023.csv")) %>%
  mutate(AreaSwept_km2 = net_width * distance_fished,
         Vessel = 1) %>% 
  select(Year = year, Catch_KG = weight, Lon = start_longitude, 
            Lat = start_latitude, Vessel, AreaSwept_km2) %>% 
    tidyr::replace_na(list(Catch_KG = 0)) %>% 
    tidyr::drop_na() %>% 
    filter(Year >= 1990)

Ndata %>% 
filter(Lon > -147) %>% 
group_by(Year) %>% 
summarise(catch = sum(Catch_KG))
saveRDS(Ndata, here::here('2024', 'R', 'vast', "Data_Geostat_Sebastes_polyspinis.RDS"))   

read.csv(here::here('2024', 'R', 'vast', "GOAThorsonGrid_Less700m.csv")) %>% 
    mutate(Area_km2 = Shape_Area / 1e6) %>% 
    select(Lat = Latitude, Lon = Longitude, Area_km2) -> input_grid

gc()

# Define settings
strata.limits <- data.frame('STRATA' = as.factor(c("Total","Western","Central")), 
                  'west_border' = c(-Inf, -Inf, -159), 
                  'east_border' = c(-147, -159, -147))

settings = make_settings( Version = "VAST_v14_0_1", #.cpp version, 
                          n_x = 750, # knots aka spatial resolution of our estimates
                          Region = "User", # Region = "gulf_of_alaska" , go to ?make_settings for other built in extrapolation grids
                          purpose = "index2", # changes default settings
                          ObsModel=c(1,1), # lognormal for BW 2nd run
                          # strata.limits = strata.limits,
                          knot_method = "grid",
                          fine_scale = TRUE,
                          bias.correct = TRUE,
                          use_anisotropy = TRUE
                      ) 

fit <- fit_model( "settings" = settings, #all of the settings we set up above
                  "Lat_i" = Ndata$Lat, #latitude of observation
                  "Lon_i" = Ndata$Lon,  #longitude of observation
                  "t_i" = Ndata$Year, #time for each observation
                  "b_i" = Ndata$Catch_KG, #in kg, raw catch or in CPUE per tow. raw catch used in 2023
                  "a_i" = Ndata$AreaSwept_km2, #sampled area for each observation
                  "v_i" = Ndata$Vessel, #set to 1 - no vessel effects
                  "input_grid" = input_grid, #only needed if you have a user input extrapolation grid, which we do for GOA
                  "optimize_args" = list("lower"=-Inf,"upper"=Inf), #TMB argument (?fit_tmb) tath can be used if you're having optimization issues, shouldn't need to change
                  "working_dir" = here::here('2024', 'R', 'vast', '1_1_2023')) #where to save the model output)
plot( fit )
saveRDS(fit, file = here::here('2024', 'R', 'vast', '1_1_2023', "VASTfit.RDS"))

# change it to lognormal 

settings = make_settings( Version = "VAST_v14_0_1", #.cpp version, 
                          n_x = 750, #knots aka spatial resolution of our estimates
                          Region = "User", #Region = "gulf_of_alaska" , go to ?make_settings for other built in extrapolation grids
                          purpose = "index2", #changes default settings
                          ObsModel=c(1,1), #lognormal for BW 2nd run
                          strata.limits = strata.limits,
                          knot_method = "grid",
                          fine_scale = TRUE,
                          bias.correct = TRUE,
                          use_anisotropy = TRUE
                      ) 

fit <- fit_model( "settings" = settings, #all of the settings we set up above
                  "Lat_i" = Ndata$Lat, #latitude of observation
                  "Lon_i" = Ndata$Lon,  #longitude of observation
                  "t_i" = Ndata$Year, #time for each observation
                  "b_i" = Ndata$Catch_KG, #in kg, raw catch or in CPUE per tow. raw catch used in 2023
                  "a_i" = Ndata$AreaSwept_km2, #sampled area for each observation
                  "v_i" = Ndata$Vessel, #set to 1 - no vessel effects
                  "input_grid" = input_grid, #only needed if you have a user input extrapolation grid, which we do for GOA
                  "optimize_args" = list("lower"=-Inf,"upper"=Inf), #TMB argument (?fit_tmb) tath can be used if you're having optimization issues, shouldn't need to change
                  "working_dir" = here::here('2024', 'R', 'vast', '1.1')) #where to save the model output)
plot( fit )
saveRDS(fit, file = here::here('2024', 'R', 'vast', "VASTfit.RDS"))

sb <- vroom::vroom(here::here('2023', "data", "output",  "goa_total_bts_biomass.csv")) %>% 
  mutate(area = 'Total', 
         type = 'db') %>% 
  select(year, area, biomass, type)

db <- vroom::vroom(here::here('2024', 'data', 'raw', 'goa_area_bts_biomass_data.csv')) %>% 
  select(year, area = regulatory_area_name, biomass = area_biomass, var = biomass_var) %>% 
  filter(year>=1990) %>% 
  group_by(year, area) %>% 
  summarise(biomass = sum(biomass), 
            var = sum(var)) %>% 
  mutate(area = case_when(area=='CENTRAL GOA' ~ 'Central',
                          area=='EASTERN GOA' ~ 'Eastern',
                          area=='WESTERN GOA' ~ 'Western'),
         type = 'db') %>% 
  select(year, area, biomass, type) %>% 
  bind_rows(sb)

vast <- read.csv(here::here('2024', 'R', 'vast', "1_1_2023_alt", "Index.csv")) %>% 
  rename_all(tolower) %>% 
  # filter(stratum!='Stratum_1') %>%
  mutate(area = case_when(stratum=='Stratum_2' ~ 'Western', 
                          stratum=='Stratum_1' ~ 'Total',
                          stratum=='Stratum_3' ~ 'Central'),
         type = 'vast') %>% 
  select(year = time, area, biomass = estimate, type) %>% 
  filter(biomass > 0)


  bind_rows(db, vast) %>% 
    mutate(area = factor(area, levels = c('Total', 'Western', 'Central', 'Eastern'))) %>% 
    ggplot(aes(year, biomass, color = type)) + 
    geom_line() +
    facet_wrap(~area, scales = 'free_y') +
    expand_limits(y=0) +
    scico::scale_color_scico_d(palette = "roma") 
    

# ba <- read.csv(here::here('2024', 'R', 'vast', "2_1_2023", "Index.csv")) %>% 
#   rename_all(tolower) %>% 
#   filter(stratum=='Stratum_1') %>% 
#   select(stratum, year = time, biomass = estimate, se = std..error.for.estimate) %>% 
#   filter(biomass>0)

read.csv(here::here('2024', 'R', 'vast', "1_1_2023_alt", "Index.csv")) %>% 
rename_all(tolower) %>% 
  filter(stratum!='Stratum_1') %>%
  mutate(area = ifelse(stratum=='Stratum_2', 'WESTERN GOA', 'CENTRAL GOA'))
select(stratum, year = time, biomass = estimate, se = std..error.for.estimate) %>% 
filter(biomass >0, stratum != 'Stratum_4') %>% 
ggplot(aes(year, biomass, color = stratum)) +
geom_line() + 
geom_line(data = sb %>% mutate(stratum = 'Stratum_1'), color = 1) +
geom_line(data = ba , color = 4) +
  expand_limits(y=0)


area_dat <- vroom::vroom(here::here('2024', 'data', 'raw', 'goa_area_bts_biomass_data.csv')) %>% 
  select(year, area = regulatory_area_name, biomass = area_biomass, var = biomass_var) %>% 
  filter(year>=1990) %>% 
  group_by(year, area) %>% 
  summarise(biomass = sum(biomass), 
            var = sum(var))

sb %>% 
  mutate(area = 'Total') %>% 
  bind_rows(area_dat) %>% 
  ggplot(aes(year, biomass, color = area)) +
  geom_line()

read.csv(here::here('2024', 'R', 'vast', "1_1_2023_alt", "Index.csv")) %>% 
  rename_all(tolower) %>% 
  # filter(stratum!='Stratum_1') %>%
  mutate(area = case_when(stratum=='Stratum_2' ~ 'Western', 
                          stratum=='Stratum_1' ~ 'Total',
                          stratum=='Stratum_3' ~ 'Central')) %>% 
  select(area, year = time, biomass = estimate, se = std..error.for.estimate) %>% 
  filter(biomass >0) %>% 
  ggplot(aes(year, biomass)) +
  geom_line() + 
  geom_point() +
  geom_line(data = area_dat, aes(year, biomass, color = area), color = 4) +
  facet_wrap(~area, scales='free_y') +
  expand_limits(y=0) +
  afscassess::theme_report()


library(rema)
# model ----
# dusky
yr = 2024
area_dat %>% 
  mutate(cv = sqrt(var)/biomass) %>% 
  select(strata = area, year, biomass, cv) %>% 
  ungroup() -> db

read.csv(here::here('2024', 'R', 'vast', "1_1_2023_alt", "Index.csv")) %>% 
  rename_all(tolower) %>% 
  filter(stratum!='Stratum_1') %>%
  mutate(area = ifelse(stratum=='Stratum_2', 'WESTERN GOA', 'CENTRAL GOA')) %>% 
  select(strata = area, year = time, biomass = estimate, se = std..error.for.estimate) %>% 
  mutate(cv = se / biomass) %>% 
  filter(biomass >0) %>% 
  select(-se) -> vast

input <- prepare_rema_input(model_name = 'db',
                            biomass_dat = db,
                            end_year = yr,
                            # how do you deal with zero biomass observations?
                            # see ?prepare_rema_input() for more options
                            zeros = list(assumption = 'NA'))
m <- fit_rema(input)
out <- tidy_rema(m)
vast %>% 
  tidyr::pivot_wider(names_from=strata, values_from=biomass, -cv)

inputv <- prepare_rema_input(model_name = 'vast',
                            biomass_dat = vast,
                            end_year = yr,
                            # how do you deal with zero biomass observations?
                            # see ?prepare_rema_input() for more options
                            zeros = list(assumption = 'NA'))
mv <- fit_rema(inputv)
outv <- tidy_rema(mv)

outv$proportion_biomass_by_strata %>% 
  vroom_write(here::here(yr, "results", "dusky_ratios.csv"), delim=",")


png(filename=here::here("figs", "dusky_re.png"), width = 6.5, height = 6.5, 
    units = "in", type ="cairo", res = 200)


out$proportion_biomass_by_strata %>% 
  tidyr::pivot_longer(-c(model_name, year)) %>% 
  bind_rows(
    outv$proportion_biomass_by_strata %>% 
      tidyr::pivot_longer(-c(model_name, year)) 
  ) %>% 
  ggplot(aes(year, value, color = name, group = interaction(name, model_name), lty= model_name)) + 
  geom_line() + 
  # facet_wrap(~name) +
    scico::scale_color_scico_d(palette = 'roma')
  
  # ggplot(aes(year, obs)) + 
  # geom_point() + 
  geom_line(aes(x=year, y=pred), color = 4) +
  # facet_wrap(~strata) +
  # geom_point(data = outv$biomass_by_strata, aes(year, obs), color = 2) +
 geom_line(data = outv$biomass_by_strata, aes(year, pred), color = 2) +
    scico::scale_color_scico_d(palette = 'roma')

out$biomass_by_strata %>% 
  select(year, strata, pred, model_name) %>% 
  bind_rows(outv$biomass_by_strata %>% 
              select(year, strata, pred, model_name)) %>% 
  mutate(strata = case_when(strata=="EASTERN GOA" ~ "Eastern",
                            strata=="CENTRAL GOA" ~ "Central",
                            strata=="WESTERN GOA" ~ "Western"),
         strata = factor(strata, levels = c("Western", "Central", "Eastern"))) %>% 
  group_by(model_name, year) %>% 
  mutate(allocation = pred / sum(pred)) -> d2

out$proportion_biomass_by_strata %>% 
  # tidyr::pivot_longer(-c(model_name, year)) %>% 
  bind_rows(
    outv$proportion_biomass_by_strata 
      # tidyr::pivot_longer(-c(model_name, year)) 
  ) %>% arrange(year, model_name) %>% 
  select(-`EASTERN GOA`) %>% 
  tail() %>% 
  tidyr::pivot_longer(-c(model_name, year)) %>% 
  tidyr::pivot_wider(names_from = model_name, 
                     values_from = value) %>% 
  arrange(year, name)

  d2 %>% 
    ggplot(aes(year, pred, color = strata, group = interaction(strata, model_name), lty= model_name)) +
  geom_line() +
  # facet_wrap(~strata) +
  scico::scale_color_scico_d(palette = "roma", end =0.6) 
  
  d2 %>% 
    ggplot(aes(year, allocation, color = model_name, group = interaction(strata, model_name), lty = strata)) +
    geom_line() +
    # facet_wrap(~strata) +
    scico::scale_color_scico_d(palette = "roma") 

  out$biomass_by_strata %>% 
    ggplot() + 
    # geom_point() + 
    geom_line(aes(x = year, y=pred), color = 4) +
    facet_wrap(~strata) +
    # geom_point(data = outv$biomass_by_strata, aes(year, obs), color = 2) +
    geom_line(data = outv$biomass_by_strata, aes(year, pred), color = 2)
  
  
  
  
  filter(year>=1990) %>% 
  mutate(strata = case_when(strata=="EASTERN GOA" ~ "Eastern",
                            strata=="CENTRAL GOA" ~ "Central",
                            strata=="WESTERN GOA" ~ "Western"),
         strata = factor(strata, levels = c("Western", "Central", "Eastern"))) %>% 
  ggplot(aes(year, pred)) + 
  geom_point(aes(y = obs), color = "darkgray") + 
  geom_errorbar(aes(ymin = obs_lci, ymax = obs_uci), color = "darkgray", width = 0.8) +
  geom_line() + 
  geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci), alpha = 0.3) +
  facet_wrap(~strata, ncol = 1, scales = "free_y") +
  ylab("Biomass (t)") +
  xlab("Year") + 
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(breaks = seq(1990,2020,5))

dev.off()
