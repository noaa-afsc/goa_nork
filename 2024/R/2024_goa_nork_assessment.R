library(afscdata)
library(afscassess)
library(ggplot2)
library(dplyr)
theme_set(theme_report())

# globals ----
year = 2024
species = "NORK"
area = "goa"
TAC <- c(4312, 5358, 5147)
rec_age = 2
plus_age = 45
norpac_species = 303
# setup 
# setup_folders(year)
# accepted_model(2022, "m22.1", 2023)

# query data ----
goa_nork(year)
c("year", "performance", "specimen_type", "join_key", "haul_join", "port_join",
  "species", "fmp_gear", "fmp_area", "fmp_subarea", 
  "age", "length", "weight",
  tolower(add_fields))
db = connect()


q_fsh_specimen(year, species = norpac_species, area, db, add_fields = c('gear', 'haul_offload_date'), save = F) -> dat

afscassess::bts_biomass(year, area=area, type='area')

dat %>% 
  tidytable::filter(age>=rec_age, !is.na(length), !is.na(performance)) %>%
  tidytable::mutate(age = ifelse(age>plus_age, plus_age, age)) %>%
  tidytable::mutate(tot = tidytable::n(), .by = year) %>%
  tidytable::filter(tot>49) %>%
  tidytable::mutate(n_h = length(unique(na.omit(haul_join))) +
                      length(unique(na.omit(port_join))),
                    .by = year) %>%
  tidytable::summarise(n_s = mean(tot),
                       n_h = mean(n_h),
                       age_tot = tidytable::n(),
                       .by = c(year, age)) %>%
  tidytable::mutate(prop = age_tot / n_s) %>%
  tidytable::left_join(expand.grid(year = unique(.$year),
                                   age = rec_age:plus_age), .) %>%
  tidytable::replace_na(list(prop = 0)) %>%
  tidytable::mutate(AA_Index = 1,
                    n_s = mean(n_s, na.rm = T),
                    n_h = mean(n_h, na.rm = T),
                    .by = year) %>%
  tidytable::select(-age_tot) %>%
  tidytable::pivot_wider(names_from = age, values_from = prop) -> fac
