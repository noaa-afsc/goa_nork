# 2023 goa northern rockfish assessment code
# ben.williams@noaa.gov

library(afscdata)
library(afscassess)
library(ggplot2)
library(dplyr)
theme_set(theme_report())

# globals ----
year = 2023
species = "NORK"
area = "goa"
TAC <- c(4312, 5358, 5147)
rec_age = 2

# setup 
# setup_folders(year)
# accepted_model(2022, "m22.1", 2023)

# query data ----
goa_dusk(year, off_yr = TRUE)
clean_catch(year, species=species, TAC=TAC)
bts_biomass(year=year, area=area, rmv_yrs=c(1984,1987)) # design-based
bts_biomass(year=year, area=area, file="vast_lognormal.csv", id="vast") # 2021 vast model

# run projection
proj_ak(year=year, last_full_assess = 2022, species="northern", region=area,
        rec_age=rec_age, folder="update", off_yr=TRUE)

# plot catch/biomass
std <- read.delim(here::here(year, "base", "nr.std"), sep="", header = TRUE)
catch <- vroom::vroom(here::here(year, "data", "output", "fsh_catch.csv"))

std %>%
  filter(name=="tot_biom") %>%
  bind_cols(filter(catch, year < max(year))) %>%
  filter(year >= 1991) %>%
  dplyr::select(year, catch, value, std.dev) %>%
  bind_rows(filter(catch, year == max(year)) %>%
              left_join(vroom::vroom(here::here(year, folder, "proj", "author_f", "bigsum.csv")) %>%
                          filter(Year == year, Alt == 2) %>%
                          mutate(value = Total_Biom * 1000) %>%
                          dplyr::select(year=Year, value))) %>%
  mutate(std.dev = ifelse(is.na(std.dev), std.dev[year==max(year)-1], std.dev)) %>%
  mutate(lci = value - std.dev * 1.96,
         uci = value + std.dev * 1.96) %>%
  mutate(perc = catch / value,
         lci = catch / lci,
         uci = catch / uci,
         mean = mean(perc)) %>%
  dplyr::select(year, value, mean, perc, lci, uci) -> df

png(filename=here::here(year, "figs", "catch_bio.png"), width = 6.5, height = 6.5,
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

# plot survey results ----
vast <- vroom::vroom(here::here(year, "data", "user_input", "2023_vast_default.csv"))
logvast <- vroom::vroom(here::here(year, "data", "user_input", "2023_vast_lognormal.csv"))
sb <- vroom::vroom(here::here(year, "data", "output",  "goa_total_bts_biomass.csv"))

png(filename=here::here(year, "figs", "bts_biomass.png"), width = 6.5, height = 6.5,
    units = "in", type ="cairo", res = 200)

vast %>%
  mutate(Model = "VAST") %>% 
  bind_rows(logvast %>% 
              mutate(Model = "VAST-lognormal")) %>% 
  dplyr::mutate(t = biomass/1000,
                lci = t - sd/1000 * 1.96,
                uci = t + sd/1000 * 1.96)  %>%
  dplyr::select(-biomass, -sd) %>%
  dplyr::bind_rows(sb %>%
                     dplyr::rename(t = biomass) %>%
                     dplyr::mutate(Model = "Design-based") %>%
                     dplyr::select(-se)) %>%
  dplyr::group_by(Model) %>%
  dplyr::mutate(mean = mean(t)) %>%
  dplyr::ungroup() %>%
  ggplot2::ggplot(ggplot2::aes(year, t, fill = Model, color = Model)) +
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lci, ymax = uci), alpha = 0.2, color = NA) +
  afscassess::scale_x_tickr(data=vast, var=year) +
  ggplot2::scale_y_continuous(labels = scales::comma) +
  ggplot2::ylab("Survey biomass (t)\n") +
  ggplot2::xlab("\nYear") +
  ggplot2::expand_limits(y = 0) +
  scico::scale_fill_scico_d(palette = "roma", begin = 0.25) +
  scico::scale_color_scico_d(palette = "roma", begin = 0.25) +
  ggplot2::geom_line(ggplot2::aes(y = mean), lty = 3) +
  afscassess::theme_report() +
  ggplot2::theme(legend.position = c(0.2, 0.8)) 

dev.off()

# setup for next year
setup_folders(year+1)
accepted_model(2022, "m22.3a", year+1)