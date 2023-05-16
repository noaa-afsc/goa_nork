library(tidyverse)
library(afscassess)
library(vroom)
library(here)
theme_set(funcr::theme_present())

# globals
year = 2022
model = "m22.3a"


# catch ----
vroom(here(year, "data", "output", "fsh_catch.csv")) %>% 
  ggplot(aes(Year, Catch)) + 
  geom_col()

# comp ----
vroom(here(year, "data", "output", "fsh_length_comp.csv")) %>% 
  mutate(id = "Fishery") %>% 
  bind_rows(vroom(here(year, "data", "output", "goa_ts_length_comp.csv")) %>% 
  mutate(id = "Survey")) %>% 
  select(-(2:4)) %>% 
  pivot_longer(-c(year, id)) %>% 
  mutate(length = as.numeric(name)) %>% 
  ggplot(aes(length, year, size = value)) + 
  geom_point(show.legend = F) + 
  scale_size_area() +
  facet_wrap(~id)

# survey ----
plot_compare_survey(year,  models = c('2022, base', '2022, m18.2b', '2022, m22', '2022, db')) +
  funcr::theme_present()
 

# inputs
models = c('2022, base', '2022, m18.2b', '2022, m22', '2022, db')
  

  dat = data.frame()
  m = list(rep(NA, length(models)))
  for(i in 1:length(models)) {
    m[[i]] = scan(text = models[i], sep = ",", what = "")
    m[[i]][2] = gsub(" ", "", m[[i]][2])
    # m[[i]][3] = gsub(" ", "", m[[i]][3])
    
    year = m[[i]][1]
    # folder = m[[i]][2]
    model = m[[i]][2]
    
    id = model
    # if(model=="db"){
    #   id = "design-based"
    # } else if(model=="m15.5a"){
    #   id = "A"
    # } else if(model=="pois_gamma_750"){
    #   id = "B"
    # } else if(model=="log_1000"){
    #   id = "C"
    # } else if(model=="log_750"){
    #   id = "D"
    # } else if(model=="pois_log_500"){
    #   id = "E"
    # } else {
    #   id = "F"
    # }
    
    dat %>%
      dplyr::bind_rows(
        read.csv(here::here(year,  model, "processed", "survey.csv")) %>%
          dplyr::rename_all(tolower) %>%
          dplyr::select(year = starts_with("y"),
                        Observed = starts_with("bio"),
                        # Predicted = pred,
                        se, lci, uci) %>%
          tidyr::pivot_longer(-c(year, se, uci, lci)) %>%
          dplyr::mutate(value = value / 1000,
                        uci = uci / 1000,
                        lci = lci / 1000,
                        model = id)) -> dat
  }
  
  dat %>%
    ggplot2::ggplot(ggplot2::aes(year, value, color = model)) +
    ggplot2::geom_point(data = dplyr::filter(dat, name == "Observed"), position=ggplot2::position_dodge(width=0.7)) +
    ggplot2::geom_errorbar(data = dplyr::filter(dat, name == "Observed"),
                           ggplot2::aes(ymin = lci, ymax = uci, color = model), width = 0.4, position=ggplot2::position_dodge(width=0.7)) +
    # ggplot2::geom_line(data = dplyr::filter(dat, name == "Predicted"),
                       # ggplot2::aes(color = model)) +
    scico::scale_color_scico_d(palette = 'batlow', begin = 0.2, end = 0.8) +
    ggplot2::scale_x_continuous(breaks = funcr::tickr(dat, year)$breaks,
                                labels = funcr::tickr(dat, year)$labels) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::xlab("Year") +
    ggplot2::ylab("Survey biomass (kt)") +
    ggplot2::expand_limits(y = 0) +
    funcr::theme_present() +
    ggplot2::theme(legend.justification=c(1,0),
                   legend.position=c(0.2,0.70))
  
# outputs
  
  models = c('2022, base', '2022, m18.2b', '2022, m22', '2022, db', '2022, m22.1', '2022, m22.1a', '2022, m22.1b' )
  
  dat = data.frame()
  m = list(rep(NA, length(models)))
  for(i in 1:length(models)) {
    m[[i]] = scan(text = models[i], sep = ",", what = "")
    m[[i]][2] = gsub(" ", "", m[[i]][2])
    # m[[i]][3] = gsub(" ", "", m[[i]][3])
    
    year = m[[i]][1]
    # folder = m[[i]][2]
    model = m[[i]][2]
    
    id = model
    # if(model=="db"){
    #   id = "design-based"
    # } else if(model=="m15.5a"){
    #   id = "A"
    # } else if(model=="pois_gamma_750"){
    #   id = "B"
    # } else if(model=="log_1000"){
    #   id = "C"
    # } else if(model=="log_750"){
    #   id = "D"
    # } else if(model=="pois_log_500"){
    #   id = "E"
    # } else {
    #   id = "F"
    # }
    
    dat %>%
      dplyr::bind_rows(
        read.csv(here::here(year,  model, "processed", "survey.csv")) %>%
          dplyr::rename_all(tolower) %>%
          dplyr::select(year = starts_with("y"),
                        Observed = starts_with("bio"),
                        Predicted = pred,
                        se, lci, uci) %>%
          tidyr::pivot_longer(-c(year, se, uci, lci)) %>%
          dplyr::mutate(value = value / 1000,
                        uci = uci / 1000,
                        lci = lci / 1000,
                        model = id)) -> dat
  }
  
  dat %>%
    ggplot2::ggplot(ggplot2::aes(year, value, color = model)) +
    # ggplot2::geom_point(data = dplyr::filter(dat, name == "Observed"), position=ggplot2::position_dodge(width=0.7)) +
    # ggplot2::geom_errorbar(data = dplyr::filter(dat, name == "Observed"),
                           # ggplot2::aes(ymin = lci, ymax = uci, color = model), width = 0.4, position=ggplot2::position_dodge(width=0.7)) +
    ggplot2::geom_line(data = dplyr::filter(dat, name == "Predicted"),
                       ggplot2::aes(color = model)) +
    scico::scale_color_scico_d(palette = 'batlow', begin = 0.2, end = 0.8) +
    ggplot2::scale_x_continuous(breaks = funcr::tickr(dat, year)$breaks,
                                labels = funcr::tickr(dat, year)$labels) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::xlab("Year") +
    ggplot2::ylab("Survey biomass (kt)") +
    ggplot2::expand_limits(y = 0) +
    funcr::theme_present() +
    ggplot2::theme(legend.justification=c(1,0),
                   legend.position=c(0.2,0.20))
  

  # biomass 
  
  models = c('2022, base', '2022, m18.2b', '2022, m22', '2022, db')
  
  models = c('2022, m22', '2022, m22.1', '2022, m22.1a', '2022, m22.1b' )
  
  
  dat = data.frame()
  m = list(rep(NA, length(models)))
  for(i in 1:length(models)) {
    m[[i]] = scan(text = models[i], sep = ",", what = "")
    m[[i]][2] = gsub(" ", "", m[[i]][2])
    # m[[i]][3] = gsub(" ", "", m[[i]][3])
    
    year = m[[i]][1]
    # folder = m[[i]][2]
    model = m[[i]][2]
    id = model
    # if(model=="db"){
    #   id = "design-based"
    # } else if(model=="m15.5a"){
    #   id = "A"
    # } else if(model=="pois_gamma_750"){
    #   id = "B"
    # } else if(model=="log_1000"){
    #   id = "C"
    # } else if(model=="log_750"){
    #   id = "D"
    # } else if(model=="pois_log_500"){
    #   id = "E"
    # } else {
    #   id = "F"
    # }
    
    yrs = read.csv(here::here(year, model, "processed", "ages_yrs.csv"))$yrs
    bio = read.csv(here::here(year, model, "processed", "bio_rec_f.csv"))
    
    dat %>%
      dplyr::bind_rows(
        read.csv(here::here(year, model, "processed", "mceval.csv"))  %>%
          dplyr::select(paste0("tot_biom_", yrs)) %>%
          dplyr::mutate(group = 1:dplyr::n()) %>%
          tidyr::pivot_longer(-group) %>%
          dplyr::mutate(year = as.numeric(gsub("tot_biom_", "", name)),
                        name = "Total biomass") %>%
          dplyr::bind_rows( read.csv(here::here(year,model, "processed", "mceval.csv")) %>%
                              dplyr::select(paste0("spawn_biom_", yrs)) %>%
                              dplyr::mutate(group = 1) %>%
                              tidyr::pivot_longer(-group) %>%
                              dplyr::mutate(year = as.numeric(gsub("spawn_biom_", "", name)),
                                            name = "Spawning biomass")) %>%
          dplyr::mutate(name = factor(name, levels = c("Total biomass", "Spawning biomass"))) %>%
          dplyr::group_by(year, name) %>%
          dplyr::summarise(median = median(value) / 1000,
                           lci = quantile(value, 0.025) / 1000,
                           uci = quantile(value, 0.975) / 1000) %>%
          dplyr::ungroup() %>%
          dplyr::left_join(data.frame(year = yrs,
                                      tot = bio$tot_biom / 1000,
                                      bio = bio$sp_biom / 1000)) %>%
          dplyr::mutate(biomass = ifelse(name == "Total biomass", tot, bio),
                        model = id) %>%
          dplyr::select(-tot, -bio)) -> dat
  }
  
  dummy = data.frame(year = rep(unique(dat$year),4),
                     name = rep(c("Total biomass", "Spawning biomass"), each = 2 * length(unique(dat$year))),
                     biomass = c(rep(0, length(unique(dat$year))), rep(160, length(unique(dat$year))),
                                 rep(0, length(unique(dat$year))), rep(60, length(unique(dat$year)))),
                     model = NA)
  
  dat %>%
    ggplot2::ggplot(ggplot2::aes(year, biomass, color = model, fill = model)) +
    ggplot2::geom_blank(data = dummy) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lci, ymax = uci), alpha = 0.1, color = NA) +
    ggplot2::facet_wrap(~name, dir = "v", scales = "free_y") +
    ggplot2::scale_y_continuous(name = "Biomass (kt)", labels = scales::comma) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::scale_x_continuous(name = "Year",
                                breaks = funcr::tickr(dat, year, 10, start = 1960)$breaks,
                                labels = funcr::tickr(dat, year, 10, start = 1960)$labels) +
    scico::scale_color_scico_d(palette = 'batlow', begin = 0.2, end = 0.8) +
    scico::scale_fill_scico_d(palette = 'batlow', begin = 0.2, end = 0.8) +
    funcr::theme_present() +
    ggplot2::theme(legend.position = c(0.2, .8))
  