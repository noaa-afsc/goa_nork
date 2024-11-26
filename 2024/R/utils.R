library(tidyverse)
theme_set(afscassess::theme_report())


# parameter plots ----
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

get_survival <- function(post, obj, reps) {
  
  mapit = function(i, post, obj) {
    surv_i = obj$report(post[i,])$Nat[,T] * obj$report(post[i,])$Sat[,T]
    list(surv_i = surv_i)
  }
  results = purrr::map(1:reps, mapit, post=post, obj=obj)
  do.call(cbind, map(results, "surv_i"))
}

# obj$report(post[i,-ncol(post)])$q

get_bio <- function(post, obj, reps) {
  
  mapit = function(i, post, obj) {
    surv_i = obj$report(post[i])$spawn_bio
    list(surv_i = surv_i)
  }
  results = purrr::map(1:reps, mapit, post=post, obj=obj)
  do.call(cbind, map(results, "surv_i")) 
  
}
# post = post[200:1000,]
# bb = get_bio(post, obj1, nrow(post))
# bb %>% 
#   as.data.frame() %>% 
#     mutate(year = years) %>% 
#   pivot_longer(-year) %>% 
#   ggplot(aes(year, value, group = name)) + 
#   geom_line(alpha = 0.2) +
#   geom_hline(yintercept = c(report1$B40, report1$B35))


proj_bio <- function(report, obj=NULL, post=NULL, reps = 500) {
  
  # values
  F40 = report$F40
  F35 = report$F35
  B40 = report$B40
  Nat = report$Nat
  Sat = report$Sat
  slx = report$slx
  ages = report$ages
  years = report$years
  waa = report$waa
  wt_mature = report$wt_mature
  spawn_frac = report$spawn_fract
  yield_ratio = report$yield_ratio
  M = report$M
  pred_rec = report$pred_rec
  stdev_rec = report$stdev_rec 
  A = nrow(Nat) # number of ages
  T = ncol(Nat) # number of years
  
  # storage
  if(!is.null(post)) {
    surv = get_survival(post=post, obj=obj, reps = reps)
    Tproj = 15
    N = Fabc_tot = Fofl_tot = replicate(reps, matrix(0, A, Tproj), simplify = FALSE)
    F40_proj = F35_proj = spawn_bio = tot_bio = replicate(reps, rep(0, Tproj), simplify = FALSE)
    for(i in 1:reps) {
      F40_proj[[i]] = rep(F40, Tproj)
      F35_proj[[i]] = rep(F35, Tproj)
      
      # total F
      Fabc_tot[[i]][,1] = report$slx[,1] * F40
      Fofl_tot[[i]][,1] = report$slx[,1] * F35
      
      # populate abundance
      N[[i]][1,] = exp(log(pred_rec) - (stdev_rec)^2 / 2 + stdev_rec + rnorm(15))
      
      for(a in 1:(A-1)) {
        N[[i]][a+1,1] = surv[,i][a]
      }
      N[[i]][A,1] = surv[,i][A-1] + surv[,i][A]
      spawn_bio[i][1] = sum(N[[i]][,1] * exp(-yield_ratio * unlist(Fabc_tot[[i]][,1]) - M)^spawn_frac * wt_mature)
      tot_bio[[i]][1] = sum(N[[i]][,1] * waa)   
      
      for(t in 1:Tproj) {
        # tier check
        if((spawn_bio[[i]][t] / B40) <= 1) {
          F40_proj[[i]][t] = F40_proj[[i]][t] * (spawn_bio[[i]][t] / B40 - 0.05) / 0.95
          F35_proj[[i]][t] = F35_proj[[i]][t] * (spawn_bio[[i]][t] / B40 - 0.05) / 0.95
        }
        # update 
        Fabc_tot[[i]][,t] = report$slx[,1] * F40_proj[[i]][t]
        Fofl_tot[[i]][,t] = report$slx[,1] * F35_proj[[i]][t]
        Z_proj = unlist(Fabc_tot) + M
        Zofl_proj = unlist(Fofl_tot) + M
        S_proj = exp(-Z_proj)
        
        # catch
        # Cat_proj[,t] = yield_ratio * N_proj[,t] * Fabc_tot_proj / Z_proj * (1 - S_proj)
        # Cat_ofl_proj[,t] = yield_ratio * N_proj[,t] * Fofl_tot_proj / Zofl_proj * (1 - exp(-Zofl_proj))
        
        if(t<Tproj) {
          for(a in 1:(A-1)){
            N[[i]][a+1,t+1] = N[[i]][a,t] * exp(-yield_ratio * unlist(Fabc_tot[i])[a] - M)
          }
          N[[i]][A,t+1] = N[[i]][A-1,t] * exp(-yield_ratio * unlist(Fabc_tot[i])[A-1] - M) +
            N[[i]][A,t] * exp(-yield_ratio * unlist(Fabc_tot[i])[A] - M)
          
          tot_bio[[i]][t+1] = sum(N[[i]][,t+1] * waa)  
          spawn_bio[[i]][t+1] = sum(N[[i]][,t+1] * exp(-yield_ratio * unlist(Fabc_tot[[i]][,t]) - M)^spawn_frac * wt_mature)
        }
      }
    }
    
    do.call(rbind, spawn_bio) %>% 
      as.data.frame() %>% 
      mutate(sim = 1:n()) %>% 
      tidyr::pivot_longer(-sim) %>% 
      mutate(year = as.numeric(gsub('V','', name)) + max(years),
             id = 'spawn_bio') %>% 
      select(-name) %>% 
      dplyr::bind_rows(
        do.call(rbind, tot_bio) %>% 
          as.data.frame() %>% 
          mutate(sim = 1:n()) %>% 
          tidyr::pivot_longer(-sim) %>% 
          mutate(year = as.numeric(gsub('V','', name)) + max(years),
                 id = 'tot_bio') %>% 
          select(-name)) 
    
  } else {
    Tproj = 2
    N = Cat = Cat_ofl= Zabc = Zofl = S = matrix(0, A, Tproj)
    tot_bio = spawn_bio = F40_proj = F35_proj= rep(0, Tproj)
    # setup
    F40_proj[1] = F40
    F35_proj[1] = F35
    
    # total F
    Fabc_tot = slx[,1] * F40_proj[1]
    Fofl_tot = slx[,1] * F35_proj[1]
    
    # first projection year
    N[1,] = pred_rec
    for(a in 1:(A-1)) {
      N[a+1,1] = Nat[a,T] * Sat[a,T]
    }
    N[A,1] = Nat[A-1,T] * Sat[A-1,T] + Nat[A,T] * Sat[A,T]
    spawn_bio[1] = sum(N[,1] * exp(-yield_ratio * Fabc_tot - M)^spawn_frac * wt_mature)
    
    for(t in 1:Tproj) {
      # tier check
      if((spawn_bio[t] / B40) > 1) {
        F40_proj[t] = F40
        F35_proj[t] = F35
      } else {
        F40_proj[t] = F40_proj[t] * (spawn_bio[t] / B40 - 0.05) / 0.95
        F35_proj[t] = F35_proj[t] * (spawn_bio[t] / B40 - 0.05) / 0.95
      }
      # update
      Fabc_tot = report$slx[,1] * F40_proj[t]
      Fofl_tot = report$slx[,1] * F35_proj[t]
      Z = Fabc_tot + M
      Zofl = Fofl_tot + M
      S = exp(-Z)
      
      # catch
      Cat[,t] = yield_ratio * N[,t] * Fabc_tot / Z * (1 - S)
      Cat_ofl[,t] = yield_ratio * N[,t] * Fofl_tot / Zofl * (1 - exp(-Zofl))
      
      if(t<Tproj) {
        for(a in 1:(A-1)){
          N[a+1,t+1] = N[a,t] * exp(-yield_ratio * Fabc_tot[a] - M)
        }
        N[A,t+1] = N[A-1,t] * exp(-yield_ratio * Fabc_tot[A-1] - M) +
          N[A,t] * exp(-yield_ratio * Fabc_tot[A] - M)
        
        tot_bio[t+1] = sum(N[,t+1] * waa)
        spawn_bio[t+1] = sum(N[,t+1] * exp(-yield_ratio * Fabc_tot - M)^spawn_frac * wt_mature)
      }
    }
    catch = colSums(Cat * waa / yield_ratio)
    catch_ofl = colSums(Cat_ofl * waa / yield_ratio)
    tot_bio = colSums(N * waa)
    
    data.frame(year = max(years)+1:Tproj,
               spawn_bio = spawn_bio,
               tot_bio = tot_bio,
               catch_abc = catch,
               catch_ofl = catch_ofl,
               F40 = F40_proj,
               F35 = F35_proj)
  }
}

single_like <- function(par_name, par_values, obj, fit) {
  n = length(par_values)
  ll <- numeric(n)
  
  for (i in 1:n) {
    params <- fit$par
    params[par_name] <- par_values[i]
    ll[i] <- obj$fn(params)
  }
  
  return(ll)
}

# Function to profile the likelihood for multiple parameters
multi_like <- function(par_names, par_ranges, obj, fit) {
  grid = expand.grid(par_ranges)
  ll = numeric(nrow(grid))
  
  for (i in 1:nrow(grid)) {
    params = fit$par
    for (j in 1:length(par_names)) {
      params[par_names[j]] = grid[i, j]
    }
    ll[i] = obj$fn(params)
  }
  
  return(data.frame(grid, log_like=ll))
}

get_spawn_bio <- function(obj, post, iters=100) {
  if(iters>nrow(post[[1]])) {iters = ncol(post[[1]])}
  
  sb = matrix(nrow=length(post)*iters, ncol=length(obj$report(post[[1]][1,])$spawn_bio))
  for (i in 1:iters) {
    for(j in 1:length(post)){
      sb[i,] = obj$report(post[[j]][i,])$spawn_bio
    }
  }
  
  sb = data.frame(sb) 
  names(sb) = obj$report(post[[1]][1,])$years
  sb %>% 
    dplyr::mutate(sim = 1:dplyr::n())  %>% 
    tidyr::pivot_longer(-sim, values_to = 'spawn_bio', names_to = 'year') %>% 
    dplyr::mutate(year = as.numeric(year))
}

get_tot_bio <- function(obj, post, iters=100) {
  if(iters>nrow(post[[1]])) {iters = ncol(post[[1]])}
  
  sb = matrix(nrow=length(post)*iters, ncol=length(obj$report(post[[1]][1,])$tot_bio))
  for (i in 1:iters) {
    for(j in 1:length(post)){
      sb[i,] = obj$report(post[[j]][i,])$tot_bio
    }
  }
  
  sb = data.frame(sb) 
  names(sb) = obj$report(post[[1]][1,])$years
  sb %>% 
    dplyr::mutate(sim = 1:dplyr::n())  %>% 
    tidyr::pivot_longer(-sim, values_to = 'tot_bio', names_to = 'year') %>% 
    dplyr::mutate(year = as.numeric(year))
}

get_rec <- function(obj, post, iters=100) {
  if(iters>nrow(post[[1]])) {iters = ncol(post[[1]])}
  
  sb = matrix(nrow=length(post)*iters, ncol=length(obj$report(post[[1]][1,])$recruits))
  for (i in 1:iters) {
    for(j in 1:length(post)){
      sb[i,] = obj$report(post[[j]][i,])$recruits
    }
  }
  
  sb = data.frame(sb) 
  names(sb) = obj$report(post[[1]][1,])$years
  sb %>% 
    dplyr::mutate(sim = 1:dplyr::n())  %>% 
    tidyr::pivot_longer(-sim, values_to = 'recruits', names_to = 'year') %>% 
    dplyr::mutate(year = as.numeric(year))
}

effn <- function(obs, pred) {
  
  colSums((1 - pred) * pred) / colSums((obs - pred)^2) 
}

sdnr <- function(obs, pred, iss) {
  n = ncol(obs) 
  sdnr = vector(length = n)
  for(i in 1:n) {
    
    sdnr[i] = sd((obs[,i] - pred[,i]) / sqrt(pred[,i] * (1 - pred[,i]) / iss[i]))
  }
  sdnr
}


osa <- function(obs, pred, iss, yrs, ind, label = 'Age', outlier=3) {
  
  obs = round(iss * obs / colSums(obs))
  pred = pred / colSums(pred)
  res = compResidual::resMulti(obs, pred)
  mat = matrix(res, nrow=nrow(res), ncol=ncol(res))
  df = as.data.frame(mat)
  names(df) <- yrs
  
  df %>%
    mutate(ind = head(ind, -1)) %>%
    tidyr::pivot_longer( -ind, names_to = 'year') %>%
    mutate(year = as.numeric(year),
           Year = factor(year),
           id = ifelse(value<0 , 'a', 'b'),
           Outlier = ifelse(abs(value) >= outlier, "Yes", "No"),
           Outlier = factor(Outlier, levels = c('No', 'Yes'))) -> df
  
  df %>%
    ggplot(aes(year, ind, color = value, size = value, shape = Outlier) ) +
    geom_point(show.legend=TRUE) +
    scale_size_area(guide=F, max_size = 3) +
    scico::scale_color_scico(limits = c(-4, 4), palette = 'vik') +
    afscassess::scale_y_tickr(data = df, var = ind, start=0) +
    afscassess::scale_x_tickr(data = df, var = year) +
    scale_shape_manual(values = c(19,8), drop = FALSE) +
    ylab(label) +
    xlab('Year') +
    ggtitle('OSA')
  
}
pearson <- function(obs, pred, iss, yrs, ind, label = 'Age', outlier=3) {
  
  as.data.frame(iss * (obs -pred) / sqrt(iss * pred)) %>%
    mutate(ind = ind) -> df
  names(df) <- c(yrs, 'ind')
  
  df %>%
    tidyr::pivot_longer( -ind, names_to = 'year') %>%
    mutate(year = as.numeric(year),
           Year = factor(year),
           id = ifelse(value<0 , 'a', 'b'),
           Outlier = ifelse(abs(value) >= outlier, "Yes", "No"),
           Outlier = factor(Outlier, levels = c("No", "Yes"))) -> df
  
  
  df %>%
    ggplot(aes(year, ind, color = value, size = value, shape = Outlier) ) +
    geom_point(show.legend=TRUE) +
    scale_size_area(guide=F, max_size = 3) +
    scico::scale_color_scico(limits = c(-4, 4), palette = 'vik') +
    afscassess::scale_y_tickr(data = df, var = ind, start=0) +
    afscassess::scale_x_tickr(data = df, var = year) +
    scale_shape_manual(values = c(19,8), drop = FALSE) +
    ylab(label) +
    xlab('Year') +
    ggtitle('Pearson')
  
}

#' aggregate residual plot for RTMB model
#'
#' @param obs observed comp
#' @param pred predictedcomp
#' @param ind vector or ages of lengths
#' @param label axis label
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obs = fish_age_obs
#' pred = report1$fish_age_pred
#' ind = ages
#' agg <- function(obs, pred, ind, label = 'Age')
#' }
agg <- function(obs, pred, ind, label = 'Age') {
  
  df = data.frame(obs = rowSums(obs)/sum(obs),
                  pred = rowSums(pred)/sum(pred),
                  ind = ind)
  
  df %>%
    ggplot(aes(ind, pred)) +
    geom_bar(aes(y=obs), stat = 'identity', alpha=0.4) +
    geom_point() +
    geom_line() +
    afscassess::scale_x_tickr(data=df, var=ind, start=0) +
    xlab(label) +
    ylab('') +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
}

annual <- function(obs, pred, ind, yrs, label = 'Age') {
  obs = as.data.frame(obs)
  pred = as.data.frame(pred)
  names(obs) <- names(pred) <- yrs
  obs %>% 
    dplyr::mutate(type = 'obs',
                  ind = ind) %>% 
    bind_rows(
      pred %>% 
        dplyr::mutate(type = 'pred',
                      ind = ind)
    ) %>% 
    tidyr::pivot_longer(-c(ind, type)) %>% 
    mutate(year = as.numeric(name)) -> df
  
  df %>% 
    filter(type=='obs') %>% 
    ggplot(aes(ind, value)) + 
    geom_col(alpha = 0.5, aes(fill = factor(ind)), show.legend = FALSE) + 
    geom_line(data = filter(df, type=='pred')) +
    facet_wrap(~year, ncol = 2, dir='v') + 
    afscassess::scale_x_tickr(data=df, var=ind, start=0) +
    scico::scale_fill_scico_d(palette = 'managua') +
    xlab(label) +
    ylab("") +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
}

sample_size <- function(obs, pred, iss, yrs){
  data.frame(year = yrs,
             ISS = iss,
             effN = effn(obs, pred),
             sdnr = sdnr(obs, pred, iss)) %>% 
    tidyr::pivot_longer(-year) %>% 
    mutate(grp = ifelse(name=='sdnr', 'SDNR', 'Sample size')) -> df
    
  df %>% 
    ggplot(aes(year, value, color = name)) +
    geom_point() +
    facet_wrap(~grp, scales = 'free_y', dir = 'h') +
    scale_color_manual("", breaks = c('effN', 'ISS'), values = c("#7E1700","#5DC0D2",1)) +
    expand_limits(y = 0) +
    theme(legend.position=c(0.1,0.8)) +
    afscassess::scale_x_tickr(data=df, var=year, to=10, start = 1960) +
    xlab('Year') +
    ylab('Value')
}
#' get all residual plots for RTMB model
#'
#' @param obs observed comp
#' @param pred predictedcomp
#' @param iss input sample size
#' @param yrs comp years
#' @param ind vector or ages of lengths
#' @param label axis label
#'
#' @export
#'
#' @examples
#' \dontrun{
#' iss = data$fish_age_iss
#' obs = fish_age_obs
#' pred = report1$fish_age_pred
#' yrs = fish_age_yrs
#' label = "Age"
#' ind = ages
#' fish_age_resids <- resids(obs = fish_age_obs, pred = report1$fish_age_pred, iss = data$fish_age_iss,
#'                          yrs = fish_age_yrs, ind = ages)
#'  fish_age_resids$osa + ggtitle('osa fishery age comp residuals')
#' }
resids <- function(obs, pred, iss, yrs, ind, label = 'Age', outlier=3) {
  list(osa = osa(obs, pred, iss, yrs, ind, label, outlier),
       pearson = pearson(obs, pred, iss, yrs, ind, label, outlier),
       agg = agg(obs, pred, ind, label),
       annual = annual(obs, pred, ind, yrs, label),
       ss = sample_size(obs, pred, iss, yrs) )
}

# out = resids(obs = data$fish_size_obs, pred = m24$fish_size_pred, yrs = data$fish_size_yrs, iss=data$fish_size_iss, ind = length_bins, label = 'Length (cm)')
# library(patchwork)
# 
# (out$pearson +
#     out$osa) /
#   (out$agg + out$ss) + plot_layout(guides = "collect") 

