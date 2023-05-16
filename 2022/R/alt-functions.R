# update selectivity
model = "m22.1"
params_table <- function(year, model, model_name){
  
  STD = read.delim(here::here(year, model, paste0(model_name, ".std")), sep="", header = TRUE)
  mceval = read.csv(here::here(year, model, "processed", "mceval.csv"))
  
  params <- function(year, STD, mceval, param){
    
    # doh - use the same names!!!
    param2 = ifelse(param == "nattymort", "natmort", param)
    
    std = STD %>%
      dplyr::filter(name == !!param)
    
    if(param == "spawn_biom_proj") {
      param2 = paste0("spawn_biom_proj_", year + 1)
      std = STD %>%
        dplyr::filter(name == !! param) %>%
        dplyr::slice_head(n=1)
    }
    
    mceval %>%
      dplyr::select(param2) %>%
      dplyr::summarise(mean = mean(!!dplyr::sym(param2)),
                       median = median(!!dplyr::sym(param2)),
                       sd = sd(!!dplyr::sym(param2)),
                       lci = quantile(!!dplyr::sym(param2), 0.025),
                       uci = quantile(!!dplyr::sym(param2), 0.975)) %>%
      dplyr::bind_cols(std) %>%
      dplyr::mutate(name := !! param) %>%
      dplyr::select(name, value, mean, median, std.dev, sd, lci, uci)
    
    
  }
  
  dplyr::bind_rows(params(year, STD, mceval, param = "q_srv1"),
                   params(year, STD, mceval, param = "nattymort"),
                   params(year, STD, mceval, param = "F40") ,
                   params(year, STD, mceval, param = "spawn_biom_proj"),
                   params(year, STD, mceval, param = "ABC")) %>%
    write.csv(here::here(year, model, "tables", "mcmc_pars.csv"), row.names = FALSE)
}

