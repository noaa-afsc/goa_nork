process_results <- function(year, model, model_name = alt, dat_name,
                            rec_age, plus_age, mcmc, mcsave, len_bins, n_retro=10){
  
  # setup
  if (!dir.exists(here::here(year, model, "processed"))){
    dir.create(here::here(year, model, "processed"), recursive=TRUE)
  }
  
  if (!dir.exists(here::here(year, model, "figs"))){
    dir.create(here::here(year, model, "figs"), recursive=TRUE)
  }
  
  if (!dir.exists(here::here(year, model, "tables"))){
    dir.create(here::here(year, model, "tables"), recursive=TRUE)
  }
  
  # helper functions
  rep_item <- function(name){
    t <- strsplit(REP[grep(name, REP)]," ")
    t <- subset(t[[1]], t[[1]]!="")
    if(t[[1]][1] == "TWL"){
      as.numeric(t[3:length(t)])
    } else {
      as.numeric(t[2:length(t)])
    }
  }
  
  
  # read in rep and ctl files
  REP <- readLines(here::here(year, model, paste0(model_name, ".rep")))
  CTL <- readLines(here::here(year, model, paste0(dat_name, ".ctl")))
  PSV <- file(here::here(year, model, paste0(model_name, ".psv")), "rb")
  STD <- read.delim(here::here(year, model, paste0(model_name, ".std")), sep="", header = TRUE)
  mceval <- read.delim(here::here(year, model, "evalout.prj"), sep="", header=FALSE)
  
  # clean rep file
  suppressWarnings(data.frame(year = unlist(strsplit(REP[grep("Year", REP)[1]]," "))) %>%
                     dplyr::mutate(year = as.numeric(year)) %>%
                     tidyr::drop_na() %>%
                     dplyr::pull(year)) -> yrs
  
  suppressWarnings(data.frame(age = unlist(strsplit(REP[grep("Age", REP)[1]]," "))) %>%
                     dplyr::mutate(age = as.numeric(age)) %>%
                     tidyr::drop_na() %>%
                     dplyr::pull(age)) -> ages
  
  styr_rec <- yrs[1] - length(ages) + 2
  
  suppressWarnings(as.data.frame(cbind(yrs = yrs, ages = ages, styr_rec = styr_rec)) %>%
                     dplyr::mutate(ages = replace(ages, duplicated(ages), NA),
                                   styr_rec = replace(styr_rec, duplicated(styr_rec), NA))) %>%
    write.csv(here::here(year, model, "processed", "ages_yrs.csv"), row.names = FALSE)
  
  # MCMC parameters ----
  
  npar = readBin(PSV, what = integer(), n=1)
  mcmcs = readBin(PSV, what = numeric(), n = (npar * mcmc / mcsave))
  close(PSV)
  mcmc_params = matrix(mcmcs, byrow=TRUE, ncol=npar)
  mcmc_params = mcmc_params[501:nrow(mcmc_params),]
  colnames(mcmc_params) = STD$name[1:ncol(mcmc_params)]
  write.csv(mcmc_params, here::here(year, model, "processed", "mcmc.csv"), row.names = FALSE)
  
  # mceval phase output ----
  
  #Curry's Change
  mceval = mceval[501:nrow(mceval),]
  
  #Length colnames = 286
  # columns mcmc_other = 271
  
  #1-8: Through objective function value
  
  colnames(mceval) = c("sigr", "q_srv1", "q_srv2", "F40", "natmort", "spawn_biom_proj",
                       "ABC", "obj_fun",
                       paste0("tot_biom_", yrs),
                       paste0("log_rec_dev_", seq(styr_rec, yrs[length(yrs)])),
                       paste0("spawn_biom_", yrs),
                       "log_mean_rec",
                       paste0("spawn_biom_proj_", max(yrs) + 1:15),
                       paste0("pred_catch_proj_", max(yrs) + 1:15),
                       paste0("rec_proj_", max(yrs) + 1:10),
                       paste0("tot_biom_proj_", max(yrs)))
  write.csv(mceval, here::here(year, model, "processed", "mceval.csv"), row.names = FALSE)
  
  # catch data ----
  
  pred = strsplit(REP[grep("Pred_Catch", REP)], " ")
  r1 = which(pred[[1]] == "Pred_Catch")
  r2 = which(pred[[1]] == "Pred_catch_later")
  r3 = which(pred[[1]] == "")
  pred = as.numeric(pred[[1]][-c(r1, r2, r3)])
  
  obs = strsplit(REP[grep("Obs_Catch", REP)], " ")
  r1 = which(obs[[1]] == "Obs_Catch")
  r2 = which(obs[[1]] == "Obs_Catch_Later")
  r3 = which(obs[[1]] == "")
  obs = as.numeric(obs[[1]][-c(r1, r2, r3)])
  
  data.frame(obs = obs, pred = pred) %>%
    write.csv(here::here(year, model, "processed", "catch.csv"))
  
  
  
  # survey data ----
  syr = REP[grep("Survey Biomass",REP)[1]:(grep("Survey Biomass",REP)[2]-2)][2]
  syr = strsplit(syr," ")
  syr = subset(syr[[1]], syr[[1]]!="")
  syr = as.numeric(syr[2:length(syr)])
  
  obs = REP[grep("Survey Biomass",REP)[1]:(grep("Survey Biomass",REP)[2]-2)][4]
  obs = strsplit(obs," ")
  obs = subset(obs[[1]], obs[[1]]!="")
  obs = as.numeric(obs[2:length(obs)])
  
  se = REP[grep("Survey Biomass",REP)[1]:(grep("Survey Biomass",REP)[2]-2)][6]
  se = strsplit(se," ")
  se = subset(se[[1]], se[[1]]!="")
  se = as.numeric(se[2:length(se)])
  
  pred = REP[grep("Survey Biomass",REP)[1]:(grep("Survey Biomass",REP)[2]-2)][3]
  pred = strsplit(pred," ")
  pred = subset(pred[[1]], pred[[1]]!="")
  pred = as.numeric(pred[2:length(pred)])
  
  
  data.frame(year = syr, biomass = obs, se = se) %>%
    dplyr::mutate(lci = biomass - 1.96 * se,
                  uci = biomass + 1.96 * se ) %>%
    dplyr::bind_cols(pred = pred) %>%
    write.csv(here::here(year, model, "processed", "survey.csv"), row.names = FALSE)
  
  
  # recruitment ----
  
  N = REP[grep("Numbers", REP):(grep("Obs_P_fish_age", REP)-2)]
  t = NA
  for(i in 1:length(yrs)){
    ts = as.numeric(strsplit(N[i+1]," ")[[1]][3])
    t = c(t, ts)}
  pred_rec = t[!is.na(t)]
  
  # biomass & F & recruits ----
  data.frame(year = yrs,
             tot_biom = rep_item("Tot_biom"),
             sp_biom = rep_item("SpBiom"),
             F = rep_item("Fully_selected_F"),
             recruits = pred_rec) %>%
    write.csv(here::here(year, model, "processed", "bio_rec_f.csv"), row.names = FALSE)
  
  
  # selectivity ----
  data.frame(age = ages,
             fish = rep_item("Fishery_Selectivity"),
             srv1 = rep_item("TWL Survey_Selectivity"),
             maturity = rep_item("Maturity")) %>%
    write.csv(here::here(year, model, "processed", "selex.csv"), row.names = FALSE)
  
  # yield ratio B40 & B35----
  
  data.frame(B40 = STD$value[which(STD$name=="B40")],
             B35 = as.numeric(REP[(grep("B_35",REP)+1):(grep("F_40",REP)[1]-1)]),
             yld_rat = as.numeric(unlist(strsplit(CTL[grep("yieldratio", CTL)], "\t"))[1])) %>%
    write.csv(here::here(year, model, "processed", "b35_b40_yld.csv"), row.names = FALSE)
  
  # size comps ----
  
  #! this will need a switch for multiple surveys
  
  obs = REP[grep("Obs_P_fish_age",REP):(grep("Pred_P_fish_age",REP)-2)]
  pred = REP[grep("Pred_P_fish_age",REP):(grep("Obs_P_fish_size",REP)-2)]
  
  obs_l = REP[grep("Obs_P_fish_size",REP):(grep("Pred_P_fish_size",REP)-2)]
  pred_l = REP[grep("Pred_P_fish_size",REP):(grep("Obs_P_srv1_age",REP)-2)]
  
  s_obs = REP[grep("Obs_P_srv1_age",REP):(grep("Pred_P_srv1_age",REP)-2)]
  s_pred = REP[grep("Pred_P_srv1_age",REP):(grep("Obs_P_srv1_size",REP)-2)]
  
  s_obs_l = REP[grep("Obs_P_srv1_size",REP):(grep("Pred_P_srv1_size",REP)-2)]
  
  rockfishr::purrit(obs, pred, rec_age, plus_age, comp = "age", lenbins = len_bins) %>%
    write.csv(here::here(year, model, "processed", "fac.csv"))
  
  rockfishr::purrit(obs_l, pred_l, rec_age, plus_age, comp = "length", lenbins = len_bins) %>%
    write.csv(here::here(year, model, "processed", "fsc.csv"))
  
  rockfishr::purrit(s_obs, s_pred, rec_age, plus_age, comp = "age", lenbins = len_bins) %>%
    write.csv(here::here(year, model, "processed", "sac.csv"))
  
  rockfishr::purrit(s_obs_l, pred = NULL, rec_age, plus_age, comp = "length", lenbins = len_bins) %>%
    write.csv(here::here(year, model, "processed", "ssc.csv"))
  
}

run_retro <- function(year, model, model_name, dat_name, mcmc, n_retro=10) {
  
  ctl <- read.delim(here::here(year, model, paste0(dat_name, ".ctl")), header=F)
  dat <- readLines(here::here(year, model, paste0(dat_name, ".dat")))
  
  
  Sec_st<-grep("#-",dat)
  Sec_end<-grep("#!",dat)
  
  st_end<-matrix(NA,nrow=length(Sec_st),ncol=2)
  st_end[,1]<-Sec_st
  st_end[,2]<-Sec_end
  
  mcmcon<-"YES"
  mcmcruns<-  mcmc# Could change these, but I like 5000 as a manageable number to deal with
  mcmcsave<-mcmcruns/500
  
  styr<-as.numeric(dat[Sec_st[2]-3]) # start of model (example 1961 for POP)
  nages<-as.numeric(dat[Sec_st[2]+3]) # number of age bins
  nlens<-as.numeric(dat[Sec_st[2]+5]) # number of length bins
  numretros<-n_retro # number of retrospective years
  
  # Set up some results files
  RES_SB <-matrix(nrow=length(seq(styr, year)),ncol=numretros)
  rownames(RES_SB) <-seq(styr, year)
  colnames(RES_SB) <-seq( year-numretros+1, year)
  RES_Rec<-matrix(nrow=length(seq(styr, year)),ncol=numretros)
  rownames(RES_Rec)<-seq(styr, year)
  colnames(RES_Rec)<-seq( year-numretros+1, year)
  
  T_start<-Sys.time() #Timer start
  
  for(y in 1:numretros){
    
    # Set endyr
    yrs_retro<-seq(year -numretros+1, year )
    endyr<-yrs_retro[y]
    nyrs<-endyr-styr+1
    DAT_retro<-c(dat[st_end[1,1]:st_end[1,2]],as.character(endyr), dat[st_end[2,1]:st_end[2,2]])
    
    # Fishery catch
    DAT_retro<-c(DAT_retro,paste(scan(text=dat [Sec_st[3]-1])[1:nyrs],collapse=" "),dat [st_end[3,1]:st_end[3,2]])
    
    
    # Trawl survey biomass
    BTSb_yrs<-length(which(scan(text=dat [Sec_st[5]-1])<=endyr))
    DAT_retro<-c(
      DAT_retro,
      as.character(BTSb_yrs),
      dat [st_end[4,1]:st_end[4,2]],
      paste(scan(text=dat [Sec_st[5]-1])[1:BTSb_yrs],collapse=" "),
      dat [st_end[5,1]:st_end[5,2]],
      paste(scan(text=dat [Sec_st[6]-1])[1:BTSb_yrs],collapse=" "),
      dat [st_end[6,1]:st_end[6,2]],
      paste(scan(text=dat [Sec_st[7]-1])[1:BTSb_yrs],collapse=" "),
      dat [st_end[7,1]:st_end[7,2]],
      paste(scan(text=dat [Sec_st[8]-1])[1:BTSb_yrs],collapse=" "),
      dat[st_end[8,1]:st_end[8,2]],
      paste(scan(text=dat[Sec_st[9]-1])[1:BTSb_yrs],collapse=" "),
      dat[st_end[9,1]:st_end[9,2]])
    
    # Fish age comp
    FAC_yrs<-length(which(scan(text=dat [Sec_st[11]-1])<(endyr-1)))
    DAT_retro<-c(DAT_retro,
                 as.character(FAC_yrs),
                 dat [st_end[10,1]:st_end[10,2]],
                 paste(scan(text=dat [Sec_st[11]-1])[1:FAC_yrs],collapse=" "),
                 dat [st_end[11,1]:st_end[11,2]],
                 paste(scan(text= dat[Sec_st[12]-1])[1:FAC_yrs],collapse=" "),
                 dat [st_end[12,1]:st_end[12,2]],
                 paste(scan(text= dat[Sec_st[13]-1])[1:FAC_yrs],collapse=" "),
                 dat [st_end[13,1]:st_end[13,2]],
                 paste(scan(text= dat[Sec_st[14]-1])[1:FAC_yrs],collapse=" "),
                 dat [st_end[14,1]:st_end[14,2]])
    for(i in 1:FAC_yrs)   DAT_retro<-c(DAT_retro,paste(scan(text=dat[Sec_st[15]-FAC_yrs-1+i]) ,collapse = " "))
    DAT_retro<-c(DAT_retro,dat[st_end[15,1]:st_end[15,2]])
    
    # Survey age comp
    SAC_yrs<-length(which(scan(text=dat[Sec_st[17]-1])<=(endyr-1)))
    DAT_retro<-c(DAT_retro,
                 as.character(SAC_yrs),
                 dat[st_end[16,1]:st_end[16,2]],
                 paste(scan(text=dat[Sec_st[17]-1])[1:SAC_yrs],collapse=" "),
                 dat[st_end[17,1]:st_end[17,2]],
                 paste(scan(text=dat[Sec_st[18]-1])[1:SAC_yrs],collapse=" "),
                 dat[st_end[18,1]:st_end[18,2]],
                 paste(scan(text=dat[Sec_st[19]-1])[1:SAC_yrs],collapse=" "),
                 dat[st_end[19,1]:st_end[19,2]],
                 paste(scan(text=dat[Sec_st[20]-1])[1:SAC_yrs],collapse=" "),
                 dat[st_end[20,1]:st_end[20,2]])
    for(i in 1:SAC_yrs)   DAT_retro<-c(DAT_retro,paste(scan(text=dat[Sec_st[21]-SAC_yrs-1+i]) ,collapse = " "))
    DAT_retro<-c(DAT_retro,dat[st_end[21,1]:st_end[21,2]])
    
    # Fish size comp
    FSC_yrs<-length(which(scan(text=dat[Sec_st[23]-1])<=(endyr-1)))
    DAT_retro<-c(DAT_retro,
                 as.character(FSC_yrs),
                 dat[st_end[22,1]:st_end[22,2]],
                 paste(scan(text=dat[Sec_st[23]-1])[1:FSC_yrs],collapse=" "),
                 dat[st_end[23,1]:st_end[23,2]],
                 paste(scan(text=dat[Sec_st[24]-1])[1:FSC_yrs],collapse=" "),
                 dat[st_end[24,1]:st_end[24,2]],
                 paste(scan(text=dat[Sec_st[25]-1])[1:FSC_yrs],collapse=" "),
                 dat[st_end[25,1]:st_end[25,2]],
                 paste(scan(text=dat[Sec_st[26]-1])[1:FSC_yrs],collapse=" "),
                 dat[st_end[26,1]:st_end[26,2]])
    for(i in 1:FSC_yrs)   DAT_retro<-c(DAT_retro,paste(scan(text=dat[Sec_st[27]-FSC_yrs-1+i]) ,collapse = " "))
    DAT_retro<-c(DAT_retro,dat[st_end[27,1]:st_end[27,2]])
    
    # Survey size comp
    SSC_yrs<-length(which(scan(text=dat[Sec_st[29]-1])<=endyr))
    DAT_retro<-c(DAT_retro,
                 as.character(SSC_yrs),
                 dat[st_end[28,1]:st_end[28,2]],
                 paste(scan(text=dat[Sec_st[29]-1])[1:SSC_yrs],collapse=" "),
                 dat[st_end[29,1]:st_end[29,2]],
                 paste(scan(text=dat[Sec_st[30]-1])[1:SSC_yrs],collapse=" "),
                 dat[st_end[30,1]:st_end[30,2]],
                 paste(scan(text=dat[Sec_st[31]-1])[1:SSC_yrs],collapse=" "),
                 dat[st_end[31,1]:st_end[31,2]],
                 paste(scan(text=dat[Sec_st[32]-1])[1:SSC_yrs],collapse=" "),
                 dat[st_end[32,1]:st_end[32,2]])
    for(i in 1:SSC_yrs)   DAT_retro<-c(DAT_retro,paste(scan(text=dat[Sec_st[33]-SSC_yrs-1+i]) ,collapse = " "))
    DAT_retro<-c(DAT_retro,dat[st_end[33,1]:st_end[33,2]])
    
    # Write data and control file
    write.table(DAT_retro, file=here::here(year, model, "retro", 'model',paste0("goa_nr_",endyr,".dat")),
                quote=FALSE,row.names=FALSE,col.names=FALSE)
    
    ctl[2,1] = paste0("goa_nr_", endyr, ".dat")
    ctl[5,1] = as.character(endyr)
    write.table(ctl, file=here::here(year, model, "retro", 'model', paste0(dat_name, ".ctl")),
                quote=FALSE,row.names=FALSE,col.names=FALSE)
    
    #Updated to account for fact that .tpl is looking for 2018.
    # write.table(CTL_retro,file=paste(pathM,"/goa_",species,"_",modelyear,".ctl",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
    
    
    #/\/\/\/\/\/\/\/\ Run model
    #/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
    
    ## set your number of MCMC runs at the top of the program... 
    setwd(here::here(year, model, "retro", "model"))
    
    #Compile the Model
    # compile_admb(MDL_name)
    
    #Determine Operating system
    
    if(mcmcon=="YES") { 
      system(paste0(model_name, '.exe',' -mcmc ',mcmcruns,' -mcsave ',mcmcsave)) 
    }else {
      # system(paste(MDL_name,'.exe ','-nox'))
      run_admb(model_name, verbose=TRUE)
    }
    
    
    #/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
    #/\/\/\/\/\/\/\/\ Get/write results
    #/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
    
    if(mcmcon=="YES") {
      system(paste0(model_name, '.exe',' -mceval'))
      file.copy(from=here::here(year, model, "retro", "model", "evalout.prj"),
                to=here::here(year, model, "retro", "results", paste0("mcmc_", endyr,".std")),
                overwrite=TRUE)
      file.copy(from=here::here(year, model, "retro", "model", paste0(model_name, ".std")),
                to=  here::here(year, model, "retro", "results", paste0("std_", endyr,".std")),
                overwrite=TRUE)
    }
    
    # Compile SSB/recruitment results
    
    STD <- read.delim(here::here(year, model, "retro", "model", paste0(model_name, ".std")), 
                      header = T, sep = "")
    RES_SB[1:nyrs,y] <- STD$value[which(STD$name=="spawn_biom")]
    RES_Rec[1:nyrs,y] <- STD$value[which(STD$name=="pred_rec")]
    
    #---------------------------------------------
    # End of retrospective model running loop
    #---------------------------------------------
  }
  
  write.csv(RES_SB, here::here(year, model, "retro", "results", "RES_SB.csv"))
  write.csv(RES_Rec, here::here(year, model, "retro", "results", "RES_Rec.csv"))
  
  T_end<-Sys.time()
  
  T_end-T_start
  
}

retro_data <- function(year, model, styr = 1961, numretros = 10, nages = 49) {
  
  pathR<-here::here(year, model, "retro", "results")
  modelyear = year
  
  setwd(pathR)
  ssb<-seq(1,(modelyear-styr+1))
  ssb_uci<-seq(1,(modelyear-styr+1))
  ssb_lci<-seq(1,(modelyear-styr+1))
  totbio<-seq(1,(modelyear-styr+1))
  i<-modelyear
  for (i in modelyear:(modelyear-numretros+1)) {
    z<-paste("std_",i,".std",sep="")
    f <- read.delim(z,sep="")
    ssb1<-subset(f$value,f$name=="spawn_biom")
    length(ssb1)<-modelyear-styr+1
    ssb<-cbind(ssb,ssb1)
    
    mcmc<-(read.table(paste("mcmc_",i,".std",sep=""),header=F,sep="")) #,nrow=5000,ncol=50+3*(i-styr+1)+nages-recage)
    mcmc<-mcmc[(0.2*length(mcmc[,1])):length(mcmc[,1]),] # remove burning
    mcmcnames<-  c("sigr","q1",  "q2",  "f40",  "M",  "ssb_next", "ABC",  "obj_fun")
    for(j in styr:i) mcmcnames<- c(mcmcnames,paste("totbio",j,sep="")) 
    for(j in (styr-nages+1):i) mcmcnames<- c(mcmcnames,paste("recdev",j,sep="")) 
    for(j in styr:i) mcmcnames<- c(mcmcnames,paste("ssb",j,sep="")) 
    mcmcnames<-c(mcmcnames,"LMR")
    for(j in (i+1):(i+15)) mcmcnames<-c(mcmcnames,paste("ssbproj",j,sep=""))
    for(j in (i+1):(i+15)) mcmcnames<-c(mcmcnames,paste("catchproj",j,sep=""))
    for(j in (i+1):(i+10)) mcmcnames<-c(mcmcnames,paste("recproj",j,sep=""))
    for(j in (i+1):(i+15)) mcmcnames<-c(mcmcnames,paste("totbioproj",j,sep=""))
    # mcmcnames<-c(mcmcnames,paste("totbioproj",i+1,sep=""))
    names(mcmc)<-mcmcnames
    mcmc_ssb<-mcmc[,which(substr(names(mcmc),1,3)=="ssb")[2]:which(substr(names(mcmc),1,3)=="ssb")[length(seq(styr,i))+1]]
    ssb_uci1<-vector(length=length(mcmc_ssb[1,]))
    ssb_lci1<-vector(length=length(mcmc_ssb[1,]))
    for(j in 1:length(mcmc_ssb[1,])) ssb_uci1[j]<-quantile(mcmc_ssb[,j],0.975)
    for(j in 1:length(mcmc_ssb[1,])) ssb_lci1[j]<-quantile(mcmc_ssb[,j],0.025)
    length(ssb_uci1)<-modelyear-styr+1
    length(ssb_lci1)<-modelyear-styr+1
    ssb_uci<-cbind(ssb_uci,ssb_uci1)
    ssb_lci<-cbind(ssb_lci,ssb_lci1)
  }
  
  
  ssb[,1]<-seq((modelyear-length(ssb[,1])+1),modelyear)
  colnames(ssb)<-c("Year",seq(modelyear,modelyear-numretros+1))
  ssb_uci[,1]<-seq((modelyear-length(ssb_uci[,1])+1),modelyear)
  colnames(ssb_uci)<-c("Year",seq(modelyear,modelyear-numretros+1))
  ssb_lci[,1]<-seq((modelyear-length(ssb_lci[,1])+1),modelyear)
  colnames(ssb_lci)<-c("Year",seq(modelyear,modelyear-numretros+1))
  
  ssb %>% 
    as.data.frame() %>% 
    mutate(id = 'ssb') %>% 
    bind_rows(as.data.frame(ssb_lci) %>% 
                mutate(id = "lci")) %>% 
    bind_rows(as.data.frame(ssb_uci) %>% 
                mutate(id = "uci")) %>% 
    tidyr::pivot_longer(-c(id, Year)) %>% 
    dplyr::rename(year=Year) 
  
}

plot_params <- function(year, model, model_name) {
  
  if (!dir.exists(here::here(year, model, "processed"))){
    stop("must run 'process_results' before creating figures")
  }
  
  read.csv(here::here(year, model, "processed", "ages_yrs.csv"))$yrs -> yrs
  
  ggplot2::theme_set(afscassess::theme_report())
  
  read.delim(here::here(year, model, paste0(model_name, ".std")), sep="", header = TRUE) %>%
    dplyr::filter(name %in% c("q_srv1", "ABC", "nattymort", "tot_biom",
                              "F40","spawn_biom")) %>%
    dplyr::group_by(name) %>%
    dplyr::slice(dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(name = dplyr::case_when(name =="q_srv1" ~ "q_srv",
                                          name =="nattymort" ~ "natmort",
                                          name == "F40" ~ "F",
                                          name == "spawn_biom" ~ "spawn_biom_",
                                          name == "tot_biom" ~ "tot_biom_",
                                          TRUE ~ name),
                  value = dplyr::case_when(name == "ABC" ~ value / 1000,
                                           name == "spawn_biom_" ~ value / 1000,
                                           name == "tot_biom_" ~ value / 1000,
                                           TRUE ~ value),
                  name = factor(name, levels = c("q_srv", "natmort", "F", "ABC", "tot_biom_", "spawn_biom_"))) -> fits
  
  read.csv(here::here(year, model, "processed", "mceval.csv"))  %>%
    dplyr::select(q_srv1, ABC, natmort, paste0("tot_biom_", yrs),
                  F40, paste0("spawn_biom_", yrs)) %>%
    dplyr::mutate(group = 1:dplyr::n()) %>%
    tidyr::pivot_longer(-group) %>%
    dplyr::mutate(years = as.numeric(gsub('\\D+','', name)),
                  name = gsub('[[:digit:]]+', '', name),
                  value = dplyr::case_when(name=="spawn_biom_" ~ value / 1000,
                                           name=="tot_biom_" ~ value / 1000,
                                           name=="ABC" ~ value / 1000,
                                           TRUE ~ value),
                  name = factor(name, levels = c("q_srv", "natmort", "F", "ABC", "tot_biom_", "spawn_biom_"))) -> dat
  
  p1 = dat %>%
    dplyr::filter(name == "q_srv") %>%
    ggplot2::ggplot(ggplot2::aes(value)) +
    # facet_wrap(~name, scales = "free", dir = "v") +
    ggplot2::geom_histogram(ggplot2::aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                            fill = "lightgray", color = "darkgray", bins = 50) +
    # scico::scale_fill_scico(palette = "grayC", direction = -1) +
    # ggplot2::scale_x_continuous(breaks = seq(0,2.5,0.5)) +
    ggplot2::geom_segment(data = dplyr::filter(fits, name == "q_srv"),
                          mapping = ggplot2::aes(x = value, xend = value, y = 0, yend = 0.05), 
                          size = 2, color = "black") +
    ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
    ggplot2::xlab(expression("Trawl survey catchability ("*italic(q)*")")) 
  
  p2 = dat %>%
    dplyr::filter(name == "natmort") %>%
    ggplot2::ggplot(ggplot2::aes(value)) +
    # facet_wrap(~name, scales = "free", dir = "v") +
    ggplot2::geom_histogram(ggplot2::aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                            fill = "lightgray", color = "darkgray", bins = 50) +
    # scico::scale_fill_scico(palette = "grayC", direction = -1) +
    # ggplot2::scale_x_continuous(breaks = seq(0,2.5,0.5)) +
    ggplot2::geom_segment(data = dplyr::filter(fits, name == "natmort"),
                          mapping = ggplot2::aes(x = value, xend = value, y = 0, yend = 0.05), size = 2, color = "black") +
    ggplot2::ylab("Probability density") +
    ggplot2::xlab(expression("Natural mortality ("*italic(M)*")")) +
    ggplot2::theme(legend.position = "none")
  
  p3 = dat %>%
    dplyr::filter(name == "F") %>%
    ggplot2::ggplot(ggplot2::aes(value)) +
    # facet_wrap(~name, scales = "free", dir = "v") +
    ggplot2::geom_histogram(ggplot2::aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                            fill = "lightgray", color = "darkgray", bins = 50) +
    # scico::scale_fill_scico(palette = "grayC", direction = -1) +
    # ggplot2::scale_x_continuous(breaks = seq(0,2.5,0.5)) +
    ggplot2::geom_segment(data = dplyr::filter(fits, name == "F"),
                          mapping = ggplot2::aes(x = value, xend = value, y = 0, yend = 0.05), 
                          size = 2, color = "black") +
    ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
    ggplot2::xlab(expression(italic(F)["40%"])) +
    ggplot2::theme(legend.position = "none")
  
  p4 = dat %>%
    dplyr::filter(name == "ABC") %>%
    ggplot2::ggplot(ggplot2::aes(value)) +
    # facet_wrap(~name, scales = "free", dir = "v") +
    ggplot2::geom_histogram(ggplot2::aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                            fill = "lightgray", color = "darkgray", bins = 50) +
    # scico::scale_fill_scico(palette = "grayC", direction = -1) +
    # ggplot2::scale_x_continuous(breaks = seq(0,2.5,0.5)) +
    ggplot2::geom_segment(data = dplyr::filter(fits, name == "ABC"),
                          mapping = ggplot2::aes(x = value, xend = value, y = 0, yend = 0.05), 
                          size = 2, color = "black") +
    ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
    ggplot2::xlab("ABC (kt)") +
    ggplot2::theme(legend.position = "none")
  
  
  p5 = dat %>%
    dplyr::filter(name == "tot_biom_", years == year) %>%
    ggplot2::ggplot(ggplot2::aes(value)) +
    # facet_wrap(~name, scales = "free", dir = "v") +
    ggplot2::geom_histogram(ggplot2::aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                            fill = "lightgray", color = "darkgray", bins = 50) +
    # scico::scale_fill_scico(palette = "grayC", direction = -1) +
    # ggplot2::scale_x_continuous(breaks = seq(0,2.5,0.5)) +
    ggplot2::geom_segment(data = dplyr::filter(fits, name == "tot_biom_"),
                          mapping = ggplot2::aes(x = value, xend = value, y = 0, yend = 0.05), 
                          size = 2, color = "black") +
    ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
    ggplot2::xlab("Current total biomass (kt)") +
    ggplot2::theme(legend.position = "none")
  
  p6 = dat %>%
    dplyr::filter(name == "spawn_biom_", years == year) %>%
    ggplot2::ggplot(ggplot2::aes(value)) +
    # facet_wrap(~name, scales = "free", dir = "v") +
    ggplot2::geom_histogram(ggplot2::aes(y=(..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]),
                            fill = "lightgray", color = "darkgray", bins = 50) +
    # scico::scale_fill_scico(palette = "grayC", direction = -1) +
    # ggplot2::scale_x_continuous(breaks = seq(0,2.5,0.5)) +
    ggplot2::geom_segment(data = dplyr::filter(fits, name == "spawn_biom_"),
                          mapping = ggplot2::aes(x = value, xend = value, y = 0, yend = 0.05), 
                          size = 2, color = "black") +
    ggplot2::theme(axis.title.y = ggplot2::element_blank()) +
    ggplot2::xlab("Current spawning biomass (kt)") +
    ggplot2::theme(legend.position = "none")
  
  png(filename=here::here(year, model, "figs", "hists.png"), width = 6.5, height = 6.5, 
      units = "in", type ="cairo", res = 200)
  
  print((cowplot::plot_grid(p1, p4, p2,  p5, p3, p6, align = "v", ncol = 2, rel_heights = c(0.5, 0.5))))
  
  dev.off()
}

plot_retro <- function(year, model, n_retro=10) {
  peels = n_retro - 1
  max_year = year
  # loop through mcmc output
  age_yr = read.csv(here::here(year, model, "processed", "ages_yrs.csv"))
  yrs = age_yr %>%
    dplyr::select(yrs) %>% 
    tidyr::drop_na() %>% 
    dplyr::pull(yrs)
  styr_rec = age_yr[1,3] 
  retro_yrs = (year - n_retro + 1):year
  
  dat = list()
  
  for(i in 1:n_retro) {
    
    read.delim(here::here(year, model, "retro", "results",
                          paste0("mcmc_", retro_yrs[i], ".std")),
               sep = "",  header = FALSE) -> df
    
    df = df[(0.2 * nrow(df)):nrow(df),] # drop burn in
    
    colnames(df) = c("sigr", "q_srv1", "q_srv2", "F40", "natmort", "spawn_biom_proj",
                     "ABC", "obj_fun",
                     paste0("tot_biom_", yrs[1]:retro_yrs[i]),
                     paste0("log_rec_dev_", styr_rec:retro_yrs[i]),
                     paste0("spawn_biom_", yrs[1]:retro_yrs[i]),
                     "log_mean_rec",
                     paste0("spawn_biom_proj_", max(retro_yrs[i]) + 1:15),
                     paste0("pred_catch_proj_", max(retro_yrs[i]) + 1:15),
                     paste0("rec_proj_", max(retro_yrs[i]) + 1:10),
                     paste0("tot_biom_proj_", max(retro_yrs[i]) + 1:15))
    
    dat[[i]] = df %>% dplyr::mutate(retro_year = retro_yrs[i])
    
  }
  
  # clean up columns so can bind all together
  col <- unique(unlist(sapply(dat, names)))
  dat <- lapply(dat, function(df) {
    df[, setdiff(col, names(df))] <- NA
    df
  })
  
  do.call(rbind, dat)  -> retro_mc
  
  # save output
  write.csv(retro_mc, here::here(year, model, "processed", "retro_mcmc.csv"), row.names = FALSE)
  
  # functions for quantiles
  q_name <- purrr::map_chr(c(.025,.975), ~ paste0("q", .x*100))
  q_fun <- purrr::map(c(.025,.975), ~ purrr::partial(quantile, probs = .x, na.rm = TRUE)) %>%
    purrr::set_names(nm = q_name)
  
  retro_mc %>%
    dplyr::select(paste0("spawn_biom_", yrs), retro_year) %>%
    tidyr::pivot_longer(c(-retro_year), values_to = "biomass") %>%
    dplyr::mutate(year = as.numeric(stringr::str_extract(name, "[[:digit:]]+")),
                  biomass = biomass / 1000) %>%
    dplyr::group_by(year, retro_year) %>%
    dplyr::summarise_at(dplyr::vars(biomass), tibble::lst(!!!q_fun, median)) %>%
    dplyr::mutate(Retro = factor(retro_year)) %>%
    dplyr::ungroup() -> dat
  

  dat %>%
    dplyr::select(year, retro_year, median) %>%
    dplyr::group_by(year) %>%
    dplyr::mutate(pdiff = (median - median[retro_year==max_year]) /
                    median[retro_year==max_year]) %>%
    tidyr::drop_na() %>%
    dplyr::filter(year %in% (max_year-peels):max_year) %>%
    dplyr::ungroup() %>%
    dplyr::filter(year == retro_year, year !=max_year) %>%
    dplyr::summarise(rho = mean(pdiff)) %>%
    dplyr::pull(rho) -> ssb_rho
  
  
  dat %>%
    # filter(retro_year==2022) %>% 
    ggplot2::ggplot(ggplot2::aes(year, median, color = Retro, group = Retro)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q2.5, ymax = q97.5, fill = Retro),
                         alpha = .05, color = NA) +
    ggplot2::ylab("Spawning biomass (kt)\n") +
    ggplot2::xlab("\nYear") +
    ggplot2::expand_limits(y = 0) +
    scico::scale_fill_scico_d(palette = "roma") +
    scico::scale_color_scico_d(palette = "roma") +
    funcr::theme_report() +
    ggplot2::scale_x_continuous(breaks = afscassess::tickr(dat, year, 10, start = 1960)$breaks,
                                labels = afscassess::tickr(dat, year, 10, start = 1960)$labels) +
    ggplot2::annotate(geom = "text", x=1963, y=Inf, hjust = -0.05, vjust = 2,
                      label = paste0("Mohn's rho = ", round(ssb_rho, 3)),
                      family = "Times") +
    ggplot2::theme(legend.position = "none") -> p1
  
  
  retro_mc %>%
    dplyr::select(paste0("spawn_biom_", yrs), retro_year) %>%
    tidyr::pivot_longer(c(-retro_year), values_to = "biomass") %>%
    dplyr::mutate(year = as.numeric(stringr::str_extract(name, "[[:digit:]]+")),
                  biomass = biomass / 1000,
                  pdiff = (biomass - biomass[retro_year==max_year]) /
                    biomass[retro_year==max_year])  %>%
    dplyr::group_by(year, retro_year) %>%
    dplyr::summarise_at(dplyr::vars(pdiff), tibble::lst(!!!q_fun, median)) %>%
    dplyr::mutate(Retro = factor(retro_year)) %>%
    dplyr::ungroup() -> dat2
  
  
  dat2 %>%
    ggplot2::ggplot(ggplot2::aes(year, median, color = Retro, group = Retro)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q2.5, ymax = q97.5, fill = Retro),
                         alpha = .05, color = NA) +
    ggplot2::ylab(paste0("Percent change from ", year, "\n")) +
    ggplot2::xlab("\nYear") +
    ggplot2::expand_limits(y = 0) +
    scico::scale_fill_scico_d(palette = "roma") +
    scico::scale_color_scico_d(palette = "roma") +
    funcr::theme_report() +
    ggplot2::scale_x_continuous(breaks = funcr::tickr(dat, year, 10, start = 1960)$breaks,
                                labels = funcr::tickr(dat, year, 10, start = 1960)$labels) +
    ggplot2::theme(legend.position = "none") -> p2
  
  png(filename=here::here(year, model, "figs", "retro.png"), width = 6.5, height = 6.5, 
      units = "in", type ="cairo", res = 200)
  
  print((cowplot::plot_grid(p1, p2, ncol = 1)))
  
  dev.off()
}

plot_compare_survey <- function(year, models = c('2022, m22')) {
  
  if (!dir.exists(here::here(year, "compare_models"))){
    dir.create(here::here(year, "compare_models"))
  }
  
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
    ggplot2::geom_point(data = dplyr::filter(dat, name == "Observed"), position=ggplot2::position_dodge(width=0.7)) +
    ggplot2::geom_errorbar(data = dplyr::filter(dat, name == "Observed"),
                           ggplot2::aes(ymin = lci, ymax = uci, color = model), width = 0.4, position=ggplot2::position_dodge(width=0.7)) +
    ggplot2::geom_line(data = dplyr::filter(dat, name == "Predicted"),
                       ggplot2::aes(color = model)) +
    scico::scale_color_scico_d(palette = 'batlow', begin = 0.2, end = 0.8) +
    ggplot2::scale_x_continuous(breaks = funcr::tickr(dat, year)$breaks,
                                labels = funcr::tickr(dat, year)$labels) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::xlab("Year") +
    ggplot2::ylab("Survey biomass (kt)") +
    ggplot2::expand_limits(y = 0) +
    funcr::theme_report() +
    ggplot2::theme(legend.justification=c(1,0),
                   legend.position=c(0.8,0.70))
  
  # ggplot2::ggsave(here::here(year, "compare_models", "figs", "compare_est_survey.png"),
  # width = 6.5, height = 6.5, units = "in", dpi = 200)
}
plot_compare_biomass <- function(year, models = c('2022, m22')) {
  
  if (!dir.exists(here::here(year, "compare_models"))){
    dir.create(here::here(year, "compare_models"))
  }
  
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
    funcr::theme_report() +
    ggplot2::theme(legend.position = c(0.2, .8))
  
  # ggplot2::ggsave(here::here(year, "compare_models", "figs", "compare_est_biomass.png"),
  # width = 6.5, height = 6.5, units = "in", dpi = 200)
}

plot_catch_bio <- function(year, model, model_name) {
  
  rep <-  readLines(here::here(here::here(year, model, paste0(model_name, ".rep"))))
  
  std <- read.delim(here::here(here::here(year, model, paste0(model_name, ".std"))), sep="", header = TRUE)
  
  # filter(catch, Year == year) %>%
  #   left_join(read.delim(here::here(here::here(year, model, "proj", "author_f", "bigsum.dat")), sep="", header = TRUE) %>%
  #               filter(Year == year, Alt == 2) %>%
  #               dplyr::select(Year, value = Total_Biom))
  
  std %>% 
    filter(name=="tot_biom") %>%
    bind_cols((catch)) %>% 
    filter(Year >= 1990, Year!=year) %>%
    dplyr::select(Year, Catch, value, std.dev) %>%
    bind_rows(filter(catch, Year == year) %>%
                left_join(read.delim(here::here(here::here(year, model, "proj", "author_f", "bigsum.dat")), sep="", header = TRUE) %>%
                            filter(Year == year, Alt == 2) %>%
                            mutate(value = Total_Biom * 1000) %>%
                            dplyr::select(Year, value))) %>%
    mutate(std.dev = ifelse(is.na(std.dev), std.dev[Year==year-1], std.dev)) %>% 
    mutate(lci = value - std.dev * 1.96,
           uci = value + std.dev * 1.96) %>% 
    mutate(perc = Catch / value,
           lci = Catch / lci,
           uci = Catch / uci,
           mean = mean(perc)) %>% 
    dplyr::select(Year, value, mean, perc, lci, uci) -> df 
  
  # df %>% 
  #   summarise(min = min(perc),
  #             max = max(perc))
  
  png(filename=here::here(year, model, "figs", "catch_bio.png"), width = 6.5, height = 6.5, 
      units = "in", type ="cairo", res = 200)
  df %>% 
    ggplot(aes(Year, perc)) + 
    geom_point() +
    geom_line() +
    geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
    geom_hline(yintercept = df$mean, lty = 3) +
    expand_limits(y = c(0, 0.08)) +
    scale_x_continuous(breaks = funcr::tickr(df, Year, start = 1990)$breaks,
                       labels = funcr::tickr(df, Year, start = 1990)$labels) +
    # theme_present() +
    xlab("Year") +
    ylab("Catch/Biomass\n") 
  
  dev.off()
}

recruit_tbl <- function(year, model, model_name, rec_age){
  
  # hard to filter year with year so change the name
  mod_year = year
  
  if (!dir.exists(here::here(year, model, "processed"))){
    stop("must run 'process_results' before creating tables")
  }
  
  # read in data
  REP <- readLines(here::here(year, model, paste0(model_name, ".rep")))
  STD <- read.delim(here::here(year, model, paste0(model_name, ".std")), sep="", header = TRUE)
  
  suppressWarnings(data.frame(year = unlist(strsplit(REP[grep("Year", REP)[1]]," "))) %>%
                     dplyr::mutate(year = as.numeric(year)) %>%
                     tidyr::drop_na() %>%
                     dplyr::pull(year)) -> yrs
  
  bio_rec = read.csv(here::here(year, model, "processed", "bio_rec_f.csv")) %>%
    dplyr::select(-F) %>%
    dplyr::mutate(recruits = round(recruits * 1000))
  
  bio_rec %>%
    dplyr::filter(year %in% (1977 + rec_age):(mod_year - rec_age)) %>%
    dplyr::summarise(recruits = round(mean(recruits))) -> pred_rec
  
  data.frame(year = mod_year + 1:2,
             tot_biom = STD$value[which(STD$name=="tot_biom_proj")][1:2],
             sp_biom = STD$value[which(STD$name=="spawn_biom_proj")][1:2],
             recruits = pred_rec$recruits) -> std_data
  
  values = dplyr::bind_rows(bio_rec, std_data) %>%
    dplyr::mutate_all(dplyr::funs(round(.)))
  
  
  
  # get mcmc data - clean it, calculate annual uci and lci
  read.csv(here::here(year, model, "processed", "mceval.csv")) %>%
    dplyr::select(dplyr::starts_with(c( "tot_biom", "spawn_biom", "log_rec_dev", "rec_proj"))) %>%
    dplyr::mutate(id = 1:dplyr::n()) %>%
    dplyr::mutate_if(is.character, dplyr::funs(as.numeric(gsub(",", "", .)))) %>%
    tidyr::pivot_longer(-id) %>%
    dplyr::mutate(value = ifelse(grepl("log", name), exp(value), value),
                  year =  as.numeric(stringr::str_extract(name, "[[:digit:]]+")),
                  name = dplyr::case_when(stringr::str_detect(name, "tot_biom") ~ "tot_biom",
                                          stringr::str_detect(name, "spawn_biom") ~ "sp_biom",
                                          stringr::str_detect(name, "log_rec") ~ "recruits",
                                          stringr::str_detect(name, "rec_proj") ~ "recruits")) %>%
    group_by(year, name) %>% 
    dplyr::summarise(lci = quantile(value, 0.025),
                     uci = quantile(value, 0.975)) %>% 
    dplyr::left_join(values, .) %>% 
    dplyr::filter(year >= 1977 & year <= mod_year + rec_age) %>%
    dplyr::mutate(tot_lci = ifelse(name == 'tot_biom', lci, NA),
                  tot_uci = ifelse(name == 'tot_biom', uci, NA),
                  sp_lci = ifelse(name == 'sp_biom', lci, NA),
                  sp_uci = ifelse(name == 'sp_biom', uci, NA),
                  rec_lci = ifelse(name == 'recruits', lci * 1000, NA),
                  rec_uci = ifelse(name == 'recruits', uci * 1000, NA)) %>% 
    group_by(year) %>% 
    dplyr::summarise(recruits = mean(recruits, na.rm = T),
                     rec_lci = mean(rec_lci, na.rm = T),
                     rec_uci = mean(rec_uci, na.rm = T),
                     tot_biom = mean(tot_biom, na.rm = T),
                     tot_lci = mean(tot_lci, na.rm = T),
                     tot_uci = mean(tot_uci, na.rm = T),
                     sp_biom = mean(sp_biom, na.rm = T),
                     sp_lci = mean(sp_lci, na.rm = T),
                     sp_uci = mean(sp_uci, na.rm = T)) %>%
    write.csv(here::here(year, model, "tables", "ts.csv"), row.names = FALSE)
  
}


plot_F <- function(year, model) {


read.csv(here::here(year, model, "processed", "bio_rec_f.csv")) %>%
  dplyr::select(year, F) -> dat

png(filename=here::here(year, model, "figs", "fF.png"), width = 6.5, height = 6.5, 
    units = "in", type ="cairo", res = 200)

dat %>%
  ggplot2::ggplot(ggplot2::aes(year, F)) +
  ggplot2::geom_line() +
  ggplot2::expand_limits(y=0) + 
  ggplot2::ylab("Fishing mortality rate (F)\n") +
  ggplot2::scale_x_continuous(name = "Year",
                              breaks = funcr::tickr(dat, year, 10, start = 0)$breaks,
                              labels = funcr::tickr(dat, year, 10, start = 0)$labels) 

dev.off()
}
