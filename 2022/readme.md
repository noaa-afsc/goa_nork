# Northern rockfish assessment order of operations

1. Code available in `northern_rockfish/year/code/year_nr_analysis.R`
2. Pull data using `gfdata` package   
   `devtools::install_github("BenWilliams-NOAA/gfdata")` 
3. Setup folder structure
4. Process data using `groundfishr` package, `afscassess` package and a couple custom functions  
   `devtools::install_github("BenWilliams-NOAA/groundfishr")`   
    `devtools::install_github("BenWilliams-NOAA/afscassess")` 
5. For each trawl survey biomass estimate create a different biomass output w/id (`ts_biomass(year, id="db")`)
6. Create a `.dat` file for each seperate model
7. Manually paste model `.tpl` and `.ctl` file in each separate model folder (also maturity data if using)  
   a. Change name of `.ctl` file in the `.tpl` file  
   b. Change name of `.dat` file in the `.ctl` file

8. Compile model via ADMB command line  
   1. set location
   2. admb modelname
   3. model_name -mcmc 10000000 -mcsave 2000
10. `process_results()`
