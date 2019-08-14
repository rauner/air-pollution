# ------------------------------------------------------------------------------
# Program Name: launch_AP.R
# Author(s): Sebastian Rauner
# Date Last Updated: Oct 13, 2018
# Program Purpose: The script launches the AP routine - soon to be a module
# Input Files: 
# Output Files:

# Info: for a detailed description check: Exposure assessment for estimation of the global burden of disease attributable to outdoor air pollution
#                                         & 10.1056/NEJMoa1506699
# To-DO:
# check the urban increment
# there is still a hardcoded year list in the downscaling and harmonization script
rm(list=ls())
closeAllConnections()


wd_dir <- "C:/Users/rauner/Documents/PIK/fasstr_sebastian"
setwd(wd_dir)


# which health impact assessment , lo is theta minus standart error...
u_mort_lvl <- 'med'     # choose between lo, med and hi 


# chose here which downscaling method should be used
GDP_downscaling          = FALSE
IIASA_aneris_downscaling = TRUE

#do the IIASA gridding
# master config gridding start year needs be adjusted
# the current pop(2015) is used
IIASA_gridding           = FALSE
grid_monthly             = FALSE
grid_yearly              = FALSE

#do the emission graphs for the unharmonized REMIND results
emission_graphs_unharmonized        = FALSE
emission_graphs_harmonized          = FALSE

#do map plots
do_map_plots                       = FALSE

#this section is not up to date, check map_plots.R for plotting possibilities
do_plots_yearly                    = FALSE
do_plots_summary_one_year          = FALSE
do_plots_summary_one_year_regional = FALSE
enavi_plots                        = FALSE
paper_plots_concentrations         = FALSE
paper_plots_deathes                = FALSE
paper_plots_deathes_relativ_to_2015_and_INDC = FALSE

#do the netcdf writeout
netcdf_write_out_yearly  <- T
netcdf_write_all_years   <- FALSE  #adds a lot of time




remind_output <- paste0(getwd()    ,'/mnt/FASST_IDL_PACK/REMIND_OUTPUT/new_ENavi')

if(emission_graphs_unharmonized == TRUE){source('emission_graphs_unharmonized.R')}




# flag that the downscaling script should use the REMIND -in as an input instead of the standart .mif IAM snapshot
# needs to be TRUE for the IIASA_aneris_downscaling to work
fasstr_connection        = TRUE

u_year     <- c(seq(2000,2015,5),seq(2020,2100,20))                  # this should be changed & harmonized in all scripts to use the REMIND output


# for the harmonization and harmonization plots
u_year     <- c(2015,seq(2020,2050,10),2100) 
# for the concentration and health impact assessment
#u_year     <- c(2015,2020, 2030,2040,  2050, 2100)
u_year     <- c(2015,2030,2050,2100)


# check what remind results are in the REMIND_OUTPUT folder
# loop over them
scenario_folders <- paste0(list.dirs(paste0(remind_output), recursive = FALSE),'/')
for( scenario_folder in scenario_folders[4]) {
  # a mix of the IIASA downscaling and the aneris harmonization scripts
  # https://github.com/iiasa/emissions_downscaling
  # https://github.com/gidden/aneris
  setwd(wd_dir)
  # reads in the REMIND gdx and does some formatting
  # there is a lot of defining and parameter setting in here
  source('quick_fasst_gridmaps_pm_gas_REMIND_read_in.R')
  
  
  # use cache when possible, for GDP and pop data
  setConfig(forcecache=T)
  
  # redundant from REMIND_read_in.R  
  downscaling_method <- 'IIASA_aneris_downscaling'   
  RUNSUFFIX <- paste0(AP_scenario)
  RUNSUFFIX <- paste0(RUNSUFFIX,'_', p_projectName,'_', ssp_scenario,'_', downscaling_method,'_')
  scenario <- sapply(strsplit(RUNSUFFIX,'_'), "[[", 2)
  
  # launches the downscaling and harmonization
  source('quick_fasst_gridmaps_pm_gas_downscaling_harmonization.R')
  
  setwd(wd_dir)
  
  #}
  closeAllConnections()
  #calculates the concentrations, times 10^9 for Mt to kg
  source('quick_fasst_gridmaps_pm_gas_concentration.R')
  
  
  # calculate the health impacts
  # currently writes the results grid, should be transferred to the IIASA gridding routine?
  source('quick_health_impacts_from_fasst_gridmaps_cluster_Seb.R')
  
  # the dalys
  source('calcMorbidity.R')
  
}



#PLOTS
if(do_map_plots == TRUE){source('map_plots.R')}

if(emission_graphs_harmonized   == TRUE){source('emission_graphs_harmonized.R')}
closeAllConnections()
