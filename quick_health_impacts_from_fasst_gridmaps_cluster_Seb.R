# ------------------------------------------------------------------------------
# Program Name: quick_health_impacts_from_fasst_gridmaps_cluster_Seb.R
# Author(s): Sebastian Rauner
#            original IDL script written by Rita van Dingenen (JRC-ISPRA)
#            rewritten in R by Jerome Hilaire (MCC & PIK) & Sebastian Rauner (PIK)
# Date Last Updated: Oct 17, 2017
# Program Purpose: This script reads the concentrations & mortality maps and calculates premature deathes 
# Input Files: 
# Output Files:
# To-Do: use the IIASA gridding routine?

#===============================================================================
# FASST-R
#===============================================================================
# original IDL script written by Rita van Dingenen (JRC-ISPRA)
# rewritten in R by Jerome Hilaire (MCC & PIK) & Sebastian Rauner (PIK)
#
# Description:
# Process project PM and O3 fields to calculate mortalities
# Needs at least 1 scenario and associated year for correct pop and base
# mortality stats. Stores mortalities per region for each scenario.
# Delta's are not calculated.
# include rrate function to calculate RR from PM2.5 using Burnett's functions
#
# WARNING:
# population and crops, but it is recommended to produce final results at the
# aggregation level at which input emissions were generated, in particular when
# looking at u_projectNameed scenarios where trends are established based on
# indicators for country groups.
#
# NOTES:
# CH4 feedback on global background O3 is included, but not longterm feedback of
# NOx and other precursors on CH4 lifetime and background O3
# CO Source-Receptors not yet included - work in progress.
#
# GLOSSARY:
# COPD  : Chronic Obstructive Pulmonary Disease (11%)
# IHD   : Ischaemic Heart Disease (40%)
# ALRI  : Acute Lower Respiratory Infections (in children) (3%)
# STROKE: Cerebrovascular disease (40%)
# RESP  :
# LC    : Lung cancer (6%)
#
# Sebastian
# added a possibility to write out results
# full write out also with pop = 0
# fraction adjusted to percentage through /100
# added /100 to all but m6m maps
# 
#
#===============================================================================
t00 <- Sys.time()
#setwd("~/fasstr_sebastian")
#===============================================================================
# USER-SECTION
#===============================================================================
# Specify scenario list
#u_year     <- c(seq(2010,2100,10))     # needed for population stats. needs to be changed according to REMIND output
# 2000 & 2005 might not be available... adjust


u_mort_lvl <- 'med'     # choose between lo, med and hi

u_VSL     <- 3600000    # cost of a premature deaths (VSL)
# see calcMonetization() for parameters such as base year and elasticity


u_zcf     <- 2   # default=2 uses Burnett built-in zcf's
u_urb_inc <- 2   # default=2 (urban increment included)

# Input data options
u_RData         <- TRUE   # Generate RData file to speed things up
u_recreateRData <- FALSE  # Re-generate RData file (other)


doMortalitiesByTM5region <- FALSE



#write these two infos into netcdf and table
u_urb_inc_mort           <- FALSE  # use urban inc for mortality calculation
u_WHO_harmonize          <- TRUE   # harmonize PM to the WHO 2014 PM25 values
# PM(all) is harmonized, this is used for the PM_mort and VSL calculation

# Debug
u_verbose       <- FALSE

# Define dimensions lon and lat @ 2880x1440
img <- 2880
jmg <- 1440
lon <- ((seq(img)-1)-img/2)/(img/360)
lat <- ((seq(jmg)-1)-jmg/2)/(jmg/180)

#===============================================================================
# INITIALISATION
#===============================================================================
t0 <- Sys.time()
if(u_verbose) cat("Initialising...\n")
# Load libraries
library(ncdf4)
library(fields)
library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)

# assigne SSP and year from REMDIN or gas script

dir_ancil                    <- paste0(dir_root,  'ANCILLARY/')                    # ancillary data
dir_spat_map                 <- paste0(dir_ancil, 'SPATIAL_MAPPING/')              # Spatial mapping

# Files
file_abase                   <- paste0(dir_ancil, 'BASE/nTM5-JRC-cy2-ipcc-v1_GLOBAL_SFC_aerosolm_2001.nc')
file_gbase                   <- paste0(dir_ancil, 'BASE/TM5-JRC-cy2-ipcc-v1_GLOBAL_SFC_tracerm_2001.nc')
file_dustsea                 <- paste0(dir_ancil, 'BASE/GBD2010_PM25_DU_SS.nc')
file_m6m                     <- paste0(dir_ancil, 'BASE/BASE_GLOBAL_M6M_boxmid.nc')
file_mf                      <- paste0(dir_ancil, 'O3_HEALTH/m3m_to_m6m_factor.nc')
dir_mort_map_grid_mask       <- paste0(dir_ancil, '/FASST_REGION_MASK')                   # mort_map_grid_mask
dir_plots                    <- paste0(getwd(),   '/plots')      
dir_harmonization            <- paste0(dir_ancil, '/CONC_HARMONIZATION')    



file_embase  <- paste0(dir_ancil, 'BASE/TOTAL_COUNTRY_EM.RCP.TXT')
load(file=file.path(dir_mort_map_grid_mask,paste0("dim_df_maps_mort_mapping_dim_id.RData"))) #full lat lon dimension



remind_emiss       <- paste0(dir_remind,     'pm_emiAPexsolve.gdx')
remind_config      <- paste0(dir_remind,     'config.Rdata')
mapping_REMIND     <- paste0(dir_spat_map,   'regionmappingREMIND.csv')
mapping_TM5        <- paste0(dir_spat_map,   'regionmappingTM5.csv')


# Load own functions
dump <- lapply(list.files("functions", pattern = "*.R", full.names = TRUE), source)

# Output NetCDF file parameters
p_progName         <- 'health_impacts_from_fasst_gridmaps.R'

# Folders
dir_root           <- paste0(getwd()    ,'/mnt/FASST_IDL_PACK/')   # main folder


dir_out    <- paste0(dir_root,  'OUTPUT/')                 # output
dir_outin  <- paste0(dir_root,  'OUTPUT/NCDF')             # output (ncdf files)
dir_outtab <- paste0(dir_root,  'OUTPUT/TABLES')           # output (tables)
dir_outcon <- paste0(getwd(),  '/concentrations')          # output concentrations
dir_out    <- dir_outin


dir_remind <- paste0(dir_root,  'REMIND_OUTPUT/')          # REMIND  output emissions
dir_who    <- file.path(dir_ancil, 'MORTALITY')               # mortality data


# Create folders (if necessary)
dir.create(dir_out, showWarnings = TRUE)



# Model parameters
H2O_FL <- 1.   # include H2O at 35% RH
m6mthr <- 33.3 # counterfactual level for M6M for O3 health impact - GBD2013
# uses between 33.3 and 41.9 for M3M

# Define mortality rate functions (IER)
# The Burnett functions assume a generic natural background threshold below
# which the risk is not increased (RR=1)
# With FASST we can actually distinguish between anthropogenic and natural
# (Dust and SS) and calculate the mortalities without the counterfactual
# threshold from actual anthropogenic PM.
# To override Burnett counterfactual concentration (threshold), set following
# flag to 1 (default=2, standard Burnett method) and provide a custom
# zcf value for the threshold (Anenberg 2010 used 5.8ug/m3 and 0ug/m3 for
# 'threshold' and 'background' runs respectively)
#zcf      <- 5.8  # when flag=2, ZCF is threshold for Anenberg method only.
zcf      <- 2.4  # when flag=2, ZCF is threshold for Anenberg method only.
# This is the values subtracted from total PM.

if(u_zcf != 2) {
  u_zcf <- 1
  #custom threshold when flag = 0. Will be applied to both Burnett and Anenberg.
  zcf <- 0
}
zcf_lab <- c('CUSTOM Zcf for Burnett & Anenberg', 'BURNETT inherent Zcf, Custom Zcf Anenberg')

zfl <- c('_ZCF-OFF','_ZCF-ON')
if (u_verbose) cat(paste('COUNTERFACTUAL VALUE FLAG:', zfl[u_zcf], '\n'))

# Urban Increment
if(u_urb_inc != 2) u_urb_inc <- 1
UFL <- c('_URBINCR_OFF','_URBINCR_ON')
if (u_verbose) cat(paste('URBAN INCREMENT FLAG :', UFL[u_urb_inc], '\n'))

# TODO: Link to our SSP population data?
SSP_POP_FLAG <- 1   # use SSP population data

SSPFL <- c('_POP-SSP')

# Initialise lists
p_rr            <- list()
p_beta          <- list()
af              <- list()

#write out list
write_out <- list()

# results list
results_list_names = c(colnames(dim_df_maps_mort_mapping_dim_id),'iso$iso.y','data_df_map_pop$pop','data_df_mort_pm$frac_05','data_df_mort_pm$frac_30',
                       'data_df_mort_pm$copd1','data_df_mort_pm$ihd1','data_df_mort_pm$lc1','data_df_mort_pm$stroke1',
                       'data_df_mort_pm$alri1',
                       'data_df_mort_pm$value','data_df_mort_pm$VSL_cost', 'data_df_map_pm$all', 'data_df_map_pm$anth',
                       'data_df_map_pm$all_no_urb','data_df_mort_pm$all_mort','data_df_map_pm$anth_no_urb', 'data_df_mort_o3$value',
                       'data_df_map_m6m$value', 'data_df_map_pop$urb_inc_density_fac',
                       'data_df_map_pop$pop_urb', 'data_df_map_pop$pop_rur')


# Initialise results array

results_array <- array(NaN, c(length(u_year),nrow(dim_df_maps_mort_mapping_dim_id), length(results_list_names)),
                       dimnames = list(u_year,
                                       row.names(dim_df_maps_mort_mapping_dim_id),
                                       results_list_names))


# Create high resolution version at same grid size as population maps
# map onto full lat range and regrid to 7.5'x7.5'
maps_dummy   <- array(0.0, c(1440,720))
maps_dummylo <- array(0.0, c(1440,720))
maps_dummyup <- array(0.0, c(1440,720))
latful  <- (0:719-360)/4.
lonful  <- (0:1439-720)/4.



# note: available population totals from SSP: 2010, 2020, ..., 2100
# available base mortalities and fraction of pop <5yr and >30yr: 2005 2010 2015
# 2030 2050.
# interpolation is carried out in the transform_data_ncdf_to_dataframe_ssp_urban...R. script
tend <- Sys.time() - t0
if(u_verbose) cat(paste0("[Time required to initialise: ", tend, "s]\n"))

#===============================================================================
# READ IN DATA
#===============================================================================
if(u_verbose) cat("Reading in data...\n")
t0 <- Sys.time()
# Old method: RRs with log-lin ER function (KREWSKI, see Anenberg et al 2010)
p_rr$cp       <- 1.129
p_rr$lc       <- 1.137
p_rr$resp_med <- 1.04
p_rr$resp_lo  <- 1.013
p_rr$resp_hi  <- 1.067

# Compute Beta coefficients. Beta = log(RR)/10
#p_beta$cp  <- log(p_rr$cp)/10.  #=0.012133229
#p_beta$lc  <- log(p_rr$lc)/10.  #=0.012839321
for (ktyp in names(p_rr)) {
  p_beta[[ktyp]] <- log(p_rr[[ktyp]])/10.
}

# Function coefficients for new Burnett IER risk rate functions. Includes the built-in counterfactuals zcf.
#rrcoef <- read.csv3(file.path(dir_who,'ER_FUNCTION_FITTINGS_ALL_ENDPOINTS.csv'), RData=u_RData, recreate_rdata=u_recreateRData, verbose=u_verbose)
rrcoef <- read.csv3(file.path(dir_who,'ER_FUNCTION_FITTINGS_ALL_ENDPOINTS_2018.csv'), RData=u_RData, recreate_rdata=u_recreateRData, verbose=u_verbose)
# Set the zcf value to zcf if non-threshold calculation is selected
if (u_zcf == 2) {
  #subsitute zcf in the structure containing the retrieved Burnett coefficients
  rrcoef[which(rrcoef$PARAM == "zcf"), -1] <- zcf
}

# Read in mortality maps
maps_mort <- load(file=file.path(dir_ancil, "MORTALITY", paste0(tolower(ssp_scenario),"_MortalityMaps_hires_df.RData")))


# Get mortality map dimensions
lons        <- dim_df_maps_mort$lon
lats        <- dim_df_maps_mort$lat
yearcatalog <- dim_df_maps_mort$time

# Read in FASST countries mapping and masks
country_mapping <- read.csv3(file.path(dir_ancil, "FASST_REGION_MASK", "country_mask_0.5x0.5_v3_mapping.csv"), RData=u_RData, recreate_rdata=u_recreateRData, verbose=u_verbose)
load(file.path(dir_ancil, "FASST_REGION_MASK", "country_mask_0.5x0.5_v3_df.RData"))

# Get FASST regions ID <> Name mapping
mapping_FASST_regions <- data.frame(
  id   = unique(country_mapping$fasst_reg_id),
  name = paste(sapply(unique(country_mapping$fasst_reg_id), function(x) {country_mapping$fasst_reg_name[which(country_mapping$fasst_reg_id == x)[1]]})),
  stringsAsFactors = FALSE
)

# Read in population data
# change
load(file.path(dir_ancil, "SSP_POP", paste0(tolower(ssp_scenario), "_pop.0.125deg_allYears_df.RData")))


#single year concentration read-in
#load concentrations and write in a list
data_df_maps_ap_years <- list()


 for (year in u_year) {
   load(file.path(dir_outcon, paste0(year,'_', RUNSUFFIX, file_output_con)))
   data_df_maps_ap_years[[paste(year)]] <- data_df_maps_ap
 
 
 }

#load(file.path(dir_outcon, paste0(RUNSUFFIX, file_output_con)))


# calculate the VSL coefficients
VSL = u_VSL * calcMonetization()


# this is the spatial explicit VSL, we can avoid the computationlly expensive joining if we only calculat the cost on a country level
# add the calcMonetization() columns to the data_df_TM5_SSP_REMIND_cell_mapping which is later joined to the mort
## load the ISO grid mapping
load(file=file.path(dir_ancil, paste0("/SPATIAL_MAPPING/mapping_iso_grid.RData")))

VSL <- as.data.frame(VSL) 
VSL <- VSL[,c('Region', 'Year', 'Value')]
VSL <- spread(VSL, 'Year', 'Value')

## join the cells with the VSL values
VSL <- left_join(data_df_TM5_SSP_REMIND_cell_mapping[,c('lon_id','lat_id','lon_id_coarse2','lat_id_coarse2', 'iso')],
                 VSL, by = c('iso' = 'Region'))

VSL <- VSL[,c('lon_id','lat_id','lon_id_coarse2','lat_id_coarse2','iso',u_year)]

tend <- Sys.time() - t0
if(u_verbose) cat(paste0("[Time required to read in data: ", tend, "s\n"))

#===============================================================================
# PROCESS AND SAVE DATA (text and NetCDF files)
#===============================================================================
if(u_verbose) cat("Processing data...\n")
t0 <- Sys.time()
# Initialise output text file containing summary information
#TODO: Add formating
txtfilepath <- file.path(dir_outtab, paste0(RUNSUFFIX, '.txt'))

cat('Change with V4: use COPD for O3\n',                           file=txtfilepath)
cat(paste0('programme: ', RUNSUFFIX, '\n'),                        file=txtfilepath, append=TRUE)
cat(paste0('date of run: ', Sys.time(), '\n'),                     file=txtfilepath, append=TRUE)
cat(paste0(p_projectName, UFL[u_urb_inc], zfl[u_zcf], '\n'),       file=txtfilepath, append=TRUE)
cat(paste0(zfl[u_zcf], '\n'),                                      file=txtfilepath, append=TRUE)
cat(paste0('Custom zcf = ', zcf, '\n'),                            file=txtfilepath, append=TRUE)
cat(paste0('H2O_FL = ', H2O_FL, '\n'),                             file=txtfilepath, append=TRUE)
cat(paste0('M6M threshold = ', m6mthr, '\n'),                      file=txtfilepath, append=TRUE)
cat('BU = BURNETT new GBD functions, AN = ANENBERG old method\n',  file=txtfilepath, append=TRUE)
cat(paste('SCENARIO','YEAR','FASST_REGION','POP_CNT','POPW_SO4','POPW_NO3','POPW_NH4','POPW_BC',
          'POPW_POM','POPW_PM2.5','POPW_DUST','POPW_SS','POPW_PMtot','POPW_M6M','BU_PM_MORT_MED',
          'BU_PM_MORT_LO','BU_PM_MORT_UP','AN_PM_MORT', 'O3_MORT_MED','O3_MORT_LO','O3_MORT_UP'),
    file=txtfilepath, append=TRUE)

# Define latitude indices where data can be found
indmin  <- which(latful == min(lats))
indmax  <- which(latful == max(lats))

#for memory reasons
#save workspace here
#save.image('health_ws.RData')


# Loop over years and compute mortality rates
# TODO: Parallelize?
# probably the who conc fac harmonization needs to be out of the loop

i <- 1
if (u_verbose) cat("Loop over years to compute mortality rates...\n")
for (year in u_year) {
  
  
  #for memory reasons
  #delete all
#  rm(list=ls())
#  closeAllConnections()
  
  #load ws from above
#  load('health_ws.RData')
  
  if (u_verbose) cat(paste0("  > Year: ", paste(year), "\n"))
  if (u_verbose) cat(paste0("    - Processing population and mortality data...\n"))
  t02 <- Sys.time()
  # Read scenario population file
  
  
  data_df_map_pop <- data_df_maps_pop[[paste0(year)]]
  
  
  cat(paste(paste(year), "population :",  paste(sum(as.numeric(na.omit(data_df_map_pop$pop))*1e-9), "billion\n")))
  
  # select correct year for mortality and pop stats
  
  maps_mort_hires <- list() 
  # Required by the Anenberg/Krewski function
  maps_mort_hires$cp <- cbind(
    data_df_maps_mort_hires$copd_map[[paste(year)]]                     %>% dplyr::rename(copd   = value),
    data_df_maps_mort_hires$alri_map[[paste(year)]]   %>% dplyr::select(value) %>% dplyr::rename(alri   = value),
    data_df_maps_mort_hires$ihd_map[[paste(year)]]    %>% dplyr::select(value) %>% dplyr::rename(ihd    = value),
    data_df_maps_mort_hires$stroke_map[[paste(year)]] %>% dplyr::select(value) %>% dplyr::rename(stroke = value)) %>%
    dplyr::mutate(value = copd + alri + ihd + stroke) %>%
    dplyr::select(-copd,-alri,-ihd,-stroke)
  
  tend2 <- Sys.time() - t02
  if(u_verbose) cat(paste0("[Time required to process population and mortality data: ", tend2, "s]\n"))
  
  # Read the scenario gridfile with air pollutant concentrations (PM & map_apO3)
  # (generated with fasst_gridmaps_pm_gas.R)
  # TODO: Save AP concentration data as data frame
  if (u_verbose) cat("    - Reading in air pollutants concentrations...\n")
  t02 <- Sys.time()
  #maps_ap <- read.ncdf(file.path(dir_outin, paste0(u_scen[isc], ".nc")), RData=u_RData, recreate_rdata=u_recreateRData, verbose=u_verbose)
  
  
  
  # Sum up aerosols
  # put into
  # Add Dust (DU) and Sea-Salt (SS) because new RR relations are non-linear
  data_df_map_pm <- cbind(
    data_df_maps_ap_years[[paste(year)]]$so4  %>% dplyr::rename(so4  = value),
    data_df_maps_ap_years[[paste(year)]]$no3  %>% dplyr::rename(no3  = value) %>% dplyr::select(no3),
    data_df_maps_ap_years[[paste(year)]]$nh4  %>% dplyr::rename(nh4  = value) %>% dplyr::select(nh4),
    data_df_maps_ap_years[[paste(year)]]$pom  %>% dplyr::rename(pom  = value) %>% dplyr::select(pom),
    data_df_maps_ap_years[[paste(year)]]$pm25 %>% dplyr::rename(pm25 = value) %>% dplyr::select(pm25),
    data_df_maps_ap_years[[paste(year)]]$du   %>% dplyr::rename(du   = value) %>% dplyr::select(du),
    data_df_maps_ap_years[[paste(year)]]$ss   %>% dplyr::rename(ss   = value) %>% dplyr::select(ss)) %>%
    dplyr::mutate(all =so4+no3+nh4+pom+pm25+du+ss + H2O_FL*(0.27*(so4+no3+nh4)+0.15*ss)) %>%
    dplyr::mutate(anth=all-du-ss) %>%
    dplyr::select(-so4,-no3,-nh4,-pom,-pm25,-du,-ss) %>% 
    filter(all >= 0)
  
  data_df_map_m6m <- data_df_maps_ap_years[[paste(year)]]$m6m
  
  # Downscale AP concentrations from 360x180 to 2880x1440
  # shouldnt we interpolate here  
  data_df_map_pm <- data_df_map_pop %>%
    left_join(data_df_map_pm, by=c("lon_id_coarse8"="lon_id", "lat_id_coarse8"="lat_id"))
  data_df_map_m6m <- data_df_map_pop %>%
    left_join(data_df_map_m6m, by=c("lon_id_coarse8"="lon_id", "lat_id_coarse8"="lat_id"))
  tend2 <- Sys.time() - t02
  if(u_verbose) cat(paste0("[Time required to read and process AP concentration data: ", tend2, "s]\n"))
  
  # correction through the urban & rural increment (urb_inc_density_fac)
  # here only for the total concentrations of PM and m6m
  
  # data_df_map_pm$all_no_urb  <- data_df_map_pm$all  
  # data_df_map_pm$anth_no_urb <- data_df_map_pm$anth
  # 
  # #data_df_map_m6m$value_no_urb <- data_df_map_m6m$value
  # 
  # data_df_map_pm$all  <- data_df_map_pm$all  * data_df_map_pop$urb_inc_density_fac
  # data_df_map_pm$anth <- data_df_map_pm$anth * data_df_map_pop$urb_inc_density_fac
  # 
  # data_df_map_m6m$value <- data_df_map_m6m$value * data_df_map_pop$urb_inc_density_fac
  
  data_df_map_pm$all_no_urb  <- data_df_map_pm$all  
  data_df_map_pm$anth_no_urb <- data_df_map_pm$anth
  
  
  data_df_map_pm$anth <- data_df_map_pm$anth * data_df_map_pop$urb_inc_density_fac
  
  # the urban increment is only applied to the anthropogenic PM, check later  
  
  if(u_urb_inc_mort == TRUE) {  data_df_map_pm$all  <- data_df_map_pm$all_no_urb  -  data_df_map_pm$anth_no_urb + data_df_map_pm$anth }
  if(u_urb_inc_mort == FALSE){  data_df_map_pm$all  <- data_df_map_pm$all_no_urb }
  
  data_df_map_m6m$value <- data_df_map_m6m$value 
  
  if(u_WHO_harmonize){
    # harmonize the PM concentrations here to historic values (WHO 2014)
    
    #harmonize to the earliest year
    if(year == min(u_year)){
      
      #read the WHO and calculate the factor on a 0.125 grid level
      #save this factor as a harmonization factor for all years
      # distribution_factor= who_concentrations_regrid/TM5_results 
      load(file=file.path(dir_harmonization, paste0("who_concentrations_regrid.RData")))
      
      who_conc <- left_join(who_conc,  data_df_map_pm, by=c("lon_id" = "lon_id", "lat_id" = "lat_id")) %>%
        dplyr::mutate(conc_fac = PM25/all) %>%
        dplyr::select(lon_id, lat_id, conc_fac)
      
      
      
    }
    
    #  left join factor and data_df_map_pm$all and multiply pm with harmonization factor
    # we only harmonize all
    data_df_map_pm <- left_join(data_df_map_pm, who_conc, by=c("lon_id" = "lon_id", "lat_id" = "lat_id")) %>%
      dplyr::mutate(all_unharm = all) %>%
      dplyr::mutate(all = all * conc_fac)
    
  }
  
  
  # Compute PM mortality maps (using Burnett's functions)
  if (u_verbose) cat("    - Computing mortality maps for PM and O3...\n")
  t02 <- Sys.time()
  # TODO: Parallelize
  # Compute for risk rate (rr) and attributable fraction (af) for each disease type
  # Calculation of attributable fraction which is later used
  for (ktyp in c("copd", "alri", "lc", "ihd", "stroke")) {
    #rr <- rrate_df(rrcoef[[paste0(toupper(ktyp), "_", toupper(u_mort_lvl))]],   data_df_map_pm)
    rr <- rrate_df_2018(rrcoef[[paste0(toupper(ktyp), "_", toupper(u_mort_lvl))]],   data_df_map_pm)
    af[[paste0("ap_",ktyp,"_",u_mort_lvl)]] <- rr %>%
      dplyr::mutate(value = (value-1)/value)
  }
  
  data_df_mort_pm <- cbind(
    data_df_map_pop,
    data_df_maps_mort_hires$frac_30_map[[paste(year)]] %>%
      dplyr::select(value) %>% dplyr::rename(frac_30 = value),
    data_df_maps_mort_hires$copd_map[[paste(year)]] %>%
      dplyr::select(value) %>% dplyr::rename(copd = value),
    data_df_maps_mort_hires$ihd_map[[paste(year)]] %>%
      dplyr::select(value) %>% dplyr::rename(ihd = value),
    data_df_maps_mort_hires$lc2_map[[paste(year)]] %>%
      dplyr::select(value) %>% dplyr::rename(lc = value),
    data_df_maps_mort_hires$stroke_map[[paste(year)]] %>%
      dplyr::select(value) %>% dplyr::rename(stroke = value),
    data_df_maps_mort_hires$frac_05_map[[paste(year)]] %>%
      dplyr::select(value) %>% dplyr::rename(frac_05 = value),
    data_df_maps_mort_hires$alri_map[[paste(year)]] %>%
      dplyr::select(value) %>% dplyr::rename(alri = value),
    af[[paste0("ap_copd_",u_mort_lvl)]] %>%
      dplyr::select(value) %>% dplyr::rename(af_copd = value),
    af[[paste0("ap_ihd_",u_mort_lvl)]] %>%
      dplyr::select(value) %>% dplyr::rename(af_ihd = value),
    af[[paste0("ap_lc_",u_mort_lvl)]] %>%
      dplyr::select(value) %>% dplyr::rename(af_lc = value),
    af[[paste0("ap_stroke_",u_mort_lvl)]] %>%
      dplyr::select(value) %>% dplyr::rename(af_stroke = value),
    af[[paste0("ap_alri_",u_mort_lvl)]] %>%
      dplyr::select(value) %>% dplyr::rename(af_alri = value)) %>%
    dplyr::mutate(copd1=pop*(frac_30/100*(copd/100*af_copd))) %>%
    dplyr::mutate(ihd1=pop*(frac_30/100*(ihd/100*af_ihd))) %>%
    dplyr::mutate(lc1=pop*(frac_30/100*(lc/100*af_lc))) %>%
    dplyr::mutate(stroke1=pop*(frac_30/100*(stroke/100*af_stroke))) %>%
    dplyr::mutate(alri1=pop*(frac_05/100*(alri/100*af_alri))) %>%
    #here the mortality is calculated
    dplyr::mutate(value=pop*(
      frac_30/100*(copd/100*af_copd + ihd/100*af_ihd + lc/100*af_lc + stroke/100*af_stroke) +
        frac_05/100*(alri/100*af_alri))) %>%
    #dplyr:select(-copd,-ihd,-lc,-stroke,-alri,-af_copd,-af_ihd,-af_lc,-af_stroke,-af_alri,-pop,-frac_30, -frac_05,-lon_id_coarse8,-lat_id_coarse8)
    dplyr::select(-copd,ihd,-lc,-stroke,-alri,-af_copd,-af_ihd,-af_lc,-af_stroke,-af_alri,-pop,-lon_id_coarse8,-lat_id_coarse8)
  
  
  
  # In GBD visualization tool O3 contributes to COPD cause of death.
  # Using COPD baseline mortalities gives much better agreement with GBD.
  # using Anenberg old function
  tmp_beta <- p_beta[[grep(u_mort_lvl, names(p_beta), value=TRUE)]]
  data_df_mort_o3 <- cbind(
    data_df_map_pop,
    data_df_maps_mort_hires$frac_30_map[[paste(year)]] %>%
      dplyr::select(value) %>% dplyr::rename(frac_30 = value),
    data_df_maps_mort_hires$copd_map[[paste(year)]] %>%
      dplyr::select(value) %>% dplyr::rename(copd = value),
    data_df_map_m6m %>%
      dplyr::select(value) %>% dplyr::rename(m6m = value)) %>%
    dplyr::mutate(value = ifelse(m6m > m6mthr, copd/100 *  frac_30/100 * pop * (1-exp(-tmp_beta*(m6m-m6mthr))), 0))
  
  
  # it should be possible to perform some of these joins out of the year loop to speed things up (at least the VSL and iso one)
  # calculate the cost of VSL here
  VSL_cost <- data_df_mort_pm[,c('lon_id', 'lat_id','copd1','ihd1','lc1','stroke1','alri1', 'value')]
  ## sum pm and o3 premature deathes
  VSL_cost$value <- (data_df_mort_pm$value+data_df_mort_o3$value)
  data_df_mort_pm$all_mort <-  data_df_mort_pm$value+data_df_mort_o3$value
  
  VSL_yearly <-               left_join(VSL,VSL, by=c("lon_id_coarse2"="lon_id", "lat_id_coarse2"="lat_id")) %>%
    dplyr::select(lon_id,lat_id, paste0(year,".y"),iso.y)
  
  iso                  <-     VSL_yearly[,c("lon_id","lat_id","iso.y")]
  
  data_df_mort_pm$VSL_cost <- as.numeric(unlist(left_join(VSL_cost, VSL_yearly[,c("lon_id","lat_id",paste0(year,".y"))], by=c("lon_id"="lon_id", "lat_id"="lat_id")) %>%
                                                  dplyr::select(value, paste0(year,".y"))    %>%
                                                  dplyr::mutate(value = Reduce(`*`, .)) %>%   # multiply deathes with VSL
                                                  dplyr::select(value)))
  
  
  
  tend2 <- Sys.time() - t02
  if(u_verbose) cat(paste0("[Time required to compute mortality maps: ", tend2, "s]\n"))
  
  # Aggregate to FASST regions (discard first region which corresponds to sea and oceans)
  # (Runtime: ~ 25 min) TODO: Parallelize
  # Setup parallel backend
  cl <- makeCluster(1)
  registerDoParallel(cl)
  if (doMortalitiesByTM5region) {
    if (u_verbose) cat("    - Looping over countries/regions and computing mortalities...\n")
    t02 <- Sys.time()
    #for (icntr in 1:56) {
    foreach(icntr=1:56, .packages=c("dplyr","tidyr")) %dopar% {
      if (u_verbose) cat(paste0("      . ", paste(icntr), ": ", mapping_FASST_regions$name[mapping_FASST_regions$id == icntr], "\n"))
      
      # Compute total number of mortalities for PM and O3
      ctot_mort_pm_med        <- as.numeric(left_join(data_df_maps_countryMask_hires %>%
                                                        dplyr::rename(reg_id = value) %>%
                                                        filter(reg_id == icntr),
                                                      data_df_mort_pm %>%
                                                        dplyr::rename(mort = value),
                                                      by=c("lon_id","lat_id")) %>%
                                              dplyr::mutate(mort=ifelse(is.na(mort),0,mort)) %>%
                                              summarise(sum=sum(mort)))
      ctot_mort_o3_2005_med   <- as.numeric(left_join(data_df_maps_countryMask_hires %>%
                                                        dplyr::rename(reg_id = value) %>%
                                                        filter(reg_id == icntr),
                                                      data_df_mort_o3 %>%
                                                        dplyr::rename(mort = value),
                                                      by=c("lon_id","lat_id")) %>%
                                              dplyr::mutate(mort=ifelse(is.na(mort),0,mort)) %>%
                                              summarise(sum=sum(mort)))
      
      
      
    } #TM5 regions
    stopCluster(cl)
  }
  
  # add a ISO country aggregation
  
  # add a REMIND region aggregation
  
  tend2 <- Sys.time() - t02
  if(u_verbose) cat(paste0("[Time required to loop over countries: ", tend2, "s]\n"))
  
  
  if(netcdf_write_out_yearly | netcdf_write_all_years){
    #for the netcdf write out we need to add the full dimensions, therefore they are joined with the full_dim
    data_df_map_pop <- left_join(dim_df_maps_mort_mapping_dim_id, data_df_map_pop, by=c("lon_id"="lon_id", "lat_id"="lat_id"))
    data_df_mort_pm <- left_join(dim_df_maps_mort_mapping_dim_id, data_df_mort_pm, by=c("lon_id"="lon_id", "lat_id"="lat_id"))
    data_df_map_pm  <- left_join(dim_df_maps_mort_mapping_dim_id, data_df_map_pm , by=c("lon_id"="lon_id", "lat_id"="lat_id"))
    data_df_mort_o3 <- left_join(dim_df_maps_mort_mapping_dim_id, data_df_mort_o3, by=c("lon_id"="lon_id", "lat_id"="lat_id"))
    data_df_map_m6m <- left_join(dim_df_maps_mort_mapping_dim_id, data_df_map_m6m, by=c("lon_id"="lon_id", "lat_id"="lat_id"))
  }
  
  if(netcdf_write_out_yearly){
    #added for checking the results, deactivate for cluster
    #-----------------------------------------------------------------------------
    # Create ncdf files of scenario PM2.5, O3 and mortalities
    #-----------------------------------------------------------------------------
    if (u_verbose) cat("    - Writing out NetCDF file...\n")
    t02 <- Sys.time()
    filenc <- file.path(dir_out, paste0(RUNSUFFIX, year  ,'.nc'))
    
    
    
    dim_lon <- ncdim_def( "lon", "degrees_east",  lon, longname='longitude gridbox bottom left corner')
    dim_lat <- ncdim_def( "lat", "degrees_north", lat, longname='latitude gridbox bottom left corner')
    
    # Define variables
    vars = list()
    vars[["iso"]]             <- ncvar_def('iso',     '',              list(dim_lon,dim_lat), longname='iso country')
    vars[["POP"]]             <- ncvar_def('POP',     '#/0.125deg cell', list(dim_lon,dim_lat), longname='# population')
    vars[["PM_MORT"]]         <- ncvar_def('PM_MORT', '#/0.125deg cell', list(dim_lon,dim_lat), longname='# of annual premature deaths (Burnett et al. 2013)')
    vars[["ALL_MORT"]]        <- ncvar_def('ALL_MORT','#/0.125deg cell', list(dim_lon,dim_lat), longname='# of annual premature deaths')
    vars[["VSL_COST"]]        <- ncvar_def('VSL_COST','$ 2005/0.125deg cell',        list(dim_lon,dim_lat), longname='annual VSL cost')
    vars[["PM"]]              <- ncvar_def('PM',      'ug/m3',         list(dim_lon,dim_lat), longname='scenario total dry PM2.5 incl dust and SS')
    #vars[["PM_UNHARM"]]       <- ncvar_def('PM_UNHARM','ug/m3',        list(dim_lon,dim_lat), longname='scenario total dry PM2.5 incl dust and SS unharmonized')
    vars[["ANTH_PM"]]         <- ncvar_def('ANTH_PM', 'ug/m3',         list(dim_lon,dim_lat), longname='scenario anthropogenic dry PM2.5')
    vars[["PM_NO_URB"]]       <- ncvar_def('PM_NO_URB','ug/m3',        list(dim_lon,dim_lat), longname='scenario total dry PM2.5 incl dust and SS, no_urb')
    vars[["ANTH_PM_NO_URB"]]  <- ncvar_def('ANTH_PM_NO_URB', 'ug/m3',  list(dim_lon,dim_lat), longname='scenario anthropogenic dry PM2.5, no_urb')
    vars[["O3_MORT"]]         <- ncvar_def('O3_MORT', '#/0.125deg cell', list(dim_lon,dim_lat), longname='# of annual premature deaths from short+longterm O3 exposure (Anenberg et al, 2010, Jerrett et. al 2009)')
    vars[["M6M"]]             <- ncvar_def('M6M',     'ppbv',          list(dim_lon,dim_lat), longname='maximal 6 months mean ozone')
    vars[["URB_INC"]]         <- ncvar_def('URB_INC', 'fac',           list(dim_lon,dim_lat), longname='urban increment factor')
    vars[["POP_URB"]]         <- ncvar_def('POP_URB', '#/0.125deg cell', list(dim_lon,dim_lat), longname='# urban population')
    vars[["POP_RUR"]]         <- ncvar_def('POP_RUR', '#/0.125deg cell', list(dim_lon,dim_lat), longname='# rural population')
    
    
    
    # Create NetCDF file
    fid <- nc_create(filenc, vars=vars)
    
    # Write out global attributes
    ncatt_put(fid, 0, 'source',  'health_impacts_from_fasst_gridmaps.R (translated from IMPACTS_FROM_FASST_GRIDDED_CONC.PRO, orginally written by Rita van Dingenen rita.van-dingenen@jrc.ec.europa.eu)')
    ncatt_put(fid, 0, 'contact', 'hilaire@pik-potsdam.de, hilaire@mcc-berlin.net')
    
    # Put values in NetCDF file
    # ncvar_put(fid, vars[["lon"]],     lon)
    # ncvar_put(fid, vars[["lat"]],     lat)
    ncvar_put(fid, vars[["iso"]],        as.factor(iso$iso.y))
    ncvar_put(fid, vars[["POP"]],        data_df_map_pop$pop)
    ncvar_put(fid, vars[["PM_MORT"]],    data_df_mort_pm$value)
    ncvar_put(fid, vars[["ALL_MORT"]],   data_df_mort_pm$all_mort)
    ncvar_put(fid, vars[["VSL_COST"]],   data_df_mort_pm$VSL_cost)
    ncvar_put(fid, vars[["PM"]],         data_df_map_pm$all)
    #ncvar_put(fid, vars[["PM_UNHARM"]],  data_df_map_pm$all_unharm)
    ncvar_put(fid, vars[["PM_NO_URB"]],  data_df_map_pm$all_no_urb)
    ncvar_put(fid, vars[["ANTH_PM"]],    data_df_map_pm$anth)
    ncvar_put(fid, vars[["ANTH_PM_NO_URB"]], data_df_map_pm$anth_no_urb)
    ncvar_put(fid, vars[["O3_MORT"]],    data_df_mort_o3$value)
    ncvar_put(fid, vars[["M6M"]],        data_df_map_m6m$value)
    ncvar_put(fid, vars[["URB_INC"]],    data_df_map_pop$urb_inc_density_fac)
    ncvar_put(fid, vars[["POP_URB"]],    data_df_map_pop$pop_urb)
    ncvar_put(fid, vars[["POP_RUR"]],    data_df_map_pop$pop_rur)
    
    
    # Close NetCDF file
    nc_close(fid)
    
    if(u_verbose) cat(paste('    - NetCDF file:',filenc,"written.\n"))
  } 
  cat(paste('mean PM all concentration'  ,round(mean(na.omit(data_df_map_pm$all)))     , 'ug/m3 in ',year,'\n'))
  cat(paste('mean PM anth concentration' ,round(mean(na.omit(data_df_map_pm$anth)))    , 'ug/m3 in ',year,'\n'))
  cat(paste('mean M6m concentration'     ,round(mean(na.omit(data_df_map_m6m$value)))  , 'ppbv in ',year,'\n'))
  cat(paste('PM premature deathes'       ,round(sum(na.omit(data_df_mort_pm$value)))   , 'in ',year,'\n'))
  cat(paste('O3 premature deathes'       ,round(sum(na.omit(data_df_mort_o3$value)))   , 'in ',year,'\n'))
  cat(paste('VSL cost'                   ,round(sum(na.omit(data_df_mort_pm$VSL_cost))), 'in $ (2005)',year,'\n'))
  cat(paste('mean PM concentration'      ,round(mean(na.omit(data_df_map_pm$all)))     , 'ug/m3 in ',year,'\n'),file=txtfilepath, append=TRUE)
  cat(paste('mean PM anth concentration' ,round(mean(na.omit(data_df_map_pm$anth)))    , 'ug/m3 in ',year,'\n'),file=txtfilepath, append=TRUE)
  cat(paste('mean M6m concentration'     ,round(mean(na.omit(data_df_map_m6m$value)))  , 'ppbv in ',year,'\n') ,file=txtfilepath, append=TRUE)
  cat(paste('PM premature deathes'       ,round(sum(na.omit(data_df_mort_pm$value)))   , 'in ',year,'\n')      ,file=txtfilepath, append=TRUE)
  cat(paste('O3 premature deathes'       ,round(sum(na.omit(data_df_mort_o3$value)))   , 'in ',year,'\n')      ,file=txtfilepath, append=TRUE)
  cat(paste('VSL cost'                   ,round(sum(na.omit(data_df_mort_pm$VSL_cost))), 'in $ (2005)',year,'\n')     ,file=txtfilepath, append=TRUE)
  cat(paste(u_WHO_harmonize                                                            , 'WHO_harmonized \n')     ,file=txtfilepath, append=TRUE)
  cat(paste(u_urb_inc_mort                                                             , 'is urban increment considered in mortality calculation \n')     ,file=txtfilepath, append=TRUE)
  
  
  
  results_array[paste(year),,] <- as.matrix(cbind( dim_df_maps_mort_mapping_dim_id, iso$iso.y, data_df_map_pop$pop,
                                                   data_df_mort_pm$frac_05,data_df_mort_pm$frac_30,
                                                   data_df_mort_pm$copd1,data_df_mort_pm$ihd1,
                                                   data_df_mort_pm$lc1,data_df_mort_pm$stroke1,
                                                   data_df_mort_pm$alri1,
                                                   data_df_mort_pm$value,data_df_mort_pm$all_mort,
                                                   data_df_mort_pm$VSL_cost, data_df_map_pm$all,
                                                   data_df_map_pm$anth, data_df_map_pm$anth_no_urb,
                                                   data_df_map_pm$all_no_urb, data_df_mort_o3$value,
                                                   data_df_map_m6m$value, data_df_map_pop$urb_inc_density_fac,
                                                   data_df_map_pop$pop_urb, data_df_map_pop$pop_rur))
  
  # write out results aggregated to the iso level
  ## preliminary - add a real write out in one table and all years
  # write.csv(aggregate(data_df_mort_pm$value  ~ iso$iso.y,
  #                     results_array[paste(year),,c('iso$iso.y', 'data_df_mort_pm$value')], sum),
  #           file.path(dir_outtab,paste0('/pm_MORT',RUNSUFFIX, year, '.csv')))
  # write.csv(aggregate(data_df_mort_o3$value  ~ iso$iso.y,
  #                     results_array[paste(year),,c('iso$iso.y', 'data_df_mort_o3$value')], sum),
  #           file.path(dir_outtab,paste0('/o3_MORT',RUNSUFFIX, year, '.csv')))
  # write.csv(aggregate(data_df_mort_pm$VSL_cost  ~ iso$iso.y,
  #                     results_array[paste(year),,c('iso$iso.y',  'data_df_mort_pm$VSL_cost')], sum),
  #           file.path(dir_outtab,paste0('/VSL',RUNSUFFIX, year, '.csv')))
  
  # write_out[[i]] <- left_join(aggregate(data_df_mort_pm$value  ~ iso$iso.y,       results_array[paste(year),,c('iso$iso.y', 'data_df_mort_pm$value')],    sum),
  #                                aggregate(data_df_mort_o3$value  ~ iso$iso.y,    results_array[paste(year),,c('iso$iso.y', 'data_df_mort_o3$value')],    sum),
  #                                aggregate(data_df_mort_pm$VSL_cost  ~ iso$iso.y, results_array[paste(year),,c('iso$iso.y', 'data_df_mort_pm$VSL_cost')], sum),
  #                                by=c("iso$iso.y"="iso$iso.y"))
  
  
  
  write_out[[i]] <- plyr::join_all(list(aggregate(data_df_mort_pm$value  ~ iso$iso.y,       results_array[paste(year),,c('iso$iso.y', 'data_df_mort_pm$value')],    sum),
                                   aggregate(data_df_mort_o3$value  ~ iso$iso.y,    results_array[paste(year),,c('iso$iso.y', 'data_df_mort_o3$value')],    sum),
                                   aggregate(data_df_mort_pm$VSL_cost  ~ iso$iso.y, results_array[paste(year),,c('iso$iso.y', 'data_df_mort_pm$VSL_cost')], sum),
                                   aggregate(data_df_mort_pm$frac_05  ~ iso$iso.y, results_array[paste(year),,c('iso$iso.y', 'data_df_mort_pm$frac_05')], mean),
                                   aggregate(data_df_mort_pm$frac_30  ~ iso$iso.y, results_array[paste(year),,c('iso$iso.y', 'data_df_mort_pm$frac_30')], mean),
                                   aggregate(data_df_mort_pm$copd1  ~ iso$iso.y, results_array[paste(year),,c('iso$iso.y', 'data_df_mort_pm$copd1')], sum),
                                   aggregate(data_df_mort_pm$ihd1  ~ iso$iso.y, results_array[paste(year),,c('iso$iso.y', 'data_df_mort_pm$ihd1')], sum),
                                   aggregate(data_df_mort_pm$lc1  ~ iso$iso.y, results_array[paste(year),,c('iso$iso.y', 'data_df_mort_pm$lc1')], sum),
                                   aggregate(data_df_mort_pm$stroke1  ~ iso$iso.y, results_array[paste(year),,c('iso$iso.y', 'data_df_mort_pm$stroke1')], sum),
                                   aggregate(data_df_mort_pm$alri1  ~ iso$iso.y, results_array[paste(year),,c('iso$iso.y', 'data_df_mort_pm$alri1')], sum),
                                   aggregate(data_df_mort_pm$all_mort  ~ iso$iso.y, results_array[paste(year),,c('iso$iso.y', 'data_df_mort_pm$all_mort')], sum)),
                                   by='iso$iso.y', type='left')
  write_out[[i]]$'year' <- year
  write_out[[i]]$'scenario' <- scenario
  
   write.csv(write_out[[i]], file.path(dir_outtab,paste0('results_',year,RUNSUFFIX, '.csv')))
  
  i<-i+1
} #year loop

names(write_out) <- u_year

#write out the scenario specific results csv aggregated to iso level
write_out_data_frame <- as.data.frame(do.call(rbind, write_out))

write_out_data_frame <- cbind(year = sapply( strsplit( rownames( write_out_data_frame), '\\.'), head, 1), write_out_data_frame)
rownames(write_out_data_frame) <- NULL

write.csv(write_out_data_frame, file.path(dir_outtab,paste0('results_',RUNSUFFIX, '.csv')))


if(netcdf_write_all_years){
  
  t02 <- Sys.time()
  filenc <- file.path(dir_out, paste0(RUNSUFFIX, '_all_years.nc'))
  
  dim_year<- ncdim_def( "year", "a",             u_year, longname='years')
  dim_lon <- ncdim_def( "lon" , "degrees_east",  lon   , longname='longitude gridbox bottom left corner')
  dim_lat <- ncdim_def( "lat" , "degrees_north", lat   , longname='latitude gridbox bottom left corner')
  
  # Define variables
  vars[["iso"]]             <- ncvar_def('iso',     '',              list(dim_year,dim_lon,dim_lat), longname='iso country')
  vars[["POP"]]             <- ncvar_def('POP',     '#/0.125deg cell', list(dim_year,dim_lon,dim_lat), longname='# population')
  vars[["PM_MORT"]]         <- ncvar_def('PM_MORT', '#/0.125deg cell', list(dim_year,dim_lon,dim_lat), longname='# of annual premature deaths (Burnett et al. 2013)')
  vars[["ALL_MORT"]]        <- ncvar_def('ALL_MORT', '#/0.125deg cell', list(dim_year,dim_lon,dim_lat), longname='# of annual premature deaths')
  vars[["VSL_COST"]]        <- ncvar_def('VSL_COST','$ 2005/0.125deg cell', list(dim_year,dim_lon,dim_lat), longname='annual VSL cost')
  vars[["PM"]]              <- ncvar_def('PM',      'ug/m3',         list(dim_year,dim_lon,dim_lat), longname='scenario total dry PM2.5 incl dust and SS')
  #vars[["PM_UNHARM"]]       <- ncvar_def('PM_UNHARM','ug/m3',        list(dim_year,dim_lon,dim_lat), longname='scenario total dry PM2.5 incl dust and SS, unharmonized')
  vars[["ANTH_PM"]]         <- ncvar_def('ANTH_PM', 'ug/m3',         list(dim_year,dim_lon,dim_lat), longname='scenario anthropogenic dry PM2.5')
  vars[["PM_NO_URB"]]       <- ncvar_def('PM_NO_URB','ug/m3',        list(dim_year,dim_lon,dim_lat), longname='scenario total dry PM2.5 incl dust and SS, no_urb')
  vars[["ANTH_PM_NO_URB"]]  <- ncvar_def('ANTH_PM_NO_URB', 'ug/m3',  list(dim_year,dim_lon,dim_lat), longname='scenario anthropogenic dry PM2.5, no_urb')
  vars[["O3_MORT"]]         <- ncvar_def('O3_MORT', '#',             list(dim_year,dim_lon,dim_lat), longname='# of annual premature deaths from short+longterm O3 exposure (Anenberg et al, 2010, Jerrett et. al 2009)')
  vars[["M6M"]]             <- ncvar_def('M6M',     'ppbv',          list(dim_year,dim_lon,dim_lat), longname='maximal 6 months mean ozone')
  vars[["URB_INC"]]         <- ncvar_def('URB_INC', 'fac',           list(dim_year,dim_lon,dim_lat), longname='urban increment factor')
  vars[["POP_URB"]]         <- ncvar_def('POP_URB', '#/0.125deg cell', list(dim_year,dim_lon,dim_lat), longname='# urban population')
  vars[["POP_RUR"]]         <- ncvar_def('POP_RUR', '#/0.125deg cell', list(dim_year,dim_lon,dim_lat), longname='# rural population')
  #
  
  
  # Create NetCDF file
  fid <- nc_create(filenc, vars=vars)
  
  # Write out global attributes
  ncatt_put(fid, 0, 'source',  'health_impacts_from_fasst_gridmaps.R (translated from IMPACTS_FROM_FASST_GRIDDED_CONC.PRO, orginally written by Rita van Dingenen rita.van-dingenen@jrc.ec.europa.eu)')
  ncatt_put(fid, 0, 'contact', 'hilaire@pik-potsdam.de, hilaire@mcc-berlin.net')
  
  # Put values in NetCDF file
  # ncvar_put(fid, vars[["lon"]],     lon)
  # ncvar_put(fid, vars[["lat"]],     lat)
  ncvar_put(fid, vars[["iso"]],            as.factor(results_array[,,'iso$iso.y']))
  ncvar_put(fid, vars[["POP"]],            results_array[,,'data_df_map_pop$pop'])
  ncvar_put(fid, vars[["PM_MORT"]],        results_array[,,'data_df_mort_pm$value'])
  ncvar_put(fid, vars[["ALL_MORT"]],        results_array[,,'data_df_mort_pm$all_mort'])
  ncvar_put(fid, vars[["VSL_COST"]],       results_array[,,'data_df_mort_pm$VSL_cost'])
  ncvar_put(fid, vars[["PM"]],             results_array[,,'data_df_map_pm$all'])
  #ncvar_put(fid, vars[["PM_UNHARM"]],      results_array[,,'data_df_map_pm$all_unharm'])
  ncvar_put(fid, vars[["PM_NO_URB"]],      results_array[,,'data_df_map_pm$all_no_urb'])
  ncvar_put(fid, vars[["ANTH_PM"]],        results_array[,,'data_df_map_pm$anth'])
  ncvar_put(fid, vars[["ANTH_PM_NO_URB"]], results_array[,,'data_df_map_pm$anth_no_urb'])
  ncvar_put(fid, vars[["O3_MORT"]],        results_array[,,'data_df_mort_o3$value'])
  ncvar_put(fid, vars[["M6M"]],            results_array[,,'data_df_map_m6m$value'])
  ncvar_put(fid, vars[["URB_INC"]],        results_array[,,'data_df_map_pop$urb_inc_density_fac'])
  ncvar_put(fid, vars[["POP_URB"]],        results_array[,,'data_df_map_pop$pop_urb'])
  ncvar_put(fid, vars[["POP_RUR"]],        results_array[,,'data_df_map_pop$pop_rur'])
  
  
  # Close NetCDF file
  nc_close(fid)  
  
  
  
  
  
  
  
}



tend <- Sys.time() - t0
if(u_verbose) cat(paste0("[Time required to process data: ", tend,']\n'))
if(u_verbose) cat("Script completed successfully.\n\n")

tendend <- Sys.time() - t00
cat(paste0("[Time required to run the whole script: ", tendend,']\n'))


