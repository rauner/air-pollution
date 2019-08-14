#===============================================================================
# FASST-R
#===============================================================================
# original IDL script written by Rita van Dingenen (JRC-ISPRA)
# rewritten in R by Jerome Hilaire (MCC & PIK)
#
# Description:
# Process project PM and O3 fields to calculate mortalities
# Needs at least 1 scenario and associated year for correct pop and base
# mortality stats. Stores mortalities per region for each scenario.
# Delta's are not calculated.
# include rrate function to calculate RR from PM2.5 using Burnett's functions
#
# WARNING:
# The grid level is provided for easy exposure calculations of
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
#===============================================================================
t00 <- Sys.time()
#===============================================================================
# USER-SECTION
#===============================================================================
# Specify scenario list
u_year     <- c('2010')    # needed for population stats.
u_scen     <- 'TEST_FASST'
u_mort_lvl <- 'med'    #choose between lo, med and hi

u_nbProc   <- 8  # Max  

u_zcf     <- 2   # default=2 uses Burnett built-in zcf's
u_urbincr <- 1   # default=1 (Urban increment not included), Set
# makes sense when the gridmaps have been created with urban
# increment on)

# Input data options
u_RData         <- TRUE   # Generate RData file to speed things up
u_recreateRData <- FALSE  # Re-generate RData file (other)

doMortalitiesByCountry <- TRUE
doMortalityMapsNCDF <- FALSE

# Debug
u_verbose       <- TRUE

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

# Load own functions
dump <- lapply(list.files("functions", pattern = "*.R", full.names = TRUE), source)

# Output NetCDF file parameters
p_projectName <- 'REMIND-FASST'
p_progName    <- 'health_impacts_from_fasst_gridmaps.R'

# Folders
dir_root   <- '/p/tmp/hilaire/REMOD/FASST-R/'              # main folder
#dir_root   <- '/mnt/FASST_IDL_PACK/'                       # main folder

dir_out    <- paste0(dir_root,  'OUTPUT/')                 # output
dir_outin  <- paste0(dir_root,  'OUTPUT/NCDF')             # output (ncdf files)
dir_outtab <- paste0(dir_root,  'OUTPUT/TABLES')           # output (tables)
dir_out    <- dir_outin

dir_ancil  <- paste0(dir_root,  'ANCILLARY/')              # ancillary data
dir_who <- file.path(dir_ancil, 'MORTALITY')               # mortality data

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
zcf      <- 5.8  # when flag=2, ZCF is threshold for Anenberg method only.
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
if(u_urbincr != 2) u_urbincr <- 1
UFL <- c('_URBINCR_OFF','_URBINCR_ON')
if (u_verbose) cat(paste('URBAN INCREMENT FLAG :', UFL[u_urbincr], '\n'))

# TODO: Link to our SSP population data?
GEA_POP_FLAG <- 1   # use UN gridded poulation and rural/urban data developed
# by IIASA for Global Energy Assessment (GEA)
GEAFL <- c('_POP-GEA')

# Initialise lists
p_rr            <- list()
p_beta          <- list()
af              <- list()

# Create high resolution version at same grid size as population maps
# map onto full lat range and regrid to 7.5'x7.5'
maps_dummy   <- array(0.0, c(1440,720))
maps_dummylo <- array(0.0, c(1440,720))
maps_dummyup <- array(0.0, c(1440,720))
latful  <- (0:719-360)/4.
lonful  <- (0:1439-720)/4.

# note: available population totals from GEA: 2000, 2005, 2010, 2020, ..., 2100
# available base mortalities and fraction of pop <5yr and >30yr: 2005 2010 2015
# 2030 2050.
# If other years are needed, an additional interpolation is carried out on the
# latter. For years > 2030 mortality stats for 2030 are used. For years < 2005,
# mortality stats for 2005 are used.

# Initialise processor array for parallelisation
cl <- makeCluster(u_nbProc)
registerDoParallel(cl)

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
rrcoef <- read.csv3(file.path(dir_who,'ER_FUNCTION_FITTINGS_ALL_ENDPOINTS.csv'), RData=u_RData, recreate_rdata=u_recreateRData, verbose=u_verbose)

# Set the zcf value to zcf if non-threshold calculation is selected
if (u_zcf == 0) {
  #subsitute zcf in the structure containing the retrieved Burnett coefficients
  rrcoef[which(rrcoef$PARAM == "zcf"), -1] <- zcf
}

# Read in mortality maps
maps_mort <- load(file.path(dir_ancil, "MORTALITY", "MortalityMaps_hires_df.RData"))

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
load(file.path(dir_ancil, "POPULATION_GEA", "total_pop.7.5min", "Population_allYears_df.RData"))

tend <- Sys.time() - t0
if(u_verbose) cat(paste0("[Time required to read in data: ", tend, "s\n"))

#===============================================================================
# PROCESS AND SAVE DATA (text and NetCDF files)
#===============================================================================
if(u_verbose) cat("Processing data...\n")
t0 <- Sys.time()
# Initialise output text file containing summary information
#TODO: Add formating
txtfilepath <- file.path(dir_outtab, paste0('FASST_', p_projectName, UFL[u_urbincr], zfl[u_zcf], GEAFL[GEA_POP_FLAG], '.txt'))

cat('Change with V4: use COPD for O3\n',                           file=txtfilepath)
cat(paste0('programme: ', p_progName, '\n'),                       file=txtfilepath, append=TRUE)
cat(paste0('date of run: ', Sys.time(), '\n'),                     file=txtfilepath, append=TRUE)
cat(paste0(p_projectName, UFL[u_urbincr], zfl[u_zcf], '\n'),       file=txtfilepath, append=TRUE)
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

# Loop over time steps (19) and compute mortality rates
# TODO: Parallelize?
if (u_verbose) cat("Loop over scenarios to compute mortality rates...\n")
for (isc in 1:length(u_scen)) {

  if (u_verbose) cat(paste0("  > Scenario: ", paste(isc), "\n"))
  if (u_verbose) cat(paste0("    - Processing population and mortality data...\n"))
  t02 <- Sys.time()
  # Read scenario population file
  data_df_map_pop <- data_df_maps_pop[[paste0(u_year[isc])]]

  cat(paste(paste(u_year[isc]), ":",  paste(sum(data_df_map_pop$pop*1e-6), "millions\n")))

  # Select correct year for mortality and pop stats
  iys  = which(yearcatalog == u_year[isc])
  #If current year is not in standard set, then interpolate or use 2030 for years > 2030
  flg_interp <- FALSE
  if (length(iys) == 0) {
    flg_interp <- TRUE
    ilow       <- max(which(yearcatalog < u_year[isc]))
    iup        <- min(which(yearcatalog > u_year[isc]))
    CL         <- length(ilow)
    CH         <- length(iup)
    if (length(ilow)) {
      iys        <- which(yearcatalog == min(yearcatalog))
      flg_interp <- FALSE
    }
    if (length(iup)) {
      iys        <- which(yearcatalog == max(yearcatalog))
      flg_interp <- FALSE
    }
  }

  # Is time interpolation required?
  if (flg_interp) {
    if (u_verbose) cat("    - Interpolating...\n")

    print("TODO: generate interpolation outside of main script")
    stop()
    # If interpolation is not required...
  } else {
    if (u_verbose) cat("    - Interpolation not required\n")

    maps_mort_hires <- data_df_maps_mort_hires[[paste(yearcatalog[iys])]]
  }

  # Required by the Anenberg/Krewski function
  maps_mort_hires$cp <- cbind(
    data_df_maps_mort_hires$copd_map[[paste(yearcatalog[iys])]]                     %>% rename(copd   = value),
    data_df_maps_mort_hires$alri_map[[paste(yearcatalog[iys])]]   %>% select(value) %>% rename(alri   = value),
    data_df_maps_mort_hires$ihd_map[[paste(yearcatalog[iys])]]    %>% select(value) %>% rename(ihd    = value),
    data_df_maps_mort_hires$stroke_map[[paste(yearcatalog[iys])]] %>% select(value) %>% rename(stroke = value)) %>%
    mutate(value = copd + alri + ihd + stroke) %>%
    select(-copd,-alri,-ihd,-stroke)

  tend2 <- Sys.time() - t02
  if(u_verbose) cat(paste0("[Time required to process population and mortality data: ", tend2, "s]\n"))

  # Read the scenario gridfile with air pollutant concentrations (PM & map_apO3)
  # (generated with fasst_gridmaps_pm_gas.R)
  # TODO: Save AP concentration data as data frame
  if (u_verbose) cat("    - Reading in air pollutants concentrations...\n")
  t02 <- Sys.time()
  maps_ap <- read.ncdf(file.path(dir_outin, paste0(u_scen[isc], ".nc")), RData=u_RData, recreate_rdata=u_recreateRData, verbose=u_verbose)

  # Convert to data frame (move to previous R script)
  data_df_maps_ap <- list()
  for (kv in names(maps_ap$var)) {
    data_df_maps_ap[[kv]] <- as.data.frame(maps_ap$var[[kv]])

   names(data_df_maps_ap[[kv]]) <- paste(seq(0,180-1))

    data_df_maps_ap[[kv]] <- data_df_maps_ap[[kv]] %>%
      cbind(data.frame(lon_id=seq(0,360-1))) %>%
      gather(lat_id,value,-lon_id) %>%
      mutate(lat_id=as.integer(lat_id))
  }

  # Sum up aerosols
  # Add Dust (DU) and Sea-Salt (SS) because new RR relations are non-linear
  data_df_map_pm <- cbind(
    data_df_maps_ap$so4  %>% rename(so4  = value),
    data_df_maps_ap$no3  %>% rename(no3  = value) %>% select(no3),
    data_df_maps_ap$nh4  %>% rename(nh4  = value) %>% select(nh4),
    data_df_maps_ap$pom  %>% rename(pom  = value) %>% select(pom),
    data_df_maps_ap$pm25 %>% rename(pm25 = value) %>% select(pm25),
    data_df_maps_ap$du   %>% rename(du   = value) %>% select(du),
    data_df_maps_ap$ss   %>% rename(ss   = value) %>% select(ss)) %>%
    mutate(all =so4+no3+nh4+pom+pm25+du+ss + H2O_FL*(0.27*(so4+no3+nh4)+0.15*ss)) %>%
    mutate(anth=all-du-ss) %>%
    select(-so4,-no3,-nh4,-pom,-pm25,-du,-ss)

  data_df_map_m6m <- data_df_maps_ap$m6m

  # Downscale AP concentrations from 360x180 to 2880x1440
  data_df_map_pm <- data_df_map_pop %>%
    left_join(data_df_map_pm, by=c("lon_id_coarse8"="lon_id", "lat_id_coarse8"="lat_id"))
  data_df_map_m6m <- data_df_map_pop %>%
    left_join(data_df_map_m6m, by=c("lon_id_coarse8"="lon_id", "lat_id_coarse8"="lat_id"))
  tend2 <- Sys.time() - t02
  if(u_verbose) cat(paste0("[Time required to read and process AP concentration data: ", tend2, "s]\n"))

  # Compute PM mortality maps (using Burnett's functions)
  if (u_verbose) cat("    - Computing mortality maps for PM and O3...\n")
  t02 <- Sys.time()
  # TODO: Parallelize
  # Compute for risk rate (rr) and attributable fraction (af) for each disease type
  for (ktyp in c("copd", "alri", "lc", "ihd", "stroke")) {
  #foreach(ktyp in c("copd", "alri", "lc", "ihd", "stroke"), .packages=c("dplyr","tidyr")) %dopar%
    rr <- rrate_df(rrcoef[[paste0(toupper(ktyp), "_", toupper(u_mort_lvl))]],   data_df_map_pm)
    af[[paste0("ap_",ktyp,"_",u_mort_lvl)]] <- rr %>%
      mutate(value = (value-1)/value)
  }

  data_df_mort_pm <- cbind(
    data_df_map_pop,
    data_df_maps_mort_hires$frac_30_map[[paste(yearcatalog[iys])]] %>%
      select(value) %>% rename(frac_30 = value),
    data_df_maps_mort_hires$copd_map[[paste(yearcatalog[iys])]] %>%
      select(value) %>% rename(copd = value),
    data_df_maps_mort_hires$ihd_map[[paste(yearcatalog[iys])]] %>%
      select(value) %>% rename(ihd = value),
    data_df_maps_mort_hires$lc2_map[[paste(yearcatalog[iys])]] %>%
      select(value) %>% rename(lc = value),
    data_df_maps_mort_hires$stroke_map[[paste(yearcatalog[iys])]] %>%
      select(value) %>% rename(stroke = value),
    data_df_maps_mort_hires$frac_05_map[[paste(yearcatalog[iys])]] %>%
      select(value) %>% rename(frac_05 = value),
    data_df_maps_mort_hires$alri_map[[paste(yearcatalog[iys])]] %>%
      select(value) %>% rename(alri = value),
    af[[paste0("ap_copd_",u_mort_lvl)]] %>%
      select(value) %>% rename(af_copd = value),
    af[[paste0("ap_ihd_",u_mort_lvl)]] %>%
      select(value) %>% rename(af_ihd = value),
    af[[paste0("ap_lc_",u_mort_lvl)]] %>%
      select(value) %>% rename(af_lc = value),
    af[[paste0("ap_stroke_",u_mort_lvl)]] %>%
      select(value) %>% rename(af_stroke = value),
    af[[paste0("ap_alri_",u_mort_lvl)]] %>%
      select(value) %>% rename(af_alri = value)) %>%
    mutate(value=1000*pop*(
      frac_30*(copd*af_copd + ihd*af_ihd + lc*af_lc + stroke*af_stroke) +
      frac_05*(alri*af_alri))) %>%
    select(-copd,-ihd,-lc,-stroke,-alri,-af_copd,-af_ihd,-af_lc,-af_stroke,-af_alri,-pop,-frac_30, -frac_05,-lon_id_coarse8,-lat_id_coarse8)

  # In GBD visualization tool O3 contributes to COPD cause of death.
  # Using COPD baseline mortalities gives much better agreement with GBD.
  # using Anenberg old function
  tmp_beta <- p_beta[[grep(u_mort_lvl, names(p_beta), value=TRUE)]]
  data_df_mort_o3 <- cbind(
    data_df_map_pop,
    data_df_maps_mort_hires$frac_30_map[[paste(yearcatalog[iys])]] %>%
      select(value) %>% rename(frac_30 = value),
    data_df_maps_mort_hires$copd_map[[paste(yearcatalog[iys])]] %>%
      select(value) %>% rename(copd = value),
    data_df_map_m6m %>%
      select(value) %>% rename(m6m = value)) %>%
  mutate(value = ifelse(m6m > m6mthr, copd * 1000 * frac_30 * pop * (1-exp(-tmp_beta*(m6m-m6mthr))), 0))

  tend2 <- Sys.time() - t02
  if(u_verbose) cat(paste0("[Time required to compute mortality maps: ", tend2, "s]\n"))

  # Aggregate to FASST regions (discard first region which corresponds to sea and oceans)
  if (doMortalitiesByCountry) {
    if (u_verbose) cat("    - Looping over countries/regions and computing mortalities...\n")
    t02 <- Sys.time()
    data_df_mort_country <- foreach(icntr=1:56, .packages=c("dplyr","tidyr"), .combine=rbind, .verbose=u_verbose) %dopar% {
      #if (u_verbose) cat(paste0("      . ", paste(icntr), ": ", mapping_FASST_regions$name[mapping_FASST_regions$id == icntr], "\n"))

      # Compute total number of mortalities for PM and O3
      ctot_mort_pm_med        <- as.numeric(left_join(data_df_maps_countryMask_hires %>%
                                                        rename(reg_id = value) %>%
                                                        filter(reg_id == icntr),
                                                      data_df_mort_pm %>%
                                                        rename(mort = value),
                                                      by=c("lon_id","lat_id")) %>%
                                              mutate(mort=ifelse(is.na(mort),0,mort)) %>%
                                              summarise(sum=sum(mort)))
      ctot_mort_o3_2005_med   <- as.numeric(left_join(data_df_maps_countryMask_hires %>%
                                                        rename(reg_id = value) %>%
                                                        filter(reg_id == icntr),
                                                      data_df_mort_o3 %>%
                                                        rename(mort = value),
                                                      by=c("lon_id","lat_id")) %>%
                                              mutate(mort=ifelse(is.na(mort),0,mort)) %>%
                                              summarise(sum=sum(mort)))

      data_df_mort_country <- data.frame(
        reg_id  = icntr,
        mort_pm = ctot_mort_pm_med,
        mort_o3 = ctot_mort_o3_2005_med 
      )

    } #cntries
  }
  tend2 <- Sys.time() - t02
  if(u_verbose) cat(paste0("[Time required to loop over countries: ", tend2, "s]\n"))

} # time step
tend <- Sys.time() - t0
if(u_verbose) cat(paste0("[Time required to process data: ", tend, "s]\n"))
if(u_verbose) cat("Script completed successfully.\n\n")

stopCluster(cl)

tendend <- Sys.time() - t00
cat(paste0("[Time required to run the whole script: ", tendend, "s]\n"))
