# ------------------------------------------------------------------------------
# Program Name: quick_fasst_gridmaps_pm_gas_concentration.R
# Author(s): Sebastian Rauner
#            original IDL script written by Rita van Dingenen (JRC-ISPRA)
#            rewritten in R by Jerome Hilaire (MCC & PIK) & Sebastian Rauner (PIK)
# Date Last Updated: Oct 17, 2017
# Program Purpose: This script calculates the concentrations from the downscaled emissions
# Input Files: 
# Output Files:
# To-Do:best would be to do the emission diff not only on the country (TM5 region) level but the grid level and use the IIASA gridding files as input
#       however not sure with the pertubation data available...


# Other
NODATA <- -9.99E4 # Missing emission data should have value -9999

# Output species
species            <- c("SO4", "NO3", "NH4", "BC", "POM", "PM25", "SS", "DU", "O3", "SO2", "NOX", "M6M")
emi_diff           <- list()
emi_diff_year_list <- list()
dgrid              <- list()
for (k in species) dgrid[[k]] <- array(0.0, dim=c(360,180))


#===============================================================================
# READING-IN DATA
#===============================================================================
if (u_verbose) cat(paste0("Reading in data...\n"))
# retrieve global base concentration fields
conc_base_aer     <- nc_readall(file_abase)
conc_base_gas     <- nc_readall(file_gbase)
conc_base_dustsea <- nc_readall(file_dustsea)
conc_base_m6m     <- nc_readall(file_m6m)
conc_base_mf      <- nc_readall(file_mf)

# CH4 perturbation fields (independant of emission location, so one normalized
# footprint per kg(CH4) emission for all source regions)
pertch4 <- nc_readall(file_pertch4)

# Get base emissions (kg/year) (AR5 base and scenarios)
emi_base  <- read.table(file_embase, header=TRUE, colClasses=c("factor", rep("numeric",10)), na.strings=c("nan"))
cntrlist <- paste(emi_base$COUNTRY)
ncntr    <- length(cntrlist)

#load emissions from previous script
load(file.path(paste0(getwd(),'/emissions/', RUNSUFFIX, 'emissions.RData')))


# Column with PM25 should contain total primary PM2.5 including BC and POM.
# The latter are subtracted later on to calculate 'other' PM2.5
# base emission data do not contain PM2.5. We set it to zero
# and define PM2.5 scenario emission data to 'net' other PM2.5 (no BC, no POM)
emi_base$pm25 <- rep(0., ncntr)

# list of emitted components. Some are currently not used (long-lived greenhouse gases)
# define header

# SHIP and AIR not in RCP, use REMIND results? if so, rematch
comp            <- c('BC','CH4','CO2','CO','N2O','NH3','NOX','OM','SO2','VOC','PM25')
names(emi_base) <- c('COUNTRY', comp)

# # list for storing the concentrations
data_df_maps_ap_years <- list()

# loop over years and calculate the emission difference and concentration 


# Setup parallel backend
# watch out there could be non flagged wrong year write out for some counties, was solved through a year specific write out
#cl <- makeCluster(detectCores())
#cl <- makeCluster(1)

#registerDoParallel(cl)


#calculate the concentrations in parallel
#data_df_maps_ap_years <-  foreach(year = u_year, .packages=c('dplyr','tidyr'), .export = c('nc_open','ncvar_get')) %dopar% { 

# one core implementation  
for(year in u_year){
  
  
  # Get scenario emissions (kg/year) (AR5 base and scenarios)
  emi_scen <- emissions_TM5 %>% filter(YEAR %in% year) %>% select(-YEAR)
  
  # transform from Mt to kg
  if(IIASA_aneris_downscaling) {emi_scen[2:ncol(emi_scen)] <- emi_scen[2:ncol(emi_scen)] *10^9}
  
  # Get country list and compute number of countries
  cntrlist <- paste(emi_scen$COUNTRY)
  ncntr    <- length(cntrlist)
  
  
  #===============================================================================
  # PROCESSING DATA
  #===============================================================================
  
  if (u_verbose) cat(paste0("Processing data...\n"))
  # Compute emission differences (normalized to 20% perturbation = factor to be
  # multiplied with Source-Receptor (SR) coefficient)
  emi_base_proc <- emi_base %>%
    dplyr::rename(country=COUNTRY) %>%
    tidyr::gather(comp, base, -country)
  emi_scen_proc <- emi_scen %>%
    tidyr::gather(comp, scen, -COUNTRY)
  
  
  colnames(emi_scen_proc)[1] <- 'country'
  
  # Special PM2.5 treatment
  # we don't have separate SR for other primary PM2.5, and base emissions
  # do not contain PM2.5 so we use the SR for BC, and we scale with base BC
  emi_scen_proc <- emi_scen_proc %>%
    tidyr::spread(comp, scen) %>%
    dplyr::mutate(PM25 = ifelse(PM25 != NODATA, PM25 - ifelse(BC != NODATA, BC, 0) -
                                  ifelse(OM != NODATA, OM, 0), PM25)) %>%    # remove BC and POM from PM2.5
    dplyr::mutate(PM25 = ifelse(PM25 != NODATA & PM25 < 0.0, 0, PM25)) %>%   # Get rid of negative values (if any)
    tidyr::gather(comp, scen, -country)
  
  # Compute scaling factors for source-receptors
  # for CH4 the SR are already normalized to delta emission of perturbation
  # factor 5 because of 20% pertubation?
  emi_diff.df <- inner_join(emi_base_proc, emi_scen_proc, by=c("country","comp")) %>%
    dplyr::mutate(value = ifelse(scen != NODATA, (scen-base)/base*5.0, 0.0)) %>%
    dplyr::mutate(value = ifelse(is.na(value), 0, value)) # Remove NAs if any
  
  #-- Compute delta grid-cell concentrations from emission differences ------
  if (u_verbose) cat(paste0("Applying Source-Receptor matrix to emission differences...\n"))
  # Loop over source regions. Multiply with delta emissions and superimpose
  # receptor grids
  for (kcntr in cntrlist) {
    if(u_verbose) cat(paste0(kcntr, "\n"))
    
    # Restore the total (for BC and POM), SO2 and NOx (for SO4, NO3_A, and NH4)
    # perturbation response fields
    # *.sav delta files now contain delta's for aerosol, gases, deposition and M6M
    if (file.exists(paste0(dir_s, 'D_0.8_', kcntr, '_ALL.NC'))) {
      F_all   <- paste0(dir_s,  'D_0.8_', kcntr, '_ALL.NC')                                    # SO2+NOx+BC+OC perturbation 20%
      F_SO2   <- list.files(dir_s, paste0('*D_0.8_', kcntr, '_SO2.NC'), full.names=TRUE)       # SO2-only perturbation 20%
      F_NOX   <- list.files(dir_s, paste0('*D_0.8_', kcntr, '_NOx.NC'), full.names=TRUE)       # NOx-only perturbation 20%
      F_NMVOC <- list.files(dir_s, paste0('*D_0.8_', kcntr, '_NMVOC_NH3.NC'), full.names=TRUE) # NMVOC-NH3 only perturbation 20%
      
      
        # Read in Source-Receptor data (Generic)
        pertall <- nc_readall(F_all)
        # Read in Source-Receptor data (Specific)
        if (length(F_NMVOC) != 0) pertvoc_nh3 <- nc_readall(F_NMVOC)
        if (length(F_NOX)   != 0) pertnox     <- nc_readall(F_NOX)
        if (length(F_SO2)   != 0) pertso2     <- nc_readall(F_SO2)
     
      # Save emission differences in 1°x1° array
      for (k in c("SO2","NOX","BC","OM","PM25","NH3","VOC","CH4")) {
        emi_diff[[k]] <- array(get_emidiff(kcntr, k),  dim=c(360,180))
      }
      emi_diff_year_list[[year]] <- emi_diff
      
      #-- SO2 perturbation ------------
      if(u_verbose) cat("  - SO2 perturbation...\n")
      # If dSO2 perturbation calc is not available...use dALL
      # (should not be the case because all regions have been run for this perturbation case)
      if (length(F_SO2)   != 0) {
        dgrid[["SO4"]] <- dgrid[["SO4"]] + emi_diff_year_list[[year]][["SO2"]]*pertso2$var$DSO4
        dgrid[["NH4"]] <- dgrid[["NH4"]] + emi_diff_year_list[[year]][["SO2"]]*pertso2$var$DNH4
        dgrid[["NO3"]] <- dgrid[["NO3"]] + emi_diff_year_list[[year]][["SO2"]]*pertso2$var$DNO3_A
        
        dgrid[["SO2"]] <- dgrid[["SO2"]] + emi_diff_year_list[[year]][["SO2"]]*pertso2$var$DSO2
        dgrid[["NOX"]] <- dgrid[["NOX"]] + emi_diff_year_list[[year]][["SO2"]]*pertso2$var$DNOX
        dgrid[["O3"]]  <- dgrid[["O3"]]  + emi_diff_year_list[[year]][["SO2"]]*pertso2$var$DO3
        dgrid[["M6M"]] <- dgrid[["M6M"]] + emi_diff_year_list[[year]][["SO2"]]*pertso2$var$DM6M
      } else {
        dgrid[["SO4"]] <- dgrid[["SO4"]] + emi_diff_year_list[[year]][["SO2"]]*pertall$var$DSO4
        dgrid[["SO2"]] <- dgrid[["SO2"]] + emi_diff_year_list[[year]][["SO2"]]*pertall$var$DSO2
      }
      
      #-- NOx perturbation ------------
      if(u_verbose) cat("  - NOx perturbation...\n")
      # if separate dNOx perturbation calc is not available...
      # use dALL-dSO2 (this is the case for most countries)
      if (length(F_NOX)   != 0) {
        dgrid[["SO4"]] <- dgrid[["SO4"]] + emi_diff_year_list[[year]][["NOX"]]*pertnox$var$DSO4
        dgrid[["NH4"]] <- dgrid[["NH4"]] + emi_diff_year_list[[year]][["NOX"]]*pertnox$var$DNH4
        dgrid[["NO3"]] <- dgrid[["NO3"]] + emi_diff_year_list[[year]][["NOX"]]*pertnox$var$DNO3_A
        
        dgrid[["SO2"]] <- dgrid[["SO2"]] + emi_diff_year_list[[year]][["NOX"]]*pertnox$var$DSO2
        dgrid[["NOX"]] <- dgrid[["NOX"]] + emi_diff_year_list[[year]][["NOX"]]*pertnox$var$DNOX
        dgrid[["O3"]]  <- dgrid[["O3"]]  + emi_diff_year_list[[year]][["NOX"]]*pertnox$var$DO3
        dgrid[["M6M"]] <- dgrid[["M6M"]] + emi_diff_year_list[[year]][["NOX"]]*pertnox$var$DM6M
      } else {
        dgrid[["SO4"]] <- dgrid[["SO4"]] + emi_diff_year_list[[year]][["NOX"]]*(pertall$var$DSO4  -pertso2$var$DSO4)
        dgrid[["O3"]]  <- dgrid[["O3"]]  + emi_diff_year_list[[year]][["NOX"]]*(pertall$var$DNO3_A-pertso2$var$DNO3_A)
        dgrid[["NH4"]] <- dgrid[["NH4"]] + emi_diff_year_list[[year]][["NOX"]]*(pertall$var$DNH4  -pertso2$var$DNH4)
        
        dgrid[["SO2"]] <- dgrid[["SO2"]] + emi_diff_year_list[[year]][["NOX"]]*(pertall$var$DSO2 -pertso2$var$DSO2)
        dgrid[["NOX"]] <- dgrid[["NOX"]] + emi_diff_year_list[[year]][["NOX"]]*(pertall$var$DNOX -pertso2$var$DNOX)
        dgrid[["O3"]]  <- dgrid[["O3"]]  + emi_diff_year_list[[year]][["NOX"]]*(pertall$var$DO3  -pertso2$var$DO3)
        dgrid[["M6M"]] <- dgrid[["M6M"]] + emi_diff_year_list[[year]][["NOX"]]*(pertall$var$DM6M -pertall$var$DSO2)
      }
      
      if(u_verbose) cat("  - NMVOC perturbation...\n")
      if (length(F_NMVOC) != 0) {
        dgrid[["SO4"]] <- dgrid[["SO4"]] + emi_diff_year_list[[year]][["NH3"]]*pertvoc_nh3$var$DSO4
        dgrid[["NH4"]] <- dgrid[["NH4"]] + emi_diff_year_list[[year]][["NH3"]]*pertvoc_nh3$var$DNH4
        dgrid[["NO3"]] <- dgrid[["NO3"]] + emi_diff_year_list[[year]][["NH3"]]*pertvoc_nh3$var$DNO3_A
        
        dgrid[["O3"]]  <- dgrid[["O3"]]  + emi_diff_year_list[[year]][["VOC"]]*pertvoc_nh3$var$DO3
        dgrid[["M6M"]] <- dgrid[["M6M"]] + emi_diff_year_list[[year]][["VOC"]]*pertvoc_nh3$var$DM6M
      }
      
      if(u_verbose) cat("  - PM2.5 perturbation...\n")
      dgrid[["BC"]]   <- dgrid[["BC"]]   + emi_diff_year_list[[year]][["BC"]]  *pertall$var$DBC
      dgrid[["POM"]]  <- dgrid[["POM"]]  + emi_diff_year_list[[year]][["OM"]]  *pertall$var$DPOM
      dgrid[["PM25"]] <- dgrid[["PM25"]] + emi_diff_year_list[[year]][["PM25"]]*pertall$var$DBC 	  # USE BC SOURCE-RECEPTORS FOR PRIMARY PM2.5 (FLY ASH)
      
      if(u_verbose) cat("  - CH4 perturbation...\n")
      # Methane: We have only one global footprint file from HTAP1 SR2.
      # Impacts on O3 and M6M
      dgrid[["O3"]]  <- dgrid[["O3"]]  + emi_diff_year_list[[year]][["CH4"]]*pertch4$var$dM24
      dgrid[["M6M"]] <- dgrid[["M6M"]] + emi_diff_year_list[[year]][["CH4"]]*pertch4$var$dM6M
      
    }
    
  }
  
  #-- Compute new grid-cell concentrations ------
  if (u_verbose) cat(paste0("Computing new concentrations...\n"))
  # Compute base concentration year-average
  SO4_BASE  <- apply(conc_base_aer$var$so4,   c(1,2), sum)/12.0
  NO3_BASE  <- apply(conc_base_aer$var$no3_a, c(1,2), sum)/12.0
  NH4_BASE  <- apply(conc_base_aer$var$nh4,   c(1,2), sum)/12.0
  BC_BASE   <- apply(conc_base_aer$var$bc,    c(1,2), sum)/12.0
  POM_BASE  <- apply(conc_base_aer$var$pom,   c(1,2), sum)/12.0
  PM25_BASE <- array(0.0, dim=c(360,180))
  
  # Calculate new grid-cell average PM concentrations
  maps_ap <- list()
  maps_ap[["SO4"]]  <- SO4_BASE  + dgrid[["SO4"]] 		#; SO4  <- ifelse(SO4  > 0, SO4,  0)
  maps_ap[["NO3"]]  <- NO3_BASE  + dgrid[["NO3"]] 		#; NO3  <- ifelse(NO3  > 0, NO3,  0)
  maps_ap[["NH4"]]  <- NH4_BASE  + dgrid[["NH4"]] 		#; NH4  <- ifelse(NH4  > 0, NH4,  0)
  maps_ap[["BC"]]   <- BC_BASE   + dgrid[["BC"]]    	#; BC   <- ifelse(BC   > 0, BC,   0)
  maps_ap[["POM"]]  <- POM_BASE  + dgrid[["POM"]] 		#; POM  <- ifelse(POM  > 0, POM,  0)
  maps_ap[["PM25"]] <- PM25_BASE + dgrid[["PM25"]]	  #; PM25 <- ifelse(PM25 > 0, PM25, 0)
  
  # make sure concentration is not negative
  #maps_ap <- ifelse(maps_ap > 0, maps_ap, 0)
  maps_ap[["SO4"]]  <- ifelse(maps_ap[["SO4"]] > 0, maps_ap[["SO4"]], 0)
  maps_ap[["NO3"]]  <- ifelse(maps_ap[["NO3"]] > 0, maps_ap[["NO3"]], 0)
  maps_ap[["NH4"]]  <- ifelse(maps_ap[["NH4"]] > 0, maps_ap[["NH4"]], 0)
  maps_ap[["BC"]]   <- ifelse(maps_ap[["BC"]]  > 0, maps_ap[["BC"]],  0)
  maps_ap[["POM"]]  <- ifelse(maps_ap[["POM"]] > 0, maps_ap[["POM"]], 0)
  maps_ap[["PM25"]] <- ifelse(maps_ap[["PM25"]]> 0, maps_ap[["PM25"]],0)
  
  # Get base concentration year-average
  SS <- conc_base_dustsea$var$ss_pm25
  DU <- conc_base_dustsea$var$du_pm25
  
  # Calculate new grid-cell average SS and DU concentrations
  maps_ap[["SS"]] <- SS + dgrid[["SS"]] #; SS <- ifelse(SS > 0, SS, 0)
  maps_ap[["DU"]] <- DU + dgrid[["DU"]] #; DU <- ifelse(DU > 0, DU, 0)
  
  maps_ap[["SS"]] <- ifelse(maps_ap[["SS"]]> 0, maps_ap[["SS"]],0)
  maps_ap[["DU"]] <- ifelse(maps_ap[["DU"]]> 0, maps_ap[["DU"]],0)
  
  # Calculate new grid-cell average trace gas concentrations (ppb)
  SO2_BASE <-  apply(conc_base_gas$var$vmr_so2  ,  c(1,2), sum) /12.0*1E9
  NOX_BASE <- (apply(conc_base_gas$var$vmr_no   ,  c(1,2), sum) +
                 apply(conc_base_gas$var$vmr_no2  ,  c(1,2), sum))/12.0*1.E9
  O3_BASE  <-  apply(conc_base_gas$var$vmr_o3   ,  c(1,2), sum) /12.0*1.E9
  
  maps_ap[["SO2"]] <- SO2_BASE + dgrid[["SO2"]]		#; SO2_PPB <- ifelse(SO2_PPB > 0, SO2_PPB, 0)
  maps_ap[["NOx"]] <- NOX_BASE + dgrid[["NOX"]]		#; NOX_PPB <- ifelse(NOX_PPB > 0, NOX_PPB, 0)
  maps_ap[["O3"]]  <- O3_BASE  + dgrid[["O3"]]	  #; O3_PPB  <- ifelse(O3_PPB  > 0, O3_PPB,  0)
  
  maps_ap[["SO2"]] <- ifelse(maps_ap[["SO2"]]> 0, maps_ap[["SO2"]],0)
  maps_ap[["NOx"]] <- ifelse(maps_ap[["NOx"]]> 0, maps_ap[["NOx"]],0)
  maps_ap[["O3"]]  <- ifelse(maps_ap[["O3"]] > 0, maps_ap[["O3"]] ,0)
  
  
  # Calculate new grid-cell average M6M concentrations (ppb)
  m6m_base <- conc_base_m6m$var$M6M
  maps_ap[["M6M"]]  <- m6m_base + dgrid[["M6M"]]  #; M6M_PPB <- ifelse(M6M_PPB > 0, M6M_PPB, 0)
  maps_ap[["M3M"]]  <- maps_ap[["M6M"]] * conc_base_mf$var$F
  maps_ap[["M6M"]]  <- ifelse(maps_ap[["M6M"]] > 0, maps_ap[["M6M"]],0)
  maps_ap[["M3M"]]  <- ifelse(maps_ap[["M3M"]] > 0, maps_ap[["M3M"]],0)
  
  
  
  # Calculate H2O in aerosol
  PM_SEC <- maps_ap[["SO4"]] + maps_ap[["NO3"]] + maps_ap[["NH4"]]
  #50% RH (relative humidity)
  maps_ap[["H2O_50"]] <- 0.43*PM_SEC + 0.27*maps_ap[["SS"]]   #see file Seasalt_hygroscopicity.xls
  #35% RH
  H2O_35 <- 0.27*PM_SEC + 0.15*maps_ap[["SS"]]
  
  # Default when no INCR applied
  BC_ENH   <- maps_ap[["BC"]]
  POM_ENH  <- maps_ap[["POM"]]
  PM25_ENH <- maps_ap[["PM25"]]
  
  # Compute total anthropogenic PM2.5
  maps_ap[["PM_TOT"]] <- maps_ap[["SO4"]] + maps_ap[["NO3"]] + maps_ap[["NH4"]] + BC_ENH + POM_ENH + PM25_ENH
  #PM_URB=PM_TOT
  #PM_RUR=PM_TOT
  
  
  #===============================================================================
  # SAVE DATA
  # #===============================================================================
  data_df_maps_ap <- list()
  for (kv in names(maps_ap)) {
    data_df_maps_ap[[kv]] <- as.data.frame(maps_ap[[kv]])
    
    names(data_df_maps_ap[[kv]]) <- paste(seq(0,180-1))
    
    data_df_maps_ap[[kv]] <- data_df_maps_ap[[kv]] %>%
      cbind(data.frame(lon_id=seq(0,360-1))) %>%
      tidyr::gather(lat_id,value,-lon_id) %>%
      dplyr::mutate(lat_id=as.integer(lat_id))
  }
  
  
  
  # saving to only one single list does not work in the parallel execution
  # instead we save them in seperate .RData for every year
  #data_df_maps_ap_years[[year]] <- data_df_maps_ap
  #names(data_df_maps_ap_years[[year]]) = tolower(names(data_df_maps_ap))
  #save(data_df_maps_ap_years, file=paste0('concentrations/', RUNSUFFIX, file_output_con))
  
  
  names(data_df_maps_ap) = tolower(names(data_df_maps_ap))
  # return(data_df_maps_ap)
  save(data_df_maps_ap, file=paste0('concentrations/', year,'_', RUNSUFFIX, file_output_con))
  
  
  
  
  
}


#stopCluster(cl)
#names(data_df_maps_ap_years) <- u_year

#save the concentrations to disk
#save(data_df_maps_ap_years, file=paste0('concentrations/', RUNSUFFIX, file_output_con))


if (u_verbose) cat(paste0("RData file: ", file_output_con, RUNSUFFIX,' written.\n'))

tendend <- Sys.time() - t00
cat(paste0("[Time required to run the whole script: ", tendend,"]\n"))
