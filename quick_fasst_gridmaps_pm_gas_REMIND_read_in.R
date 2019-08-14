# ------------------------------------------------------------------------------
# Program Name: quick_fasst_gridmaps_pm_gas_cluster_Seb_REMIND.R
# Author(s): Sebastian Rauner
#            original IDL script written by Rita van Dingenen (JRC-ISPRA)
#            rewritten in R by Jerome Hilaire (MCC & PIK) & Sebastian Rauner (PIK)
# Date Last Updated: Oct 17, 2017
# Program Purpose: The script reads the REMIND AP emission gdx output and does some formatting
# Input Files: 
# Output Files:
# To-Do: 

#===============================================================================
# FASST-R
#===============================================================================

# Produces gridmaps of all components (SO4, NO3a, NH4, BC, POM, other PM2.5,
# dust, seasalt, SO2, NOx, O3, M6M, M3M) from the input emission scenarios.
# The resulting fields are calculated as total (from delta emissions with base
# emissions, and respective SR, superimposed on the AR5 base run).
# The analysis of the fields (e.g. delta concentrations between scenarios,
# associated impacts) has to be done separately from produced output.

# WARNING: the grid level is provided for easy exposure calculations of
# population and crops, but it is recommended to produce final results at the
# aggregation level at which input emissions were generated, in particular when
# looking at u_projectNameed scenarios where trends are established based on
# indicators for country groups.

# CH4 feedback on global background O3 is included, but not longterm feedback of
# NOx and other precursors on CH4 lifetime and background O3
# CO Source-Receptors not yet included - work in progress.
# M6M Maximal 6-month mean of daily maximal hourly ozone

# use the new SSP emission data when available
# add a more sophisticated reporting with all aerosols?
#===============================================================================

t00 <- Sys.time()
u_verbose      = TRUE




#===============================================================================
# INITIALISATION
#===============================================================================
if (u_verbose) cat(paste0("Initialising...\n"))
# Load libraries
library(ncdf4)
#library(xlsx)
library(luscale) # rename_dimnames
library(remind)
library(quitte)  #read.gdx
library(gdx)
library(doParallel) # cluster
library(moinput) # calcGDP
library(tidyr)
library(dplyr)
library(stringr)


# Own functions
nc_readall <- function (ncfile) {
  fid = nc_open(ncfile)
  out = list()
  out[["dim"]] = list()
  out[["var"]] = list()
  out[["units"]] = list()
  out[["long_name"]] = list()
  out[["standard_name"]] = list()
  #out[["global_attributes"]] = list()

  for (kdim in names(fid$dim)) {
    out[["dim"]][[kdim]] = ncvar_get(fid, kdim)
  }

  for (kvar in names(fid$var)) {
    out[["var"]][[kvar]] = ncvar_get(fid, kvar)
    out[["units"]][[kvar]] = fid$var[[kvar]]$units
    out[["long_name"]][[kvar]] = fid$var[[kvar]]$long_name
    out[["standard_name"]][[kvar]] = fid$var[[kvar]]$standard_name
  }
  
  return(out)
  nc_close(fid)

}

get_emidiff <- function(country, comp) {
  data = emi_diff.df$value[which(emi_diff.df$country == country & emi_diff.df$comp == comp)]
  return(ifelse(is.finite(data), data, 0.0))
}

# is not used anymore, instead the gdx file is read
#Function to read all mifs in folder
readAllReportingMIFinFolder <- function (dir) {
  filelist <- list.files(dir, "*.mif", full.names = TRUE)
  if (length(filelist) == 0) {
    filelist <- list.files(dir, "*.csv", full.names = TRUE)
  }
  if (length(filelist) == 0) {
    stop(paste("No mif or csv files found in directory", 
               dir))
  }
  
  else {
    cat("Reading data from mif files...\n")
    print(filelist,sep="\n")
    data <- do.call("mbind", lapply(filelist, function(x) read.report(x, as.list = FALSE)))
  }
  return(data)
}

# Function to rename variables
# from: SSP2-Ref-SPA0-V13.REMIND-MAGPIE.Emissions|VOC|Energy Demand|Industry (Mt VOC/yr)
#   to: SSP2-Ref-SPA0-V13.REMIND-MAGPIE.VOC.indst
name_mif2magpie    <- function(x,str_emi=NULL,str_sector=NULL) {
  
  if (is.null(str_emi))    str_emi    <- "(CO2|CO|NOX|(Sulfur|SO2)|VOC|NH3|BC|OC|N2O|CH4|NMHC|NO2)"
  if (is.null(str_sector)) str_sector <- "(\\(|Total) \\(|Energy Supply|Electricity) \\(|Residential and Commercial|ResCom|Ground (Transportation|Transport|Trans)|Energy Demand\\|Industry|Solvents|Waste \\(|Aviation|Shipping|Land Use \\(|Land Use\\|(Agriculture |Forest Burning|Savannah Burning|Agricultural (waste|Waste) Burning))"
  
  # change |emissions to |Total
  getNames(x) <- gsub(paste0(str_emi," \\("),paste0("\\1\\|Total \\("),getNames(x))
  vars <- getNames(x)[grepl(paste0(".*(Emissions|Emi)\\|",str_emi,"\\|.*",str_sector),getNames(x))]
  y    <- x[,,sort(vars)]
  
  # remove "Emissions|"
  getNames(y) <- gsub(paste0("Emissions\\|"),"",getNames(y)) # for output of calcOutput (MOINPUT) without scenario and model
  # remove "Emi|"
  getNames(y) <- gsub(paste0("Emi\\|"),"",getNames(y))       # for output of calcOutput (MOINPUT) without scenario and model
  # remove unit
  getNames(y) <- gsub(" \\(Mt( |/).*","",getNames(y))
  getNames(y) <- gsub(" \\(kt( |/).*","",getNames(y))
  getNames(y) <- gsub(" \\(MCOO( |/).*","",getNames(y))
  # rename Sulfur
  getNames(y) <- gsub("Sulfur","SO2",getNames(y))
  # replace dash after emission name with dot
  getNames(y) <- gsub(paste0(str_emi,"\\|"),"\\1.",getNames(y))
  
  
  # reverse order of emi and sector to the same order used in CEDS16
  # SSP2-19-SPA0-V16.REMIND-MAGPIE.SO2.Aircraft -> REMIND-MAGPIE.SSP2-19-SPA0-V16.SO2.Aircraft
  getNames(y) <- gsub("^([^\\.]*)\\.([^\\.]*)\\.([^\\.]*)\\.(.*$)","\\2.\\1.\\3.\\4",getNames(y))
  getSets(y) <- c("Region","Year","Model","Scenario","Species","Sector")
  return(y)
}

# Function to rename variables from reporting to match IIASA & aneris format
gdx_to_IIASA_aneris <- function(gdx){
  
  
  ######### initialisation  ###########
  airpollutants <- c("ch4","so2","bc","oc","CO","VOC","NOx","NH3")
 
  
  ######### internal function  ###########
  generateReportingEmiAirPol <- function(pollutant,i_emiAPexsolve=pm_emiAPexsolve,i_emiAPexo=pm_emiAPexo){
    poll_rep <- toupper(pollutant)
    tmp <- NULL
    
    # reduce to the pollutant
    emiAPexsolve      <- collapseNames(i_emiAPexsolve[,,pollutant])
    emiAPexo          <- collapseNames(i_emiAPexo[,,pollutant])
    getSets(emiAPexo) <- getSets(emiAPexsolve)
    
    
    #add ch4 manually
    # double counting, no problem because of harmonization
    if(pollutant == 'ch4'){
    emiAPexsolve[,,"power"] <- vm_emiMacSector[,,"ch4coal"] + vm_emiMacSector[,,"ch4oil"] + vm_emiMacSector[,,"ch4gas"]
    emiAPexsolve[,,"indst"] <- vm_emiMacSector[,,"ch4coal"] + vm_emiMacSector[,,"ch4oil"] + vm_emiMacSector[,,"ch4gas"]
    emiAPexsolve[,,"res"]   <- vm_emiTe

    emiAPexo[,,"AgWasteBurning"]  <- vm_emiMacSector[,,"ch4agwaste"]
    emiAPexo[,,"ForestBurning"]   <- vm_emiMacSector[,,"ch4forest"]
    emiAPexo[,,"GrasslandBurning"]<- vm_emiMacSector[,,"ch4savan"]
    }
    
    if(pollutant != 'ch4'){
    # add indprocess to indst
    emiAPexsolve[,,"indst"] <- emiAPexsolve[,,"indst"] + emiAPexsolve[,,"indprocess"]
    }
    # Replace REMIND sector names by reporting ones
    mapping = data.frame(
      remind = c("power", "indst", "res", "trans", "solvents", "extraction"),
      reporting = c(paste0("Emissions|", poll_rep, "|Energy Supply|Electricity (Mt ", poll_rep, "/yr)"),
                    paste0("Emissions|", poll_rep, "|Energy Demand|Industry (Mt ", poll_rep, "/yr)"), 
                    paste0("Emissions|", poll_rep, "|Energy Demand|Residential and Commercial (Mt ", poll_rep, "/yr)"),
                    paste0("Emissions|", poll_rep, "|Energy Demand|Transportation|Ground Transportation (Mt ", poll_rep, "/yr)"),
                    paste0("Emissions|", poll_rep, "|Solvents (Mt ", poll_rep, "/yr)"), 
                    paste0("Emissions|", poll_rep, "|Energy Supply|Extraction (Mt ", poll_rep, "/yr)")))
    
    emiAPexsolve <- setNames(emiAPexsolve[,,mapping$remind],as.character(mapping$reporting)) 
    
    tmp <- 
      mbind(emiAPexsolve,
            setNames(emiAPexo[,,"AgWasteBurning"],  paste0("Emissions|",poll_rep,"|Land Use|Agricultural Waste Burning (Mt ",poll_rep,"/yr)")),
            setNames(emiAPexo[,,"Agriculture"],     paste0("Emissions|",poll_rep,"|Land Use|Agriculture (Mt ",poll_rep,"/yr)")),
            setNames(emiAPexo[,,"ForestBurning"],   paste0("Emissions|",poll_rep,"|Land Use|Forest Burning (Mt ",poll_rep,"/yr)")),
            setNames(emiAPexo[,,"GrasslandBurning"],paste0("Emissions|",poll_rep,"|Land Use|Savannah Burning (Mt ",poll_rep,"/yr)")),
            setNames(emiAPexo[,,"Waste"],           paste0("Emissions|",poll_rep,"|Waste (Mt ",poll_rep,"/yr)")))
    
    
    # Set NAs to 0
    tmp[is.na(tmp)] <- 0
    
    return(tmp)
  }
  
  ####### conversion factors ##########
  pm_conv_TWa_EJ    <- 31.536
  conv_MtSO2_to_MtS <- 1/2     # 32/(32+2*16)
  
  ####### read in needed data #########
  ## sets
  ttot  <-  as.numeric(readGDX(gdx, name=c("ttot"), format="first_found"))
  ## parameter
  pm_emiAPexsolve   <- readGDX(gdx, name=c("pm_emiAPexsolve"), field="l", format="first_found")[,ttot,]
  pm_emiAPexo       <- readGDX(gdx, name=c("pm_emiAPexo"), field="l", format="first_found")[,ttot,airpollutants]
  pm_emiAPexoGlob   <- readGDX(gdx, name=c("pm_emiAPexoGlob"), field="l", format="first_found")[,ttot,airpollutants]
  
  #ch4
  vm_emiTe          <- readGDX(gdx, name=c("vm_emiTe"), field="l", format="first_found")[,ttot,'ch4']
  vm_emiMacSector   <- readGDX(gdx, name=c("vm_emiMacSector"), field="l", format="first_found")[,ttot,]
  
  ####### prepare parameter ########################
  magclass::getNames(pm_emiAPexsolve) <- gsub("SOx","so2",magclass::getNames(pm_emiAPexsolve))
  magclass::getNames(pm_emiAPexsolve) <- gsub("NMVOC","VOC",magclass::getNames(pm_emiAPexsolve))
  
  ####### calculate reporting parameters ############
  # Loop over air pollutants and call reporting generating function
  out <- do.call("mbind", lapply(airpollutants, generateReportingEmiAirPol))
  
  # Add global values
  out   <- mbind(out, dimSums(out,dim=1))
  
  # Loop over air pollutants and add some variables
  for (pollutant in airpollutants) {
    poll_rep <- toupper(pollutant)
    tmp <- NULL
    # Add Aviation and Int. Shipping emissions
    tmp <- mbind(tmp,setNames(pm_emiAPexoGlob["GLO",,"InternationalShipping"][,,pollutant],paste0("Emissions|",poll_rep,"|Energy Demand|Transportation|International Shipping (Mt ",poll_rep,"/yr)")),
                 setNames(pm_emiAPexoGlob["GLO",,"Aviation"][,,pollutant],             paste0("Emissions|",poll_rep,"|Energy Demand|Transportation|Aviation (Mt ",poll_rep,"/yr)"))
    )
    tmp1 <- new.magpie(getRegions(out),getYears(out),magclass::getNames(tmp),fill=0)
    tmp1["GLO",,] <- tmp["GLO",,]
    out  <- mbind(out,tmp1)
    # Aggregation: Transportation and Energy Supply
    out <- mbind(out,
                 setNames(dimSums(out[,,
                                      c(paste0("Emissions|",poll_rep,"|Energy Demand|Transportation|Ground Transportation (Mt ",poll_rep,"/yr)"),
                                        paste0("Emissions|",poll_rep,"|Energy Demand|Transportation|International Shipping (Mt ",poll_rep,"/yr)"),
                                        paste0("Emissions|",poll_rep,"|Energy Demand|Transportation|Aviation (Mt ",poll_rep,"/yr)"))],dim = 3),
                          paste0("Emissions|",poll_rep,"|Energy Demand|Transportation (Mt ",poll_rep,"/yr)")),
                 setNames(dimSums(out[,,
                                      c(paste0("Emissions|",poll_rep,"|Energy Supply|Electricity (Mt ",poll_rep,"/yr)"),
                                        paste0("Emissions|",poll_rep,"|Energy Supply|Extraction (Mt ",poll_rep,"/yr)"))],dim = 3),
                          paste0("Emissions|",poll_rep,"|Energy Supply (Mt ",poll_rep,"/yr)"))
    )
    # Aggregation: Energy Demand + Energy Supply, Land Use
    out <- mbind(out, 
                 setNames(dimSums(out[,,c(paste0("Emissions|",poll_rep,"|Energy Demand|Industry (Mt ",poll_rep,"/yr)"),
                                          paste0("Emissions|",poll_rep,"|Energy Demand|Residential and Commercial (Mt ",poll_rep,"/yr)"),
                                          paste0("Emissions|",poll_rep,"|Energy Demand|Transportation (Mt ",poll_rep,"/yr)"),
                                          paste0("Emissions|",poll_rep,"|Energy Supply (Mt ",poll_rep,"/yr)"))],dim = 3),
                          paste0("Emissions|",poll_rep,"|Energy Supply and Demand (Mt ",poll_rep,"/yr)")),
                 setNames(dimSums(out[,,c(paste0("Emissions|",poll_rep,"|Land Use|Savannah Burning (Mt ",poll_rep,"/yr)"),
                                          paste0("Emissions|",poll_rep,"|Land Use|Forest Burning (Mt ",poll_rep,"/yr)"),
                                          paste0("Emissions|",poll_rep,"|Land Use|Agriculture (Mt ",poll_rep,"/yr)"),
                                          paste0("Emissions|",poll_rep,"|Land Use|Agricultural Waste Burning (Mt ",poll_rep,"/yr)"))],dim = 3),
                          paste0("Emissions|",poll_rep,"|Land Use (Mt ",poll_rep,"/yr)"))
    )
    # Compute total
    out <- mbind(out,
                 setNames(dimSums(out[,,c(paste0("Emissions|",poll_rep,"|Energy Supply and Demand (Mt ",poll_rep,"/yr)"),
                                          paste0("Emissions|",poll_rep,"|Solvents (Mt ",poll_rep,"/yr)"),
                                          paste0("Emissions|",poll_rep,"|Land Use (Mt ",poll_rep,"/yr)"),
                                          paste0("Emissions|",poll_rep,"|Waste (Mt ",poll_rep,"/yr)"))],dim=3),
                          paste0("Emissions|",poll_rep," (Mt ",poll_rep,"/yr)"))
    )
  }
  
  #adjust the format and naming according to what the IIASA & aneris script needs
  out <- as.quitte(out)
  
  out$variable <- out$all_sectorEmi 
  out$all_sectorEmi <- NULL
  
  # put the units in the units column and cut them from the variables column
  re <- "\\(([^()]+)\\)"
  out$unit <- gsub(re, "\\1", str_extract_all( out$variable, re))
  
  out$variable <- gsub('\\ \\(.*?\\)', '', out$variable)
  
  # get the model (hardcoded) and scenario from the remind config
  out$model    <- paste(model)
  out$scenario <- ssp_scenario
  
  #append the harmonization status Unharmonized
  out$variable <- paste0(out$variable,'|Unharmonized')
  
  
  #capitalize first character of colnames
  colnames(out) <- str_to_title(colnames(out))
  
  
  #spread to wide format
  out <-  spread(out, Period, Value)
  
  return(out)
  }


# Output NetCDF file parameters
u_projectName <- 'REMIND-FASST'
u_progName    <- 'quick_fasst_gridmaps_pm_gas_cluster_Seb_REMIND.R'

# Folders
dir_root      <- paste0(getwd()    ,'/mnt/FASST_IDL_PACK/')               # main folder

dir_out       <- paste0(dir_root,  'OUTPUT/')                             # output
dir_ancil     <- paste0(dir_root,  'ANCILLARY/')                          # ancillary data ()
dir_reg       <- paste0(dir_ancil, 'FASST_REGION_MASK/')                  # region data
dir_s         <- paste0(dir_root,  'CNTRY_2_GRID_V4/NCDF/')               # Source-Receptor data
dir_emi       <- paste0(dir_root,  'RCP_EMISS')                           # RCP input emissions
dir_remind    <- scenario_folder                                          # REMIND  output emissions
dir_spat_map  <- paste0(dir_ancil, 'SPATIAL_MAPPING/')                    # Spatial mapping

# Files
file_abase   <- paste0(dir_ancil, 'BASE/nTM5-JRC-cy2-ipcc-v1_GLOBAL_SFC_aerosolm_2001.nc')
file_gbase   <- paste0(dir_ancil, 'BASE/TM5-JRC-cy2-ipcc-v1_GLOBAL_SFC_tracerm_2001.nc')
file_dustsea <- paste0(dir_ancil, 'BASE/GBD2010_PM25_DU_SS.nc')
file_m6m     <- paste0(dir_ancil, 'BASE/BASE_GLOBAL_M6M_boxmid.nc')
file_mf      <- paste0(dir_ancil, 'O3_HEALTH/m3m_to_m6m_factor.nc')

file_pertch4 <- paste0(dir_s,     'HTAP_CH4/D_NORMALIZED_GLOBAL_CH4.nc')

file_embase  <- paste0(dir_ancil, 'BASE/TOTAL_COUNTRY_EM.RCP.TXT')
#load(file=file.path(dir_emi,  paste0('RCP_EMISS_TM5_FAAST_RCP',u_rcp,'.RData')))


remind_emiss       <- paste0(dir_remind,     'fulldata.gdx')
remind_config      <- paste0(dir_remind,     'config.Rdata')
mapping_REMIND     <- paste0(dir_spat_map,   'regionmappingREMIND.csv')
mapping_TM5        <- paste0(dir_spat_map,   'regionmappingTM5.csv')



file_output_con  <- 'concentrations.RData'

# read REMIND config file to get SSP scenario
load(remind_config)
ssp_scenario <- cfg$gms$c_LU_emi_scen
#for testing
model        <- 'REMIND-MAGPIE'
p_projectName <- cfg$gms$c_expname

AP_scenario <- cfg$gms$cm_APscen



#if rcp scenario not set use the SSP to rcp mapping
if(cfg$gms$cm_rcp_scen!= 'none'){u_rcp <- cfg$gms$cm_rcp_scen} else{
  
  # rcp to ssp mapping, only necessary untill the SSP data is available
  if(ssp_scenario == "SSP1"){u_rcp <- '26'}
  if(ssp_scenario == "SSP2"){u_rcp <- '26'}
  if(ssp_scenario == "SSP3"){u_rcp <- '45'}
  if(ssp_scenario == "SSP4"){u_rcp <- '60'}
  if(ssp_scenario == "SSP5"){u_rcp <- '85'}
}

if(GDP_downscaling == T)
{

  downscaling_method <- 'GDP_downscaling'

## read in REMIND emission data from fulldata.gdx, rename, reorder, append and filter for u_year

#pm_emiAPexsolve
remind_emissions_gdx_pm_emiAPexsolve <-  read.gdx(remind_emiss, 'pm_emiAPexsolve')
remind_emissions_gdx_pm_emiAPexsolve <-  remind_emissions_gdx_pm_emiAPexsolve %>% dplyr::rename(  t = tall
                                                                                          ,regi = all_regi
                                                                                          ,sector = all_sectorEmi
                                                                                          ,emi = emiRCP)
                                                                              
#pm_emiAPexo
remind_emissions_gdx_pm_emiAPexo <-  read.gdx(remind_emiss, 'pm_emiAPexo')
remind_emissions_gdx_pm_emiAPexo <-  remind_emissions_gdx_pm_emiAPexo %>% dplyr::rename(  t = ttot
                                                                                  ,regi = all_regi
                                                                                  ,emi = all_enty
                                                                                  ,sector = all_exogEmi)
remind_emissions_gdx_pm_emiAPexo < - remind_emissions_gdx_pm_emiAPexo[c("t", "regi", "sector","emi","value")]

#pm_emiAPexoGlob
remind_emissions_gdx_pm_emiAPexoGlob <-  read.gdx(remind_emiss, 'pm_emiAPexoGlob')
remind_emissions_gdx_pm_emiAPexoGlob <-  remind_emissions_gdx_pm_emiAPexoGlob %>% dplyr::rename(  t = ttot
                                                                                          ,emi = all_enty
                                                                                          ,sector = all_exogEmi) %>%   
                                                                                  dplyr::mutate(regi = 'GLO')
remind_emissions_gdx_pm_emiAPexoGlob <- remind_emissions_gdx_pm_emiAPexoGlob[c("t", "regi", "sector","emi","value")]

#remind_emissions_gdx <- the three gdx appended, filter for t in u_year
remind_emissions_gdx <- rbind(remind_emissions_gdx_pm_emiAPexsolve,
                              remind_emissions_gdx_pm_emiAPexo) %>%
                        filter(t %in% u_year)


remind_emissions_gdx_avi_intship <- remind_emissions_gdx_pm_emiAPexoGlob %>%
                                    filter(t %in% u_year)

# adjust units to TM5 ones used #here
# REMIND uses Mt, TM5 kg? therefore 10^9
# N2O reported in kt, adjust if added

remind_emissions_gdx$value <- remind_emissions_gdx$value * 10^9
remind_emissions_gdx_avi_intship$value <- remind_emissions_gdx_avi_intship$value * 10^9


# sum up all industries
remind_emissions_gdx <- remind_emissions_gdx %>%
                        dplyr::group_by(t, regi, emi) %>%
                        dplyr::summarize(value = sum(value))
# change global sector names to the ones used in TM5 (AIR, SHIP)

levels(remind_emissions_gdx_avi_intship$sector) <- c(levels(remind_emissions_gdx_avi_intship$sector), "SHIP", "AIR")

remind_emissions_gdx_avi_intship[remind_emissions_gdx_avi_intship$sector == 'Aviation', 'sector'] <- 'AIR'
remind_emissions_gdx_avi_intship[remind_emissions_gdx_avi_intship$sector == 'InternationalShipping','sector'] <- 'SHIP'

droplevels(remind_emissions_gdx_avi_intship$sector)


remind_emissions_gdx_avi_intship <- remind_emissions_gdx_avi_intship %>%
                                    dplyr::select(-regi)%>%
                                    dplyr::rename(period = t, region = sector)

remind_emissions_gdx <-  as.magpie(remind_emissions_gdx)

}


if(IIASA_aneris_downscaling == T)
{
  downscaling_method <- 'IIASA_aneris_downscaling'
  remind_emissions_gdx <- gdx_to_IIASA_aneris(gdx = remind_emiss) #filter for years in the gdx_to_IIASA_aneris function?

# # all the rest is for extracting shipping and aviation and assigne it to regi=GLO
# # is that really needed, or is it double counted later?
#   best would be to do the emission diff not only on the country (TM5 region) level but the grid level and use the IIASA gridding files as input
# 
# #pm_emiAPexoGlob
# remind_emissions_gdx_pm_emiAPexoGlob <-  read.gdx(remind_emiss, 'pm_emiAPexoGlob')
# remind_emissions_gdx_pm_emiAPexoGlob <-  remind_emissions_gdx_pm_emiAPexoGlob %>% dplyr::rename(  t = ttot
#                                                                                                   ,emi = all_enty
#                                                                                                   ,sector = all_exogEmi) %>%   
#   dplyr::mutate(regi = 'GLO')
# remind_emissions_gdx_pm_emiAPexoGlob <- remind_emissions_gdx_pm_emiAPexoGlob[c("t", "regi", "sector","emi","value")]
# 
# 
# remind_emissions_gdx_avi_intship <- remind_emissions_gdx_pm_emiAPexoGlob %>%
#   filter(t %in% u_year)
# 
# # adjust units to TM5 ones used #here
# # REMIND uses Mt, TM5 kg? therefore 10^9
# # N2O reported in kt, adjust if added
# 
# remind_emissions_gdx_avi_intship$value <- remind_emissions_gdx_avi_intship$value * 10^9
# 
# 
# # change global sector names to the ones used in TM5 (AIR, SHIP)
# 
# levels(remind_emissions_gdx_avi_intship$sector) <- c(levels(remind_emissions_gdx_avi_intship$sector), "SHIP", "AIR")
# 
# remind_emissions_gdx_avi_intship[remind_emissions_gdx_avi_intship$sector == 'Aviation', 'sector'] <- 'AIR'
# remind_emissions_gdx_avi_intship[remind_emissions_gdx_avi_intship$sector == 'InternationalShipping','sector'] <- 'SHIP'
# 
# droplevels(remind_emissions_gdx_avi_intship$sector)
# 
# 
# remind_emissions_gdx_avi_intship <- remind_emissions_gdx_avi_intship %>%
#   dplyr::select(-regi)%>%
#   dplyr::rename(period = t, region = sector)

}
