# ------------------------------------------------------------------------------
# Program Name: quick_fasst_gridmaps_pm_gas_downscaling_harmonization.R
# Author(s): Sebastian Rauner
#            original IDL script written by Rita van Dingenen (JRC-ISPRA)
#            rewritten in R by Jerome Hilaire (MCC & PIK) & Sebastian Rauner (PIK)
# Date Last Updated: Oct 17, 2017
# Program Purpose: This script downscales and harmonizes the emissions according to the chosen method of the launchpad
#                  You can chose between a simple GDP based downscaling without a harmonization and a IIASA downscaling with aneris CEDS based harmonization
# Input Files: 
# Output Files:
# To-Do: 

harm_scaling_launch_dir <- paste0(getwd(),'/downscaling/emissions_downscaling/exe/launchpad/')
harm_scaling_result_dir <- paste0(getwd(),'/downscaling/emissions_downscaling/final-output/module-B')

# read the region mappings
regionmapping_REMIND     <- read.csv2(mapping_REMIND,check.names=FALSE)
regionmapping_TM5        <- read.csv2(mapping_TM5,check.names=FALSE)


#downscaling with simply the GDP
if(GDP_downscaling == T)
{
  ################
  # start the downscaling and harmoniztaion substituing the first GDP based approach which follows
  ################
  
  # read the GDP data and filter for the SSP scenario
  gdp <- calcOutput("GDPppp",     years=u_year, aggregate = 'reg') 
  gdp_ssp <- gdp[,,paste0('gdp_',ssp_scenario)]
  
  # rename variable name to value (required by speed_aggregate)
  dimnames(gdp_ssp)$variable <- 'value'
  
  #TO-DO: scale REMIND emissions to CEDS
  
  
  #TO-DO: SCALE to country with more sophisticated weights than GDP
  
  # the weighting works, tested for ABW(LAM) solvents SOx 2010
  emissions_ISO <- speed_aggregate(remind_emissions_gdx, rel=regionmapping_REMIND, weight=gdp_ssp,
                                   from = 'RegionCode', to = 'CountryCode')
  
  # aggregate to TM5 regions
  emissions_TM5 <- speed_aggregate(emissions_ISO, rel=regionmapping_TM5,
                                   from = 'CountryCode', to = 'RegionCode')
}


#downscaling IIASA and aneris scripts
if(IIASA_aneris_downscaling == T)
{
  
  
  # harmonized year lists?
  
  source(paste0(harm_scaling_launch_dir,'launch_downscaling.R'))
  if(IIASA_gridding){source(paste0(harm_scaling_launch_dir,'launch_gridding.R'))}
  
  # read the IIASA downscaling and aneris harmonization results .csv
  emissions_ISO <-  read.csv(file = paste0(harm_scaling_result_dir,'/B.', iam, '_Harmonized_emissions_downscaled_for_gridding', '_', RUNSUFFIX, '.csv'))
  
  # delete first column which is empty
  emissions_ISO[1] <- NULL
  
  # what are the years columns:
  years_colnames <- colnames(emissions_ISO[grepl( "X" , names( emissions_ISO ) )])
  
  
  
  # filter here for sector '', these are the totals which should not be summed up
  emissions_ISO      <- emissions_ISO %>% filter(sector != '')
  
  emissions_AIR  <- emissions_ISO %>% filter(sector == 'Aircraft')
  emissions_SHIP <- emissions_ISO %>% filter(sector == 'International Shipping')
  
  # select iso, em and year columns
  # summarize over iso and em
  emissions_ISO <-    emissions_ISO                            %>%
    dplyr::select(model, scenario,iso, em, years_colnames)  %>%
    dplyr::group_by(model, scenario, iso, em)                %>%
    dplyr::summarise_all(.funs =  sum)      %>%
    dplyr::ungroup(iso, em)                 %>%
    dplyr::group_by(model, scenario, iso, em)
  
  emissions_AIR <-    emissions_AIR                            %>%
    dplyr::select(model, scenario,iso, em, years_colnames)  %>%
    dplyr::group_by(model, scenario, iso, em)                %>%
    dplyr::summarise_all(.funs =  sum)      %>%
    dplyr::ungroup(iso, em)                 %>%
    dplyr::group_by(model, scenario, iso, em)%>%
    dplyr::rename(FASST = iso)
  emissions_AIR$FASST <- 'AIR'
  
  emissions_SHIP<-    emissions_SHIP                            %>%
    dplyr::select(model, scenario,iso, em, years_colnames)  %>%
    dplyr::group_by(model, scenario, iso, em)                %>%
    dplyr::summarise_all(.funs =  sum)      %>%
    dplyr::ungroup(iso, em)                 %>%
    dplyr::group_by(model, scenario, iso, em)%>%
    dplyr::rename(FASST = iso)
  emissions_SHIP$FASST <- 'SHIP'
  
  # map TM5 regions and sum up over them 
  mapping_REMIND_ISO_FASST       <- read.csv2(file.path(dir_ancil, "OTHER", "mapping_ISO-FASST-REMIND.csv"))
  mapping_REMIND_ISO_FASST$ISO.3 <- tolower(mapping_REMIND_ISO_FASST$ISO.3)
  
  # join mapping here  
  emissions_ISO <- left_join(
    mapping_REMIND_ISO_FASST %>%
      dplyr::select(FASST, ISO.3),
    emissions_ISO,
    by=c("ISO.3"="iso"))
  
  emissions_ISO <- na.omit(emissions_ISO)
  
  

  emissions_TM5 <-    emissions_ISO                            %>%
    dplyr::select(model, scenario, FASST, em, years_colnames)  %>%
    dplyr::group_by(model, scenario, FASST, em)                %>%
    dplyr::summarise_all(.funs =  sum)                         %>%
    dplyr::ungroup()                             
  
  emissions_TM5 <- dplyr::bind_rows(emissions_TM5, emissions_AIR, emissions_SHIP)
  #  wide format
  emissions_TM5 <-   emissions_TM5                     %>%
    dplyr::select(-model, -scenario, FASST, em, years_colnames) %>%
    tidyr::gather(data =  emissions_TM5, key = 'period', years_colnames) 
  # some renaming of colnames
  names(emissions_TM5)[names(emissions_TM5) == '.']     <- 'value'
  names(emissions_TM5)[names(emissions_TM5) == 'em']    <- 'emi'
  names(emissions_TM5)[names(emissions_TM5) == 'FASST'] <- 'region'
  
  emissions_TM5$period <- gsub('X', '', emissions_TM5$period)
  
  remind_emissions_gdx_avi_intship <- NULL
}




emissions_TM5 <- as.quitte(emissions_TM5) %>% dplyr::select(period, region, emi, value) 

# Add global shipping and aviation as SHIP and AIR
# REMIND 

emissions_TM5 <- rbind(emissions_TM5, remind_emissions_gdx_avi_intship)

# rename emi to upper
levels(emissions_TM5$emi) <-toupper(levels(emissions_TM5$emi))

# is OC (REMIND) the same as OM (used in TM5)
levels(emissions_TM5$emi) <- c(levels(emissions_TM5$emi), "OM")
levels(emissions_TM5$emi) <- c(levels(emissions_TM5$emi), "SO2")
levels(emissions_TM5$emi) <- c(levels(emissions_TM5$emi), "NOX")

emissions_TM5[emissions_TM5$emi == 'OC', 'emi']     <- 'OM'
emissions_TM5[emissions_TM5$emi == 'Sulfur', 'emi'] <- 'SO2'
emissions_TM5[emissions_TM5$emi == 'SULFUR', 'emi'] <- 'SO2'
emissions_TM5[emissions_TM5$emi == 'NOx', 'emi']    <- 'NOX'
droplevels(emissions_TM5$emi)

# these emis are used by the model c("SO2","NOX","BC","OM","PM25","NH3","VOC","CH4")
# REMIND gives us BC OC CH4 VOC CO NOX SO2 NMVOC SOX


# convert to a wide format according to emi [COMP]
emissions_TM5 <-
  emissions_TM5 %>%
  tidyr::spread(emi,value)   %>%
  dplyr::mutate("PM25" = VOC * 0) %>%
  dplyr::mutate("NH3" = VOC * 0) %>%
  dplyr::select(period, region,'BC','CH4','CO','NH3','NOX','OM','SO2','VOC','PM25') %>%
  dplyr::rename('YEAR' = period, 'COUNTRY' = region)

setwd(wd_dir)
#save the emissions to disk
save(emissions_TM5, file=paste0(getwd(),'/emissions/', RUNSUFFIX, 'emissions.RData'))

