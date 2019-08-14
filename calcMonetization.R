# calculate the VSL as a function of time and region
# this script implements the methodology proposed by the OECD and used by the lancet commission
# http://www.oecd-ilibrary.org/environment/mortality-risk-valuation-in-environment-health-and-transport-policies_9789264130807-en
# 10.1016/S0140-6736(17)32345-0


# see http://www.thelancet.com/cms/attachment/2114605792/2084555537/mmc1.pdf pp 25
# we calculate the VSL_coeff, we groupe the countries in high, middle and low and
# assigne correspondinfg eslasticities

library(dplyr)
library(tidyr)
library(moinput)

#setwd("C:/Users/rauner/Documents/PIK/fasstr_sebastian")

calcMonetization <- function(ssp_scenario                  = 'ssp2',
                             u_set_ref_regi                = 'EUR' ,
                             u_set_ref_per                 = 2005  ,
                             u_set_eleasticity_interval    = 0.2   ,  
                             u_set_eleasticity_high_income = 0.8   ,
                             u_set_eleasticity_low_income  = 1.2    ) {

ssp_scenario   <- toupper(ssp_scenario)


if(!(u_set_ref_per %in% u_year)){year <- c(u_set_ref_per, u_year)} else(year <-  u_year)




# the source takes per capita income, we use per capita GDP
# GDPppp
gdp <- calcOutput("GDPppp",     years=year, aggregate = F) 
gdp_ssp <- gdp[,,paste0('gdp_',ssp_scenario)]


# population
pop <- calcOutput("Population",     years=year, aggregate = F) 
pop_ssp <- pop[,,paste0('pop_',ssp_scenario)]

# GDPppp per capita
gdp_capita <- gdp_ssp/pop_ssp

# get the reference GDPppp/capita for the reference region for which we have a WTP for VSL
gdp <- calcOutput("GDPppp",     years=u_set_ref_per, aggregate = TRUE) 
gdp_ssp_ref <- gdp[u_set_ref_regi,,paste0('gdp_',ssp_scenario)]

pop <- calcOutput("Population",     years=u_set_ref_per, aggregate = TRUE) 
pop_ssp_ref <- pop[u_set_ref_regi,,paste0('pop_',ssp_scenario)]

gdp_capita_ref <- gdp_ssp_ref/pop_ssp_ref

# assgine the elasticity according to the relation to the reference gdp_capita_ref in the reference year
# keep the assignment static, e.g. it does not matter how the EUR to other regions develop
# we use 1.2 for low income and 0.8 for high income countries, < and > then EUR

gdp_capita <- add_columns(gdp_capita, addnm = "eleasticity", dim = 3.1)

gdp_capita[,u_set_ref_per,'eleasticity'] <- case_when(gdp_capita[,u_set_ref_per,paste0('gdp_',ssp_scenario,'.','pop_',ssp_scenario)] * (1 + u_set_eleasticity_interval) > as.numeric(gdp_capita_ref)  ~ 0.8,
                                         gdp_capita[,u_set_ref_per,paste0('gdp_',ssp_scenario,'.','pop_',ssp_scenario)] * (1 - u_set_eleasticity_interval) < as.numeric(gdp_capita_ref)  ~ 1.2,
                                         TRUE ~ 1)


VSL_coeff  <- (gdp_capita[,,paste0('gdp_',ssp_scenario,'.','pop_',ssp_scenario)]/as.numeric(gdp_capita_ref)) ** as.numeric(gdp_capita[,u_set_ref_per,'eleasticity'])

# only return the u_years and not the reference year if that is not in u_years

# input write out
  # VSL_coeff <- as.data.frame(VSL_coeff)
  # VSL_coeff$Cell <-NULL
  # VSL_coeff$Data1<-NULL
  # VSL_coeff$Data2<-'vsl_coeff'
  # write.csv(as.data.frame(VSL_coeff), 'vsl_coeff.csv')

return(VSL_coeff[,u_year,])
}

