#calculate DALYS from the mortality results

library(dplyr)
library(tidyr)

#setwd("C:/Users/rauner/Documents/PIK/fasstr_sebastian")

wd_dir <- "C:/Users/rauner/Documents/PIK/fasstr_sebastian"
setwd(wd_dir)

dir_root   <- paste0(getwd(),   '/mnt/FASST_IDL_PACK/')            # main folder
dir_ancil  <- paste0(dir_root,  'ANCILLARY')                       # ancillary data
dir_DALYS  <- paste0(dir_ancil, '/DALYS')                          # DALYS
dir_outtab <- paste0(dir_root,  'OUTPUT/TABLES')           # output (tables)
#VSL from OECD (2012) Mortality risk valuation in environment, health and transport policies. OECD Publishing
#divded by the 2017 ration of mortality to DALYS, all cause all ages
u_DALY <- 3600000 / 30.4

#read the premature deaths to DALY ratio
daly_ratio           <- read.csv(file=file.path(dir_DALYS, 'daly_ratio.csv'))

daly_ratio[is.na(daly_ratio)] <- 0

#read the results table
write_out_data_frame <- read.csv(file.path(dir_outtab,paste0('results_',RUNSUFFIX, '.csv')),sep=',')


#calcualte the DALYS and write them in the results table
#alri are <5, the rest >30

#join results and DALYS
results <- inner_join(write_out_data_frame,daly_ratio, by=c("iso.iso.y"="ISO")) %>% mutate(daly_copd = data_df_mort_pm.copd1* copd) %>%
                                                                         mutate(daly_alri = data_df_mort_pm.alri1* alri) %>%
                                                                         mutate(daly_ihd  = data_df_mort_pm.ihd1 * ihd)  %>%
                                                                         mutate(daly_stroke = data_df_mort_pm.stroke1* stroke)%>%
                                                                         mutate(daly_lc = data_df_mort_pm.lc1* lc)%>%
                                                                         mutate(daly_total = daly_copd+daly_alri+
                                                                                             daly_ihd +daly_lc+
                                                                                             daly_stroke)

#calculate the new social cost of DALYS
source(paste0(getwd(),   '/functions/calcMonetization.R')   )
# calculate the DALY coefficients
DALY = u_DALY * calcMonetization()
DALY[is.na(DALY)] <- 0
DALY <- as.data.frame(DALY) %>% dplyr::select('Region', 'Year', 'Value')
names(DALY)[3] <- 'cost_per_daly'

DALY$Year <- as.character(DALY$Year)

DALY$Year <- as.integer(DALY$Year)


#multiply dalys with results dalys
results <- inner_join(results, DALY, by=c("iso.iso.y"="Region", 'year'='Year')) %>% mutate(daly_cost = daly_total * cost_per_daly)
results$'scenario' <- scenario
write.table(results, file.path(dir_outtab,paste0('results_dalys_',RUNSUFFIX, '.csv')),sep=';')


