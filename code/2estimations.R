# Package library ############################################################
library(pacman)
#Basic packages
p_load("tidyverse","here","viridis","patchwork","kableExtra",install=T)

#DiD packages
p_load("did2s","did","didimputation")

#Other packages
p_load("ggthemes","knitr","trackdown",install=T)

# Directories ####

projectfolder<-here::here()
source(here("code","0afunctionSheet.R"))
source(here("code","0bfunctionSheet.R"))
data <- here(projectfolder,"data")
figures <- here(projectfolder,"figures")
scrap <- here(projectfolder,"scrap")
tables <- here(projectfolder,"tables")
   
# Estimate effects ####

#Clusters outside of forest areas (to be dropped)
load(file=here(data,"o2DefBuff1000_5000_60.Rdata"))
drop_clusters <- c(unique(o2DefBuff$cluster[o2DefBuff$fcover2000<0.3]),178,243)

# 1 Rings - Main specifications ####
    
load(here(dataOut,"o1DefRingCtrls1000_1000x10_60.Rdata"))
  o1DefRing <- o1DefRing %>% filter(!(cluster %in% c(drop_clusters))) %>% 
    mutate(defAll_per_mine_ha = (defAll_haCum/size_mine)/area_ring*100,
           defAgr_per_mine_ha = (defAgr_haCum/size_mine)/area_ring*100,
           defSettl_per_mine_ha = (defSettl_haCum/size_mine)/area_ring*100,
           defMining_per_mine_ha = (defMining_haCum/size_mine)/area_ring*100)
  
  # Rings
  r1DefRing <- DynamicEst(
    outcomes = "defAll_share",
    data = o1DefRing,
    estimator = c("bjs")
  )
  save(r1DefRing,file=here(dataOut,"r1DefRing30km_60FC.Rdata"))

# 2 Buffer 5 km - Main specifications ####
  load(file=here(dataOut,"o2DefBuff1000_5000_60.Rdata"))

  o2DefBuff <- o2DefBuff %>% filter(fcover2000>0.3,!(cluster %in% c(178,243)),!is.na(size_mine))
  
  r2DefBuff <- DynamicEst(
    outcomes = c(paste(c("defAll","defAgr","defSettl","defMining"),sep="_","share"),paste(c("defAll","defAgr","defSettl","defMining"),sep="_","haCum")),
    data = o2DefBuff,
    estimator = c("csa","gar")
  )
  save(r2DefBuff,file=here(dataOut,"r2DefBuff1000_5000_60FC.Rdata"))
  
# 5 Sensitivity analysis (Appendix) ####
    
  ## A5.1 Spillover robust estimator ####
    
    load(here(dataOut,"o2DefBuff1000_5000_60.Rdata"))
    o2DefBuff <- o2DefBuff %>% filter(!is.na(size_mine),fcover2000>0.3,!(cluster %in% c(178,243)))
    r3SpillRob <- SpilloverEst(outcomes = c(paste(c("defAll","defAgr","defSettl","defOther"),sep="_","share"),paste(c("defAll","defAgr","defSettl","defOther"),sep="_","haCum")),
                               data = o2DefBuff,
                                firststage=NULL,
                                pre_periods = -5,
                                dist=5000
                                ) 
    
    save(r3SpillRob,file=here(dataOut,"r3SpillRob1000_5000_sp10000.Rdata"))

  ## A5.2 Compare different cluster thresholds ####
    load(here(dataOut,"oADefBuff500_5000_60.Rdata"))
    oADefBuff <- oADefBuff %>% filter(fcover2000 > 0.3,!is.na(size_mine))
          
          rADefBuff500 <- DynamicEst(
            outcomes = c("defAll_share"),
            data = oADefBuff,
            estimator = c("csa")
          )
          
    save(rADefBuff500,file=here(dataOut,"rADefBuff500_5000.Rdata"))
    
    load(here(dataOut,"oADefBuff1500_5000_60.Rdata"))
    oADefBuff <- oADefBuff %>% filter(fcover2000 > 0.3,!is.na(size_mine))
    
    rADefBuff1500 <- DynamicEst(
      outcomes = c("defAll_share"),
      data = oADefBuff,
      estimator = c("csa")
    )
    
    save(rADefBuff1500,file=here(dataOut,"rADefBuff1500_5000.Rdata"))
    
    load(here(dataOut,"oADefBuff3000_5000.Rdata"))
    oADefBuff <- oADefBuff %>% filter(fcover2000 > 0.3,!is.na(size_mine))
    
    rADefBuff3000 <- DynamicEst(
      outcomes = c("defAll_share"),
      data = oADefBuff,
      estimator = c("csa")
    )
    
    save(rADefBuff3000,file=here(dataOut,"rADefBuff3000_5000.Rdata"))
    
    load(here(dataOut,"rADefBuff500_5000.Rdata"))
    rADefBuff500 <- rADefBuff500 %>% filter(estimator == "Callaway & Sant'Anna (2021)", period %in% c(-5:10),outcome=="defAll_share") %>%  mutate(across(where(is.numeric), round, 3)) %>% mutate(att = att*100,se = se*100,CI = paste0("[",as.character(lci),",",as.character(uci),"]")) %>% select(period,att,se) 
    load(here(dataOut,"r2DefBuff1000_5000_60FC.Rdata"))
    rADefBuff1000 <- r2DefBuff %>% filter(estimator == "Callaway & Sant'Anna (2021)", period %in% c(-5:10),outcome=="defAll_share")%>%  mutate(across(where(is.numeric), round, 3)) %>% mutate(att = att*100,se = se*100,CI = paste0("[",as.character(lci),",",as.character(uci),"]")) %>% select(att1000 = att,se1000 = se)
    load(here(dataOut,"rADefBuff1500_5000_60.Rdata"))
    rADefBuff1500 <- rADefBuff1500 %>% filter(estimator == "Callaway & Sant'Anna (2021)", period %in% c(-5:10),outcome=="defAll_share")%>%  mutate(across(where(is.numeric), round, 3)) %>% mutate(att = att*100,se = se*100,CI = paste0("[",as.character(lci),",",as.character(uci),"]")) %>% select(att1500 = att,se1500 = se)
    load(here(dataOut,"rADefBuff3000_5000.Rdata"))
    rADefBuff3000 <- rADefBuff3000 %>% filter(estimator == "Callaway & Sant'Anna (2021)", period %in% c(-5:10),outcome=="defAll_share")%>%  mutate(across(where(is.numeric), round, 3)) %>% mutate(att = att*100,se = se*100,CI = paste0("[",as.character(lci),",",as.character(uci),"]")) %>% select(att3000 = att,se3000 = se)
    
    Circsa <- rADefBuff500 %>% bind_cols(rADefBuff1000,rADefBuff1500,rADefBuff3000) %>%
      kbl(format = "latex",booktabs=T) %>%
      kable_styling(latex_options = c("repeat_header")) %>%
      add_header_above(c("period" = 1, "500m clustering" = 2, "1000m clustering" = 2, "1500m clustering" = 2,"3000m clustering" = 2)) %>%
      save_kable(file=here(tables,"CircCsa5000.tex"))

  ## A5.4 GFC data ####
    load(here(dataOut,"oAGFC1000_5000.Rdata"))
    oAGFC <- oAGFC %>% filter(fcover2000>0.3,!(cluster %in% c(178,243)))
    
    
    rAGFC <- DynamicEst(
      data = oAGFC,
      estimator = c("csa"),
      outcomes = c("defRate_acc","def_ha_acc")
    )
    
    save(rAGFC,file=here(dataOut,"rAGFC1000_5000.Rdata"))
    
    ## A5.6 TMF data ####
    
    load(here(dataOut,"oATMF1000_5000.Rdata"))
    oATMF <- oATMF %>% filter(fcover2000>0.3,!(cluster %in% c(178,243)))
    
    rATMF <- DynamicEst(
      outcomes = c("defRate_acc","def_ha_acc","degrRate_acc","degr_ha_acc"),
      data = oATMF,
      estimator = c("csa")
    )
    
    save(rATMF, file=here(dataOut,"rATMF1000_5000.Rdata"))
    
    
    ## A6.1 Alternative specifications #### 
    load(here(dataOut,"o2DefBuff1000_5000.Rdata"))
    drop_clusters <- c(unique(o2DefBuff[o2DefBuff$fcover2000<0.3,"cluster"]),178,243)
    
    ### A6.1.1 Not-yet-treated + Covariates ####
      
    # Rings
    
      load(here(dataOut,"o1DefRingCovs1000_1000x10.Rdata"))
      o1DefRing <- o1DefRing %>% filter(!(cluster %in% c(drop_clusters))) %>% mutate(defAll_per_mine_ha = (defAll_haCum/size_fmine)/area_ring*100,defAgr_per_mine_ha = (defAgr_haCum/size_fmine)/area_ring*100,defSettl_per_mine_ha = (defSettl_haCum/size_fmine)/area_ring*100,defOther_per_mine_ha = (defOther_haCum/size_fmine)/area_ring*100)
      
      r1DefRingCov <- DynamicEst(
        outcomes = c(paste(c("defAll","defAgr","defSettl"),sep="_","share")),
        data = o1DefRing,
        xformla = ~ fcover2000 + access + roads_other_dist + altitude + slope + river_dist + loggroads_old_dist + border_dist + temperature + annual_precipitation,
        estimator = c("csa")
      )
      
    #Circle
      
      load(file=here(dataOut,"o2DefBuff1000_5000_60.Rdata"))
      o2DefBuff <- o2DefBuff %>% filter(fcover2000>0.3,!(cluster %in% c(178,243)),!is.na(size_mine))
      
      r2DefBuffCov <- DynamicEst(
        outcomes = c(paste(c("defAll","defAgr","defSettl","defMining"),sep="_","share")),
        data = o2DefBuff,
        xformla = ~ fcover2000 + access + roads_other_dist + altitude + slope + river_dist + loggroads_old_dist + border_dist + temperature + annual_precipitation + s3aPAs,
        estimator = c("csa")
      )
      
      save(r2DefBuffCov,file=here(dataOut,"r2DefBuffCov1000_5000_60FC.Rdata"))
    
    ### A6.1.2 Random controls + matching ####

      # Never treated Rings
      
      load(here(dataOut,"o1DefRingCtrls1000_1000x10_60.Rdata"))
      
      o1DefRing <- o1DefRing %>% mutate(treat = if_else(treat<0,0,treat),time = as.numeric(time)) %>%
        rename(altitude = altmed30s) %>%
        filter(!(cluster %in% c(drop_clusters)))
      
      r1DefRingNeverTreat <- DynamicEst(
        outcomes = c(paste(c("defAll"),sep="_","share")),
        data = o1DefRing,
        control = "nevertreated",
        xformla = ~ fcover2000 + access + roads_other_dist + altitude + slope + river_dist + loggroads_old_dist + border_dist + temperature + annual_precipitation,
        estimator = c("csa")
      )
      
      # Nevertreated Buffers
      
      load(file=here(dataOut,"o2DefBuff1000_5000_60.Rdata"))
      o2DefBuff <- o2DefBuff %>% filter(fcover2000>0.3,!(cluster %in% c(178,243)))
      
      r2DefBuffNeverTreat <- DynamicEst(
        outcomes = c(paste(c("defAll","defAgr","defSettl","defMining"),sep="_","share")),
        data = o2DefBuff,
        xformla = ~ fcover2000 + access + roads_other_dist + altitude + slope + river_dist  + border_dist + s3aPAs + temperature + annual_precipitation,
        estimator = c("csa")
      )
      
      save(r2DefBuffNeverTreat,file=here(dataOut,"r2DefBuffNeverTreat1000_5000_60FC.Rdata"))
    