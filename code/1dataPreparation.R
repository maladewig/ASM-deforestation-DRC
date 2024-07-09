#1 Undisturbed tropical moist forest
#2 Degraded tropical moist forest
#3 deforested land
#4 tropical moist forest regrowth
#5 Permanent and seasonal water
#6 Other land cover

# 0 Load directories ####

#Prepare data for artisanal mining analysis
#Rerun without preprojecting everything to increase accuracy!
#Packages
library(terra)
library(here)
library(dplyr)
library(tidyr)
library(sf)

#Directories
projectfolder<-here::here()
dataPrep <- here(projectfolder,"dataPrep")
dataInt <- here(projectfolder,"dataInt")
dataOut <- here(projectfolder,"dataOut")
figures <- here(projectfolder,"figures")
scrap <- here(projectfolder,"scrap")
setwd(dataPrep)
sourcedata<-("C:/Users/maltela/OneDrive - Norwegian University of Life Sciences/Documents/Data/DRC")

#Standard parameters
source(here("code","0bfunctionSheet.R"))
#targproj<-"+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m no_defs"

# Spatial extent for analysis
s0extent <- vect(here(dataPrep,"s0extent_v2.shp"))
#provinces<-vect(here(sourcedata,"Forest_atlas_provinces","provinces.shp"))
#provinces<-project(provinces,targproj)
#s0extent<-provinces[provinces$nom_prov %in% c("Nord-Kivu","Sud-Kivu","Maniema","Ituri","Haut-Uele")]
#s0extent <- buffer(s0extent,20000)

# Mining polygons
s2aMinesPolyMerged <- vect(here(dataPrep,"s2aMinesPolyMerged1000.shp"))
s0extent<-project(s0extent,s2aMinesPolyMerged)
s2aMinesPolyMerged <- crop(s2aMinesPolyMerged,s0extent)

#system('shutdown -s')

#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~

#1) Clustering of mines ####

s2Mines<-rast(here(dataPrep,"s0rastExtEmpt.tif"))

## 1.1) Single rasters to polygon layer ####
filelist<-list.files(here(sourcedata,"Mining_Robert"))
for (year in c(1:21)) {
  Mines<-rast(filelist[year])
  Mines<-as.polygons(Mines,dissolve=T,values=T)
  Mines<-project(Mines,targproj)
  #Mines<-as.polygons(Mines)
  Mines<-crop(Mines,s0eastDRC)
  if (nrow(Mines)>1) {
    Mines<-Mines[2,]
    Mines<-disagg(Mines)
    Mines$treat<-year
    Mines<-Mines[,"treat"]
    if (exists("s2aMinesPoly")) {
      s2aMinesPoly<-rbind(s2aMinesPoly,Mines)
    } else {
      s2aMinesPoly<-Mines
    }
  }
}
writeVector(s2aMinesPoly,here(dataPrep,"s2aMinesPoly.shp"),overwrite=T)
#s2aMinesPoly <- vect(here(dataPrep,"s2aMinesPoly.shp"))

## 1.2) Cluster mines with density based clustering ####
s2aMinesPolyMerged<-calibrateMinesPoly(MinesPoly=s2aMinesPoly,clusterradius=c(3000))
s2aMinesPolyMerged$cluster <- c(1:nrow(s2aMinesPolyMerged))
writeVector(s2aMinesPolyMerged,here(dataPrep,"s2aMinesPolyMerged3000.shp"),overwrite=T)
s2aMinesPolyMerged <- vect(here(dataPrep,"s2aMinesPolyMerged1000.shp"))
s0extent <- project(s0extent,s2aMinesPolyMerged)
s2aMinesPolyMerged <- crop(s2aMinesPolyMerged,s0extent)

#2) GFC data for forest cover calculation ####

s1GFC <-rast(here(dataPrep,"s1GFC.tif"))
fcover <- s1GFC[[1]]
#fcover <- ifel(s1GFC[[1]] >=60,1,0) # Count pixels with more than 60% forest cover in canopy cover layer
#fcover <- ifel(s1GFC[[1]] >=30,1,0)
#writeRast(fcover,here(dataInt,"s1GFC_fcover60.tif"))
#writeRast(fcover,here(dataInt,"s1GFC_fcover30.tif"))
#fcover <- rast(here(dataInt,"s1GFC_fcover30.tif"))
fcover <- rast(here(dataInt,"s1GFC_fcover60.tif"))

#s1GFC[[1]] <- fcover

# 3) Load covariates ####

covariates_vect = c("s3aPAs","s3bRivers","s3c1RoadsPrim","s3c2RoadsOther","s3eLoggRoadsOld","DRCborders")
covariates_rast = c("s3hAccess","s3iSuitability","slpmed30s","altmed30s","annual_precipitation","annual_temp","growing_period_days","P_PET_ratio","modified_fournier_index","popdens2000")

# 4) Extract panel data sets ####
#load(here(dataOut,"o1aRingGFC1000_1000x20.Rdata"))

# 4.1) Rings ####

# Ring polygons (pre-buffered in GIS software)
i1aMinesBuffed <- vect(here(dataInt,"i1aMines1000Buff1000x30.shp"))

covariates_vect = c("s3aPAs","s3bRivers","s3c1RoadsPrim","s3c2RoadsOther","s3eLoggRoadsOld","DRCborders")
covariates_rast = c("s3hAccess","s3iSuitability","slpmed30s","altmed30s","annual_precipitation","annual_temp")

  ## 4.1.1) All loss ####
s2hAllLandUse <- rast(here(dataPrep,"s2hAllLandUse_60.tif"))
o1DefRing <-PolyToPanel(flossRast = s2hAllLandUse, 
                           fcoverRast = fcover,
                           mineVect = i1aMinesBuffed,
                          #covariates_vect = covariates_vect,
                          #covariates_rast = covariates_rast,
                           circle=F)
  
 # size_fmine <- expanse(s2aMinesPolyMerged,unit="ha")
#  size_fmine <- data.frame(cluster=c(1:length(size_fmine)),size_fmine)
  o1DefRing <- o1DefRing %>% 
    #left_join(size_fmine,by=c("cluster")) %>%
    rename(defAll_share = defRate_acc,defAll_annualRate = defRate,defAll_haCum = def_ha_acc,defAll_haAnnual= def_ha)

  #save(o1DefRing,file = here(scrap,"o1DefRing1.Rdata"))
    
  ## 4.1.2) Small-scale farming ####
  s2dAgrRast <- rast(here(dataPrep,"s2dAgrRast_60.tif"))
  o1DefRing <- AddOutcomeRings(rast = s2dAgrRast, 
                            data=o1DefRing,
                            circle = F,
                            mineVect = i1aMinesBuffed)

  o1DefRing <- o1DefRing %>% 
    rename(defAgr_share = defRate_acc,defAgr_annualRate = defRate,defAgr_haCum = def_ha_acc,defAgr_haAnnual= def_ha) 
  
  #save(o1DefRing,file=here(dataInt,"i1DefRingAllFarm1000x10.Rdata"))
  #load(here(dataInt,"i1DefRingAllFarm1000x10.Rdata"))
  
  ## 4.1.3) Settlement expansion ####
  s2eSettlRast <- rast(here(dataPrep,"s2eSettlRast_60.tif"))
  o1DefRing <- AddOutcomeRings(rast = s2eSettlRast, 
                            data=o1DefRing,
                            circle = F,
                            mineVect = i1aMinesBuffed)
  
  o1DefRing <- o1DefRing %>% 
    rename(defSettl_share = defRate_acc,defSettl_annualRate = defRate,defSettl_haCum = def_ha_acc,defSettl_haAnnual= def_ha)

  #save(o1DefRing,file = here(dataOut,"o1DefRing1000_1000x30.Rdata"))
  
  #save(o1DefRing,file=here(scrap,"i1DefRingAllFarmSettl1000x10.Rdata"))
  
  ## 4.1.4) Mining ####
  s2cMinesRast <- rast(here(dataPrep,"s2cMinesRast_60.tif"))
  o1DefRing <- AddOutcomeRings(rast = s2cMinesRast,
                               data=o1DefRing,
                               circle = F,
                               mineVect = i1aMinesBuffed)
  
  o1DefRing <- o1DefRing %>% 
    rename(defMining_share = defRate_acc,defMining_annualRate = defRate,defMining_haCum = def_ha_acc,defMining_haAnnual= def_ha)
  
  save(o1DefRing,file=here(scrap,"i1DefRingAllFarmSettlMine1000x10.Rdata"))
  save(o1DefRing,file=here(dataOut,"o1DefRingCtrls1000_1000x10_60.Rdata"))
  
  ## 4.1.5) Other ####
  s2fOtherRast <- rast(here(dataPrep,"s2fOtherRast.tif"))
  o1DefRing <- AddOutcomeRings(rast = s2fOtherRast, 
                              data=o1DefRing,
                              circle = F,
                              mineVect = i1aMinesBuffed)
  
  o1DefRing <- o1DefRing %>% 
  rename(defOther_share = defRate_acc,defOther_annualRate = defRate,defOther_haCum = def_ha_acc,defOther_haAnnual= def_ha)
  
  save(o1DefRing,file=here(dataOut,"o1DefRingCtrls1000_1000x10.Rdata"))
  
# 4.2) Circle ####

 i2MinesBuffed <- vect(here(dataInt,"i2CircleFullSample1000_5000.shp"))
#i1bMinesBuffed <- vect(here(dataInt,"i1bMines1000Buff5000.shp"))
#i1bMinesBuffed <- vect(here(dataInt,"iAMines500Buff5000.shp"))
#i1bMinesBuffed <- vect(here(dataInt,"iAMines1500Buff6000.shp"))
#i1bMinesBuffed <- vect(here(dataInt,"iAMines3000Buff6000.shp"))

  covariates_vect = c("s3aPAs","s3bRivers","s3c1RoadsPrim","s3c2RoadsOther","s3eLoggRoadsOld","DRCborders")
  covariates_rast = c("s3hAccess","s3iSuitability","slpmed30s","altmed30s","annual_precipitation","annual_temp","growing_period_days","P_PET_ratio","modified_fournier_index","popdens2000","FFI2000")
  
  ## 4.2.1) All loss ####
  
o2DefBuff <-PolyToPanel(flossRast = s2hAllLandUse, 
                               fcoverRast = fcover,
                               mineVect = i2MinesBuffed,
                               covariates_vect = covariates_vect, 
                               covariates_rast = covariates_rast,
                               )
  
  #covariates <- covariates %>% group_by(cluster) %>% summarise(across(everything(), list(mean)))
  #covariates <- o2DefBuff %>% select(cluster,all_of(c(covariates_rast,covariates_vect)))
  #save(covariates,file=here(dataInt,"covariates5km.Rdata"))
  
  
  o2DefBuff <- o2DefBuff %>% rename(defAll_share = defRate_acc,defAll_annualRate = defRate,defAll_haCum = def_ha_acc,defAll_haAnnual= def_ha)#%>%
   # rename(PA_dist = s3aPAs,river_dist = s3bRivers,roads_prim_dist = s3c1RoadsPrim, roads_other_dist = s3c2RoadsPrim)
  #save(o2DefBuff,file=here(dataOut,"oADefBuff500_5000.Rdata"))
  #save(o2DefBuff,file=here(dataOut,"oADefBuff1500_5000.Rdata"))
  #save(o2DefBuff,file=here(dataOut,"oADefBuff3000_5000.Rdata"))
  
  
  
  ## 4.2.2) Small-scale farming ####
  
  o2DefBuff <- AddOutcomeRings(rast = s2dAgrRast, 
                               data=o2DefBuff,
                               mineVect = i2MinesBuffed)
  
  o2DefBuff <- o2DefBuff %>% 
    rename(defAgr_share = defRate_acc,defAgr_annualRate = defRate,defAgr_haCum = def_ha_acc,defAgr_haAnnual= def_ha) 
  
  ## 4.2.3) Settlement expansion ####
  o2DefBuff <- AddOutcomeRings(rast = s2eSettlRast, 
                               data=o2DefBuff,
                               mineVect = i2MinesBuffed)
  
  o2DefBuff <- o2DefBuff %>% 
    rename(defSettl_share = defRate_acc,defSettl_annualRate = defRate,defSettl_haCum = def_ha_acc,defSettl_haAnnual= def_ha)
  
  ## 4.2.4) Mining ####
  o2DefBuff <- AddOutcomeRings(rast = s2cMinesRast,
                               data=o2DefBuff,
                               mineVect = i2MinesBuffed)
  
  o2DefBuff <- o2DefBuff %>% 
    rename(defMining_share = defRate_acc,defMining_annualRate = defRate,defMining_haCum = def_ha_acc,defMining_haAnnual= def_ha)
  
  ## 4.2.5) Other ####
  o2DefBuff <- AddOutcomeRings(rast = s2fOtherRast, 
                               data=o2DefBuff,
                               mineVect = i2MinesBuffed)
  
  o2DefBuff <- o2DefBuff %>% 
    rename(defOther_share = defRate_acc,defOther_annualRate = defRate,defOther_haCum = def_ha_acc,defOther_haAnnual= def_ha)
  
  o2DefBuff <- o2DefBuff %>% 
    rename(access = s3hAccess, 
           agrSuitability = s3iSuitability,
           slope = slpmed30s,
           temperature=annual_temp,
           growingperiod_length = growing_period_days, 
           river_dist = s3bRivers,
           road_prim_dist = s3c1RoadsPrim,
           roads_other_dist = s3c2RoadsOther,
           #city_town_dist = s3dCityTown,
           #loggroads_dist = s3eLoggRoads,
           loggroads_old_dist = s3eLoggRoadsOld,
           #loggroads_old_open_dist = s3eLoggRoadsOldOpen,
           border_dist = DRCborders,
           altitude = altmed30s)
  
MinDist<- DistHist(MinesPolys=s2aMinesPolyMerged)
o2DefBuff <- merge(o2DefBuff,MinDist,by="cluster",all.x=T)

s3conflict <- vect(here(dataInt,"i3conflict20km_1000.shp"))
cdf <- as.data.frame(s3conflict)
cdf[,"time"] <- as.numeric(strftime(as.Date(cdf[,"event_date"],format ="%d %b %Y"),"%y"))

o2DefBuff <- cdf %>% group_by(cluster,time) %>% 
  summarise(conflict = n()) %>%
  filter(time>0,time<=20) %>%
  select(cluster,time,conflict) %>%
  right_join(o2DefBuff,by=c("cluster","time")) %>%
  mutate(conflict = replace_na(conflict,0)) %>%
  mutate(conflict = if_else(is.na(size_mine),NA,conflict)) %>%
  arrange(cluster,time,distance)

save(o2DefBuff,file=here(dataOut,"o2DefBuff1000_5000_60.Rdata"))

### 4.5) Add spillover exposure variables ####
#load(here(dataOut,"o2DefBuff1000_5000.Rdata"))
s2aMinesPolyMerged <- vect(here(dataPrep,"s2aMinesPolyMerged1000.shp"))

rings = c(10000)
o2DefBuff <- o2DefBuff %>% mutate(treat_ind = if_else((time>=treat) & treat!=0,1,0),rel_treat= if_else(treat!=0,time - treat,NA))
dist_treat <- DistToTreat(buffers = rings,polys = s2aMinesPolyMerged)
o2DefBuff<-merge(o2DefBuff,dist_treat,by=c("cluster"),all.x=T)
spillovers <- paste("ring",rings,sep="")
for (dist in spillovers) {
  o2DefBuff[,paste("rel",dist,sep="_")]<-o2DefBuff[,"time"] - o2DefBuff[,dist]
  o2DefBuff[,paste(dist,"ind",sep="_")] <- ifelse(o2DefBuff[,paste("rel",dist,sep="_")]>=0,1,0)
}
o2DefBuff[,"treat_or_spill"] <- ifelse(rowSums(o2DefBuff[,c("treat_ind",paste(spillovers,"ind",sep="_"))],na.rm=T)>0,1,0)

save(o2DefBuff,file=here(dataOut,"o2DefBuff1000_5000_60.Rdata"))

# APPENDIX STUFF  ----------------------------------------------------------------------------

# A1) Different control units ####
  ## A1.1) Randomly selected controls ####
    ### A1.1.1) Rings ####
    
    # Ring polygons (pre-buffered in GIS software)
    i2aRandomControls <- vect(here(dataInt,"i1cRandomControls1000_1000x10000.shp"))
    #i2aRandomControls <- project(i2aRandomControls,i1aMinesBuffed)
    #i2aRandomControls[,"cluster"] <- 1:nrow(i2aRandomControls)
    #i2aRandomControls$treat <- NA
    #writeVector(i2aRandomControls,here(dataInt,"i1cRandomControls1000_1000x10000.shp"),overwrite=T)
    
      #### A1.1.1.1) All loss ####
    s2hAllLandUse <- rast(here(dataPrep,"s2hAllLandUse.tif"))
    fcover <- rast(here(dataInt,"s1GFC_fcover60.tif"))
    
    AllRingCtrl <-PolyToPanel(flossRast = s2hAllLandUse, 
                              fcoverRast = fcover,
                              mineVect = i2aRandomControls,
                              covariates_vect = c("s3aPAs","s3bRivers","s3c1RoadsPrim","s3c2RoadsOther","s3dCityTownVillage","s3eLoggRoads","s3eLoggRoadsOld","s3eLoggRoadsOldOpen","DRCborders"),
                              covariates_rast = c("s3hAccess","s3iSuitability","slpmed30s","altmed30s","annual_precipitation","annual_temp","growing_period_days","P_PET_ratio","modified_fournier_index","popdens2000"),
                              circle=F)
      
    AllRingCtrl <- AllRingCtrl %>% 
        rename(defAll_share = defRate_acc,defAll_annualRate = defRate,defAll_haCum = def_ha_acc,defAll_haAnnual= def_ha)
      
    save(AllRingCtrl,file=here(dataInt,"i2RandCtrlDefRingAll1000x10.Rdata"))
      
      
      #### A1.1.1.2) Small-scale farming ####
      s2dAgrRast <- rast(here(dataPrep,"s2dAgrRast.tif"))
      o2RingCtrl <- AddOutcomeRings(rast = s2dAgrRast, 
                                   data=AllRingCtrl,
                                   circle = F,
                                   mineVect = i2aRandomControls)
      
      o2RingCtrl <- o2RingCtrl %>% 
        rename(defAgr_share = defRate_acc,defAgr_annualRate = defRate,defAgr_haCum = def_ha_acc,defAgr_haAnnual= def_ha) 
      
    #### A1.1.1.3) Settlement expansion ####
      s2eSettlRast <- rast(here(dataPrep,"s2eSettlRast.tif"))
      o2RingCtrl <- AddOutcomeRings(rast = s2eSettlRast, 
                                   data=o2RingCtrl,
                                   circle = F,
                                   mineVect = i2aRandomControls)
      
      o2RingCtrl <- o2RingCtrl %>% 
        rename(defSettl_share = defRate_acc,defSettl_annualRate = defRate,defSettl_haCum = def_ha_acc,defSettl_haAnnual= def_ha)
      
      #### A1.1.1.4) Other ####
      s2fOtherRast <- rast(here(dataPrep,"s2fOtherRast.tif"))
      o2RingCtrl <- AddOutcomeRings(rast = s2fOtherRast, 
                                   data=o2RingCtrl,
                                   circle = F,
                                   mineVect = i2aRandomControls)
      
      o2RingCtrl <- o2RingCtrl %>% 
        rename(defOther_share = defRate_acc,defOther_annualRate = defRate,defOther_haCum = def_ha_acc,defOther_haAnnual= def_ha)
      
      save(o2RingCtrl,file=here(dataInt,"o2RingCtrl1000_1000x10.Rdata"))

      ### A1.1.2) Circles ####

      i2MinesBuffed <- vect(here(dataInt,"i2CircleFullSample1000_5000.shp"))
      
      covariates_vect = c("s3aPAs","s3bRivers","s3c1RoadsPrim","s3c2RoadsOther","s3eLoggRoadsOld","DRCborders")
      covariates_rast = c("s3hAccess","s3iSuitability","slpmed30s","altmed30s","annual_precipitation","annual_temp","growing_period_days","P_PET_ratio","modified_fournier_index","popdens2000","FFI2000")
      
      ## 4.2.1) All loss ####
      
      o2DefBuff <-PolyToPanel(flossRast = s2hAllLandUse, 
                              fcoverRast = fcover,
                              mineVect = i2MinesBuffed,
                              covariates_vect = covariates_vect, 
                              covariates_rast = covariates_rast,
      )
      
      #covariates <- covariates %>% group_by(cluster) %>% summarise(across(everything(), list(mean)))
      #covariates <- o2DefBuff %>% select(cluster,all_of(c(covariates_rast,covariates_vect)))
      #save(covariates,file=here(dataInt,"covariates5km.Rdata"))
      
      
      o2DefBuff <- o2DefBuff %>% rename(defAll_share = defRate_acc,defAll_annualRate = defRate,defAll_haCum = def_ha_acc,defAll_haAnnual= def_ha)#%>%
      # rename(PA_dist = s3aPAs,river_dist = s3bRivers,roads_prim_dist = s3c1RoadsPrim, roads_other_dist = s3c2RoadsPrim)
      #save(o2DefBuff,file=here(dataOut,"oADefBuff500_5000.Rdata"))
      #save(o2DefBuff,file=here(dataOut,"oADefBuff1500_5000.Rdata"))
      #save(o2DefBuff,file=here(dataOut,"oADefBuff3000_5000.Rdata"))
      
      
      
      ## 4.2.2) Small-scale farming ####
      
      o2DefBuff <- AddOutcomeRings(rast = s2dAgrRast, 
                                   data=o2DefBuff,
                                   mineVect = i2MinesBuffed)
      
      o2DefBuff <- o2DefBuff %>% 
        rename(defAgr_share = defRate_acc,defAgr_annualRate = defRate,defAgr_haCum = def_ha_acc,defAgr_haAnnual= def_ha) 
      
      ## 4.2.3) Settlement expansion ####
      o2DefBuff <- AddOutcomeRings(rast = s2eSettlRast, 
                                   data=o2DefBuff,
                                   mineVect = i2MinesBuffed)
      
      o2DefBuff <- o2DefBuff %>% 
        rename(defSettl_share = defRate_acc,defSettl_annualRate = defRate,defSettl_haCum = def_ha_acc,defSettl_haAnnual= def_ha)
      
      ## 4.2.4) Mining ####
      o2DefBuff <- AddOutcomeRings(rast = s2cMinesRast,
                                   data=o2DefBuff,
                                   mineVect = i2MinesBuffed)
      
      o2DefBuff <- o2DefBuff %>% 
        rename(defMining_share = defRate_acc,defMining_annualRate = defRate,defMining_haCum = def_ha_acc,defMining_haAnnual= def_ha)
      
      ## 4.2.5) Other ####
      o2DefBuff <- AddOutcomeRings(rast = s2fOtherRast, 
                                   data=o2DefBuff,
                                   mineVect = i2MinesBuffed)
      
      o2DefBuff <- o2DefBuff %>% 
        rename(defOther_share = defRate_acc,defOther_annualRate = defRate,defOther_haCum = def_ha_acc,defOther_haAnnual= def_ha)
      
      o2DefBuff <- o2DefBuff %>% 
        rename(access = s3hAccess, 
               agrSuitability = s3iSuitability,
               slope = slpmed30s,
               temperature=annual_temp,
               growingperiod_length = growing_period_days, 
               river_dist = s3bRivers,
               road_prim_dist = s3c1RoadsPrim,
               roads_other_dist = s3c2RoadsOther,
               #city_town_dist = s3dCityTown,
               #loggroads_dist = s3eLoggRoads,
               loggroads_old_dist = s3eLoggRoadsOld,
               #loggroads_old_open_dist = s3eLoggRoadsOldOpen,
               border_dist = DRCborders,
               altitude = altmed30s)
      
      MinDist<- DistHist(MinesPolys=s2aMinesPolyMerged)
      o2DefBuff <- merge(o2DefBuff,MinDist,by="cluster",all.x=T)
      
      s3conflict <- vect(here(dataInt,"i3conflict20km_1000.shp"))
      cdf <- as.data.frame(s3conflict)
      cdf[,"time"] <- as.numeric(strftime(as.Date(cdf[,"event_date"],format ="%d %b %Y"),"%y"))
      
      o2DefBuff <- cdf %>% group_by(cluster,time) %>% 
        summarise(conflict = n()) %>%
        filter(time>0,time<=20) %>%
        select(cluster,time,conflict) %>%
        right_join(o2DefBuff,by=c("cluster","time")) %>%
        mutate(conflict = replace_na(conflict,0)) %>%
        mutate(conflict = if_else(is.na(size_mine),NA,conflict)) %>%
        arrange(cluster,time,distance)
      
      save(o2DefBuff,file=here(dataOut,"oADefBuff1000_5000_60.Rdata"))
      
  ## A2) Mining hubs ####
      iAHubs <- vect(here(dataInt,'iAHubs5000.shp'))
      iAHubs$cluster <- iAHubs$n
      iAHubs$distance <- 5000
      VillRoad <- iAHubs %>% group_by(cluster) %>% filter(time==1) %>% select(cluster,distance)
      
      ## A2.1) All loss ####
      
      oADefVillBuff <-PolyToPanel(flossRast = s2hAllLandUse, 
                              fcoverRast = fcover,
                              mineVect = VillRoad,
                              #covariates_vect = covariates_vect, 
                              #covariates_rast = covariates_rast
      )
      
      #covariates <- covariates %>% group_by(cluster) %>% summarise(across(everything(), list(mean)))
      #covariates <- o2DefBuff %>% select(cluster,all_of(c(covariates_rast,covariates_vect)))
      #save(covariates,file=here(dataInt,"covariates5km.Rdata"))
      
      
      oADefVillBuff <- oADefVillBuff %>% rename(defAll_share = defRate_acc,defAll_annualRate = defRate,defAll_haCum = def_ha_acc,defAll_haAnnual= def_ha)#%>%
      # rename(PA_dist = s3aPAs,river_dist = s3bRivers,roads_prim_dist = s3c1RoadsPrim, roads_other_dist = s3c2RoadsPrim)
      save(oADefVillBuff,file=here(dataOut,"oADefMinHubBuff_5000.Rdata"))
      
      
      
      ## A2.2) Small-scale farming ####
      
      oADefVillBuff <- AddOutcomeRings(rast = s2dAgrRast, 
                                   data=oADefVillBuff,
                                   mineVect = VillRoad)
      
      oADefVillBuff <- oADefVillBuff %>% 
        rename(defAgr_share = defRate_acc,defAgr_annualRate = defRate,defAgr_haCum = def_ha_acc,defAgr_haAnnual= def_ha) 
      
      ## A2.3) Settlement expansion ####
      s2eSettlRast <- rast(here(dataPrep,"s2eSettlRast.tif"))
      oADefVillBuff <- AddOutcomeRings(rast = s2eSettlRast, 
                                   data=oADefVillBuff,
                                   mineVect = VillRoad)
      
      oADefVillBuff <- oADefVillBuff %>% 
        rename(defSettl_share = defRate_acc,defSettl_annualRate = defRate,defSettl_haCum = def_ha_acc,defSettl_haAnnual= def_ha)
      
      ## A2.4) Mining ####
      s2cMinesRast <- rast(here(dataPrep,"s2cMinesRast.tif"))
      oADefVillBuff <- AddOutcomeRings(rast = s2cMinesRast,
                                   data=oADefVillBuff,
                                   mineVect = VillRoad)
      
      oADefVillBuff <- oADefVillBuff %>% 
        rename(defMining_share = defRate_acc,defMining_annualRate = defRate,defMining_haCum = def_ha_acc,defMining_haAnnual= def_ha)
      
      ## A2.5) Other ####
      s2fOtherRast <- rast(here(dataPrep,"s2fOtherRast.tif"))
      oADefVillBuff <- AddOutcomeRings(rast = s2fOtherRast, 
                                   data=oADefVillBuff,
                                   mineVect = VillRoad)
      
      oADefVillBuff <- oADefVillBuff %>% 
        rename(defOther_share = defRate_acc,defOther_annualRate = defRate,defOther_haCum = def_ha_acc,defOther_haAnnual= def_ha)
      
      save(oADefVillBuff,file=here(dataOut,"oADefMinHubBuff_5000.Rdata"))
      #load(here(dataOut,"oADefVillBuff_5000.Rdata"))
      
      #o3AllLandUseCirc <- AddOutcomeRings(data = o2aCircGFC,rast=s2hAllLandUse,buffers=6000)
      MinDist<- DistHist(MinesPolys=s2aMinesPolyMerged)
      oADefVillBuff <- merge(oADefVillBuff,MinDist,by="cluster",all.x=T)
      
      s3conflict <- vect(here(dataInt,"i3conflict20km_1000.shp"))
      cdf <- as.data.frame(s3conflict)
      cdf[,"time"] <- as.numeric(strftime(as.Date(cdf[,"event_date"],format ="%d %b %Y"),"%y"))
      
      oADefVillBuff <- cdf %>% group_by(cluster,time) %>% 
        summarise(conflict = n()) %>%
        filter(time>0,time<=20) %>%
        select(cluster,time,conflict) %>%
        right_join(oADefVillBuff,by=c("cluster","time")) %>%
        mutate(conflict = replace_na(conflict,0)) %>%
        arrange(cluster,time,distance)
      
      save(oADefVillBuff,file=here(dataOut,"oADefMinHubBuff_5000.Rdata"))
      
  ## A3) Within 1km of loggin roads ####
      
      # Ring polygons (pre-buffered in GIS software)
      i2aRandomControls <- vect(here(dataInt,"i1cRandomControls1000_1000x10000.shp"))
      #i2aRandomControls <- project(i2aRandomControls,i1aMinesBuffed)
      #i2aRandomControls[,"cluster"] <- 1:nrow(i2aRandomControls)
      #i2aRandomControls$treat <- NA
      #writeVector(i2aRandomControls,here(dataInt,"i1cRandomControls1000_1000x10000.shp"),overwrite=T)
      
      #### A3.1) All loss ####
      s2hAllLandUse <- rast(here(dataPrep,"s2hAllLandUse.tif"))
      fcover <- rast(here(dataInt,"s1GFC_fcover60.tif"))
      
      AllRingCtrl <-PolyToPanel(flossRast = s2hAllLandUse, 
                                fcoverRast = fcover,
                                mineVect = i2aRandomControls,
                                covariates_vect = c("s3aPAs","s3bRivers","s3c1RoadsPrim","s3c2RoadsOther","s3dCityTown","s3eLoggRoads","s3eLoggRoadsOld","s3eLoggRoadsOldOpen","DRCborders"),
                                covariates_rast = c("s3hAccess","s3iSuitability","slpmed30s","altmed30s","annual_precipitation","annual_temp","growing_period_days","P_PET_ratio","modified_fournier_index","popdens2000"),
                                circle=F)
      
      AllRingCtrl <- AllRingCtrl %>% 
        rename(defAll_share = defRate_acc,defAll_annualRate = defRate,defAll_haCum = def_ha_acc,defAll_haAnnual= def_ha)
      
      save(AllRingCtrl,file=here(dataInt,"i2RandCtrlDefRingAll1000x10.Rdata"))
      
      
      #### A3.2) Small-scale farming ####
      s2dAgrRast <- rast(here(dataPrep,"s2dAgrRast.tif"))
      o2RingCtrl <- AddOutcomeRings(rast = s2dAgrRast, 
                                    data=AllRingCtrl,
                                    circle = F,
                                    mineVect = i2aRandomControls)
      
      o2RingCtrl <- o2RingCtrl %>% 
        rename(defAgr_share = defRate_acc,defAgr_annualRate = defRate,defAgr_haCum = def_ha_acc,defAgr_haAnnual= def_ha) 
      
      #### A3.3) Settlement expansion ####
      s2eSettlRast <- rast(here(dataPrep,"s2eSettlRast.tif"))
      o2RingCtrl <- AddOutcomeRings(rast = s2eSettlRast, 
                                    data=o2RingCtrl,
                                    circle = F,
                                    mineVect = i2aRandomControls)
      
      o2RingCtrl <- o2RingCtrl %>% 
        rename(defSettl_share = defRate_acc,defSettl_annualRate = defRate,defSettl_haCum = def_ha_acc,defSettl_haAnnual= def_ha)
      
      #### A3.4) Other ####
      s2fOtherRast <- rast(here(dataPrep,"s2fOtherRast.tif"))
      o2RingCtrl <- AddOutcomeRings(rast = s2fOtherRast, 
                                    data=o2RingCtrl,
                                    circle = F,
                                    mineVect = i2aRandomControls)
      
      o2RingCtrl <- o2RingCtrl %>% 
        rename(defOther_share = defRate_acc,defOther_annualRate = defRate,defOther_haCum = def_ha_acc,defOther_haAnnual= def_ha)
      
      save(o2RingCtrl,file=here(dataInt,"o2RingCtrl1000_1000x10.Rdata"))
      ### A3.5) Circles ####
      i1cRandCtrl <- vect(here(dataInt,"i1cRandomControls1000_5000m.shp"))
      i1cRandCtrl <- i1cRandCtrl %>% tidyterra::mutate(cluster = c((max(i1bMinesBuffed$cluster)+1):(nrow(i1cRandCtrl)+max(i1bMinesBuffed$cluster))),treat = 0,size_mine = NA,distance=5000) %>% select(-id)
      i2MinesBuffed <- rbind(i1bMinesBuffed,i1cRandCtrl)
      writeVector(i2MinesBuffed,here(dataInt,"i2CircleFullSample1000_5000.shp"))

    ## A4) Other forest loss data ####
      ### A4.1) GFC data - 5000m ####
    oAGFC <-PolyToPanel(flossRast = s1GFC[[2]],
                             fcoverRast = fcover,
                             mineVect = i1bMinesBuffed)
    
    save(oAGFC,file=here(dataOut,"oAGFC1000_5000.Rdata"),overwrite=T)
    
    
    ### A4.2) TMF data ####
    s1aTMFdegr <-rast(here(dataPrep,"s1aTMFdegradation.tif"))
    s1cTMFfcover <-rast(here(dataPrep,"s1cTMFfcover2000.tif"))
    
    oATMF <- PolyToPanel(flossRast = s1aTMFdegr,
                         fcoverRast = s1cTMFfcover,
                         mineVect = i1bMinesBuffed,
                         TMF = T)
    
    oATMF <- oATMF %>% 
      rename(degrRate_acc = defRate_acc,degrRate = defRate,degr_ha_acc = def_ha_acc,degr_ha= def_ha) 
    
    s1bTMFdef <-rast(here(dataPrep,"s1bTMFdef.tif"))
    oATMF <- AddOutcomeRings(rast = s1bTMFdef, 
                             data=oATMF,
                             mineVect = i1bMinesBuffed,
                             TMF = T)

    save(oATMF,file=here(dataOut,"oATMF1000_5000.Rdata"),overwrite=T)
    
    #### 3.1.1.ii) Outcome as 5km circle 
    #i) Deforestation
    s1TMF <- c(s1cTMFfcover,s1bTMFdef)
    o2bCircTMFdef <-PolyToPanel(rastStack = s1TMF, mineVect = s2aMinesPolyMerged,covariates_vect = c("s3aPAs","s3bRivers"),circle=T,buffers = c(6000))
    
    #ii) Degradation
    s1TMF <- c(s1cTMFfcover,s1aTMFdegr)
    o2bCircTMFdeg <-PolyToPanel(rastStack = s1TMF, mineVect = s2aMinesPolyMerged,covariates_vect = c("s3aPAs","s3bRivers"),circle=T,buffers = c(6000))
    colnames(o2bCircTMFdegr)[colnames(o2bCircTMFdegr) %in% c("def","def_ha","defRate","def_ha_acc","defRate_acc")] <- c("degr","degr_ha","degrRate","degr_ha_acc","degrRate_acc")
    o2bCircTMF <- base::merge(o2bCircTMFdef,o2bCircTMFdegr[,c("cluster","time","degr","degr_ha","degrRate","degr_ha_acc","degrRate_acc")],by=c("cluster","time"),all.x=T)
    
    o2bCircTMF$treat <- o2bCircTMF$treat + 2000
    o2bCircTMF$treat_ind <- ifelse(o2bCircTMF$time>=o2bCircTMF$treat,1,0)
    o2bCircTMF$rel_treat <- o2bCircTMF$time - o2bCircTMF$treat

    s1TMF <- c(s1cTMFfcover,s1aTMFdegr)
    o2bCircTMFdegr <-PolyToPanel(rastStack = s1TMF, mineVect = s2aMinesPolyMerged,circle=T,buffers = c(5000))
    
    save(o2bCircTMF,file=here(dataOut,"o2bCircTMF500.Rdata"))
    
    rings = c(10000)
    dist_treat <- DistToTreat(buffers = rings,polys = s2aMinesPolyMerged)
    dist_treat$ring10000 <- dist_treat$ring10000 + 2000
    o2bCircTMF<-merge(o2bCircTMF,dist_treat,by=c("cluster"),all.x=T)
    spillovers <- paste("ring",rings,sep="")
    for (dist in spillovers) {
      o2bCircTMF[,paste("rel",dist,sep="_")]<-o2bCircTMF[,"time"] - o2bCircTMF[,dist]
      o2bCircTMF[,paste(dist,"ind",sep="_")] <- ifelse(o2bCircTMF[,paste("rel",dist,sep="_")]>=0,1,0)
    }
    o2bCircTMF[,"treat_or_spill"] <- ifelse(rowSums(o2bCircTMF[,c("treat_ind",paste(spillovers,"ind",sep="_"))],na.rm=T)>0,1,0)
    
    save(o2bCircTMF,file=here(dataOut,"o2bCircTMF1000_6000.Rdata"))
    
    s1aTMFdegr <- rast(here(dataPrep,"s1aTMFdegradation.tif"))
    s1bTMFdef <- rast(here(dataPrep,"s1bTMFdef.tif"))
    s1cTMFfcover <- rast(here(dataPrep,"s1cTMFfcover2000.tif"))
    
## A5) Cluster threshold variations ####
    
    s2hAllLandUse <- rast(here(dataPrep,"s2hAllLandUse_60.tif"))
    s2dAgrRast <- rast(here(dataPrep,"s2dAgrRast_60.tif"))
    s2eSettlRast <- rast(here(dataPrep,"s2eSettlRast_60.tif"))
    s2cMinesRast <- rast(here(dataPrep,"s2cMinesRast_60.tif"))
    fcover <- rast(here(dataInt,"s1GFC_fcover60.tif"))
    
    for (threshold in c(500,1500,3000)) {
      vect(here(dataInt,glue::glue("iAMines{threshold}Buff5000.shp")))
      
      oADefBuff <-PolyToPanel(flossRast = s2hAllLandUse, 
                              fcoverRast = fcover,
                              mineVect = i2MinesBuffed,
      )
      oADefBuff <- oADefBuff %>% rename(defAll_share = defRate_acc,defAll_annualRate = defRate,defAll_haCum = def_ha_acc,defAll_haAnnual= def_ha)#%>%
      save(oADefBuff,file=here(dataOut,glue::glue("oADefBuff{threshold}_5000_60.Rdata")))
    }
    
#system('shutdown -s')

