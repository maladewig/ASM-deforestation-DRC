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
dataPrep <- here(projectfolder,"dataPrep")
dataInt <- here(projectfolder,"dataInt")
dataOut <- here(projectfolder,"dataOut")
figures <- here(projectfolder,"figures")
scrap <- here(projectfolder,"scrap")
tables <- here(projectfolder,"tables")
   
# Estimate effects ####

#Clusters outside of forest areas (to be dropped)
#load(file=here(dataOut,"o2DefBuffCtrls1000_5000.Rdata"))
load(file=here(dataOut,"o2DefBuff1000_5000_60.Rdata"))
drop_clusters <- c(unique(o2DefBuff$cluster[o2DefBuff$fcover2000<0.3]),178,243)

# 1 Rings - Main specifications ####
    
#load(here(dataOut,"o1DefRingCovs1000_1000x20.Rdata"))
load(here(dataOut,"o1DefRingCtrls1000_1000x10_60.Rdata"))
  o1DefRing <- o1DefRing %>% filter(!(cluster %in% c(drop_clusters))) %>% 
    mutate(defAll_per_mine_ha = (defAll_haCum/size_mine)/area_ring*100,
           defAgr_per_mine_ha = (defAgr_haCum/size_mine)/area_ring*100,
           defSettl_per_mine_ha = (defSettl_haCum/size_mine)/area_ring*100,
           defMining_per_mine_ha = (defMining_haCum/size_mine)/area_ring*100)
  
  # Rings
  r1DefRing <- DynamicEst(
    outcomes = "defAll_share",
    #outcomes = c(paste(c("defAll","defAgr","defSettl"),sep="_","share"),paste(c("defAll","defAgr","defSettl"),sep="_","haCum"),paste(c("defAll","defAgr","defSettl","defMining"),sep="_","per_mine_ha")),
    data = o1DefRing,
    estimator = c("bjs")
  )
  save(r1DefRing,file=here(dataOut,"r1DefRing30km_60FC.Rdata"))
  #load(file=here(dataOut,"r1DefRing30km.Rdata"))

# 2 Buffer 5 km - Main specifications ####
  #load(here(dataOut,"o2DefBuff1000_1000-6000.Rdata"))
  load(file=here(dataOut,"o2DefBuff1000_5000_60.Rdata"))
  #save(o2DefBuff,file=here(dataOut,"o2DefBuffCtrls1000_5000_30FC.Rdata"))
  
  o2DefBuff <- o2DefBuff %>% filter(fcover2000>0.3,!(cluster %in% c(178,243)),!is.na(size_mine))
  
  r2DefBuff <- DynamicEst(
    #outcomes = c(paste(c("defAll","defAgr","defSettl","defMining","defOther"),sep="_","share")),
    #outcomes = c("defAll_haCum","defAgr_haCum","defSettl_haCum"),
    outcomes = c(paste(c("defAll","defAgr","defSettl","defMining"),sep="_","share"),paste(c("defAll","defAgr","defSettl","defMining"),sep="_","haCum")),
    data = o2DefBuff,
    #xformla = ~ access + roads_other_dist + altitude + slope + river_dist + loggroads_old_dist + border_dist + growingperiod_length + fcover2000,
    estimator = c("csa","gar")
  )
  save(r2DefBuff,file=here(dataOut,"r2DefBuff1000_5000_60FC.Rdata"))
  
# 3 Heterogeneity analysis ####
  
  load(file=here(dataOut,"o2DefBuff1000_5000_60.Rdata"))
  o2DefBuff <- o2DefBuff %>% filter(fcover2000>0.3,!(cluster %in% c(178,243)),!is.na(size_mine))
  
  outcome_hetero <- "defAll_share"
  #outcome_hetero <- "defSettl_share"
  outcome_hetero <- "defAgr_share"
  
  r3Hetero <- ATThetero(outcomes = c("defAll_share","defSettl_share","defAgr_share"), data=o2DefBuff)
  r3Hetero <- r3Hetero %>% 
    arrange(cluster, time) %>%
    group_by(cluster) %>% 
    mutate(year_of_mining = treat+2000,
           access=access/60,
           PA_dist = s3aPAs/1000,
           road_prim_dist = road_prim_dist/1000,
           roads_other_dist = roads_other_dist/1000,
           loggroads_old_dist = loggroads_old_dist/1000,
           confAft5 = max(ifelse(rel_treat %in% c(0:5),cumsum(conflict),0)), 
           confBef5 = max(ifelse(rel_treat %in% c(15:-1),cumsum(conflict),0))) %>%
    filter(!is.na(adj_out))
  
  plots <- list()
  all_vars <- c("PA_dist","access","agrSuitability","loggroads_old_open_dist","river_dist","road_prim_dist","roads_other_dist","year_of_mining","fcover2000","confAft5","popdens2000","city_town_village")
  #all_vars <- c("confBef5")
  for (covariate in all_vars) {
    for (t in c(5,-1)) {
      ResidOut <- r3Hetero %>% filter(rel_treat == t)
      xlab_title <- c("Distance to protected areas (km)","Travel time to city (h)","Agricultural Suitability score","Distance to logging road (km)","Distance to river (m)", "Distance to primary road (km)","Distance to road (km)","Year of mining start","Share of forest cover in 2000","Conflicts since mining start","Population density","city_town_village_dist")[which(all_vars==covariate)]
      sp = 0.75
      step = 1
      #seq.sl <- seq(min(cov),max(cov),(max(cov)-min(cov))/100)
      drop_outliers = T
      if (t==5) {
        if (drop_outliers ==T) {
          quartiles <- quantile(as.matrix(ResidOut[,covariate]), probs=c(.25, .75), na.rm = T)
          IQR <- IQR(as.matrix(ResidOut[,covariate]))
          Lower <- quartiles[1] - 1.5*IQR
          Upper <- quartiles[2] + 1.5*IQR 
          ResidOut <- subset(ResidOut, ResidOut[,covariate] > Lower & ResidOut[,covariate] < Upper)
        } 
        cov <- as.matrix(ResidOut[,covariate])
        seq.sl <- seq(min(cov), max(cov), length.out=30)
      }
      # Drop outliers
      
      cov <- as.matrix(ResidOut[,covariate])
      low3 <- loess(ResidOut$adj_out~cov, span=sp)
      low3p <- predict(low3, seq.sl, se=TRUE)
      cboots <- lo.boot(y="adj_out",x=covariate,seq=seq.sl,data=as.data.frame(ResidOut),cil=0.05, cih=0.95)
      
      cild <- cboots$cil
      cihd <- cboots$cih
      #cild <- low3p$fit - low3p$se.fit * 1.96 
      #cihd <- low3p$fit + low3p$se.fit * 1.96 
      
      fit3 <- low3p$fit
      g <- data.frame(cild, cihd, fit3,seq.sl,rel_treat=as.factor(t))
      
      if (t==5) {
        plots[[length(plots)+1]] <- ggplot(g,aes(x=seq.sl,color = rel_treat,fill=rel_treat)) +
          geom_ribbon(aes(ymin=cild, ymax=cihd,fill = rel_treat),alpha = 0.4,color = NA) +
          geom_hline(yintercept=0, linetype="dashed", color = "black") +
          theme_bw() +
          geom_line(aes(x=seq.sl, y=fit3,fill = rel_treat,color=rel_treat)) +
          scale_color_manual(name="Years to mining",values = c("#221c2c","#1c351e","#6b4716","#6a5e18","#7d6756")) +
          scale_fill_manual(name="Years to mining",values = c("#766297","#64ad6b","#c7842a","#d8c656","#7d6756")) +
          xlab(xlab_title) +
          ylab("Share of forest lost")
        
      } else {
        plots[[length(plots)]] <- plots[[length(plots)]] + 
          geom_ribbon(data = g,aes(ymin=cild, ymax=cihd,fill = rel_treat),alpha = 0.4,color = NA) +
          geom_line(data = g,aes(x=seq.sl, y=fit3,fill = rel_treat,color=rel_treat)) 
      }
    }
  }
  
  plots[[6]] + plots[[7]] + plot_layout(guides = "collect",ncol = 2,widths = c(0.5,0.5)) + plot_annotation(tag_levels = "a") & theme(legend.position='bottom')
  plots[[2]]  + plots[[3]] + plots[[4]]  + plots[[5]] +plots[[10]] + plots[[8]] + plot_layout(guides = "collect",ncol = 2,widths = c(0.5,0.5)) + plot_annotation(tag_levels = "a") & theme(legend.position='bottom')
  plots[[1]]  + plots[[2]] + plots[[6]]  + plots[[7]] +plots[[9]] + plots[[11]] + plot_layout(guides = "collect",ncol = 2,widths = c(0.5,0.5)) + plot_annotation(tag_levels = "a") & theme(legend.position='bottom')

  #NEW HETEROGENEITY: ---------------------------------------------------
  load(file=here(dataOut,"o2DefBuffCtrls1000_5000.Rdata"))
  o2DefBuff <- o2DefBuff %>% filter(fcover2000>0.3,!(cluster %in% c(178,243)),!is.na(size_mine))
  
  
  o2DefBuff <- o2DefBuff %>% 
    arrange(cluster, time) %>%
    group_by(cluster) %>% 
    mutate(year_of_mining = treat+2000,
           access=access/60,
           PA_dist = s3aPAs/1000,
           road_prim_dist = road_prim_dist/1000,
           roads_other_dist = roads_other_dist/1000,
           loggroads_old_dist = loggroads_old_dist/1000,
           confAft5 = max(ifelse(rel_treat %in% c(0:5),cumsum(conflict),0)), 
           confBef5 = max(ifelse(rel_treat %in% c(15:-1),cumsum(conflict),0))) %>%
    ungroup()
  
  plots <- list()
  
  all_vars <- c("access","agrSuitability","loggroads_dist","river_dist","road_prim_dist","roads_other_dist","year_of_mining","fcover2000","confAft5","city_town_village_dist","city_town_village_road")
  #all_vars <- c("road_prim_dist")
  r3Hetero <- ATThetero(outcomes = c("defAll_share","defSettl_share","defAgr_share"), data=o2DefBuff)
  clusters <- intersect(r3Hetero$cluster[r3Hetero$rel_treat == -1 & !is.na(r3Hetero$adj_out_defAll_share)],r3Hetero$cluster[r3Hetero$rel_treat == 5 & !is.na(r3Hetero$adj_out_defAll_share)])
  pre_out <- r3Hetero %>% filter(cluster %in% clusters, rel_treat ==-1)
  post_out <- r3Hetero %>% filter(cluster %in% clusters, rel_treat ==5)
  ResultsDf <- post_out %>% mutate(diff_All = post_out$adj_out_defAll_share - pre_out$adj_out_defAll_share,
                                   diff_Agr = post_out$adj_out_defAgr_share - pre_out$adj_out_defAgr_share,
                                   diff_Settl = post_out$adj_out_defSettl_share - pre_out$adj_out_defSettl_share)
  
  for (covariate in all_vars) {
    print(covariate)
    #if (covariate == "road_prim_dist") {
    #  ResultsDf <- ResultsDf %>% filter(road_prim_dist <=79)
    #}
    xlab_title <- c("Travel time to city (h)","Agricultural Suitability score","Distance to logging road (km)","Distance to river (m)", "Distance to primary road (km)","Distance to road (km)","Year of mining start","Share of forest cover in 2000","Conflicts since mining start","Distance to nearest village","city_town_village_main_dist","city_town_village_other_dist")[which(all_vars==covariate)]
    sp = 0.75
    step = 1
    #seq.sl <- seq(min(cov),max(cov),(max(cov)-min(cov))/100)
    drop_outliers = T
    if (drop_outliers ==T) {
      quartiles <- quantile(as.matrix(ResultsDf[,covariate]), probs=c(.25, .75), na.rm = T)
      IQR <- IQR(as.matrix(ResultsDf[,covariate]))
      Lower <- quartiles[1] - 1.5*IQR
      Upper <- quartiles[2] + 1.5*IQR 
      ResidOut <- subset(ResultsDf, ResultsDf[,covariate] > Lower & ResultsDf[,covariate] < Upper)
    } else {
      ResidOut <- ResultsDf
    }
    cov <- as.matrix(ResidOut[,covariate])
    seq.sl <- seq(min(cov), max(cov), length.out=30)
    cov <- as.matrix(ResidOut[,covariate])
    
    # All
    low3 <- loess(ResidOut$diff_All~cov, span=sp)
    low3p <- predict(low3, seq.sl, se=TRUE)
    cboots <- lo.boot(y="diff_All",x=covariate,seq=seq.sl,data=as.data.frame(ResidOut),cil=0.05, cih=0.95)
    
    cild <- cboots$cil
    cihd <- cboots$cih
    fit3 <- low3p$fit
    g_all <- data.frame(cild, cihd, fit3,seq.sl,Driver ="All")
    
    # Agriculture
    low3 <- loess(ResidOut$diff_Agr~cov, span=sp)
    low3p <- predict(low3, seq.sl, se=TRUE)
    cboots <- lo.boot(y="diff_Agr",x=covariate,seq=seq.sl,data=as.data.frame(ResidOut),cil=0.05, cih=0.95)
    
    cild <- cboots$cil
    cihd <- cboots$cih
    fit3 <- low3p$fit
    g_agr <- data.frame(cild, cihd, fit3,seq.sl,Driver ="Farming")
    
    # Settlement
    low3 <- loess(ResidOut$diff_Settl~cov, span=sp)
    low3p <- predict(low3, seq.sl, se=TRUE)
    cboots <- lo.boot(y="diff_Settl",x=covariate,seq=seq.sl,data=as.data.frame(ResidOut),cil=0.05, cih=0.95)
    
    cild <- cboots$cil
    cihd <- cboots$cih
    fit3 <- low3p$fit
    g_settl <- data.frame(cild, cihd, fit3,seq.sl,Driver ="Settlements")
    
    g <- rbind(g_all,g_agr,g_settl)
    
    plots[[length(plots)+1]] <- ggplot(g,aes(x=seq.sl,group = Driver,color=Driver)) +
      geom_ribbon(aes(ymin=cild, ymax=cihd,fill=Driver),alpha = 0.4,color=NA) +
      geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
      geom_line(aes(x=seq.sl, y=fit3,linetype=Driver),linewidth = .8) +
      scale_color_manual(values = c("#221c2c","#1c351e","#6b4716","#6a5e18","#7d6756")) +
      scale_fill_manual(values = c("#766297","#64ad6b","#c7842a","#d8c656","#7d6756")) +
      theme_bw() +
      xlab(xlab_title) +
      ylab("Triggered forest loss (pp)")
  }
  
  plots[[11]]
  
  plots[[1]] + plots[[6]] + plot_layout(guides = "collect",ncol = 2,widths = c(0.5,0.5)) + plot_annotation(tag_levels = "a") & theme(legend.position='bottom')
  plots[[1]]  + plots[[2]] + plots[[3]]  + plots[[4]] +plots[[9]] + plots[[5]] + plot_layout(guides = "collect",ncol = 2,widths = c(0.5,0.5)) + plot_annotation(tag_levels = "a") & theme(legend.position='bottom')
  plots[[7]] + plots[[5]]  + plots[[6]] +plots[[8]] + plots[[10]] + plots[[11]] + plot_layout(guides = "collect",ncol = 2,widths = c(0.5,0.5)) + plot_annotation(tag_levels = "a") & theme(legend.position='bottom')
  
  
  # 4 Table 1 ####    
  # Hectares loss per hectare mine (distance rings)
  agr <- r1DefRing %>% filter(outcome =="defAgr_per_mine_ha",estimator == "Callaway & Sant'Anna (2021)",period == 10,buffers <= 10000) %>% select(buffers,att) %>% mutate(col1 = "Smallholder agriculture",buffers = paste(buffers/1000,"km",sep = " ")) %>% pivot_wider(names_from = buffers,values_from = c(att))
  settl <- r1DefRing %>% filter(outcome =="defSettl_per_mine_ha",estimator == "Callaway & Sant'Anna (2021)",period == 10,buffers <= 10000) %>% select(buffers,att) %>% mutate(col1 = "Settlement expansion",buffers = paste(buffers/1000,"km",sep = " ")) %>% pivot_wider(names_from = buffers,values_from = c(att))
  #other <- r1DefRing %>% filter(outcome =="defOther_per_mine_ha",estimator == "Callaway & Sant'Anna (2021)",period == 10,buffers <= 10000) %>% select(buffers,att) %>% mutate(col1 = "Other land-use",buffers = paste(buffers/1000,"km",sep = " ")) %>% pivot_wider(names_from = buffers,values_from = c(att))
  all <- r1DefRing %>% filter(outcome =="defAll_per_mine_ha",estimator == "Callaway & Sant'Anna (2021)",period == 10,buffers <= 10000) %>% select(buffers,att) %>% mutate(col1 = "All land-use",buffers = paste(buffers/1000,"km",sep = " ")) %>% pivot_wider(names_from = buffers,values_from = c(att))
  
  
  tabl1 <- agr %>% bind_rows(settl,other,all) %>% column_to_rownames(var="col1") %>%
    kableExtra::kbl(digits=2,booktabs = T,caption = "Deforestation for each hectare mine",label = "def_by_mine_ha",format = "latex") %>%
    kableExtra::kable_styling(latex_options = "scale_down") %>%
    kableExtra::save_kable(tabl1,file=here(tables,"def_mine_ha.tex"))
  
# 5 Sensitivity analysis (Appendix) ####
    
  ## A5.1 Spillover robust estimator ####
    
    load(here(dataOut,"o2DefBuff1000_5000_60.Rdata"))
    #load(here(dataOut,"o2aCircGFC1000_6000.Rdata"))
    o2DefBuff <- o2DefBuff %>% filter(!is.na(size_mine),fcover2000>0.3,!(cluster %in% c(178,243)))
    #second_stage = stats::as.formula(glue::glue("~ i(rel_treat, ref = c(-1), keep = (-5:20)) + i({rel_spill}, ref = c(-Inf, -1),keep = (-5:20))"))
    # second_stage = ~i(rel_treat, ref = c(-1))
    # second_stage = ~ treat_ind + ring10000_ind
    # second_stage = ~i(rel_treat, ref = c(-1)) + i(rel_ring10000, ref = c(-Inf, -1))
    
    r3SpillRob <- SpilloverEst(outcomes = c(paste(c("defAll","defAgr","defSettl","defOther"),sep="_","share"),paste(c("defAll","defAgr","defSettl","defOther"),sep="_","haCum")),
                               #outcomes = c("defAll_share"),
                               data = o2DefBuff,
                                firststage=NULL,
                                pre_periods = -5,
                                dist=5000
                                ) #second_stage = ~ treat_ind + ring1000_ind + ring2000_ind + ring3000_ind + ring4000_ind + ring5000_ind + ring6000_ind + ring7000_ind + ring8000_ind + ring9000_ind + ring10000_ind,
    
    save(r3SpillRob,file=here(dataOut,"r3SpillRob1000_5000_sp10000.Rdata"))

  ## A5.2 Compare different cluster thresholds ####
    load(here(dataOut,"oADefBuff500_5000_60.Rdata"))
    oADefBuff <- oADefBuff %>% filter(fcover2000 > 0.3,!is.na(size_mine))
          
          rADefBuff500 <- DynamicEst(
            #outcomes = c("defAll_share","defAll_haCum"),
            outcomes = c("defAll_share"),
            data = oADefBuff,
            estimator = c("csa")
          )
          
    save(rADefBuff500,file=here(dataOut,"rADefBuff500_5000.Rdata"))
    
    load(here(dataOut,"oADefBuff1500_5000_60.Rdata"))
    oADefBuff <- oADefBuff %>% filter(fcover2000 > 0.3,!is.na(size_mine))
    
    rADefBuff1500 <- DynamicEst(
      #outcomes = c("defAll_share","defAll_haCum"),
      outcomes = c("defAll_share"),
      data = oADefBuff,
      estimator = c("csa")
    )
    
    save(rADefBuff1500,file=here(dataOut,"rADefBuff1500_5000.Rdata"))
    
    load(here(dataOut,"oADefBuff3000_5000.Rdata"))
    oADefBuff <- oADefBuff %>% filter(fcover2000 > 0.3,!is.na(size_mine))
    
    rADefBuff3000 <- DynamicEst(
      #outcomes = c("defAll_share","defAll_haCum"),
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
    
  ## A5.3 Subsample analysis ####
    
    # mean(MinDist$minDist)
    # d <- density(MinDist)
    # hist(MinDist,breaks=c(seq(0,100000,by=5000)),xlab="Distance bins (5000m)",ylab="Frequency",main="Distance to nearest neighbouring mine")
    # dplot <- plot(d, main="Kernel Density of minium distance to other mines")
    
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
    
    #Calculate share of triggered agriculture in overall agr forest loss
    
    mod_csa<-did::att_gt(
      yname = "defAgr_haCum",
      tname="time",
      idname="cluster",
      gname = "treat",
      data = o2DefBuff,
      clustervars = "cluster",
      control_group = "notyettreated",
      bstrap=T,
      cband=T
    )
    
    est_agg<-did::aggte(mod_csa, type = "calendar",na.rm=T)
    
    s2dAgrRast <- ifel(s2dAgrRast %in% c(0:2,19,20),NA,s2dAgrRast)
    i1bMinesBuffed <- vect(here(dataInt,"i1bMines1000Buff5000.shp"))
    buff <- project(i1bMinesBuffed,s2dAgrRast)
    t1 <- mask(s2dAgrRast,buff)
    expanse(t1,unit="ha")
    
    ## A6.1 Alternative specifications #### 
    load(here(dataOut,"o2DefBuff1000_5000.Rdata"))
    drop_clusters <- c(unique(o2DefBuff[o2DefBuff$fcover2000<0.3,"cluster"]),178,243)
    
    ### A6.1.1 Not-yet-treated + Covariates ####
      
    # Rings
    
      load(here(dataOut,"o1DefRingCovs1000_1000x10.Rdata"))
      o1DefRing <- o1DefRing %>% filter(!(cluster %in% c(drop_clusters))) %>% mutate(defAll_per_mine_ha = (defAll_haCum/size_fmine)/area_ring*100,defAgr_per_mine_ha = (defAgr_haCum/size_fmine)/area_ring*100,defSettl_per_mine_ha = (defSettl_haCum/size_fmine)/area_ring*100,defOther_per_mine_ha = (defOther_haCum/size_fmine)/area_ring*100)
      
      r1DefRingCov <- DynamicEst(
        outcomes = c(paste(c("defAll","defAgr","defSettl"),sep="_","share")),
        #outcomes = c(paste(c("defAll","defAgr","defSettl","defOther"),sep="_","share"),paste(c("defAll","defAgr","defSettl","defOther"),sep="_","haCum"),paste(c("defAll","defAgr","defSettl","defOther"),sep="_","per_mine_ha")),
        data = o1DefRing,
        xformla = ~ fcover2000 + access + roads_other_dist + altitude + slope + river_dist + loggroads_old_dist + border_dist + temperature + annual_precipitation,
        estimator = c("csa")
      )
      
    #Circle
      
      load(file=here(dataOut,"o2DefBuff1000_5000_60.Rdata"))
      o2DefBuff <- o2DefBuff %>% filter(fcover2000>0.3,!(cluster %in% c(178,243)),!is.na(size_mine))
      
      r2DefBuffCov <- DynamicEst(
        outcomes = c(paste(c("defAll","defAgr","defSettl","defMining"),sep="_","share")),
        #outcomes = c("defAll_share","defAll_haCum"),
        #outcomes = c(paste(c("defAll","defAgr","defSettl","defOther","defMining"),sep="_","share"),paste(c("defAll","defAgr","defSettl","defOther","defMining"),sep="_","haCum")),
        data = o2DefBuff,
        #pretrend.base ="universal",
        #xformla = ~ fcover2000 + access + roads_other_dist + altitude + slope + river_dist + loggroads_old_dist + border_dist + growingperiod_length,
        xformla = ~ fcover2000 + access + roads_other_dist + altitude + slope + river_dist + loggroads_old_dist + border_dist + temperature + annual_precipitation + s3aPAs,
        estimator = c("csa")
      )
      
      save(r2DefBuffCov,file=here(dataOut,"r2DefBuffCov1000_5000_60FC.Rdata"))
    
    ### A6.1.2 Random controls + matching ####

      # Never treated Rings
      
      load(here(dataOut,"o1DefRingCtrls1000_1000x10.Rdata"))
      
      o1DefRing <- o1DefRing %>% mutate(treat = if_else(treat<0,0,treat),time = as.numeric(time)) %>%
        #mutate(treat = if_else(treat>0,treat+2000,treat),time = time+2000) %>%
        rename(altitude = altmed30s) %>%
        filter(!(cluster %in% c(drop_clusters)))
      
      r1DefRingNeverTreat <- DynamicEst(
        outcomes = c(paste(c("defAll"),sep="_","share")),
        #outcomes = c(paste(c("defAll","defAgr","defSettl","defOther"),sep="_","share"),paste(c("defAll","defAgr","defSettl","defOther"),sep="_","haCum"),paste(c("defAll","defAgr","defSettl","defOther"),sep="_","per_mine_ha")),
        data = o1DefRing,
        #pretrend.base = "universal",
        control = "nevertreated",
        xformla = ~ fcover2000 + access + roads_other_dist + altitude + slope + river_dist + loggroads_old_dist + border_dist + temperature + annual_precipitation,
        estimator = c("csa")
      )
      
      # Nevertreated Buffers
      
      load(file=here(dataOut,"o2DefBuff1000_5000_60.Rdata"))
      o2DefBuff <- o2DefBuff %>% filter(fcover2000>0.3,!(cluster %in% c(178,243)))
      
      r2DefBuffNeverTreat <- DynamicEst(
        outcomes = c(paste(c("defAll","defAgr","defSettl","defMining"),sep="_","share")),
        #outcomes = c("defAll_share","defAll_haCum"),
        #outcomes = c(paste(c("defAll","defAgr","defSettl","defOther","defMining"),sep="_","share"),paste(c("defAll","defAgr","defSettl","defOther","defMining"),sep="_","haCum")),
        data = o2DefBuff,
        #control = "nevertreated",
        #pretrend.base ="universal",
        xformla = ~ fcover2000 + access + roads_other_dist + altitude + slope + river_dist  + border_dist + s3aPAs + temperature + annual_precipitation,
        estimator = c("csa")
      )
      
      save(r2DefBuffNeverTreat,file=here(dataOut,"r2DefBuffNeverTreat1000_5000_60FC.Rdata"))
    
    ### A6.1.3 Village controls ####
    
    
    
  ## A7 Mine characteristics ####
      load(file=here(dataOut,"o2DefBuff1000_5000_60.Rdata"))
      #save(o2DefBuff,file=here(dataOut,"o2DefBuffCtrls1000_5000_30FC.Rdata"))
      
      o2DefBuff <- o2DefBuff %>% filter(fcover2000>0.3,!(cluster %in% c(178,243)),!is.na(size_mine))
    t2 <- plotBoxplots(MinesDF = o2DefBuff,covariates = c("agrSuitability","access","river_dist","road_prim_dist","roads_other_dist","loggroads_old_dist","fcover2000"))
  
  # Scrap code ----------------------------------------------------------------------
    #load(here(dataInt,"o2RingCtrl1000_1000x10.Rdata"))
    o2RingCtrl <- o2RingCtrl %>% 
      rename(access = s3hAccess, 
      agrSuitability = s3iSuitability,
      slope = slpmed30s,
      temperature=annual_temp,
      growingperiod_length = growing_period_days, 
      river_dist = s3bRivers,
      road_prim_dist = s3c1RoadsPrim,
      roads_other_dist = s3c2RoadsOther,
      city_town_dist = s3dCityTown,
      loggroads_dist = s3eLoggRoads,
      loggroads_old_dist = s3eLoggRoadsOld,
      loggroads_old_open_dist = s3eLoggRoadsOldOpen,
      border_dist = DRCborders,
      altitude = altmed30s)
    
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
    
      

    #o2RingNeverTreat <- o2RingNeverTreat %>% select(intersect(names(o2RingTreat),names(o2RingCtrl))) %>% arrange(names(o2RingTreat))
    # Rings

    data <- o3AllLandUseCirc
    mod_csa<-did::att_gt(
      yname = "def_ha_acc",
      tname="time",
      idname="cluster",
      gname = "treat",
      data = data,
      clustervars = "cluster",
      control_group = "notyettreated",
      bstrap=T,
      #est_method=est_type,
      cband=T
    )
    
    did::aggte(mod_csa, type = "calendar",na.rm=T)
    
    sum(data %>% filter(treat %in% c(2:19),time == 18) %>% select(def_ha_acc))
  
    # Compare datasets
    load(here(dataOut,"r2aCircGFC1000_5.Rdata"))
    r2aCircGFC <- r2aCircGFC %>% filter(estimator == "Callaway & Sant'Anna (2021)", period %in% c(-5:10),outcome=="defRate_acc")%>%  mutate(across(where(is.numeric), round, 3)) %>% mutate(CI = paste0("[",as.character(lci),",",as.character(uci),"]")) %>% select(period,attGFC = att,seGFC = se)
    load(here(dataOut,"r2bCircTMF1000_5.Rdata"))
    r2bCircTMFDef <- r2bCircTMF %>% filter(estimator == "Callaway & Sant'Anna (2021)", period %in% c(-5:10),outcome=="defRate_acc")%>%  mutate(across(where(is.numeric), round, 3)) %>% mutate(CI = paste0("[",as.character(lci),",",as.character(uci),"]")) %>% select(attTMFdef = att,seTMFdef = se)
    r2bCircTMFDegr <- r2bCircTMF %>% filter(estimator == "Callaway & Sant'Anna (2021)", period %in% c(-5:10),outcome=="degrRate_acc")%>%  mutate(across(where(is.numeric), round, 3)) %>% mutate(CI = paste0("[",as.character(lci),",",as.character(uci),"]")) %>% select(attTMFdegr = att,seTMFdegr = se)
    
    CircDatComp <- r2aCircGFC %>% bind_cols(r2bCircTMFDef,r2bCircTMFDegr) %>%
      kbl(format = "latex",booktabs=T) %>%
      kable_styling(latex_options = c("repeat_header")) %>%
      add_header_above(c("period" = 1, "GFC" = 2, "TMF deforestation" = 2, "TMF degradation" = 2)) %>%
      save_kable(file=here(tables,"CircDatComp.tex"))
    
    
    load(here(dataOut,"r2aCircGFC1000_5.Rdata"))
    r2aCircGFC <- r2aCircGFC %>% mutate(type = "All")
    load(file=here(dataOut,"r3MinesOnlyCirc1000_5.Rdata"))
    r3MinesOnlyCirc <- r3MinesOnlyCirc %>% mutate(type = "Mining")
    load(here(dataOut,"o3AgrOnlyCirc1000_5000.Rdata"))
    r3AgrOnlyCirc <- r3AgrOnlyCirc %>% mutate(type = "Farming")
    
   r2aCircGFC %>% bind_rows(r2aCircGFC,r3AgrOnlyCirc,r3MinesOnlyCirc) %>%
      mutate(type = factor(type)) %>%
      filter(outcome == "defRate_acc",estimator == "Callaway & Sant'Anna (2021)",period>=-5,period <= 8) %>%
      ggplot(aes(x=period, y=att,group = type,shape=type,color=type)) + 
      geom_errorbar(aes(ymin=lci, ymax=uci), width=.2,position=position_dodge(width=0.3)) + 
      geom_line(position=position_dodge(width=0.3),linewidth=.8) +
      geom_point(size=2.3,position=position_dodge(width=0.3)) +
      geom_vline(xintercept=-0.5, linetype="dashed", color = "black") +
      geom_hline(yintercept=0, linetype="solid", color = "black") +
      theme_bw() +
      #theme(legend.position="bottom") +
      scale_color_viridis(option = "cividis",discrete = T) +
      #scale_fill_continuous(guide = "colourbar") +
      xlab("Years to mining") +
      ylab("Additional share of forest lost") 
    #labs(color=legend_title) 


   
    
    
    
    load(here(dataOut,"r2bCircTMF500_5.Rdata"))
    load(here(dataOut,"r2aCircGFC500_5.Rdata"))
    load(here(dataOut,"r3AgrOnlyCirc500_5.Rdata"))
    load(file=here(dataOut,"r3MinesOnlyCirc500_5.Rdata"))
    load(here(dataOut,"r3aSpillGFC500.Rdata"))
    load(here(dataOut,"r3bSpillTMF500.Rdata"))

    t2Circ <- r2aCircGFC %>%
      filter(period == 10,outcome %in% c("defRate_acc","def_ha_acc")) %>% 
      arrange(outcome,estimator) %>%
      mutate(across(where(is.numeric), round, 2)) %>%
      mutate(data = "GFC",'CI' = glue::glue("[{lci},{uci}]")) %>%
      dplyr::select(data,outcome,estimator,att,CI) %>%
      #rename_with(~paste(cnames,"GFC",sep="_",recycle0 = TRUE)) 
      rename("ATT (Main)" = att,"CI (Main)" = CI)
      
    t2Circ <- r2bCircTMF %>%
      filter(period == 10,outcome %in% c("defRate_acc","def_ha_acc","degrRate_acc","degr_ha_acc")) %>% 
      arrange(outcome,estimator) %>%
      mutate(across(where(is.numeric), round, 2)) %>%
      mutate(data = "TMF",'CI' = glue::glue("[{lci},{uci}]")) %>%
      dplyr::select(data,outcome,estimator,att,CI) %>%
      rename("ATT (Main)" = att,"CI (Main)" = CI) %>%
      bind_rows(t2Circ)
    
    t3Spill <- r3aSpillGFC %>%
      filter(period == 10,outcome %in% c("defRate_acc","def_ha_acc"),estimator == "Treatment effect") %>% 
      arrange(outcome,estimator) %>%
      mutate(across(where(is.numeric), round, 2)) %>%
      mutate('CI' = glue::glue("[{lci},{uci}]"),estimator = "Gardner (2022)",data="GFC") %>%
      dplyr::select(data,outcome,estimator,att,CI) %>%
      rename("ATT (spill. robust)" = att,"CI (spill. robust)" = CI)
      #right_join(t2Circ,by = join_by(data,Outcome,Estimator))
    
      t2Circ <- r3bSpillTMF %>%
      filter(period == 10,outcome %in% c("defRate_acc","def_ha_acc","degrRate_acc","degr_ha_acc"),estimator == "Treatment effect") %>% 
      arrange(outcome,estimator) %>%
      mutate(across(where(is.numeric), round, 2)) %>%
      mutate('CI' = glue::glue("[{lci},{uci}]"),estimator = "Gardner (2022)",data="TMF") %>%
      dplyr::select(data,outcome,estimator,att,CI) %>%
      rename("ATT (spill. robust)" = att,"CI (spill. robust)" = CI) %>%
      bind_rows(t3Spill) %>%
      right_join(t2Circ,by = join_by(data,outcome,estimator))
    
      t2Circ <- r3AgrOnlyCirc %>%
        filter(period == 10,outcome %in% c("defRate_acc","def_ha_acc")) %>% 
        arrange(outcome,estimator) %>%
        mutate(across(where(is.numeric), round, 2)) %>%
        mutate(data = "GFC",'CI' = glue::glue("[{lci},{uci}]")) %>%
        dplyr::select(data,outcome,estimator,att,CI) %>%
        rename("ATT (sm. farming)" = att,"CI (sm. farming)" = CI) %>%
        right_join(t2Circ,by = join_by(data,outcome,estimator))
      
      t2Circ <- r3MinesOnlyCirc %>%
        filter(period == 10,outcome %in% c("defRate_acc","def_ha_acc")) %>% 
        arrange(outcome,estimator) %>%
        mutate(across(where(is.numeric), round, 2)) %>%
        mutate(data = "GFC",'CI' = glue::glue("[{lci},{uci}]")) %>%
        dplyr::select(data,outcome,estimator,att,CI) %>%
        rename("ATT (mines only)" = att,"CI (mines only)" = CI) %>%
        right_join(t2Circ,by = join_by(data,outcome,estimator)) %>%
        dplyr::select(Data = data,Outcome = outcome,Estimator=estimator,"ATT (Main)","CI (Main)","ATT (spill. robust)","CI (spill. robust)","ATT (sm. farming)","CI (sm. farming)","ATT (mines only)","CI (mines only)") %>%
        arrange(Data,Outcome,Estimator)
      
      save(t2Circ,file=here(tables,"t2Circ1000.Rdata"))
      load(here(tables,"t2Circ1000.Rdata"))
      t2Circ  %>% dplyr::select(-Data) %>%
        kbl(format = "latex") %>%
        kable_styling(latex_options = c("repeat_header")) %>%
        add_header_above(c(" " = 2, "Main" = 2, "Spillover robust" = 2,"Smallholder farming" = 2,"Mining only" = 2)) %>%
        pack_rows("GFC", 1, 4) %>%
        pack_rows("TMF", 5, 12)   %>%
        save_kable(file=here(tables,"t1.tex"))

    
    load(here(dataOut,"r1aRingGFC500_1000x20.Rdata"))
    colnames(r1aRingGFC) <- paste(colnames(r1aRingGFC),"GFCmain")
    load(here(dataOut,"r1bRingTMF500_1000x20.Rdata"))
    r1bDegr <- r1bRingTMF %>% 
      filter(outcome %in% c("degrRate","degrRate_acc","degr_ha","degr_ha_acc")) %>%
      rename(paste(colnames(r1bRingTMF),"TMFmainDegr"))
    r1bDef <- r1RingTMF %>% 
      filter(outcome %in% c("defRate","defRate_acc","def_ha","def_ha_acc")) %>%
      rename(paste(colnames(r1bRingTMF),"TMFmainDef"))
    

    