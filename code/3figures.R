library(pacman)
p_load("tidyverse","here","viridis","patchwork","kableExtra","extrafont",install=T)

# Directories
projectfolder<-here::here()
source(here("code","0afunctionSheet.R"))
source(here("code","0bfunctionSheet.R"))
dataInt <- here(projectfolder,"dataInt")
dataOut <- here(projectfolder,"dataOut")
figures <- here(projectfolder,"figures")
tables <- here(projectfolder,"tables")

# Results DFs
#load(here(dataOut,"r1aRingGFC1000_1000x10.Rdata"))
#load(here(dataOut,"r1bRingTMF1000_1000x10.Rdata"))



# 1) Plot distance ring estimates (Figure 2) ------------------------------------

  load(here(dataOut,"r1DefRing30km_60FC.Rdata"))

  # Redefine color gradients around treatment timing
  precol_magma <- colorRampPalette(c("#C6507EFF", "#0D0887FF"))
  postcol_magma <- colorRampPalette(c("#EF7F4FFF", "#F0F921FF"))
  precol_viridis <- colorRampPalette(c("#39568Cff","#440154FF"))
  postcol_viridis <- colorRampPalette(c("#29AF7FFF","#FDE725FF"))
  
  p1AllLandUse <- plot_ATT_dist_yr(resultsDF = r1DefRing[r1DefRing$buffers<=10000,],
                                   outcomes = c("defAll_share"),
                                   events = c(-5:-1,1,5,8,10),
                                   maxdist = 10,
                                   plot_type="bars")
  p1 <- p1AllLandUse[[1]] + 
    scale_color_gradientn(colors=c(precol_magma(5),postcol_magma(11)),name="Years since mining",breaks = c(-5,0,10)) + 
    ylab("Converted forest (pp)")
  
  p2SpatTempAgr <- plot_ATT_dist_yr(resultsDF = r1DefRing[r1DefRing$buffers<=6000,],
                                    outcomes = c("defAgr_share"),
                                    estimator = c("csa"),
                                    maxdist = 6,
                                    events = c(-5:-1,1,5,8,10),
                                    plot_type="bars")
  
  p2 <- p2SpatTempAgr[[1]]+ 
    scale_color_gradientn(colors=c(precol_magma(5),postcol_magma(11)),name="Years since mining",breaks = c(-5,0,10)) + 
    ylim(-.01*100, 0.126*100) + 
    ylab("Farmland expansion (pp)")
  
  p3SpatTempSettl <- plot_ATT_dist_yr(resultsDF = r1DefRing[r1DefRing$buffers<=4000,],
                                      outcomes = c("defSettl_share"),
                                      maxdist = 4,
                                      estimator = c("csa"),
                                      events = c(-5:-1,1,5,8,10),
                                      plot_type="bars")
  
  p3 <- p3SpatTempSettl[[1]]  + 
    scale_color_gradientn(colors=c(precol_magma(5),postcol_magma(11)),name="Years since mining",breaks = c(-5,0,10)) + 
    ylim(-.01*100, 0.126*100)  + 
    ylab("Settlement expansion (pp)")
  
  p1 / (p2 + p3 + plot_layout(ncol = 2,widths = c(0.6,0.4))) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a")  & theme(legend.position='bottom')
  
  ggsave(filename = here(figures,"fig2.pdf"),width = 180,height = 160, units = "mm")
  
#2) Plot 5-km buffer estimates (Figure 3) ------------------------------------
  
  load(here(dataOut,"r2DefBuff1000_5000_60FC.Rdata"))
  
  r2All <- r2DefBuff %>% 
    filter(outcome == "defAll_share") %>%
    mutate(Driver = "All land uses")
  r2MiningOnly <- r2DefBuff %>% 
    filter(outcome == "defMining_share") %>%
    mutate(Driver = "Mining")
  r2AgrOnly <- r2DefBuff %>% 
    filter(outcome == "defAgr_share") %>%
    mutate(Driver = "Farming")
  r2SettlOnly <- r2DefBuff %>% 
    filter(outcome == "defSettl_share") %>%
    mutate(Driver = "Settlements")
  r2All %>% bind_rows(r2AgrOnly,r2MiningOnly,r2SettlOnly) %>%
    mutate(Driver = factor(Driver,levels = c("All land uses","Farming","Settlements","Mining")),att = att*100,lci = lci*100,uci = uci*100) %>%
    filter(estimator == "Callaway & Sant'Anna (2021)",period>=-5,period <= 10) %>%
    ggplot(aes(x=period, y=att,group = Driver,shape=Driver,color=Driver,linetype = Driver)) + 
    geom_ribbon(aes(ymin=lci, ymax=uci,fill=Driver),alpha=0.18,position=position_dodge(width=0.5),linewidth=0.1) + 
    geom_line(position=position_dodge(width=0.5),linewidth=.3,linetype = "solid") +
    geom_point(size=.8,position=position_dodge(width=0.5)) +
    geom_vline(xintercept=-0.5, linetype="dashed", color = "black",,linewidth = 0.4) +
    geom_hline(yintercept=0, linetype="solid", color = "black",linewidth = 0.4) +
    theme_bw(base_size = 7) +
   # theme(legend.position='bottom') +
    scale_color_manual(values = c("#766297","#64ad6b","#c7842a","#7d6756","#d8c656")) +
    scale_fill_manual(values = c("#766297","#64ad6b","#c7842a","#7d6756","#d8c656")) +
    scale_shape_manual(values = c(15,16,17,8,18)) +
    scale_linetype_manual(values = c("solid", "dashed", "dotdash", "dotted")) +
    xlab("Years since mining") +
    ylab("Additionally converted forest (pp)") +
    guides(color = guide_legend(keywidth = 0.7,keyheight = 0.7))
  
  ggsave(filename = here(figures,"fig3.pdf"),width = 88,height = 65, units = "mm")
  
# 3) Plot heterogeneities (Figure 4) ------------------------------------
  
  load(file=here(dataOut,"o2DefBuff1000_5000_60.Rdata"))
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
  all_vars <- c("access","agrSuitability","river_dist","road_prim_dist","roads_other_dist","year_of_mining","fcover2000","confAft5")
  r3Hetero <- ATThetero(outcomes = c("defAll_share","defSettl_share","defAgr_share"), data=o2DefBuff)
  clusters <- intersect(r3Hetero$cluster[r3Hetero$rel_treat == -1 & !is.na(r3Hetero$adj_out_defAll_share)],r3Hetero$cluster[r3Hetero$rel_treat == 5 & !is.na(r3Hetero$adj_out_defAll_share)])
  pre_out <- r3Hetero %>% filter(cluster %in% clusters, rel_treat ==-1)
  post_out <- r3Hetero %>% filter(cluster %in% clusters, rel_treat ==5)
  ResultsDf <- post_out %>% mutate(diff_All = post_out$adj_out_defAll_share - pre_out$adj_out_defAll_share,
                                   diff_Agr = post_out$adj_out_defAgr_share - pre_out$adj_out_defAgr_share,
                                   diff_Settl = post_out$adj_out_defSettl_share - pre_out$adj_out_defSettl_share)
  
  for (covariate in all_vars) {
    xlab_title <- c("Travel time to city (h)","Agricultural Suitability score","Distance to river (m)", "Distance to primary road (km)","Distance to road (km)","Year of mining start","Share of forest cover in 2000","Conflicts since mining start")[which(all_vars==covariate)]
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
      geom_line(aes(x=seq.sl, y=fit3,linetype=Driver),linewidth = .4) +
      scale_color_manual(values = c("#221c2c","#1c351e","#6b4716","#6a5e18","#7d6756")) +
      scale_fill_manual(values = c("#766297","#64ad6b","#c7842a","#d8c656","#7d6756")) +
      theme_bw(base_size = 7) +
      xlab(xlab_title) +
      ylab("Converted forest (pp)") +
      guides(color = guide_legend(keywidth = 0.7,keyheight = 0.7))
  }
  
  plots[[1]]  + plots[[5]] + plots[[2]]  + plots[[3]] +plots[[7]] + plots[[8]] + plot_layout(guides = "collect",ncol = 2,widths = c(0.5,0.5)) + plot_annotation(tag_levels = "a") & theme(legend.position='bottom')
  
  ggsave(filename = here(figures,"fig4.pdf"),width = 180,height = 180, units = "mm")
  