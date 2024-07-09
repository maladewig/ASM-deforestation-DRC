############################# Data preparation #################################

calibrateMinesPoly<-function(MinesPoly,clusterradius) {
  #MinesPoly:  shapefile with mining polygons
  #clusterradius: Maximum distance between mines for density based clustering
  
  #Cluster mines that are density reachable seperately for each year
  MinesPoly$ID <- c(1:nrow(MinesPoly))
  MinesPoly$cluster <- 0
  MinesDist<-distance(MinesPoly,MinesPoly)
  #MinesDist<-ifelse(MinesDist==0,99999,MinesDist)
  clustercount <- 0
  for (ID in 1:nrow(MinesPoly)) {
    clustercount<-clustercount+1
    neighbours<-which(MinesDist[,ID]<clusterradius) #Find other clusters within defined clusterradius
    if (length(neighbours)>0) {
      clusters <- unique(MinesPoly$cluster[neighbours])
      clusters <- clusters[clusters!=0]
      MinesPoly$cluster[MinesPoly$cluster %in% clusters]<-clustercount
    }
    MinesPoly$cluster[MinesPoly$ID==ID] <- clustercount
    MinesPoly$cluster[MinesPoly$ID %in% neighbours] <- clustercount
    
  }
  
  #Keep only mining pixels within cluster with smallest year of establishment
  MinesPoly$size_mine <- expanse(MinesPoly,unit="ha")
  
  size_ha <- MinesPoly %>% 
    tidyterra::group_by(cluster) %>% 
    tidyterra::summarise(size_mine=sum(size_mine))
  size_ha <- as.data.frame(size_ha)
  
  MinesPoly <- MinesPoly %>% 
    tidyterra::group_by(cluster) %>% 
    tidyterra::filter(Y == min(Y)) %>% 
    tidyterra::summarise(treat = mean(Y))
  
  MinesPoly <- merge(MinesPoly,size_ha,by="cluster")
  MinesPoly$cluster <- c(1:nrow(MinesPoly))
  
  # MinesDF<-as.data.frame(MinesPoly)
  # MinesPolyShape<-convHull(MinesPoly,"cluster")
  # MinesPoly <- merge(MinesPolyShape,MinesDF,by="cluster")
  #MinesPoly$size_hull <- expanse(MinesPoly,unit="ha")
  return(MinesPoly)
}

PolyToPanel<-function(buffers= NULL, covariates_vect = NULL, circle=F,covariates_rast = NULL, fcoverRast = fcover_lay,flossRast,mineVect = i1aMinesBuffed,TMF=F) {
#Function to generate panel data table with buffersize-dependent deforestation values as time variable
  
  # 1a) If circles around mines: Calculate buffers
  if (circle == T) {
    print("Calculate circles around mines")
    for (b in buffers) {
      if (b == buffers[[1]]) {
        bufflayer <- buffer(mineVect,b)
        bufflayer$distance <- b
      } else {
        templayer <- buffer(mineVect,b)
        templayer$distance <- b
        bufflayer <- rbind(bufflayer,templayer)
      }
    }
    
  } else {
    # 1b) If rings instead of circles:
    #Create buffer rings IF NOT given as input vector
    if (length(buffers) != 0) {
      print("Calculate buffer rings around mines")
      print(buffers)
      
      buffercreation<-function(buffersize) {
        minebuff <- buffer(mine,buffersize)
        return(minebuff)
      }
      
      for (m in 1:nrow(mineVect)) {
        mine <- mineVect[m,]
        minebuff <- lapply(buffers,buffercreation)
        for(b in length(buffers):1){
          if(b>1) {
            minebuff[[b]]<-minebuff[[b]]-minebuff[[b-1]]
          }
        }
        if (m == 1) {
          bufflayer<-do.call(rbind,minebuff)
        } else {
          templayer<-do.call(rbind,minebuff)  
          bufflayer<-rbind(bufflayer,templayer)
        }
      }
      
      bufflayer$distance <- rep(buffers,times=nrow(mineVect))
    } else {
      bufflayer <- mineVect
    }
  }
  bufflayer$area_ring <- expanse(bufflayer,unit="ha")
  # minesDef <-data.frame(rep(1:length(unique(bufflayer$cluster)),each=length(unique(bufflayer$distance))*21),rep(1:21,times=length(unique(bufflayer$cluster))*length(unique(bufflayer$distance))),rep(unique(bufflayer$distance),each=21,times=length(unique(bufflayer$cluster))))
  # colnames(minesDef) <- c("cluster","time","distance")
  
  # 2) Calculate statistics within buffers
  aoi<-sf::st_as_sf(bufflayer)
  
  # 2.i) Forest cover
  fcover <- exactextractr::exact_extract(fcoverRast,aoi,"mean",append_cols=T)
  colnames(fcover)[ncol(fcover)]<-c("fcover2000")
  #fcover <- fcover %>% select(cluster,time,fcover2000)
  

  # 2.ii) Forest loss
  lossyear <- exactextractr::exact_extract(flossRast, aoi, append_cols=T,function(value, coverage_fraction) {
    data.frame(value = value,
               frac = coverage_fraction / sum(coverage_fraction)) %>%
      group_by(value) %>%
      summarize(freq = sum(frac), .groups = 'drop') %>%
      pivot_wider(names_from = 'value',
                  names_prefix = '',
                  values_from = 'freq')
  }) %>% 
    mutate(across(!names(bufflayer), replace_na, 0)) %>%
    gather("time","def",!names(bufflayer),convert=T) %>%
    filter(time>0) %>%
    left_join(fcover) %>% 
    arrange(cluster,distance,time) %>%
    group_by(cluster,distance) %>%
    mutate(def_ha = def * area_ring,fcover2000_ha = fcover2000*area_ring) %>%
    mutate(defRate = def_ha/fcover2000_ha,def_ha_acc = cumsum(def_ha)) %>%
    #mutate(defRate_acc = def_ha_acc/fcover2000_ha,treat_ind = ifelse(time>=treat,1,0),rel_treat = time-treat) %>%
    mutate(defRate_acc = def_ha_acc/fcover2000_ha) %>%
    select(!def)
  
  if (TMF == T) {
    lossyear <- lossyear %>% mutate(treat = treat + 2000) %>% mutate(treat_ind = ifelse(time>=treat,1,0),rel_treat = time-treat)
  }
  # 2.iv) Data mask values (no data = 0, land = 1, water = 2)
  #Potentially to be added
  
  # 3) Calculate covariate values for each mine
  aoi_cov<-aoi %>% 
    dplyr::filter(distance==min(distance)) %>% 
    dplyr::select(cluster)
  
  # 3.i) Raster Covariates
  for (rast in covariates_rast) {
    print(rast)
    tempRast<-rast(here(dataPrep,paste(rast,".tif",sep="")))
    if (crs(tempRast) != crs(bufflayer)) {
      tempRast <- project(tempRast,bufflayer)
    }
    tempRastDf<-exactextractr::exact_extract(tempRast,aoi_cov,"mean",append_cols=T)
    colnames(tempRastDf)[ncol(tempRastDf)]<-c(rast)
    lossyear<-merge(lossyear,tempRastDf,by=c("cluster"),all.x=T)
  }
  
  # 3.ii) Vector covariates
  aoi_cov <- vect(aoi_cov)
  aoi_cov <- centroids(aoi_cov)
  for (vect in covariates_vect) {
    print(vect)
    if (vect=="s3aPAs") {
      tempVect<-vect(here(dataPrep,paste(vect,".shp",sep="")))
      tempVect <- project(tempVect,aoi_cov)
      tempVectDf<-as.data.frame(distance(aoi_cov,tempVect))
      colnames(tempVectDf)<-vect
      tempVectDf$cluster<-c(1:nrow(tempVectDf))
      invPAs<-vect(here(dataPrep,"s3aInvPAs.shp"))
      invPAsDf<-as.data.frame(distance(aoi_cov,invPAs))
      tempVectDf[tempVectDf$s3aPAs==0,"s3aPAs"]<-invPAsDf[tempVectDf$s3aPAs==0,1]*-1
      lossyear<-merge(lossyear,tempVectDf,by=c("cluster"),all.x=T)
    } else {
      tempVect<-vect(here(dataPrep,paste(vect, ".shp", sep="")))
      tempVect <- project(tempVect,aoi_cov)
      tempVectDf<-as.data.frame(distance(aoi_cov,tempVect))
      colnames(tempVectDf)<-vect
      tempVectDf$cluster<-c(1:nrow(tempVectDf))
      lossyear<-merge(lossyear,tempVectDf,by="cluster",all.x=T)
    }
  }
  #s3gConflict<-vect("s3gConflict.shp")
  
  return(lossyear)
}

AddOutcomeRings <-function(rast,mineVect = s2aMinesPolyMerged,data,circle = F, buffers = 5000,TMF=F) {
  if (circle == T) {
    print("Calculate circles around mines")
    for (b in buffers) {
      if (b == buffers[[1]]) {
        bufflayer <- buffer(mineVect,b)
        bufflayer$distance <- b
      } else {
        templayer <- buffer(mineVect,b)
        templayer$distance <- b
        bufflayer <- rbind(bufflayer,templayer)
      }
    }
  } else {
      bufflayer <- mineVect # if pre-buffered layer
    }
    
  aoi<-sf::st_as_sf(bufflayer)
ResDf <- exactextractr::exact_extract(rast, aoi, append_cols=T,function(value, coverage_fraction) {
  data.frame(value = value,
             frac = coverage_fraction / sum(coverage_fraction)) %>%
    group_by(value) %>%
    summarize(freq = sum(frac), .groups = 'drop') %>%
    pivot_wider(names_from = 'value',
                names_prefix = '',
                values_from = 'freq')
}) %>%
  mutate(across(!names(bufflayer), replace_na, 0)) %>%
  gather("time","def",!names(bufflayer),convert=T) %>%
  filter(time>0) 

if (TMF == T) {
  ResDf$treat <- ResDf$treat + 2000
}

ResDf <- ResDf %>%
  dplyr::right_join(data) %>%
  arrange(cluster,distance,time) %>%
  group_by(cluster,distance) %>%
  mutate(def_ha = def * area_ring,fcover2000_ha = fcover2000*area_ring) %>%
  mutate(defRate = def_ha/fcover2000_ha,def_ha_acc = cumsum(def_ha)) %>%
  mutate(defRate_acc = def_ha_acc/fcover2000_ha) %>%
  select(!def)

return(ResDf)
}

DistToTreat<-function(polys=s2aMinesPolyMerged,years=c(1:21),buffers) {
  #Function to calculate the distance to the closest treated unit
  #treat_thres is a parameter to indicate the distance from a mining coordinate to not be defined as treated. Can e.g. be the radius used for clustering mines.
  
  mindisttreat<-as.data.frame(polys[,c("cluster")])
  for (max_b in 1:length(buffers)) {
    print(buffers[max_b])
    mindisttreat <- mindisttreat %>%
      add_column(newring = Inf) # value for observations that do not experience spillover exposure
    
    if(max_b==1) {
      min_b <- 0 
    } else {
      min_b <- buffers[max_b-1]
    }
    
    for (y in 20:1) {
      dist<- distance(polys,polys[polys$treat==y,])
      dist<-ifelse(dist==0,Inf,dist)
      dist <- apply(dist, 1, FUN = min)
      with_neighbours <- which(dist<=buffers[max_b])
      mindisttreat[with_neighbours,ncol(mindisttreat)] <- y
    }
    colnames(mindisttreat)[ncol(mindisttreat)]<- paste("ring",buffers[max_b],sep="")
  }
  return(mindisttreat)
}

DistHist <- function(MinesPolys) {
  #Function to plot histogram by minimum distance to other mines
  
  dist<-distance(MinesPolys,MinesPolys)
  dist<-ifelse(dist==0,99999,dist)
  minDistDf <- data.frame()
  clusters <- MinesPolys$cluster
  for (cluster in clusters) {
    row <- which(MinesPolys$cluster == cluster)
    NN <- which.min(dist[row,])
    minDist <- dist[cluster,NN]
    minDistDf <- rbind(minDistDf,data.frame(row,minDist))
  }
  colnames(minDistDf) <- c("cluster","minDist")
  return(minDistDf)
}

plot_ATT_dist <- function(resultsDF_GFC,resultsDF_TMF,outcomes) {
  
  est_plot<-ggplot(ES_Df, aes(x=buffers, y=att, colour=outcome, group=estimator)) + 
  geom_errorbar(aes(ymin=lci, ymax=uci), width=.5,position=position_dodge(width=150)) + 
  geom_line(position=position_dodge(width=150)) +
  geom_point(size=2,position=position_dodge(width=150)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme_minimal() +
  scale_color_viridis(discrete=TRUE) +
  xlab("Distance to mine (m)") +
  ylab("ATT estimate") +
  labs(title = "years after treatment") +
  ggtitle("ATE by distance to mine")
  
  return(est_plot)
}
  
plot_ATT_dist_yr <- function(resultsDF,outcomes,events = c(5,10,15),data_source, estimator ="csa", plot_color="viridis",plot_type="lines",maxdist=12000) {
  plots <- list()
  all_est <- c("csa","gar","bjs","sunab","twfe")
  ref <- c("Callaway & Sant'Anna (2021)","Gardner (2022)","Borusyak et al. (2022)","Sun & Abraham (2021)","TWFE")[which(all_est==estimator)]
 # if (data_source == "TMF") {
    for (out in outcomes) {
      all_outcomes <- c("defRate","defRate_acc","def_ha","def_ha_acc","degrRate","degrRate_acc","degr_ha","degr_ha_acc")
      ylab_title <- c("Additional annual deforestation rate","Additional share of forest lost"," Additional annually deforested hectares","Additionally deforested hectares","Additional annual degradation rate","Additional share of forest degraded","Additional annually degraded hectares","Additionally degraded hectares")[which(all_outcomes==out)]
      
    if (plot_type =="lines")  {
      #if (out %in% c("degrRate","degrRate_acc","degr_ha","degr_ha_acc")) plot_color = "magma"
    p <- resultsDF %>% filter(outcome==out, period %in% events,estimator == ref) %>% mutate(period = as.factor(period)) %>%
      ggplot(aes(x=buffers, y=att, colour=period, group=period)) + 
      geom_errorbar(aes(ymin=lci, ymax=uci), width=.5,position=position_dodge(width=250)) +
      geom_line(position=position_dodge(width=250)) +
      geom_point(size=2,position=position_dodge(width=250)) +
      geom_hline(yintercept=0, linetype="dashed", color = "black") +
      theme_bw() +
      scale_color_viridis(option=plot_color,discrete=TRUE) +
      xlab("Distance from mine (m)") +
      ylab(ylab_title)
    
    plots[[length(plots)+1]] <- p
    }
      
    if (plot_type == "bars") {
      #axlabels=c("[0,1]","[1,2]","2 - 3","3 - 4","4 - 5","5 - 6","6 - 7","7 - 8","8 - 9","9 - 10","10 - 11","11 - 12","12 - 13","13 - 14","14 - 15","15 - 16","16 - 17","17 - 18","18 - 19","19 - 20")
      axlabels=c("[0-1]","[1-2]","[2-3]","[3-4]","[4-5]","[5-6]","[6-7]","[7-8]","[8-9]","[9-10]","[10-11]","[11-12]","[12-13]","[13-14]","[14-15]","[15-16]","[16-17]","[17-18]","[18-19]","[19-20]","[20-21]","[21-22]","[22-23]","[23-24]","[24-25]")
      p <- resultsDF %>% filter(outcome==out, period >= min(events),period <= max(events),estimator == ref) %>%
        group_by(buffers) %>%
        #mutate(att = if_else(period>0,att - lag(att,order_by = buffers),att),buffers = buffers-500) %>%
        mutate(buffers = (buffers-500)/1000,att = att*100,lci = lci*100,uci = uci*100) %>%
        arrange(period) %>%
        #mutate(period=factor(period))%>%
        ggplot(aes(x=buffers, y=att,group=buffers)) + 
        geom_vline(xintercept=seq(from =0, to=maxdist,by=1), linetype="dashed", color = "darkgrey", linewidth = .5) +
        #geom_bar(position="stack",stat = "identity") +
        geom_point(aes(color=period,group=period),size=1,position=position_dodge(width=.9)) +
        geom_errorbar(aes(ymin=lci, ymax=uci,group=period,color=period), width=0.05,position=position_dodge(width=.9)) +
        geom_hline(yintercept=0, linetype="dashed", color = "black") +
        theme_bw() +
        theme(legend.position="bottom",axis.ticks = element_blank()) +
        #theme(legend.position="bottom",axis.ticks = element_blank(),text=element_text(family="Arial")) +
        #scale_fill_viridis(option=plot_color) +
        scale_color_viridis(option="magma",direction=-1,name="Years to mining",breaks = c(min(events),0,max(events))) +
        scale_x_continuous(breaks = seq(from =0.5, to=maxdist-0.5,by=1),labels=axlabels[1:length(seq(from =0.5, to=maxdist-0.5,by=1))]) +
        xlab("Distance from mine (km)") +
        ylab(ylab_title)
      
      plots[[length(plots)+1]] <- p
    }
    }

return(plots)
}

plotEvents <- function(ES_Df,outcomes,pre_min = c(-5), post_max = c(15),estimator="csa",distances=NULL,ylab_title = NULL,plot_color = "viridis",legend_title = NULL) {
  e_plots <- list()
  if (is.null(distances)) distances = unique(ES_Df$buffers)
  if (is.null(legend_title)) legend_title = "Distance in m"
  all_est <- c("csa","gar","bjs","sunab","twfe")
  ref <- c("Callaway & Sant'Anna (2021)","Gardner (2022)","Borusyak et al. (2022)","Sun & Abraham (2021)","TWFE")[which(all_est==estimator)]
  for (out in outcomes) {
    if (is.null(ylab_title)) {
    all_outcomes <- c("defRate","defRate_acc","def_ha","def_ha_acc","degrRate","degrRate_acc","degr_ha","degr_ha_acc")
    ylab_title <- c("Additional deforestation p.a. in %","Additional share of forest lost"," Additionally deforested hectares p.a. in %","ATT (deforested hectares)","ATT (annual degradation rate)","ATT (share of forest degraded)","ATT (annually degraded hectares)","ATT (degraded hectares)")[which(all_outcomes==out)]
    }
  p <- ES_Df %>%
    filter(outcome == out,estimator == ref,buffers %in% distances,period>=pre_min,period <= post_max) %>%
    mutate(buffers = buffers/1000,att = att*100,lci = lci*100, uci = uci*100) %>%
    ggplot(aes(x=period, y=att,color=buffers,group=buffers)) + 
    geom_errorbar(aes(ymin=lci, ymax=uci), width=.2,position=position_dodge(width=0.5)) + 
    geom_line(position=position_dodge(width=0.5)) +
    geom_point(size=1.5,position=position_dodge(width=0.5)) +
    geom_vline(xintercept=-0.5, linetype="dashed", color = "black") +
    geom_hline(yintercept=0, linetype="solid", color = "black") +
    theme_bw() +
    #theme(legend.position="bottom") +
    scale_color_viridis(option = plot_color) +
    #scale_fill_continuous(guide = "colourbar") +
    xlab("Years since mining") +
    ylab(ylab_title) 
    #labs(color=legend_title) 
  
  e_plots[[length(e_plots)+1]] <- p
  }
  return(e_plots)
}

plotSpillRob <- function(ResDf,outcomes,events = c(-5:10)) {
  
  spill_plots <- list()
  
  for (out in outcomes) {
  all_outcomes <- c("defRate","defRate_acc","def_ha","def_ha_acc","degrRate","degrRate_acc","degr_ha","degr_ha_acc")
  ylab_title <- c("ATT (annual deforestation rate)","ATT (share of forest lost)"," ATT (annually deforested hectare)","ATT (deforested hectares)","ATT (annual degradation rate)","ATT (share of forest degraded)","ATT (annually degraded hectares)","ATT (degraded hectares)")[which(all_outcomes==out)]
  p <- ResDf %>% 
    filter(outcome == out, period %in% events) %>%
    #filter(group %in% c("Treatment effect","Spillover effect")) %>%
    ggplot() +
    geom_vline(xintercept=-0.5, linetype="dashed", color = "black") +
    geom_hline(yintercept = 0, color = "black") +
    geom_point(aes(x = period, y = att, color = estimator),position=position_dodge(width=0.6)) +
    geom_errorbar(aes(x = period, ymin = lci, ymax = uci, color = estimator),position=position_dodge(width=0.6)) +
    theme_bw() +
    theme(legend.position="bottom") +
    scale_color_viridis(discrete=TRUE) +
    #scale_x_continuous(minor_breaks = seq(-7, 14, 1)) + 
    #scale_y_continuous(minor_breaks = seq(-75, 45, 50), limits = c(-75, 50)) + 
    labs(y = ylab_title, x = "Years since mining", color = NULL)
  
  spill_plots[[length(spill_plots)+1]] <- p
  }
  return(spill_plots)
}

plotBoxplots <- function(MinesDF,covariates) {
  plots <- list()
  p <- 1
  for (cov in covariates) {
    tempDF <- MinesDF[,c("cluster","treat",cov)]
    colnames(tempDF)[[3]] <- "cov"
    plots[[p]] <- tempDF %>% filter(treat==1) %>%
      group_by(cluster) %>%
      summarise(covariate = mean(cov)) %>%
      ggplot(aes(x=covariate)) +
      geom_boxplot() +
      xlab(cov)
    p <- p+1
  }
  return(plots)
}