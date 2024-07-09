##################### Functions for Model estimation ###########################

ES_estimation<-function(buffers=NULL, outcomes="defRate", data, estimator, xformla=NULL,pre_min = c(-19), post_max = c(19),output_folder) {
  e_plots <- list()
  # estimator:    -twfe    Two way fixed effects estimator
  #               -csa     Callaway & Sant'Anna (2020)
  #               -sunab   Sun & Abraham (2020)  
  #               -dch     de Chaisemartin & d'haultfoueille (2020) SLOW!!!
  #               -gar     Gardner (2022)
  
  #  timeperiods<-c(leadslags[1]:leadslags[2])
  #keepvars<-sapply(timeperiods, function(x) paste("rel_year_",x,sep=""))
  if (is.null(buffers)) {
    buffers <- unique(data$distance)
  }
 # data$rel_year <- data$time - data$treat # Time to treatment
  for (out in outcomes) {
  for (dist in buffers) {
    if ("twfe" %in% estimator) {
      #data$tdummy<-1 # In this case all treated
      if (is.null(xformla)) {
        xformla_null = "0 "
      } else {
        xformla_null = paste0("0 + ", as.character(xformla)[[2]])
      }
      twfe_formla = stats::as.formula(glue::glue("{out} ~ 1 + {xformla_null} + i(rel_treat, ref = -1) | time + cluster"))
      mod_twfe<-fixest::feols(twfe_formla, data = data[data$distance==dist,])
      ES_twfe<-broom::tidy(mod_twfe)
      ES_twfe<-ES_twfe %>% separate(term,into=c(NA,NA,"rel_time",NA),":")
      ES_twfe<-ES_twfe[,c(1:3)]
      colnames(ES_twfe)<-c("rel_time","ATT","SE")
      ES_twfe$rel_time<-as.numeric(ES_twfe$rel_time)
      ES_twfe$lci <- ES_twfe$ATT - 1.96*ES_twfe$SE
      ES_twfe$uci <- ES_twfe$ATT + 1.96*ES_twfe$SE
      ES_twfe$estimator <- "TWFE"
      ES_twfe$outcome <- out
      ES_twfe$distance <- dist
      
      if (exists("ES_Df")) {
        ES_Df<-rbind(ES_Df,ES_twfe)
      } else {
        ES_Df <- ES_twfe
      }
    }
    
    if ("sunab" %in% estimator) {
      if (is.null(xformla)) {
        sunab_xformla = "1"
      }
      else {
        sunab_xformla = paste0("1 + ", as.character(xformla)[[2]])
      }
      sunab_formla = stats::as.formula(glue::glue("{out} ~ {sunab_xformla} + sunab(treat, rel_treat, ref.c =0, ref.p = -1) | cluster + time"))
      ES_sunab <- fixest::feols(sunab_formla, data = data[data$distance==dist,])
      ES_sunab<-broom::tidy(ES_sunab)
      ES_sunab <- ES_sunab %>% separate(term, into=(c(NA,"rel_time")),"::")
      ES_sunab<-ES_sunab[,c(1:3)]
      colnames(ES_sunab)<-c("rel_time","ATT","SE")
      ES_sunab$rel_time<-as.numeric(ES_sunab$rel_time)
      ES_sunab$lci <- ES_sunab$ATT - 1.96*ES_sunab$SE
      ES_sunab$uci <- ES_sunab$ATT + 1.96*ES_sunab$SE
      ES_sunab$estimator <- "Sun & Abraham (2021)"
      ES_sunab$outcome <- out
      ES_sunab$distance <- dist
      
      
      if (exists("ES_Df")) {
        ES_Df<-rbind(ES_Df,ES_sunab)
      } else {
        ES_Df <- ES_sunab
      }
    }
    
    if ("csa" %in% estimator) {
      mod_csa<-did::att_gt(
        yname = out,
        tname="time",
        idname="cluster",
        gname = "treat",
        data = data[data$distance==dist,],
        control_group = "notyettreated",
        xformla = xformla,
        clustervars = "cluster",
        bstrap=T,
        cband=T)
      
      ES_csa <- did::aggte(mod_csa, type = "dynamic",na.rm=T)
      ES_csa <- broom::tidy(ES_csa)
      ES_csa <- ES_csa[,c(3:5,8,9)]
      colnames(ES_csa)<-c("rel_time","ATT","SE","lci","uci")
      ES_csa$estimator <- "Callaway & Sant'Anna (2021)"
      ES_csa$outcome <- out
      ES_csa$distance <- dist
      
      
      if (exists("ES_Df")) {
        ES_Df<-rbind(ES_Df,ES_csa)
      } else {
        ES_Df <- ES_csa
      }
    }
    
    if ("dch" %in% estimator) {
      data$tdummy<-ifelse(data$treat <= data$time,1,0)
      mod_dch <- DIDmultiplegt::did_multiplegt(
        df = data[data$distance==dist,],
        Y = outcome,
        G = "cluster",
        T = "time",
        D = "tdummy", 
        controls = c(),
        placebo = 19,
        brep = 20,
        #dynamic = 0,
        dynamic = 19,
        average_effect = NULL)
      ests = mod_dch[grepl("^placebo_|^effect|^dynamic_", names(mod_dch))]
      ses = as.numeric(mod_dch[grepl("^se_placebo|^se_effect|^se_dynamic", names(mod_dch))])
      
      ES_dch = data.frame(
        term      = names(ests),
        ATT       = as.numeric(ests),
        SE        = as.numeric(mod_dch[grepl("^se_placebo|^se_effect|^se_dynamic", names(mod_dch))])
      ) |>
        # For CIs we'll assume standard normal distribution
        within({
          lci  = ATT - SE*(qnorm(1-(1-0.95)/2))
          uci = ATT + SE*(qnorm(1-(1-0.95)/2))
        })
      ES_dch <- ES_dch[3:37,-1]
      ES_dch$rel_time <- c(-17:17)
      ES_dch$estimator <- "de Chaisemartin & D'Haultfoeuille (2020)"
      
      if (exists("ES_Df")) {
        ES_Df<-rbind(ES_Df,ES_dch)
      } else {
        ES_Df <- ES_dch
      }
    }  
    
    if ("gar" %in% estimator) {
      if (is.null(xformla)) {
        first_stage = ~ 0 | cluster + time 
      } else {
        first_stage <- stats::as.formula(glue::glue("~ 0 + {xformla} | cluster + time")[[2]])
      }
      
      data$treat_ind <- 0
      data[data$time>=data$treat,"treat_ind"]<-1
      data$rel_time <- data$time - data$treat
      
      mod_gar<-did2s::did2s(
        data = data[data$distance==dist & data$rel_time>=pre_min,],
        yname = out,
        first_stage = first_stage,
        second_stage = ~ i(rel_time, ref=c(-1, Inf)),
        treatment = "treat_ind",
        cluster_var = "cluster"
      )
      
      ES_gar<-broom::tidy(mod_gar)
      ES_gar$term <- sub("rel_time::", "", ES_gar$term)
      ES_gar <- ES_gar[ES_gar$term %in% c(-20:20),1:3]
      ES_gar$term <- as.numeric(ES_gar$term)
      colnames(ES_gar) <- c("rel_time","ATT","SE")
      ES_gar$lci <- ES_gar$ATT - 1.96*ES_gar$SE
      ES_gar$uci <- ES_gar$ATT + 1.96*ES_gar$SE
      ES_gar <- bind_rows(ES_gar,
                                tibble(rel_time = -1, ATT = 0, lci = 0, uci = 0, SE = 0))
      
      colnames(ES_gar) <- c("rel_time","ATT","SE","lci","uci")
      ES_gar$estimator <- "Gardner (2022)"
      ES_gar$outcome <- out
      ES_gar$distance <- dist
      
      if (exists("ES_Df")) {
        ES_Df<-rbind(ES_Df,ES_gar)
      } else {
        ES_Df <- ES_gar
      }
    }
    
    if ("bjs" %in% estimator) {
      if (is.null(xformla)) {
        bjs_first_stage = ~1 | time + cluster
      } else {
        xformla_null = paste0("0 + ", as.character(xformla)[[2]])
        bjs_first_stage = stats::as.formula(glue::glue("~ 1 + {xformla_null} | time"))
      }
      mod_bjs<-didimputation::did_imputation(
        data = data[data$distance==dist & data$rel_time>=pre_min,],
        yname = out,
        gname = "treat",
        tname = "time",
        idname = "cluster",
        first_stage = bjs_first_stage,
        horizon = T,
        pretrends = c(pre_min:-1)
      )  
      ES_bjs <- mod_bjs[,-1]
      ES_bjs$term <- as.numeric(ES_bjs$term)
      colnames(ES_bjs) <- c("rel_time","ATT","SE","lci","uci")
      ES_bjs$estimator <- "Borusyak et al. (2022)"
      ES_bjs$outcome <- out
      ES_bjs$distance <- dist
      
      if (exists("ES_Df")) {
        ES_Df<-rbind(ES_Df,ES_bjs)
      } else {
        ES_Df <- ES_bjs
      }
    }
    
    ES_Df <- ES_Df[ES_Df$rel_time>=pre_min,]
    ES_Df <- ES_Df[ES_Df$rel_time<=post_max,]
    
    #Plot different event studies  
   
    # plot_title <- paste(estimator,collapse="_")
    # ggsave(filename = here(figures,output_folder,glue::glue("ES_{plot_title}_{out}_{dist}.jpg")), plot = p,width = 6,height = 5)
    }
  }
  return(ES_Df)
}

DynamicEst <- function(data, outcomes, buffers=NULL, xformla=NULL,estimator = "csa",anticipation = 0,pretrend.base="varying",control = "notyettreated") {
  # Function to calculate ATT by distance to mine and event time.
  #Currently only for CSA and BJS estimator
  # max_e:     defines the post-treatment periods for which to aggregate the ATT. 
  # buffers:   buffer distances to be included
  if (is.null(buffers)) buffers <- unique(data$distance)
  ES_Df <- data.frame()
  for (out in outcomes) {
    print(out)
    for (buffer in buffers) {
      if ("gar" %in% estimator) {
        if (is.null(xformla)) {
          first_stage = ~ 0 | cluster + time 
        } else {
          first_stage <- stats::as.formula(glue::glue("~ 0 + {xformla} | cluster + time")[[2]])
        }
        
        data$treat_ind <- 0
        data[data$time>=data$treat,"treat_ind"]<-1
        data$rel_time <- data$time - data$treat
        
        mod_gar<-did2s::did2s(
          data = data[data$distance==buffer & data$rel_time>=-10,],
          yname = out,
          first_stage = first_stage,
          second_stage = ~ i(rel_time, ref=c(-1, Inf)),
          treatment = "treat_ind",
          cluster_var = "cluster"
        )
        
        ES_gar<-broom::tidy(mod_gar)
        ES_gar$term <- sub("rel_time::", "", ES_gar$term)
        ES_gar <- ES_gar[ES_gar$term %in% c(-20:20),1:3]
        ES_gar$term <- as.numeric(ES_gar$term)
        colnames(ES_gar) <- c("period","att","se")
        ES_gar$lci <- ES_gar$att - 1.96*ES_gar$se
        ES_gar$uci <- ES_gar$att + 1.96*ES_gar$se
        ES_gar <- bind_rows(ES_gar,
                            tibble(period = -1, att = 0, se = 0, lci = 0, uci = 0))
        ES_gar$estimator <- "Gardner (2022)"
        ES_gar$outcome <- out
        ES_gar$buffers <- buffer
        
        if (exists("ES_Df")) {
          ES_Df<-rbind(ES_Df,ES_gar)
        } else {
          ES_Df <- ES_gar
        }
        }
        
        if ("bjs" %in% estimator) {
          if (is.null(xformla)) {
            bjs_first_stage = ~1 | time + cluster
          } else {
            xformla_null = paste0("0 + ", as.character(xformla)[[2]])
            bjs_first_stage = stats::as.formula(glue::glue("~ 1 + {xformla_null} | time + cluster"))
          }
          data$treat_ind <- 0
          data[data$time>=data$treat,"treat_ind"]<-1
          data$rel_time <- data$time - data$treat
          
          mod_bjs<-didimputation::did_imputation(
            data = data[data$distance==buffer & data$rel_treat>=5,],
            yname = out,
            gname = "treat",
            tname = "time",
            idname = "cluster",
            first_stage = bjs_first_stage,
            #horizon = T,
            #pretrends = T
          )  
          ES_bjs <- mod_bjs[,-c(1:2)]
            ES_bjs$buffers <- buffer
            ES_bjs$period <- paste("+",event,sep="")
            colnames(ES_bjs) <- c("att","se","lci","uci","buffers","period")
            ES_bjs$estimator <- "Borusyak et al. (2022)"
            ES_bjs$outcome <- out
            
          if (exists("ES_Df")) {
            ES_Df<-rbind(ES_Df,ES_bjs)
          } else {
            ES_Df <- ES_bjs
          }
        }
        
        if ("csa" %in% estimator) {
          mod_csa<-did::att_gt(
            yname = out,
            tname= "time",
            idname="cluster",
            gname = "treat",
            data = data[data$distance==buffer,],
            xformla = xformla,
            base_period = pretrend.base,
            #anticipation = anticipation,
            clustervars = "cluster",
            control_group = control,
            bstrap=T,
            #est_method=est_type,
            cband=T
          )
          
          est_agg <- did::aggte(mod_csa, type = "dynamic",na.rm=T)
          ES_csa <- data.frame(period = est_agg$egt, 
                               att = est_agg$att.egt, 
                               se = est_agg$se.egt,
                               lci = est_agg$att.egt - est_agg$crit.val.egt*est_agg$se.egt,
                               uci = est_agg$att.egt + est_agg$crit.val.egt*est_agg$se.egt)
          
          ES_csa <- ES_csa %>% 
            mutate(buffers = buffer, estimator = "Callaway & Sant'Anna (2021)", outcome = out)

          if (!exists("ES_Df")) {
            ES_Df <- ES_csa
          } else {
            ES_Df<-rbind(ES_Df,ES_csa)
          }
        }
      }
    }
  return(ES_Df) 
  }


ATT_by_dist <- function(data, outcomes = "defRate_acc", buffers, max_e = Inf, xformla=NULL,estimator = "csa",plot_title=NULL,ring=F) {
  #Plot aggregate ATTs by distance for different estimators
  # max_e:     defines the post-treatment periods for which to aggregate the ATT. 
  # buffers:   buffer distances to be included
  ES_Df <-data.frame()
  for (out in outcomes) {
    for (buffer in buffers) {
      if ("gar" %in% estimator) {
        if (is.null(xformla)) {
          first_stage = ~ 0 | cluster + time 
        } else {
          first_stage <- stats::as.formula(glue::glue("~ 0 + {xformla} | cluster + time")[[2]])
        }
        
        data$treat_ind <- 0
        data[data$time>=data$treat,"treat_ind"]<-1
        data$rel_time <- data$time - data$treat
        
        mod_gar<-did2s::did2s(
          data = data[data$distance==buffer,],
          yname = outcome,
          first_stage = first_stage,
          second_stage = ~ treat_ind,
          treatment = "treat_ind",
          cluster_var = "cluster"
        )
        ES_gar<-broom::tidy(mod_gar)
        ES_gar$term <- sub("rel_time::", "", ES_gar$term)
        ES_gar <- ES_gar[ES_gar$term %in% c(-20:20),1:3]
        ES_gar$term <- as.numeric(ES_gar$term)
        colnames(ES_gar) <- c("rel_time","ATT","SE")
        ES_gar$lci <- ES_gar$ATT - 1.96*ES_gar$SE
        ES_gar$uci <- ES_gar$ATT + 1.96*ES_gar$SE
        
        colnames(ES_gar) <- c("rel_time","ATT","SE","lci","uci")
        ES_gar$estimator <- "Gardner (2022)"
      }
      
      if ("bjs" %in% estimator) {
        if (is.null(xformla)) {
          bjs_first_stage = ~1 | time + cluster
        } else {
          xformla_null = paste0("0 + ", as.character(xformla)[[2]])
          bjs_first_stage = stats::as.formula(glue::glue("~ 1 + {xformla_null} | time + cluster"))
        }
        mod_bjs<-didimputation::did_imputation(
          data = data[data$distance==buffer,],
          yname = outcome,
          gname = "treat",
          tname = "time",
          idname = "cluster",
          first_stage = bjs_first_stage,
          #horizon = T,
          #pretrends = T
        )
        ES_bjs <- mod_bjs[,-c(1:2)]
        ES_bjs$buffers <- buffer
        ES_bjs$period <- paste("+",event,sep="")
        colnames(ES_bjs) <- c("att","se","lci","uci","buffers","period")
        ES_bjs$estimator <- "Borusyak et al. (2022)"
        
        if (exists("ES_Df")) {
          ES_Df<-rbind(ES_Df,ES_bjs)
        } else {
          ES_Df <- ES_bjs
        }
      }
      
      if ("csa" %in% estimator) {  
        mod_csa<-did::att_gt(
          yname = out,
          tname="time",
          idname="cluster",
          gname = "treat",
          data = data[data$distance==buffer,],
          xformla = xformla,
          clustervars = "cluster",
          control_group = "notyettreated",
          bstrap=T,
          cband=T
        )
        
        est_agg<-did::aggte(mod_csa, type = "simple",na.rm=T)
        estcofs <- data.frame(est_agg$overall.att,est_agg$overall.se,est_agg$overall.att - 1.96*est_agg$overall.se,est_agg$overall.att + 1.96*est_agg$overall.se,buffer)
        colnames(estcofs) <- c("att","se","lci","uci","buffers")
        estcofs$outcome <- out
        estcofs$estimator <- "Callaway & Sant'Anna (2021)"
        
        if (!exists("ES_Df")) {
          ES_Df <- estcofs
        } else {
          ES_Df<-rbind(ES_Df,estcofs)
        }
      }
    }
}
return(ES_Df)
}
# 

SpilloverEst<-function(outcomes,data,firststage=NULL,dist,pre_periods) {
  #Function to plot area of common support to test PSM
  ResDF <- data.frame()
  if (is.null(firststage)) {
    first_stage = ~ 0 | cluster + time 
  } else {
    first_stage <- stats::as.formula(glue::glue("~ 0 + {xformla} | cluster + time")[[2]])
  }
  
  data$treat_ind <- 0
  data[data$time>=data$treat,"treat_ind"]<-1
  data$rel_time <- data$time - data$treat
  # if (type=="static") {
  #   second_stage_spill = ~ treat_ind + ring10000_ind
  #   second_stage_nospill = ~ treat_ind
  # } else {
    second_stage_spill = ~ i(rel_treat, ref = c(-1)) + i(rel_ring10000, ref = c(-Inf, -1))
    second_stage_nospill = ~ i(rel_treat, ref = c(-1))
  # }
   for (out in outcomes) {
     print(out)
  # Estimate model with spillover exposure
  mod_spill<-did2s::did2s(
    data = data[data$distance==dist & data$rel_treat>= pre_periods,],
    yname = out,
    first_stage = first_stage,
    second_stage = second_stage_spill,
    treatment = "treat_or_spill",
    cluster_var = "cluster"
  )
  
# Results as dynamic event study
  
  # With Spill 
  
  ES_spill <- broom::tidy(mod_spill)
  ES_spill <- ES_spill %>% separate(term, into = c("group", "period"),sep="::") %>% 
    mutate(period = as.numeric(period),
           lci = estimate - 1.96*std.error,
           uci = estimate + 1.96*std.error) %>%
    dplyr::select(period,att=estimate,se=std.error,lci,uci,group) %>%
    bind_rows(tibble(period = -1, att = 0, se = 0, lci = 0, uci = 0, group = "rel_treat"),
              tibble(period = -1, att = 0, se = 0, lci = 0, uci = 0, group = "rel_ring10000")) %>%
              mutate(estimator = case_when(
                group == "rel_treat" ~ "Treatment effect",
                group == "rel_ring10000" ~ "Spillover effect"),
                outcome = out,
                buffers = dist) %>%
    dplyr::select(!group)
    
  
  ResDF <- rbind(ResDF,ES_spill)
  
  # Without Spill 
  
  mod_nospill<-did2s::did2s(
    data = data[data$distance==dist & data$rel_treat>= pre_periods,],
    yname = out,
    first_stage = first_stage,
    second_stage = second_stage_nospill,
    treatment = "treat_ind",
    cluster_var = "cluster"
  )
  
  ES_nospill <- broom::tidy(mod_nospill)
  ES_nospill <- ES_nospill %>% 
    mutate(
    period = as.numeric(sub(".*::", "", term)),
    lci = estimate - 1.96*std.error,
    uci = estimate + 1.96*std.error) %>%
    dplyr::select(period,att=estimate,se=std.error,lci,uci) %>%
    bind_rows(tibble(period = -1, att = 0, se = 0, lci = 0, uci = 0)) %>%
    mutate(estimator = "Gardner (2022)",
      outcome = out,
      buffers = dist)
  
    ResDF <- rbind(ResDF,ES_nospill)
}
  return(ResDF)
  }
  
  # Results as static coefficients
  # if (type=="static") {
  #   plot_spill_static <- jtools::plot_summs(mod_spill,mod_nospill,model.names=c("with spillover","without spillover"),colors=c("viridis"),coefs=c( "Treatment effect"= "treat_ind", "Spillover effect (10km)" = "ring10000_ind")) + theme(legend.position="bottom")
  #   plot_title=glue::glue("static_spillover_{outcome}.jpg")
  #   ggsave(filename = here(figures,output_folder,plot_title), plot = plot_spill_static,width = 5.5,height = 3)
  # }
  
  #effect_plot(mod_gar,pred=s3aPAs,data=data[data$distance==dist,])


ATThetero <- function(data,distance=5000, first_stage= NULL,outcomes="def_ha_acc", gname = "treat",wname = NULL, wtr = NULL, horizon = NULL,pretrends = NULL, cluster_var = NULL, second_stage = NULL) {
  #Calculate individual treatment effects for heterogeneity analyisis
  
  zz000weight <- zz000adj <- NULL
  lhs <- term <- NULL
  
  # Extract first stage vars from formula
  if (is.null(first_stage)) {
    first_stage <- glue::glue("0 | cluster + time")
  } else if (inherits(first_stage, "formula")) {
    first_stage <- as.character(first_stage)[[2]]
  }
  
  for (yname in outcomes) {
  
  # Formula for fitting the first stage
  formula <- stats::as.formula(glue::glue("{yname} ~ {first_stage}"))
  
  # Get list of event_time
  event_time <- unique(data$rel_treat[is.finite(data$rel_treat)])
  
  # horizon/allhorizon options
  if (is.null(wtr)) {
    
    # event-study
    if (!is.null(horizon)) {
      # create event time weights
      
      # allhorizons
      if (all(horizon == TRUE)) horizon <- event_time
      
      # Create wtr of horizons
      wtr <- paste0("zz000wtr", event_time[event_time >= 0])
      
      # Generate one wtr dummy variable for each event time
     # data <- copy(as.data.frame(data)) %>% setDT()
      #  setDT(data)
      e = 0
      for (var in wtr) {
        data[,var] <- if_else(is.na(rel_treat), 0, 1 * (rel_treat == e))
      }

    } else {
      wtr <- "zz000wtrtreat"
      data[, wtr] <- 1 * data$treat_ind == 1
    }
  }
  
  # Weights specified or not
  if (is.null(wname)) {
    data[, "zz000weight"] <- 1 #All get weight of 1 assigned
  } else {
    data[, "zz000weight"] <- data[,wname]
  }
  
  # Estimate Y(0) using untreated observations
  first_stage_est <- fixest::feols(formula,
                                   se = "standard",
                                   data = data[data$treat_ind == 0, ],
                                   weights = ~zz000weight,
                                   warn = FALSE, notes = FALSE
  )
  
  # Residualize outcome variable(s)
    data[, paste("adj_out",yname,sep="_")] <- data[[yname]] - stats::predict(first_stage_est, newdata = data)
  #colnames(data)[colnames(data) %in% c("size_mine","s3hAccess","s3iSuitability","s3aPAs","s3bRivers")] <- c("Size of Mine","Travel time to cities","Agricultural suitability","Distance to PAs","Distance to river")
  }
  return(data) 
  }

lo.boot <- function(y="string.y", x="string.x", data, N=500, seq=NULL, cil=0.025, cih=0.975, span=0.75,
                    plot.boot=FALSE, plot.ci=FALSE, plot.o=FALSE, ...){ 
  
  #if seq is NULL then use equal interval along range of X data
  if(is.null(seq)){
    seq <- seq(min(data[x]), max(data[x]), length.out=30)
  }
  
  #store Y and X data
  data <- data[c(y, x)]
  
  #number of rows in data
  datasize <- dim(data)[1]
  
  #plot original data
  if(plot.o){
    plot(data[,2], data[,1])
  }#end if(plot.o)
  
  #create matrix to hold loess predictions
  predmat <- matrix(0, N, length(seq))
  
  #loop for bootstraps
  for(i in 1:N){
    
    #take sample of rows of length(X)
    xind <- sample(1:datasize, datasize, replace=T)
    #create X,Y data based on random sample
    x=data[xind,]
    #loess on random sample
    low <- loess(x[,1] ~ x[,2], span=span, ...)
    predmat[i,] <- predict(low, seq)     
    
    #plot bootstraps if called
    if(plot.boot){
      
      lines(seq, predmat[i,])     #Plot a sample bootstrap curve.   
      
    }#end if
    
  }#end for
  
  #dim confidence interval vectors
  cih.lo <- vector()
  cil.lo <- vector()
  warn <-vector()
  j <- 1
  
  #store CI data foreach prediction of seq
  for(i in 1:length(seq)){
    
    #check for missing values in the prdicted data
    if(any(is.na(predmat[,i]))){
      warn[j] <- i
      j <- j + 1
    }#end if
    
    cih.lo[i] <- quantile(predmat[,i], cih, na.rm=T)
    cil.lo[i] <- quantile(predmat[,i], cil, na.rm=T)
  }#end for
  
  #give warning if there were any NAs in the predicted data
  if(j > 1){
    warning("There were missing predicted values in column(s) ", warn, "\n", 
            "You may want to reduce the range of your sequence to avoid predictions at the tails", call.=F)
  }
  
  #plot pointwise bootstrap confidence intervals if called
  if(plot.ci){		
    lines(seq, cih.lo, col="red")
    lines(seq, cil.lo, col="red")
  }
  
  #store output
  out <- list(predmat=predmat, cil=cil.lo, cih=cih.lo)
  return(out)
}


#   for (covar in c("Size of Mine","Travel time to cities","Agricultural suitability","Distance to PAs","Distance to river")) {
#   pdf(file=here(figures,"Heterogeneity",glue::glue("hetero_{yname}_{covar}.pdf")))
#   est <- binsreg::binsreg(
#     y = data$adj_out,x = data[,covar],w=data$treat_ind,
#     line = T, ci = T,cb=T
#   )
#   est$bins_plot + theme_minimal() + xlab(covar) + ylab("ATT in hectar")
#   dev.off()
#   #het_plot <- est$bins_plot + theme_minimal() + xlab("var") + ylab("ATT in hectar")
#   #ggsave(file=here(figures,"Heterogeneity",glue::glue("hetero_{yname}_{covar}.jpg")),plot = het_plot)
#   }
#   return(est)
# }

#Compare coefficients across models

# modComp<- function(outcomes = c("defRatecirc","defRate_acccirc","def_hacirc","def_ha_acccirc")) {
# 
# ES_Df <- c()
# for (outcome in outcomes) {
#   print(outcome)
#   mod_csa<-did::att_gt(
#     yname = outcome,
#     gname = "treat",
#     tname = "time",
#     idname = "cluster",
#     data = o1aDfPanelRing[o1aDfPanelRing$distance==5000,],
#     xformla = NULL,
#     clustervars = "cluster",
#     control_group = "notyettreated",
#     bstrap=T,
#     cband=T
#   )
#   
#   est_agg<-did::aggte(mod_csa, type = "simple",na.rm=T)
#   ES_csa <- data.frame(est_agg$overall.att,est_agg$overall.se,est_agg$overall.att - 1.96*est_agg$overall.se,est_agg$overall.att + 1.96*est_agg$overall.se)
#   colnames(ES_csa) <- c("ATT","SE","lci","uci")
#   ES_csa$estimator <- "Callaway & Sant'Anna (2021)"
#   
#   mod_gar <- did2s::did2s(
#     data = o1aDfPanelRing[o1aDfPanelRing$distance==5000,],
#     yname = outcome,
#     first_stage = ~ 0 | cluster + time,
#     second_stage =  ~ treat_ind,
#     treatment = "treat_ind",
#     cluster_var = "cluster"
#   )
#   
#   ES_gar<-broom::tidy(mod_gar)
#   ES_gar <- ES_gar[,2:3]
#   colnames(ES_gar) <- c("ATT","SE")
#   ES_gar$lci <- ES_gar$ATT - 1.96*ES_gar$SE
#   ES_gar$uci <- ES_gar$ATT + 1.96*ES_gar$SE
#   ES_gar$estimator <- "Gardner (2022); Borusyak et al. (2022)"
#   
#   mod_but<-did2s::did2s(
#     data = o1aDfPanelRing[o1aDfPanelRing$distance==5000,],
#     yname = outcome,
#     first_stage = ~ 0 | cluster + time,
#     second_stage =  ~ treat_ind + ring1000_ind + ring2000_ind + ring3000_ind + ring4000_ind + ring5000_ind + ring6000_ind + ring7000_ind + ring8000_ind + ring9000_ind + ring10000_ind,
#     treatment = "treat_or_spill",
#     cluster_var = "cluster"
#   )
#   
#   ES_but<-broom::tidy(mod_but)
#   ES_but <- ES_but[ES_but$term=="treat_ind",2:3]
#   colnames(ES_but) <- c("ATT","SE")
#   ES_but$lci <- ES_but$ATT - 1.96*ES_but$SE
#   ES_but$uci <- ES_but$ATT + 1.96*ES_but$SE
#   ES_but$estimator <- "Butts (2023), spillover-robust"
#   
#   
#   
#   mod_bjs <- didimputation::did_imputation(
#     data = o1aDfPanelRing[o1aDfPanelRing$distance==5000,],
#     yname = outcome,
#     gname = "treat",
#     tname = "time",
#     idname = "cluster",
#     first_stage = ~1 | time + cluster,
#   )  
#   
#   
#   ES_bjs <- mod_bjs[,3:6]
#   colnames(ES_bjs) <- c("ATT","SE","lci","uci")
#   ES_bjs$estimator <- "Borusyak et al. (2022)"
#   
#   ES_Dfout <- rbind(ES_csa,ES_bjs,ES_gar,ES_but)
#   ES_Dfout$outcome <- outcome
#   ES_Df <- rbind(ES_Df,ES_Dfout)
#   
#   
#   ggplot(data = ES_Df[ES_Df$outcome %in% c("defRatecirc"),], 
#          aes(x = ATT, y = outcome, xmin = lci, xmax = uci, 
#              color = estimator, shape = estimator)) +
#     geom_pointrange(position = position_dodge(width = 0.5)) +
#     labs(title = "Model Estimates of Brain and Body Weight on REM Sleep",
#          x = "Coefficient Estimate",
#          y = "Predictor",
#          caption = "abcde") +
#     scale_y_discrete(labels = c("Annual share of forest loss")) +
#     ggpubr::theme_pubclean(flip = TRUE)
# }
# 
# 
# ggplot(data = ES_Df[ES_Df$outcome %in% c("defRatecirc","defRate_acccirc"),], 
#        aes(x = ATT, y = outcome, xmin = lci, xmax = uci, 
#            color = estimator, shape = estimator)) +
#   geom_pointrange(position = position_dodge(width = 0.5)) +
#   labs(title = "Estimates for different models and outcomes",
#        x = "ATT estimate",
#        y = "",
#        caption = "abcde") +
#   scale_y_discrete(labels = c("Annual share of forest loss","Accumulated share of forest loss")) +
#   ggpubr::theme_pubclean(flip = TRUE)
# 
# # 
# # huxtable::quick_pdf(
# #   ES_bjs,
# #   file = here(scrap,"huxtable-output.pdf"),
# #   borders = 0.4,
# #   open = interactive()
# # )
# }
