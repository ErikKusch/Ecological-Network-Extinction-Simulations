#' ####################################################################### #
#' PROJECT: [Biodiversity Simplification & Ecological Network Topology] 
#' CONTENTS: 
#'  - Gbif Data Retrieval
#'  DEPENDENCIES:
#'  - 
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# GBIF Occurrences =========================================================
Gbif_Species <- function(species = NULL, year_vec = 2000:2020){
  n_species <- length(unique(species))
  for(k in year_vec){
    message(paste("Year", k, "/" , max(year_vec)))
    Gbif <- NA
    while(is.na(Gbif)){
      try(Gbif <- occ_data(scientificName = unique(species),
                           hasGeospatialIssue = FALSE,
                           hasCoordinate = TRUE,
                           year = k,
                           limit = 1e5
      ))
    }
    if(n_species > 1){
      Data_ls <- lapply(Gbif, '[[', 2) # extract downloaded gbif data into list object, each element corresponds to a genuskey
      if(k == year_vec[[1]]){
        Data <- Data_ls
      }else{
        Data <- pblapply(X = 1:n_species, FUN = function(x){
          data.frame(
            rbindlist(list(Data[[x]],
                           Data_ls[[x]]),
                      fill = TRUE)
          )
        }
        )
      }
      names(Data) <- names(Data_ls)
    }else{
      if(k == year_vec[[1]]){
        Data <- as.data.frame(Gbif$data)
      }else{
        Data <- data.frame(rbindlist(list(Data, as.data.frame(Gbif$data)), fill = TRUE))
      }
    }
  }
  return(Data)
}

# GBIF Outliers ============================================================
Gbif_Outliers <- function(x = NULL, Enviro_ras = Enviro_ras, Centroids = Shapes_ct){
  coordinates(x) <- ~ decimalLongitude + decimalLatitude
  ## Environmental outliers, classified as beyond mean +/- 4*sd
  enviro_vals <- data.frame(raster::extract(Enviro_ras, x))
  enviro_vals$key <- x$key
  enviro_vals <- na.omit(enviro_vals)
  x <- x[x$key %in% enviro_vals$key, ]
  sd_vec <- apply(enviro_vals[,-3], MARGIN = 2, FUN = sd, na.rm = TRUE)
  mean_vec <- apply(enviro_vals[,-3], MARGIN = 2, FUN = mean, na.rm = TRUE)
  x$Out_Enviro <- FALSE
  Out_pos <- which(enviro_vals[,1] > mean_vec[1]+4*sd_vec[1] | enviro_vals[,1] < mean_vec[1]-4*sd_vec[1] |
                     enviro_vals[,2] > mean_vec[2]+4*sd_vec[2] | enviro_vals[,2] < mean_vec[2]-4*sd_vec[2])
  if(length(Out_pos)>0){x$Out_Enviro[Out_pos] <- TRUE}
  ## Centroids, matched to centroids defined above and coordinates rounded to fourth digit
  x$Out_Centroid <- FALSE
  Out_pos <- which(round(coordinates(x), 4)[,1] %in% round(coordinates(Centroids), 4)[,1] &
                     round(coordinates(x), 4)[,2] %in% round(coordinates(Centroids), 4)[,2])
  if(length(Out_pos)>0){x$Out_Centroid[Out_pos] <- TRUE}
  return(x)
}

# Climate Preferences ======================================================
Clim_Preferences <- function(data = occ_ls, Enviro_ras = Enviro_ras, Outliers = TRUE, Boot = 1e3){
  occ_spd <- do.call(rbind, data)
  occ_spd$spec <- rep(names(data), unlist(lapply(data, nrow)))
  spec_vec <- unique(occ_spd$spec)
  preferences_df <- data.frame(spec = NA,
                               Temp_median = NA, 
                               Temp_sd = NA,
                               Water_median = NA,
                               Water_sd = NA) 
  set.seed(42)
  counter <- 1
  pb <- txtProgressBar(min = 0, max = length(spec_vec), style = 3)
  
  for(spec in spec_vec){
    occ_iter <- occ_spd[occ_spd$spec == spec, ]
    if(Outliers){
      NonOut_pos <- rowSums(as.data.frame(occ_iter)[, 4:5]) == 0 # 4:5 are the outlier colums here
      occ_iter <- occ_iter[NonOut_pos,]
    }
    extract_df <- na.omit(raster::extract(Enviro_ras, occ_iter))
    values_df <- data.frame(
      X1 = NA,
      X2 = NA
    )
    for(i in 1:Boot){
      rows <- sample(x = 1:nrow(extract_df), size = round(nrow(extract_df)*0.7), replace = TRUE)
      values_df <- rbind(values_df, extract_df[rows, ])
    }
    values_df <- na.omit(values_df)
    
    preferences_df <- rbind(preferences_df, 
                            c(spec,
                              apply(values_df, 2, median)[1], apply(values_df, 2, sd)[1],
                              apply(values_df, 2, median)[2], apply(values_df, 2, sd)[2])
    )
    setTxtProgressBar(pb, counter)
    counter <- counter + 1
  }
  preferences_df <- na.omit(preferences_df)
  return(preferences_df)
}

# Network Topology =========================================================
FUN_Topo <- function(plot_df){
  Nes <- try(bipartite::networklevel(web = plot_df, index = "weighted nestedness"), silent = TRUE)
  Mod <- try(bipartite::NOS(web = plot_df)$mod, silent = TRUE)
  round(
    data.frame(
      n_species = ncol(plot_df),
      n_animals = sum(colnames(plot_df) %in% rownames(animals_gowdis)),
      n_plants = sum(colnames(plot_df) %in% rownames(plants_gowdis)),
      n_links = sum(plot_df>0),
      Nestedness = as.numeric(ifelse(class(Nes) == "try-error", NA, Nes)),
      Modularity = as.numeric(ifelse(class(Mod) == "try-error", NA, Mod))
    ), 4)
}

# Extinction Simulation ====================================================
FUN_SimComp <- function(PlantAnim = NULL, # should be set either to a vector of plant species names or animal species names
                        RunName = "ALL", # for file naming
                        IS = 0.3,
                        CutOffs,
                        WHICH = c("Strength", "Climate", "IUCN")
){
  
  print(paste0(RunName, "; IS = ", IS))
  if(file.exists(file.path(Dir.Exports, paste0(RunName, "SimulationNets_", 
                                               IS, 
                                               "_CutOffs_", paste(unlist(CutOffs), collapse = "-"), 
                                               ".RData")))){
    print("Extinctions already simulated")
    load(file.path(Dir.Exports, paste0(RunName, "SimulationNets_", 
                                       IS, 
                                       "_CutOffs_", paste(unlist(CutOffs), collapse = "-"), 
                                       ".RData")))
    return(Sim_ls)
  }else{
    Sim_ls <- pblapply(names(AnalysisData_ls), 
                       cl = cl,
                       function(x){
                         # print(x)
                         # x <- names(AnalysisData_ls)[1]
                         x <- AnalysisData_ls[[x]]
                         net <- as.network(x$Adjacency, matrix.type = "adjacency", 
                                           ignore.eval=FALSE, names.eval='weight')
                         
                         ## Subsetting proxies for bottom-up/top-down simulations
                         if(!is.null(PlantAnim)){
                           x$prox_centrality <- x$prox_centrality[names(x$prox_centrality) %in% PlantAnim]
                           x$prox_climate <- x$prox_climate[names(x$prox_climate) %in% PlantAnim]
                           x$prox_IUCN <- x$prox_IUCN[names(x$prox_IUCN) %in% PlantAnim]
                         }
                         
                         ## Centrality-Driven -------------------------------------------------------
                         # print("Extinction of Keystone Species (Centrality)")
                         proxcen <- x$prox_centrality[x$prox_centrality > quantile(x$prox_centrality, CutOffs$Strength)] # just eliminate upper 25% quantile
                         primext_namesS <- names(proxcen)
                         primext_order <- match(primext_namesS, rownames(x$Adjacency))
                         CustOrder_ExtS <- ExtS_Rand <- as.list(c(NA, NA))
                         if("Strength" %in% WHICH){
                           CustOrder_ExtS <- SimulateExtinctions(Network = net, Method = "Ordered", Order = primext_order, IS = IS)
                           ExtS_Rand <- RandomExtinctions(Network = net, nsim = 100, 
                                                            parallel = FALSE, ncores = parallel::detectCores(), 
                                                            SimNum = length(proxcen),
                                                            IS = IS)
                           # CompareExtinctions(Nullmodel = ExtS_Rand[[1]], Hypothesis = CustOrder_ExtS[[1]]) 
                         }
                         
                         ## Climate-Driven ----------------------------------------------------------
                         # print("Extinction of Threatened Species (Climate Projections)")
                         primext_namesC <- names(x$prox_climate)[x$prox_climate > CutOffs$Climate] # random cutoff of climate risk severity selected here
                         primext_order <- match(primext_namesC, rownames(x$Adjacency))
                         CustOrder_ExtC <- ExtC_Rand <- as.list(c(NA, NA))
                         if("Climate" %in% WHICH){
                           if(length(primext_namesC) != 0){
                             CustOrder_ExtC <- SimulateExtinctions(Network = net, Method = "Ordered", Order = primext_order,
                                                                 IS = IS)
                             ExtC_Rand <- RandomExtinctions(Network = net, nsim = 100, 
                                                              parallel = FALSE, ncores = parallel::detectCores(), 
                                                              SimNum = length(primext_namesC),
                                                              IS = IS)
                             # CompareExtinctions(Nullmodel = Rando_Ext, Hypothesis = CustOrder_ExtC)
                           } 
                         }
                         
                         ## IUCN-Driven -------------------------------------------------------------
                         # print("Extinction of Threatened Species (IUCN Categories)")
                         primext_namesI <- names(x$prox_IUCN)[x$prox_IUCN > CutOffs$Climate] # random cutoff of climate risk severity selected here
                         primext_order <- match(primext_namesI, rownames(x$Adjacency))
                         CustOrder_ExtI <- ExtI_Rand <- as.list(c(NA, NA))
                         if("IUCN" %in% WHICH){
                           if(length(primext_namesI) != 0){
                             CustOrder_ExtI <- SimulateExtinctions(Network = net, Method = "Ordered", Order = primext_order,
                                                                              IS = IS)
                             ExtI_Rand <- RandomExtinctions(Network = net, nsim = 100, 
                                                              parallel = FALSE, ncores = parallel::detectCores(), 
                                                              SimNum = length(primext_namesI),
                                                              IS = IS)
                             # CompareExtinctions(Nullmodel = Rando_Ext, Hypothesis = CustOrder_ExtI)
                           } 
                         }
                         
                         ## Export ------------------------------------------------------------------
                         list(Strength = list(Removed = primext_namesS,
                                              Prediction = as.matrix(CustOrder_ExtS[[2]]),
                                              Random = lapply(ExtS_Rand$nets, as.matrix)
                         ),
                         Climate = list(Removed = primext_namesC,
                                        Prediction = as.matrix(CustOrder_ExtC[[2]]),
                                        Random = lapply(ExtC_Rand$nets, as.matrix)
                         ),
                         IUCN = list(Removed = primext_namesI,
                                     Prediction = as.matrix(CustOrder_ExtI[[2]]),
                                     Random = lapply(ExtI_Rand$nets, as.matrix)
                         )
                         )
                       })
    names(Sim_ls) <- names(AnalysisData_ls)
    save(Sim_ls, file = file.path(Dir.Exports, paste0(RunName, "SimulationNets_", 
                                                      IS, 
                                                      "_CutOffs_", paste(unlist(CutOffs), collapse = "-"), 
                                                      ".RData"))) 
    return(Sim_ls)
  }
}

# Network Topology Comparison ==============================================
FUN_TopoComp <- function(Sim_ls = NULL, RunName = "ALL", IS, CutOffs){
  
  if(file.exists(file.path(Dir.Exports, paste0(RunName, "SimulationTopo_", IS, 
                                               "_CutOffs_", paste(unlist(CutOffs), collapse = "-"), 
                                               ".RData")))){
    print("Topology already extracted")
    load(file.path(Dir.Exports, paste0(RunName, "SimulationTopo_", IS, 
                                       "_CutOffs_", paste(unlist(CutOffs), collapse = "-"), 
                                       ".RData")))
    return(Save_ls)
  }else{
    ## Topology Calculation
    PostExt_ls <- pblapply(names(Sim_ls), 
                           cl = cl,
                           FUN = function(netID){
                             Storage_ls <- list(Strength = list(Prediction = NA, Random = NA),
                                                Climate = list(Prediction = NA, Random = NA),
                                                IUCN = list(Prediction = NA, Random = NA))
                             for(i in names(Storage_ls)){
                               Storage_ls[[i]][["Prediction"]] <- FUN_Topo(as.matrix(Sim_ls[[netID]][[i]][["Prediction"]]))
                               Rand_ls <- lapply(Sim_ls[[netID]][[i]][["Random"]], FUN_Topo)
                               Storage_ls[[i]][["Random"]] <- do.call(rbind, Rand_ls)
                             }
                             Storage_ls
                           })
    # parallel::stopCluster(cl)
    names(PostExt_ls) <- names(Sim_ls)
    
    
    ## Topology Extraction
    Topo_ls <- list(Strength = NA,
                    Climate = NA,
                    IUCN = NA)
    for(k in c("Strength", "Climate", "IUCN")){
      # print(k)
      Pred_df <- do.call(rbind, lapply(lapply(PostExt_ls, "[[", k), "[[", "Prediction"))
      Pred_df$netID <- names(Sim_ls)
      Pred_df$Proxy <- k
      Pred_df$Simulation <- "Prediction"
      Rand_df <- do.call(rbind, lapply(lapply(PostExt_ls, "[[", k), "[[", "Random"))
      if(is.null(Rand_df)){ # this happens when no simulation could be run for the entire list of networks for this proxy
        next()
      }else{
        Rand_df$netID <- rep(names(unlist(lapply(lapply(lapply(PostExt_ls, "[[", k), "[[", "Random"), nrow))), # these are the names for which random sims could be run
                             each = 1e2)
        Rand_df$Proxy <- k
        Rand_df$Simulation <- "Random"
        # print(head(Pred_df))
        # print(head(Rand_df))
        Topo_ls[[k]] <- rbind(Pred_df, Rand_df) 
      }
    }
    Topo_df <- do.call(rbind, Topo_ls)
    ## Effect size calculation
    MeanPred <- aggregate(.~netID+Proxy+Simulation, FUN = mean, data = Topo_df[Topo_df$Simulation == "Prediction",])
    MeanRand <- aggregate(.~netID+Proxy+Simulation, FUN = mean, data = Topo_df[Topo_df$Simulation == "Random",])
    SDRand <- aggregate(.~netID+Proxy+Simulation, FUN = sd, data = Topo_df[Topo_df$Simulation == "Random",])
    Merge1_df  <- merge(MeanPred, MeanRand, by = c("netID", "Proxy"), all = TRUE)
    Merge1_df  <- merge(Merge1_df, SDRand, by = c("netID", "Proxy"), all = TRUE)
    Eff_df <-  (
      Merge1_df[,4:ncol(MeanPred)] - # mean predictions
        Merge1_df[,(ncol(MeanPred)+2):(ncol(MeanPred)+2+ncol(MeanPred)-4)]  # mean randoms
    )/  Merge1_df[,(ncol(Merge1_df)-(ncol(MeanPred)-4)):ncol(Merge1_df)]
    Eff_df$netID <- Merge1_df$netID
    Eff_df$Proxy <- Merge1_df$Proxy
    
    
    Save_ls <- list(Topo_ls = PostExt_ls, Topo_df = Topo_df, Eff_df = Eff_df)
    
    ## Data Return
    save(Save_ls, file = file.path(Dir.Exports, paste0(RunName, "SimulationTopo_", IS, 
                                                       "_CutOffs_", paste(unlist(CutOffs), collapse = "-"), 
                                                       ".RData")))
    return(Save_ls)
  }
}

# Loading topology data for each extinction cascade ========================
loadTopo <- function(RunName = "ALL", CutOffs){
  fs <- list.files(path = Dir.Exports, pattern = paste0(RunName, "SimulationTopo"))
  fs <- fs[grep(pattern = paste(unlist(CutOffs), collapse = "-"), fs)]
  IS_vec <- as.numeric(unlist(
    lapply(
      regmatches(fs, gregexpr("[[:digit:]]+\\.*[[:digit:]]*",fs)), 
      "[[", 1
    )
  ))
  fs <- fs[order(IS_vec)]
  
  for(i in 1:length(fs)){
    Eff2_df <- loadRData(file.path(Dir.Exports, fs[i]))$Eff_df
    Topo2_df <- loadRData(file.path(Dir.Exports, fs[i]))$Topo_df
    Eff2_df$IS <- Topo2_df$IS <- sort(IS_vec)[i]
    if(i == 1){
      Eff_df <- Eff2_df
      Topo_df <- Topo2_df
    }else{
      Eff_df <- rbind(Eff_df, Eff2_df)
      Topo_df <- rbind(Topo_df, Topo2_df)
    }
  }
  colnames(Eff_df) <- gsub(pattern = ".x", replacement = "", colnames(Eff_df))
  return(list(Topo_df, Eff_df))
}
