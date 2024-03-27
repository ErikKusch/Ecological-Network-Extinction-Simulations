#' ####################################################################### #
#' PROJECT: [Biodiversity Simplification & Ecological Network Topology] 
#' CONTENTS: 
#'  - Gbif Data Retrieval
#'  DEPENDENCIES:
#'  - 
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# SAVING & LOADING with Progress Bars ======================================
## shamelessly plagiarised from: https://stackoverflow.com/questions/6165077/how-do-i-create-a-progress-bar-for-data-loading-in-r
saveObj <- function(object, file.name){
  outfile <- file(file.name, "wb")
  serialize(object, outfile)
  close(outfile)
}

loadObj <- function(file.name){
  library(foreach)
  filesize <- file.info(file.name)$size
  chunksize <- ceiling(filesize / 100)
  pb <- txtProgressBar(min = 0, max = 100, style=3)
  infile <- file(file.name, "rb")
  data <- foreach(it = icount(100), .combine = c) %do% {
    setTxtProgressBar(pb, it)
    readBin(infile, "raw", chunksize)
  }
  close(infile)
  close(pb)
  return(unserialize(data))
}

# GBIF Occurrences =========================================================
Gbif_Species <- function(species = NULL, year_vec = 2000:2020){
  taxonKeys <- pbsapply(species, FUN = function(x){
    keys <- na.omit(name_backbone(name = x, rank = "species", verbose = TRUE)$speciesKey)
    keys[1]
  })
  res <- occ_download(
    pred_in("taxonKey", as.numeric(unlist(taxonKeys))),
    pred_not(pred("basisOfRecord", "FOSSIL_SPECIMEN")),
    pred("occurrenceStatus", "PRESENT"),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    pred_gte("year", min(year_vec)),
    pred_lte("year", max(year_vec))
  )
  print(occ_download_meta(res)$doi)
  # Downloading and Loading
  res_meta <- occ_download_wait(res, status_ping = 600, curlopts = list(), quiet = FALSE)
  res_get <- occ_download_get(res, overwrite = TRUE)
  Gbif <- occ_download_import(res_get)
  
  Data <- split(Gbif, f = Gbif$speciesKey)
  names(Data) <- names(taxonKeys)[match(names(Data), as.character(taxonKeys))]
  
  return(list(Data = Data, NonOcc = species[species %nin% names(Data)]))
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
      NonOut_pos <- rowSums(as.data.frame(occ_iter)[, c("Out_Enviro", "Out_Centroid")]) == 0 
      occ_iter <- occ_iter[NonOut_pos,]
    }
    extract_df <- na.omit(raster::extract(Enviro_ras, occ_iter))
    print(paste(spec, nrow(extract_df), sep = " - "))
    values_df <- data.frame(
      X1 = NA,
      X2 = NA
    )
    if(nrow(extract_df)<1e4){ #bootstrap only when there are less than 1e4 useable locations - this saves RAM
      for(i in 1:Boot){
        rows <- sample(x = 1:nrow(extract_df), size = round(nrow(extract_df)*0.7), replace = TRUE)
        values_df <- rbind(values_df, extract_df[rows, ])
      } 
    }else{
      values_df <- extract_df
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
FUN_Topo <- function(plot_df, 
                     animals = animals_sp, 
                     plants = plants_sp){
  if(sum(plot_df, na.rm = TRUE) == 0){
    data.frame(
      n_species = 0,
      n_animals = 0,
      n_plants = 0,
      n_links = 0,
      mean_interac = 0,
      connectedness = 0,
      Nestedness = NA,
      Modularity = NA
    )
  }else{
    bothtris <- as.matrix(
      as_adjacency_matrix(
        as.undirected(
          graph_from_adjacency_matrix(plot_df, mode = "directed")
        )
      )
    )
    
    adj_mat <- bothtris[rownames(bothtris) %in% plants_sp,
                        colnames(bothtris) %in% animals_sp]
    Nes <- try(bipartite::networklevel(web = adj_mat, index = "weighted nestedness"), silent = TRUE)
    Mod <- try(bipartite::NOS(web = adj_mat)$mod, silent = TRUE)
    round(
      data.frame(
        n_species = ncol(plot_df),
        n_animals = sum(colnames(plot_df) %in% rownames(animals_gowdis)),
        n_plants = sum(colnames(plot_df) %in% rownames(plants_gowdis)),
        n_links = sum(plot_df>0),
        mean_interac = mean(plot_df[plot_df>0]),
        connectedness = sum(plot_df>0)/length(plot_df),
        Nestedness = as.numeric(ifelse(class(Nes) == "try-error", NA, Nes)),
        Modularity = as.numeric(ifelse(class(Mod) == "try-error", NA, Mod))
      ), 4)
  }
  
  
}

# Extinction Simulation ====================================================
FUN_SimComp <- function(PlantAnim = NULL, # should be set either to a vector of plant species names or animal species names
                        RunName = "ALL", # for file naming
                        IS = 0.3,
                        Rewiring = FALSE,
                        CutOffs,
                        PotPartners, 
                        Traits,
                        WHICH = c("SSP245") # or SSP585
){
  
  writeLines(paste0("## ", RunName, " ## \nIS (% of IS required for continued existence) = ", IS, " \nRewiring Cutoff (probability of rewiring required) = ", Rewiring))
  if(file.exists(file.path(Dir.Exports, paste0(RunName, "SimulationNets_", 
                                               IS, "_", Rewiring,
                                               "_CutOffs_", paste(unlist(CutOffs), collapse = "-"), 
                                               ".RData")))){
    print("Extinctions already simulated")
    load(file.path(Dir.Exports, paste0(RunName, "SimulationNets_", 
                                       IS, "_", Rewiring,
                                       "_CutOffs_", paste(unlist(CutOffs), collapse = "-"), 
                                       ".RData")))
    return(Sim_ls)
  }else{
    Sim_ls <- pblapply(names(AnalysisData_ls), 
                       cl = cl,
                       function(y){
                         # y <- names(AnalysisData_ls)[16]
                         message(y)
                         x <- AnalysisData_ls[[y]]
                         
                         ## network object creation
                         net <- as.network(x$Adjacency, matrix.type = "adjacency", 
                                           ignore.eval=FALSE, names.eval='weight')
                         
                         ## distance matrix merging
                         NA_UR <- matrix(NA, ncol = ncol(x$gow_animals), nrow = nrow(x$gow_plants))
                         dimnames(NA_UR) <- list(rownames(x$gow_plants), colnames(x$gow_animals))
                         NA_LL <- matrix(NA, ncol = ncol(x$gow_plants), nrow = nrow(x$gow_animals))
                         dimnames(NA_LL) <- list(rownames(x$gow_animals), colnames(x$gow_plants))
                         dist_mat <- cbind(rbind(x$gow_plants, NA_LL), rbind(NA_UR, x$gow_animals))
                         
                         ## Rewiring function
                         if(Rewiring == 0){
                           RewiringFun <- FALSE
                           decay <- FALSE
                         }else{
                           # print(dim(dist_mat))
                           diag(dist_mat) <- NA
                           dist <- quantile(dist_mat, na.rm = TRUE, Rewiring) # rewiring distance distance quantile
                           decay <- 1/as.numeric(dist) # assuming mean rewiring capability lies at dist_10
                           # print(decay)
                           RewiringFun <- function(x){1-pexp(x, rate = decay)} 
                         }
                         
                         ## Rewiring Probabilities according to PotPartners
                         ### Matrix for storing rewiring probabilities
                         prob_mat <- dist_mat 
                         prob_mat[!is.na(prob_mat)] <- NA
                         ### identify relevant RF models for this network
                         specs <- get.vertex.attribute(net, "vertex.names") 
                         RewPartners_ls <- PotPartners[specs]
                         ### identify relevant trait subset
                         yTraits <- Traits[Traits$net.id == y, ]
                         anim_trait <- yTraits[!duplicated(yTraits$animal.phylo.id), ]
                         anim_trait$sp <- anim_trait$animal.phylo.id
                         plant_trait <- yTraits[!duplicated(yTraits$plant.phylo.id), ]
                         plant_trait$sp <- plant_trait$plant.phylo.id
                         ### populate probability matrix
                         for(v in colnames(prob_mat)){
                           if(v %in% colnames(animals_gowdis)){
                             calc_df <- plant_trait
                           }else{
                             calc_df <- anim_trait
                           }
                           if(class(RewPartners_ls[[v]]) == "randomForest"){
                             Probs <- predict(RewPartners_ls[[v]], calc_df)
                             names(Probs) <- calc_df$sp
                             Probs <- Probs[names(Probs) %in% rownames(prob_mat)]
                             prob_mat[names(Probs),v] <- Probs    
                           }else{
                             if(RewPartners_ls[[v]] == 0){
                               Probs <- rep(0, nrow(calc_df))
                               names(Probs) <- calc_df$sp
                               Probs <- Probs[names(Probs) %in% rownames(prob_mat)]
                               prob_mat[names(Probs),v] <- Probs    
                             }else{
                               Probs <- rep(1, nrow(calc_df))
                               names(Probs) <- calc_df$sp
                               Probs <- Probs[names(Probs) %in% rownames(prob_mat)]
                               prob_mat[names(Probs),v] <- Probs    
                             } 
                           }
                         }
                         
                         ## Subsetting proxies for bottom-up/top-down simulations
                         if(!is.null(PlantAnim)){
                           x$prox_centrality <- x$prox_centrality[names(x$prox_centrality) %in% PlantAnim]
                           x$prox_climate <- x$prox_climate[names(x$prox_climate) %in% PlantAnim]
                           x$prox_IUCN <- x$prox_IUCN[names(x$prox_IUCN) %in% PlantAnim]
                         }
                         
                         ## Climate-Driven ----------------------------------------------------------
                         # print("Extinction of Threatened Species (Climate Projections)")
                         
                         if(WHICH == "SSP245"){
                           primext_namesC <- names(x$prox_climate)[x$prox_climate > CutOffs$Climate] # random cutoff of climate risk severity selected here
                           primext_order <- match(primext_namesC, rownames(x$Adjacency))
                           RMNum <- length(primext_namesC) 
                         }
                         if(WHICH == "SSP585"){
                           primext_namesC <- names(x$prox_climateSSP585)[x$prox_climateSSP585 > CutOffs$Climate] # random cutoff of climate risk severity selected here
                           primext_order <- match(primext_namesC, rownames(x$Adjacency))
                           RMNum <- length(primext_namesC) 
                         }

                         CustOrder_ExtC <- ExtC_Rand <- as.list(c(NA, NA, NA, NA))
                         names(CustOrder_ExtC) <- names(ExtC_Rand) <- c("sims", "R50", "R100", "Network")
                         # if("SSP245" %in% WHICH){
                           if(length(primext_namesC) != 0){
                             CustOrder_ExtC <- SimulateExtinctions(Network = net, 
                                                                   Method = "Ordered", 
                                                                   Order = primext_order,
                                                                   IS = IS,
                                                                   NetworkType = "Mutualistic",
                                                                   ## PDF-driven rewiring block
                                                                   Rewiring = function(x){x},
                                                                   # decay = Rewiring
                                                                   # RewiringDist = dist_mat, #
                                                                   ### Probability matrix-driven block
                                                                   RewiringDist = prob_mat, # 
                                                                   RewiringProb = Rewiring,
                                                                   forceFULL = TRUE
                             )
                             ExtC_Rand <- RandomExtinctions(Network = net, nsim = 500,
                                                            parallel = FALSE, ncores = parallel::detectCores(),
                                                            SimNum = RMNum,
                                                            IS = IS,
                                                            NetworkType = "Mutualistic",
                                                            ## PDF-driven rewiring block
                                                            Rewiring = function(x){x},
                                                            # decay = Rewiring
                                                            # RewiringDist = dist_mat, #
                                                            ### Probability matrix-driven block
                                                            RewiringDist = prob_mat, #
                                                            RewiringProb = Rewiring,
                                                            forceFULL = TRUE
                             )
                             ExtC_Rand <- list(sims = ExtC_Rand$sims,
                                               R50 = ExtC_Rand$R50result,
                                               R50 = ExtC_Rand$R100result,
                                               Network = ExtC_Rand$nets)
                           } 
                         
                         ## Centrality-Driven -------------------------------------------------------
                         # print("Extinction of Keystone Species (Centrality)")
                         namesML <- namesLM <- primext_namesC
                         CustOrder_ExtSLM <- CustOrder_ExtSML <- as.list(c(NA, NA, NA, NA))
                         names(CustOrder_ExtSLM) <- names(CustOrder_ExtSML) <- c("sims", "R50", "R100", "Network")
                         if(length(primext_namesC) != 0){
                           ### Most Connected ----
                           print("Strength - Most to Least")
                           proxcen <- x$prox_centrality[1:RMNum]
                           primext_namesS <- names(proxcen)
                           namesML <- primext_namesS
                           primext_order <- match(primext_namesS, rownames(x$Adjacency))
                           CustOrder_ExtSML <- SimulateExtinctions(Network = net, 
                                                                   Method = "Ordered",
                                                                   Order = primext_order, 
                                                                   IS = IS,
                                                                   NetworkType = "Mutualistic",
                                                                   ## PDF-driven rewiring block
                                                                   Rewiring = function(x){x},
                                                                   # decay = Rewiring
                                                                   # RewiringDist = dist_mat, #
                                                                   ### Probability matrix-driven block
                                                                   RewiringDist = prob_mat, # 
                                                                   RewiringProb = Rewiring,
                                                                   forceFULL = TRUE)
                           
                           
                           ### Least Connected ----
                           print("Strength - Least to Most")
                           proxcen <- x$prox_centrality[(length(x$prox_centrality)-RMNum+1):length(x$prox_centrality)]
                           primext_namesS <- names(proxcen)
                           namesLM <- primext_namesS
                           primext_order <- match(primext_namesS, rownames(x$Adjacency))
                           CustOrder_ExtSLM <- SimulateExtinctions(Network = net, 
                                                                   Method = "Ordered",
                                                                   Order = primext_order, 
                                                                   IS = IS,
                                                                   NetworkType = "Mutualistic",
                                                                   ## PDF-driven rewiring block
                                                                   Rewiring = function(x){x},
                                                                   # decay = Rewiring
                                                                   # RewiringDist = dist_mat, #
                                                                   ### Probability matrix-driven block
                                                                   RewiringDist = prob_mat, # 
                                                                   RewiringProb = Rewiring,
                                                                   forceFULL = TRUE)
                         }
                         
                         ## Export ------------------------------------------------------------------
                         Fun.Save <- function(x = CustOrder_ExtC, 
                                              y = ExtC_Rand
                                              ){
                           
                           ## This happens when no simulation was run
                           if(sum(is.na(as.matrix(x$Network))) != 0){
                             Pred <- as.matrix(NA)
                             Rand <- as.matrix(NA)
                           }else{
                             ## this is complete network annihilation, if this happens, we manually assign a zero-matrix
                             if(sum(as.matrix.network.adjacency(x$Network, attrname = "weight")) == 0 | 
                                nrow(as.matrix(x$Network)) == 1){
                               Pred <- as.matrix(0)
                             }else{
                               Pred <- as.matrix.network.adjacency(x$Network, attrname = "weight")
                             }
                             Rand <- lapply(y$Network, as.matrix.network.adjacency, attrname = "weight")
                             Rand[unlist(lapply(Rand, sum)) == 0 |
                                    unlist(lapply(Rand, FUN = function(k){nrow(as.matrix(k))})) == 1] <- as.matrix(0)
                           }
                           if(nrow(as.matrix(x$Network)) != 1){ # a matrix with more than one row and column (i.e., non-total annihilation)
                             Pred <- as.matrix.network.adjacency(x$Network,
                                                                 attrname = "weight")
                             Rand <- lapply(y$Network, as.matrix.network.adjacency,
                                            attrname = "weight")
                           }else{
                             if(class(x$Network) == "network"){
                               if(!is.na(get.vertex.attribute(x$Network, "vertex.names "))){
                                 Pred <- as.matrix.network.adjacency(x$Network,
                                                                     attrname = "weight")
                                 Rand <- lapply(y$Network, as.matrix.network.adjacency,
                                                attrname = "weight")
                               }else{
                                 Pred <- as.matrix(x$Network)
                                 Rand <- lapply(y$Network, as.matrix)
                               }
                             }else{
                               if(!is.na(x$Network)){
                                 Pred <- as.matrix.network.adjacency(x$Network,
                                                                     attrname = "weight")
                                 Rand <- lapply(y$Network, as.matrix.network.adjacency,
                                                attrname = "weight")
                               }else{
                                 Pred <- as.matrix(x$Network)
                                 Rand <- lapply(y$Network, as.matrix)
                               }
                           }
                           }
                           return(list(Pred = Pred, Rand = Rand))
                         }
                         
                         list(
                         Climate = list(Removed = primext_namesC,
                                        Prediction = Fun.Save(x = CustOrder_ExtC, 
                                                              y = ExtC_Rand
                                                              )$Pred,
                                        Random = Fun.Save(x = CustOrder_ExtC, y = ExtC_Rand)$Rand
                         ),
                         MostToLeast = list(Removed = namesML,
                              Prediction = Fun.Save(x = CustOrder_ExtSML, 
                                                    y = ExtC_Rand
                              )$Pred,
                              Random = Fun.Save(x = CustOrder_ExtSML, y = ExtC_Rand)$Rand
                         ),
                         LeastToMost = list(Removed = namesLM,
                              Prediction = Fun.Save(x = CustOrder_ExtSLM, 
                                                    y = ExtC_Rand
                              )$Pred,
                              Random = Fun.Save(x = CustOrder_ExtSLM, y = ExtC_Rand)$Rand
                         )
                         )
                       })
    names(Sim_ls) <- names(AnalysisData_ls)
    save(Sim_ls, file = file.path(Dir.Exports, paste0(RunName, "SimulationNets_", 
                                                      IS, "_", Rewiring,
                                                      "_CutOffs_", paste(unlist(CutOffs), collapse = "-"), 
                                                      ".RData"))) 
    return(Sim_ls)
  }
}

# Network Topology Comparison ==============================================
FUN_TopoComp <- function(Sim_ls = NULL, RunName = "ALL", IS, Rewiring, CutOffs, Pre = PreExt_df){
  
  if(file.exists(file.path(Dir.Exports, paste0(RunName, "SimulationTopo_", 
                                               IS, "_", Rewiring,
                                               "_CutOffs_", paste(unlist(CutOffs), collapse = "-"), 
                                               ".RData")))){
    print("Topology already extracted")
    load(file.path(Dir.Exports, paste0(RunName, "SimulationTopo_", 
                                       IS, "_", Rewiring,
                                       "_CutOffs_", paste(unlist(CutOffs), collapse = "-"), 
                                       ".RData")))
    return(Save_ls)
  }else{
    ## Topology Calculation
    PostExt_ls <- pblapply(names(Sim_ls),
                           cl = cl,
                           FUN = function(netID){
                             # netID <- names(Sim_ls)[1]
                             print(netID)
                             Storage_ls <- list(Strength = list(Removed = NA, Prediction = NA, Random = NA),
                                                Climate = list(Removed = NA, Prediction = NA, Random = NA),
                                                IUCN = list(Removed = NA, Prediction = NA, Random = NA)
                                                # ,
                                                # IUCN_Climate = list(Removed = NA, Prediction = NA, Random = NA)
                             )
                             for(i in names(Storage_ls)){
                               if(length(Sim_ls[[netID]][[i]][["Removed"]]) != 0){
                                 Storage_ls[[i]][["Removed"]] <- length(Sim_ls[[netID]][[i]][["Removed"]])
                                 Storage_ls[[i]][["Prediction"]] <- FUN_Topo(as.matrix(Sim_ls[[netID]][[i]][["Prediction"]]))
                                 # Rand_ls <- lapply(Sim_ls[[netID]][[i]][["Random"]], FUN_Topo)
                                 # Storage_ls[[i]][["Random"]] <- do.call(rbind, Rand_ls) 
                               }else{
                                 Storage_ls[[i]][["Removed"]] <- 0
                                 Storage_ls[[i]][["Prediction"]] <- PreExt_df[PreExt_df$netID == netID, 1:8]
                               }
                             }
                             Storage_ls
                           })
    # parallel::stopCluster(cl)
    names(PostExt_ls) <- names(Sim_ls)
    
    
    ## Topology Extraction
    if(RunName == "SSP585"){
      Topo_ls <- list(Climate = NA
      )
    }else{
      Topo_ls <- list(Strength = NA,
                      Climate = NA,
                      IUCN = NA
                      # ,
                      # IUCN_Climate = NA
      )
    }
    
    for(k in names(Topo_ls)){
      # k <- names(Topo_ls)[3]
      # print(k)
      Rem_df <- do.call(rbind, lapply(lapply(PostExt_ls, "[[", k), "[[", "Removed"))
      Pred_df <- do.call(rbind, lapply(lapply(PostExt_ls, "[[", k), "[[", "Prediction"))
      Pred_df$netID <- names(Sim_ls)
      Pred_df$Proxy <- k
      Pred_df$Simulation <- "Prediction"
      Pred_df$Removed <- Rem_df[,1]
      # Rand_df <- do.call(rbind, lapply(lapply(PostExt_ls, "[[", k), "[[", "Random"))
      # if(is.null(Rand_df)){ # this happens when no simulation could be run for the entire list of networks for this proxy
      #   Topo_ls[[k]] <- Pred_df
      #   next()
      # }else{
      #   Rand_df <- Rand_df[unlist(lapply(strsplit(rownames(Rand_df), split = "[.]"), "[[", 1)) %in% rownames(Rem_df)[Rem_df != 0], ]
      #   if(is.logical(Rand_df)){
      #     Topo_ls[[k]] <- Pred_df
      #     next()
      #   }
      #   Rand_df$netID <- rep(names(unlist(lapply(lapply(lapply(PostExt_ls, "[[", k), "[[", "Random"), nrow))), # these are the names for which random sims could be run
      #                        each = 1e2)
      #   Rand_df$Proxy <- k
      #   Rand_df$Simulation <- "Random"
      #   Rand_df$Removed <- 9999
        # print(head(Pred_df))
        # print(head(Rand_df))
        # Topo_ls[[k]] <- rbind(Pred_df, Rand_df) 
        Topo_ls[[k]] <- Pred_df
      }
    }
    Topo_df <- do.call(rbind, Topo_ls)
    ## Effect size calculation
    MeanPred <- aggregate(.~netID+Proxy+Simulation, data = Topo_df[Topo_df$Simulation == "Prediction",],
                          FUN = mean, na.rm = TRUE, na.action=NULL)
    # MeanRand <- aggregate(.~netID+Proxy+Simulation, data = Topo_df[Topo_df$Simulation == "Random",],
    #                       FUN = mean, na.rm = TRUE, na.action=NULL)
    # SDRand <- aggregate(.~netID+Proxy+Simulation, data = Topo_df[Topo_df$Simulation == "Random",],
    #                     FUN = sd, na.rm = TRUE, na.action=NULL)
    # Merge1_df  <- merge(MeanPred, MeanRand, by = c("netID", "Proxy"), all = TRUE)
    # Merge1_df  <- merge(Merge1_df, SDRand, by = c("netID", "Proxy"), all = TRUE)
    # Eff_df <-  Merge1_df[,4:(ncol(MeanPred)-1)] - # mean predictions
    #   Merge1_df[,(ncol(MeanPred)+2):(ncol(MeanPred)+2+ncol(MeanPred)-5)]  # mean randoms
    # Eff_df$netID <- Merge1_df$netID
    # Eff_df$Proxy <- Merge1_df$Proxy
    # 
    # EffSD_df <- Merge1_df[,(ncol(Merge1_df)-(ncol(MeanPred)-4)):(ncol(Merge1_df)-1)]
    # EffSD_df$netID <- Merge1_df$netID
    # EffSD_df$Proxy <- Merge1_df$Proxy
    # colnames(Eff_df) <- colnames(EffSD_df)
    
    Save_ls <- list(Topo_ls = PostExt_ls, Topo_df = Topo_df
                    # , Eff_df = Eff_df, EffSD_df = EffSD_df
                    # , Eff_df = Eff_df, EffSD_df = EffSD_df
                    )
    
    ## Data Return
    save(Save_ls, file = file.path(Dir.Exports, paste0(RunName, "SimulationTopo_", 
                                                       IS, "_", Rewiring,
                                                       "_CutOffs_", paste(unlist(CutOffs), collapse = "-"), 
                                                       ".RData")))
    return(Save_ls)
  }
# }

# Loading topology data for each extinction cascade ========================
loadTopo <- function(RunName = "ALL", CutOffs, Pre){
  Topos_vec <- c("n_species", "n_plants", "n_animals", "n_links", "Nestedness", "Modularity")
  fs <- list.files(path = Dir.Exports, pattern = paste0(RunName, "SimulationTopo"))
  fs <- fs[grep(pattern = paste(unlist(CutOffs), collapse = "-"), fs)]
  # IS_vec <- as.numeric(unlist(
  #   lapply(
  #     regmatches(fs, gregexpr("[[:digit:]]+\\.*[[:digit:]]*",fs)), 
  #     "[[", 1
  #   )
  # ))
  # Rew_vec <- as.numeric(unlist(
  #   lapply(
  #     regmatches(fs, gregexpr("[[:digit:]]+\\.*[[:digit:]]*",fs)), 
  #     "[[", 2
  #   )
  # ))
  # fs <- fs[order(IS_vec)]
  pb <- txtProgressBar(min = 0, max = length(fs), style = 3)
  for(i in 1:length(fs)){
    ## data extraction
    # Eff2_df <- loadRData(file.path(Dir.Exports, fs[i]))$Eff_df
    # EffSD2_df <- loadRData(file.path(Dir.Exports, fs[i]))$EffSD_df
    Topo2_df <- loadRData(file.path(Dir.Exports, fs[i]))$Topo_df
    IS_Iter <- as.numeric(unlist(lapply(
      regmatches(fs[i], gregexpr("[[:digit:]]+\\.*[[:digit:]]*", fs[i])), 
      "[[", 1
    )))
    add <- ifelse(IS_Iter == 585, 1, 0)
    # EffSD2_df$IS <- Eff2_df$IS <- 
    Topo2_df$IS <- as.numeric(unlist(lapply(
      regmatches(fs[i], gregexpr("[[:digit:]]+\\.*[[:digit:]]*", fs[i])),
      "[[", 1+add
    )))
    # EffSD2_df$RE <- Eff2_df$RE <- 
    Topo2_df$RE <- as.numeric(unlist(
      lapply(
        regmatches(fs[i], gregexpr("[[:digit:]]+\\.*[[:digit:]]*",fs[i])),
        "[[", 2+add
      )
    ))
    
    ## difference to pre-extinction
    ## Relative effect
    RelTopo_df <- Topo2_df[Topo2_df$Simulation == "Prediction", ]
    PlotCombin_df <- merge(Pre, RelTopo_df, by = c("netID"))
    Rel_ls <- lapply(Topos_vec, function(x){
      Change_vec <- apply(PlotCombin_df[ , grepl(x,colnames(PlotCombin_df))], 1, diff)
      Plot_df <- data.frame(
        netID = PlotCombin_df$netID,
        Proxy = PlotCombin_df$Proxy.y,
        IS = PlotCombin_df$IS,
        RE = PlotCombin_df$RE,
        Topology = x,
        # Pre = PlotCombin_df[,grep(x,colnames(PlotCombin_df))[1]],
        Post = PlotCombin_df[,grep(x,colnames(PlotCombin_df))[2]],
        AbsChange = Change_vec
      )
      Plot_df$RelChange <- abs(Plot_df$AbsChange)/PlotCombin_df[,grep(x,colnames(PlotCombin_df))[1]]
      Plot_df
    })
    Change2_df <- do.call(rbind, Rel_ls)
    if(i == 1){
      # Eff_df <- Eff2_df
      # EffSD_df <- EffSD2_df
      Topo_df <- Topo2_df
      Change_df <- Change2_df
    }else{
      # Eff_df <- rbind(Eff_df, Eff2_df)
      # EffSD_df <- rbind(EffSD_df, EffSD2_df)
      Topo_df <- rbind(Topo_df, Topo2_df)
      Change_df <- rbind(Change_df, Change2_df)
    }
    setTxtProgressBar(pb, i)
  }
  # colnames(Eff_df) <- gsub(pattern = ".x", replacement = "", colnames(Eff_df))
  return(list(Topology = Topo_df, 
              EffectSize = NULL, #Eff_df, 
              RandomSD = NULL, #EffSD_df, 
              Change = Change_df))
}
