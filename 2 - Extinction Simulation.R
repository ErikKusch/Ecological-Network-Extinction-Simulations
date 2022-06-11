#' ####################################################################### #
#' PROJECT: [Biodiversity Simplification & Ecological Network Topology] 
#' CONTENTS: 
#'  - Extinction Simulation
#'  DEPENDENCIES:
#'  - "0 - Preamble.R"
#'  - "1 - DataRetrieval.R" has to have been run and produced "AnalysesData.RData" in "Dir.Data" directory
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE =================================================================
rm(list=ls())
set.seed(42)

## Sourcing ----------------------------------------------------------------
source("0 - Preamble.R")
source("X - NetworkExtinctionFuns.R")
source("0 - Data_Functions.R")

message("########### STARTING ANALYSIS AND EXTINCTION SIMULATION ###########")

# DATA LOADING & MANIPULATING ==============================================
message("### DATA PREPARATION ###")
## Data Loading ------------------------------------------------------------
print("Loading Data")
load(file.path(Dir.Data, "AnalysesData.RData")) # load the data needed for the analysis
## data now loaded:
#' - List_ls: List object of which each element is a frugivory adjacency matrix with plants as rows and animals as columns
#' - networks_df: Metadata for all networks expressed as adjacency matrices in List_ls (column "net.id" in networks_df can be matched to List_ls names)
#' - Prox.Centrality_ls: List object of which each element is a named vector that sorts the species for each network in our study from most to least central (as measured by connection strength)
#' - Prox.Climate_ls: List object of which each element a list containing a named vector that sorts the species for each network in our study from most to least at risk from a climate-standpoint and a dataframe showing the climate risk for each species for temperature and soil moisture
#' - Prox.IUCN_df: data frame containing IUCN criteria for all species in our analyses.
#' - traits_df: data frame of trait expressions per species in our analysis
#' - animals_gowdis: square matrix of animal species dissimilarity in trait space
#' - plants_gowdis: square matrix of plant species dissimilarity in trait space
#' 

## Data Manipulation -------------------------------------------------------
print("Reformatting Data")
### IUCN Categories as numbers
Prox.IUCN_df$Category <- toupper(substrRight(Prox.IUCN_df$Category, 2)) # reduce categories to just two letters in uppercase
## assigning levels to categories and storing them as numbers in the data frame
IUCN_levels <- factor(Prox.IUCN_df$Category,
                      levels = c(
                        "DD", # data deficient
                        "LC", # least concern
                        "CD", # conservation dependant
                        "NT", # near threatened
                        "VU", # vulnerable
                        "EN", # endangered
                        "CR" # critically endangered
                      )
)
Prox.IUCN_df$CategoryNum <- as.numeric(IUCN_levels)

### Readying objects in List_ls for analyses
AnalysisData_ls <- pblapply(names(List_ls), function(y){
  # extracting adjacency matrix
  x <- List_ls[[y]]
  
  # creating square adjacency matrix (needed for package NetworkExtinction)
  x_df <- expand.grid(colnames(x), rownames(x))
  x_df1 <- data.frame(x)
  x_df$Freq <- apply(X = x_df, MARGIN = 1, FUN = function(z, x){
    x[rownames(x) == as.character(z[2]), colnames(x) == as.character(z[1])]
  }, x = x)
  g <- graph_from_data_frame(x_df, directed = TRUE)
  E(g)$weight <- x_df$Freq
  x_sq <- t(as_adjacency_matrix(g, attr = "weight", sparse = FALSE))
  
  # IUCN sorting
  IUCN_vec <- Sort.DF(Data = Prox.IUCN_df[Prox.IUCN_df$Species %in% unlist(dimnames(x)), ], Column = "CategoryNum", decreasing = TRUE)$CategoryNum
  names(IUCN_vec) <- Sort.DF(Data = Prox.IUCN_df[Prox.IUCN_df$Species %in% unlist(dimnames(x)), ], Column = "CategoryNum", decreasing = TRUE)$Species
  
  list(
    # square adjacency matrix
    Adjacency = x_sq,
    # plant gower distance of traits for the plant species in x
    gow_plants = plants_gowdis[rownames(plants_gowdis) %in% rownames(x), rownames(plants_gowdis) %in% rownames(x)],
    # animal gower distance of traits for the animal species in x
    gow_animals = animals_gowdis[rownames(animals_gowdis) %in% colnames(x), rownames(animals_gowdis) %in% colnames(x)],
    # Centrality
    prox_centrality = Prox.Centrality_ls[[y]],
    # Climate
    prox_climate = Prox.Climate_ls[[y]]$Order,
    # IUCN
    prox_IUCN = IUCN_vec
  )
})
names(AnalysisData_ls) <- names(List_ls)

### Removing networks which have only one row/column of realised interactions
ErrorCheck <- rowSums(
  cbind(
    unlist(pblapply(AnalysisData_ls, function(x){sum(rowSums(x$Adjacency) > 0) > 1})),
    unlist(pblapply(AnalysisData_ls, function(x){sum(colSums(x$Adjacency) > 0) > 1}))
  )
)
AnalysisData_ls <- AnalysisData_ls[ErrorCheck == 2]

# EXTINCTION SIMULATION(S) =================================================
message("### EXTINCTION SIMULATION(S) ###")
plants_sp <- rownames(plants_gowdis)
animals_sp <- rownames(animals_gowdis)

# install.packages("NetworkExtinction")
library(NetworkExtinction)

FUN_SimComp <- function(PlantAnim = NULL, # should be set either to a vector of plant species names or animal species names
                        RunName = "ALL" # for file naming
){
  
  # cl <- parallel::makeCluster(parallel::detectCores()) # for parallel pbapply functions
  # parallel::clusterExport(cl,
  #                         varlist = c('FUN_Topo', "animals_gowdis", "plants_gowdis", "AnalysisData_ls"), 
  #                         envir = environment()
  # )
  
  Sim_ls <- pblapply(names(AnalysisData_ls), 
                     # cl = cl, 
                     function(x){
    # print(x)
    
    x <- AnalysisData_ls[[x]]
    net <- as.network(x$Adjacency)
    
    ## Subsetting proxies for bottom-up/top-down simulations
    if(!is.null(PlantAnim)){
      x$prox_centrality <- x$prox_centrality[names(x$prox_centrality) %in% PlantAnim]
      x$prox_climate <- x$prox_climate[names(x$prox_climate) %in% PlantAnim]
      x$prox_IUCN <- x$prox_IUCN[x$prox_IUCN %in% PlantAnim]
    }
    
    ## Centrality-Driven -------------------------------------------------------
    # print("Extinction of Keystone Species (Centrality)")
    proxcen <- x$prox_centrality[x$prox_centrality > quantile(x$prox_centrality, 0.75)] # just eleiminate upper 25% quantile
    primext_names <- names(proxcen)
    primext_order <- match(primext_names, rownames(x$Adjacency))
    CustOrder_ExtS <- ExtinctionOrderEK(Network = net, Order = primext_order)
    ExtS_Rand <- RandomExtinctionsEK(Network = net, nsim = 100, 
                                     parallel = TRUE, ncores = parallel::detectCores(), 
                                     SimExt = length(proxcen))
    closeAllConnections()
    # CompareExtinctions(Nullmodel = ExtS_Rand[[1]], Hypothesis = CustOrder_ExtS[[1]])
    
    ## Climate-Driven ----------------------------------------------------------
    # print("Extinction of Threatened Species (Climate Projections)")
    primext_names <- names(x$prox_climate)[x$prox_climate > 2] # random cutoff of climate risk severity selected here
    CustOrder_ExtC <- ExtC_Rand <- as.list(c(NA, NA))
    if(length(primext_names) != 0){
      primext_order <- match(primext_names, rownames(x$Adjacency))
      CustOrder_ExtC <- ExtinctionOrderEK(Network = net, Order = primext_order)
      ExtC_Rand <- RandomExtinctionsEK(Network = net, nsim = 100, 
                                       parallel = TRUE, ncores = parallel::detectCores(), 
                                       SimExt = length(primext_names))
      # CompareExtinctions(Nullmodel = Rando_Ext, Hypothesis = CustOrder_ExtC)
      closeAllConnections()
    }
    
    ## IUCN-Driven -------------------------------------------------------------
    # print("Extinction of Threatened Species (IUCN Categories)")
    primext_names <- names(x$prox_IUCN)[x$prox_IUCN > 3] # random cutoff of climate risk severity selected here
    CustOrder_ExtI <- ExtI_Rand <- as.list(c(NA, NA))
    if(length(primext_names) != 0){
      primext_order <- match(primext_names, rownames(x$Adjacency))
      CustOrder_ExtI <- ExtI_Rand <- ExtinctionOrderEK(Network = net, Order = primext_order)
      ExtI_Rand <- RandomExtinctionsEK(Network = net, nsim = 100, 
                                       parallel = TRUE, ncores = parallel::detectCores(), 
                                       SimExt = length(primext_names))
      # CompareExtinctions(Nullmodel = Rando_Ext, Hypothesis = CustOrder_ExtI)
      closeAllConnections()
    }
    # stopCluster(cl)
    
    ## Export ------------------------------------------------------------------
    list(Strength = list(Prediction = as.matrix(CustOrder_ExtS[[2]]),
                         Random = lapply(ExtS_Rand$nets, as.matrix)
    ),
    Climate = list(Prediction = as.matrix(CustOrder_ExtC[[2]]),
                   Random = lapply(ExtC_Rand$nets, as.matrix)
    ),
    IUCN = list(Prediction = as.matrix(CustOrder_ExtI[[2]]),
                Random = lapply(ExtI_Rand$nets, as.matrix)
    )
    )
  })
  names(Sim_ls) <- names(AnalysisData_ls)
  save(Sim_ls, file = file.path(Dir.Exports, paste0(RunName, "SimulationNets.RData")))
  return(Sim_ls)
}

SimComp_ALL <- if(file.exists(file.path(Dir.Exports, "ALLSimulationNets.RData"))){
                  loadRData(file.path(Dir.Exports, "ALLSimulationNets.RData"))
                }else{
                  FUN_SimComp(PlantAnim = NULL, RunName = "ALL")
                }
SimComp_Plants <- if(file.exists(file.path(Dir.Exports, "PlantsSimulationNets.RData"))){
  loadRData(file.path(Dir.Exports, "PlantsSimulationNets.RData"))
}else{
  FUN_SimComp(PlantAnim = plants_sp, RunName = "ALL")
}
SimComp_Animals <- if(file.exists(file.path(Dir.Exports, "AnimalsSimulationNets.RData"))){
  loadRData(file.path(Dir.Exports, "AnimalsSimulationNets.RData"))
}else{
  FUN_SimComp(PlantAnim = animals_sp, RunName = "Animals")
}

# NETWORK TOPOLOGY =========================================================
message("### NETWORK TOPOLOGIES ###")

## Pre-Extinction ----------------------------------------------------------
print("Pre-Extinction")
## Setting up parlallel cluster
PreExt_df <- pblapply(lapply(AnalysisData_ls, "[[", "Adjacency"), 
                      FUN = FUN_Topo)
PreExt_df <- do.call(rbind, PreExt_df)
PreExt_df$netID <- names(AnalysisData_ls)
PreExt_df$Proxy <- "Pre-Extinction"
PreExt_df$Simulation <- "Pre-Extinction"

## Post-Extinction ---------------------------------------------------------
print("Post-Extinction")

FUN_TopoComp <- function(Sim_ls = NULL, RunName = "ALL"){
  # 
  # cl <- parallel::makeCluster(parallel::detectCores()) # for parallel pbapply functions
  # parallel::clusterExport(cl,
  #                         varlist = c('FUN_Topo', "animals_gowdis", "plants_gowdis"), 
  #                         envir = environment()
  # )
  # 
  ## Topology Calculation
  PostExt_ls <- pblapply(names(Sim_ls), 
                         # cl = cl,
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
  Save_ls <- list(Topo_ls = PostExt_ls, Topo_df = do.call(rbind, Topo_ls))

  ## Data Return
  save(Save_ls,
       file = file.path(Dir.Exports, paste0(RunName, "SimulationTopo.RData")))
  return(Save_ls)
}

TopoComp_ALL <- if(file.exists(file.path(Dir.Exports, "ALLSimulationTopo.RData"))){
  loadRData(file.path(Dir.Exports, "ALLSimulationTopo.RData"))
}else{
  FUN_TopoComp(Sim_ls = SimComp_ALL, RunName = "ALL")
}
TopoComp_Plants <- if(file.exists(file.path(Dir.Exports, "PlantsSimulationTopo.RData"))){
  loadRData(file.path(Dir.Exports, "PlantsSimulationTopo.RData"))
}else{
  FUN_TopoComp(Sim_ls = SimComp_Plants, RunName = "Plants")
}
TopoComp_Animals <- if(file.exists(file.path(Dir.Exports, "AnimalsSimulationTopo.RData"))){
  loadRData(file.path(Dir.Exports, "AnimalsSimulationTopo.RData"))
}else{
  FUN_TopoComp(Sim_ls = SimComp_Animals, RunName = "Animals")
}


# VISUALISATION ============================================================
message("### RESULT VISUALISATION ###")

FUN_PlotMod <- function(Pre_df, Post_df){
  Plot_df <- rbind(Pre_df, Post_df)
  Plot_df<- reshape(data = Plot_df,
                    idvar = "netID",
                    varying = colnames(Plot_df)[1:6],
                    v.name = "Value",
                    timevar = "Topology",
                    times = colnames(Plot_df)[1:6],
                    new.row.names = 1:(nrow(Plot_df)*6),
                    direction = "long"
  )
  
  ## PROXY COMPARISON ----
  ProxyComp_ls <- as.list(rep(NA, length(unique(Plot_df$Topology))))
  names(ProxyComp_ls) <- unique(Plot_df$Topology)
  Titles <- c("Number of Species in Network following Extinction",
              "Number of Plant Species in Network following Extinction",
              "Number of Animal Species in Network following Extinction",
              "Number of Realised Links in Network following Extinction",
              "Nestedness of Network following Extinction", 
              "Modularity of Network following Extinction")
  for(i in 1:length(unique(Plot_df$Topology))){
    ## Comparison of random vs. predicted
    comps <- list(c("Random", "Prediction"))
    RandPred_gg <- ggplot(Plot_df[Plot_df$Simulation != "Pre-Extinction" & 
                                    Plot_df$Topology == unique(Plot_df$Topology)[i], ], 
                          aes(x = factor(Simulation, levels = c("Random", "Prediction")), y = Value)) +
      geom_boxplot() +
      facet_wrap(~factor(Proxy, levels = c("Strength", "Climate", "IUCN"))) + 
      stat_compare_means(comparisons = comps, method = 't.test', label = 'p.signif') +
      theme_bw() + labs(x = "Simulation", y = "Network Topology Metric")
    ## Comparison of proxy effects in predictions
    comps <- list(c("Strength", "Climate"),
                  c("Strength", "IUCN"),
                  c("Climate", "IUCN")
    )
    PredComp_gg <- ggplot(Plot_df[Plot_df$Simulation != "Pre-Extinction" & 
                                    Plot_df$Simulation != "Random" & 
                                    Plot_df$Topology == unique(Plot_df$Topology)[i], ], 
                          aes(x = factor(Proxy, levels = c("Strength", "Climate", "IUCN")), y = Value)) +
      geom_boxplot() +
      stat_compare_means(comparisons = comps, method = 't.test', label = 'p.signif') +
      theme_bw() + labs(x = "Proxy", y = "Network Topology Metric")
    ## Fusing of Plots
    title <- ggdraw() + draw_label(Titles[i], fontface='bold')
    ProxyComp_ls[[i]] <- cowplot::plot_grid(title, RandPred_gg, PredComp_gg, ncol = 1, rel_heights = c(0.05, 0.4, 0.5)) 
  }
  ProxyComp_ls
  
  ## NETWORK TOPOLOGIES PRE AND POST EXTINCTION ----
  comps <- list(c("Pre-Extinction", "Strength"),
                c("Pre-Extinction", "Climate"),
                c("Pre-Extinction", "IUCN")
                # ,
                # c("Strength", "Climate"),
                # c("Strength", "IUCN"),
                # c("Climate", "IUCN")
  )
  
  PrePost_gg <- ggplot(Plot_df[Plot_df$Simulation != "Random", ], 
                       aes(x = factor(Proxy, levels = c("Pre-Extinction", "Strength", "Climate", "IUCN")), y = Value)) +
    geom_boxplot() +
    facet_wrap(~factor(Topology, levels = c("n_species", "n_plants", "n_animals",
                                            "n_links", "Nestedness", "Modularity")), 
               scales = "free") + 
    stat_compare_means(comparisons = comps, method = 't.test', label = 'p.signif') +
    theme_bw() + labs(x = "Proxy", y = "Network Topology Metric")
  
  ## RETURN PLOT LIST ----
  return(list(ProxyComp = ProxyComp_ls,
              PrePost = PrePost_gg))
  
}

Plots_ALL <- FUN_PlotMod(Pre_df = PreExt_df, Post_df = TopoComp_ALL$Topo_df)
Plots_Plants <- FUN_PlotMod(Pre_df = PreExt_df, Post_df = TopoComp_Plants$Topo_df)
Plots_Animals <- FUN_PlotMod(Pre_df = PreExt_df, Post_df = TopoComp_Animals$Topo_df)




