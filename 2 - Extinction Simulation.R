#' ######################################################################
#' PROJECT: [Biodiversity Simplification & Ecological Network Topology] 
#' CONTENTS: 
#'  - Extinction Simulation
#'  DEPENDENCIES:
#'  - "0 - Preamble.R"
#'  - "1 - DataRetrieval.R" has to have been run and produced "AnalysesData.RData" in "Dir.Data" directory
#' AUTHOR: [Erik Kusch]
#' ######################################################################

# PREAMBLE ==============================================================
rm(list=ls())
set.seed(42)

## Sourcing -------------------------------------------------------------
source("0 - Preamble.R")
# source("X - NetworkExtinctionFunsRewiring.R")
source("0 - Data_Functions.R")

message("########### STARTING ANALYSIS AND EXTINCTION SIMULATION ###########")

# DATA LOADING & MANIPULATING ===========================================
message("### DATA PREPARATION ###")
## Data Loading ---------------------------------------------------------
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

## Data Manipulation ----------------------------------------------------
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
    prox_climate = ProxClim_ls[["ssp245"]][[y]]$Order,
    prox_climateSSP585 = ProxClim_ls[["ssp585"]][[y]]$Order,
    # IUCN
    prox_IUCN = IUCN_vec
  )
})
names(AnalysisData_ls) <- names(List_ls)

### Removing networks which have only one row/column of realised interactions
print("Checking for characteristics that make networks unusable")
ErrorCheck <- rowSums(
  cbind(
    unlist(pblapply(AnalysisData_ls, function(x){sum(rowSums(x$Adjacency) > 0) > 1})),
    unlist(pblapply(AnalysisData_ls, function(x){sum(colSums(x$Adjacency) > 0) > 1}))
  )
)
AnalysisData_ls <- AnalysisData_ls[ErrorCheck == 2]

# POTENTIAL ASSOCIATIONS ================================================
message("### IDENTIFYING POTENTIAL REWIRING PARTNERS ###")
plants_sp <- rownames(plants_gowdis)
animals_sp <- rownames(animals_gowdis)

## Metaweb --------------------------------------------------------------
# metaweb_mat <- add_matrices(lapply(lapply(AnalysisData_ls, "[[", "Adjacency"), function(x){x > 0}))
# metaweb_mat <- metaweb_mat[ , colnames(metaweb_mat) %in% colnames(animals_gowdis)]
# metaweb_mat <- metaweb_mat[rownames(metaweb_mat) %in% rownames(plants_gowdis), ]
meta_df <- traits_df
meta_df$value <- meta_df$value>0
meta_df <- meta_df[meta_df$animal.phylo.id %in% animals_sp | meta_df$plant.phylo.id %in% plants_sp, ]

## Potential Partner Identification -------------------------------------
options(warn=-1)
### Animals ----------------------------
print("Animals")
if(!file.exists(file.path(Dir.Exports, "RewiringAnimals.RData"))){
  RewClass_Animals <- pblapply(animals_sp, function(x){
    # message(x)
    sub_df <- meta_df[meta_df$animal.phylo.id == x, ]
    reg_df <- sub_df[,c(4, 22:29)] # value, and plant traits
    if(sum(reg_df$value) == 0 | nrow(reg_df) == 1){
      0
    }else{
      if(sum(reg_df$value) == nrow(reg_df)){
        1
      }else{
        invisible(capture.output(
          modRF <- tuneRF(y = reg_df[,1], x = reg_df[,-1], plot = FALSE, doBest = TRUE, trace = FALSE, do.trace = FALSE)
        ))
        modRF
      }
    }
  })
  names(RewClass_Animals) <- animals_sp
  save(RewClass_Animals, file = file.path(Dir.Exports, "RewiringAnimals.RData"))
}else{
  load(file.path(Dir.Exports, "RewiringAnimals.RData"))
}

### Plants ----------------------------
print("Plants")
if(!file.exists(file.path(Dir.Exports, "RewiringPlants.RData"))){
  RewClass_Plants <- pblapply(plants_sp, function(x){
    sub_df <- meta_df[meta_df$plant.phylo.id == x, ]
    reg_df <- sub_df[,c(4, 2, 5:15)] # value, and plant traits
    if(sum(reg_df$value) == 0 | nrow(reg_df) == 1){
      0
    }else{
      if(sum(reg_df$value) == nrow(reg_df)){
        1
      }else{
        invisible(capture.output(
          modRF <- tuneRF(y = reg_df[,1], x = reg_df[,-1], plot = FALSE, doBest = TRUE, trace = FALSE, do.trace = FALSE)
        ))
        modRF
      }
    }
  })
  names(RewClass_Plants) <- plants_sp
  save(RewClass_Plants, file = file.path(Dir.Exports, "RewiringPlants.RData"))
}else{
  load(file.path(Dir.Exports, "RewiringPlants.RData"))
}

### Full List of Models ---------------
RewClass_ls <- c(RewClass_Animals, RewClass_Plants)
options(warn=0)

# PARALLEL EXECTIONS ====================================================
CutOffs <- list(Strength = 0.75,
                Climate = 2,
                IUCN = 5)

message("### REGISTERING CLUSTER")
nCores <- ifelse(parallel::detectCores()>(length(AnalysisData_ls)/2),
                 (length(AnalysisData_ls)/2), parallel::detectCores())
cl <- parallel::makeCluster(as.integer(nCores)+1) # for parallel pbapply functions
parallel::clusterExport(cl,
                        varlist = c('FUN_Topo', "CutOffs",  "AnalysisData_ls", 
                                    "animals_gowdis", "plants_gowdis", "plants_sp", "animals_sp", 
                                    "RewClass_ls", "meta_df",
                                    "install.load.package", "package_vec"
                                    # ,
                                    # ".DataInit", "SimulateExtinctions", "RandomExtinctions", "ExtinctionOrder"
                                    ),
                        envir = environment()
)
clusterpacks <- clusterCall(cl, function() sapply(package_vec, install.load.package))

# PRE-EXCTINCTION =======================================================
message("### PRE-EXTINCTION ###")
print("Network Topologies")
PreExt_df <- pblapply(lapply(AnalysisData_ls, "[[", "Adjacency"),
                      cl = cl,
                      FUN = FUN_Topo)
PreExt_df <- do.call(rbind, PreExt_df)
PreExt_df$netID <- names(AnalysisData_ls)
PreExt_df$Proxy <- "Pre-Extinction"
PreExt_df$Simulation <- "Pre-Extinction"
parallel::clusterExport(cl,
                        varlist = c("PreExt_df"),
                        envir = environment()
)

# POST-EXCTINCTION ======================================================
message("### EXTINCTION SIMULATION(S) ###")

if(all(unlist(
  lapply(list(file.path(Dir.Exports, "PlotTopoAll_ls.RData"),
              file.path(Dir.Exports, "PlotTopoPlants_ls.RData"),
              file.path(Dir.Exports, "PlotTopoAnimals_ls.RData"),
              file.path(Dir.Exports, "PlotTopoClimSSP585_ls.RData")
  ),
  file.exists)
))){
  print("Simulations and topologies already calculated - Loading 4 files from hard drive")
  PlotTopoAll_ls <- loadObj(file.path(Dir.Exports, "PlotTopoAll_ls.RData"))
  PlotTopoPlants_ls <- loadObj(file.path(Dir.Exports, "PlotTopoPlants_ls.RData"))
  PlotTopoAnimals_ls <- loadObj(file.path(Dir.Exports, "PlotTopoAnimals_ls.RData"))
  PlotTopoClimSSP585_ls <- loadObj(file.path(Dir.Exports, "PlotTopoClimSSP585_ls.RData"))
}else{
  ## Extinction Simulations ---------------------------------------------
  for(Rewiring_Iter in seq(0,  1, 0.05)){
    for(IS_iter in seq(0, 1, 0.05)){
      Sim_ls <- FUN_SimComp(PlantAnim = NULL, RunName = "ALL",
                            IS = IS_iter, Rewiring = Rewiring_Iter,
                            CutOffs = CutOffs, PotPartners = RewClass_ls, Traits = meta_df)
      TopoComp_ls <- FUN_TopoComp(Sim_ls = Sim_ls, RunName = "ALL",
                                  IS = IS_iter, Rewiring = Rewiring_Iter,
                                  CutOffs = CutOffs, Pre = PreExt_df)
      Sim_ls <- FUN_SimComp(PlantAnim = plants_sp, RunName = "Plants",
                            IS = IS_iter, Rewiring = Rewiring_Iter,
                            CutOffs = CutOffs, PotPartners = RewClass_ls, Traits = meta_df)
      TopoComp_ls <- FUN_TopoComp(Sim_ls = Sim_ls, RunName = "Plants",
                                  IS = IS_iter, Rewiring = Rewiring_Iter,
                                  CutOffs = CutOffs)
      Sim_ls <- FUN_SimComp(PlantAnim = animals_sp, RunName = "Animals",
                            IS = IS_iter, Rewiring = Rewiring_Iter,
                            CutOffs = CutOffs, PotPartners = RewClass_ls, Traits = meta_df)
      TopoComp_ls <- FUN_TopoComp(Sim_ls = Sim_ls, RunName = "Animals",
                                  IS = IS_iter, Rewiring = Rewiring_Iter,
                                  CutOffs = CutOffs)
      Sim_ls <- FUN_SimComp(PlantAnim = NULL, RunName = "SSP585", WHICH = "SSP585",
                            IS = IS_iter, Rewiring = Rewiring_Iter,
                            CutOffs = CutOffs, PotPartners = RewClass_ls, Traits = meta_df)
      TopoComp_ls <- FUN_TopoComp(Sim_ls = Sim_ls, RunName = "SSP585",
                                  IS = IS_iter, Rewiring = Rewiring_Iter, CutOffs = CutOffs)
    }
  }

  ## Topology Loading and Storing as one object -------------------------
  ## while loading in the topologies, we also compute absolute and relative change of each simulation to the pre-extinction network topologies
  PlotTopoAll_ls <- loadTopo(RunName = "ALL", CutOffs = CutOffs, Pre = PreExt_df)
  saveObj(PlotTopoAll_ls, file.name = file.path(Dir.Exports, "PlotTopoAll_ls.RData"))
  PlotTopoPlants_ls <- loadTopo(RunName = "Plants", CutOffs = CutOffs, Pre = PreExt_df)
  saveObj(PlotTopoPlants_ls, file.name = file.path(Dir.Exports, "PlotTopoPlants_ls.RData"))
  PlotTopoAnimals_ls <- loadTopo(RunName = "Animals", CutOffs = CutOffs, Pre = PreExt_df)
  saveObj(PlotTopoAnimals_ls, file.name = file.path(Dir.Exports, "PlotTopoAnimals_ls.RData"))
  PlotTopoAnimals_ls <- loadTopo(RunName = "SSP585", CutOffs = CutOffs, Pre = PreExt_df)
  saveObj(PlotTopoAnimals_ls, file.name = file.path(Dir.Exports, "PlotTopoClimSSP585_ls.RData"))
}
