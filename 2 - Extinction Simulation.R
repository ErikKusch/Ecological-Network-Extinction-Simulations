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

# EXTINCTION SIMULATION(S) =================================================
message("### EXTINCTION SIMULATION(S) ###")

animals_sp <- rownames(animals_gowdis)
plants_sp <- rownames(plants_gowdis)

PlantAnim <- NULL # should be set either to a vector of plant species names or animal species names
RunName <- "ALL" # for file naming

# install.packages("NetworkExtinction")
library(NetworkExtinction)

Sim_ls <- pblapply(names(AnalysisData_ls), function(x){
  print(x)
  
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
                                   parallel = TRUE, ncores = 4, 
                                   SimExt = length(proxcen))
  # CompareExtinctions(Nullmodel = ExtS_Rand[[1]], Hypothesis = CustOrder_ExtS[[1]])
  
  ## Climate-Driven ----------------------------------------------------------
  # print("Extinction of Threatened Species (Climate Projections)")
  primext_names <- names(x$prox_climate)[x$prox_climate > 2] # random cutoff of climate risk severity selected here
  CustOrder_ExtC <- as.list(c(NA, NA))
  if(length(primext_names) != 0){
    primext_order <- match(primext_names, rownames(x$Adjacency))
    CustOrder_ExtC <- ExtinctionOrderEK(Network = net, Order = primext_order)
    # CompareExtinctions(Nullmodel = Rando_Ext, Hypothesis = CustOrder_ExtC)
  }
  ExtC_Rand <- RandomExtinctionsEK(Network = net, nsim = 100, 
                                   parallel = TRUE, ncores = 4, 
                                   SimExt = length(primext_names))
  
  
  ## IUCN-Driven -------------------------------------------------------------
  # print("Extinction of Threatened Species (IUCN Categories)")
  primext_names <- names(x$prox_IUCN)[x$prox_IUCN > 3] # random cutoff of climate risk severity selected here
  CustOrder_ExtI <- as.list(c(NA, NA))
  if(length(primext_names) != 0){
    primext_order <- match(primext_names, rownames(x$Adjacency))
    CustOrder_ExtI <- ExtinctionOrderEK(Network = net, Order = primext_order)
    # CompareExtinctions(Nullmodel = Rando_Ext, Hypothesis = CustOrder_ExtI)
  }
  ExtI_Rand <- RandomExtinctionsEK(Network = net, nsim = 100, 
                                   parallel = TRUE, ncores = 4, 
                                   SimExt = length(primext_names))
  
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

# NETWORK TOPOLOGY =========================================================
message("### NETWORK TOPOLOGIES ###")

FUN_Topo <- function(plot_df){
  Nes <- try(networklevel(web = plot_df, index = "weighted nestedness"), silent = TRUE)
  Mod <- try(NOS(web = plot_df)$mod, silent = TRUE)
  round(
    data.frame(
      n_species = ncol(plot_df),
      n_animals = sum(colnames(plot_df) %in% rownames(animals_gowdis)),
      n_plants = sum(colnames(plot_df) %in% rownames(plants_gowdis)),
      n_links = sum(plot_df>0),
      Nestedness = as.numeric(ifelse(class(Nes) == "try-error", NA, Nes)),
      Modularity = as.numeric(ifelse(class(Mod) == "try-error", NA, Mod))
    ), 2)
}

## Pre-Extinction ----------------------------------------------------------
print("Pre-Extinction")
PreExt_df <- pblapply(lapply(AnalysisData_ls, "[[", "Adjacency"), FUN = FUN_Topo)
PreExt_df <- do.call(rbind, PreExt_df)
rownames(PreExt_df) <- names(AnalysisData_ls)

## Post-Extinction ---------------------------------------------------------
print("Post-Extinction")
PostExt_ls <- pblapply(names(Sim_ls), FUN = function(netID){
  Storage_ls <- list(Random = NA,
                     Strength = NA,
                     Climate = NA,
                     IUCN = NA
  )
  for(i in 1:4){
    if(i == 1){
      plot_ls <- as.matrix(Sim_ls[[netID]][[i]])
      Rand_ls <- lapply(plot_ls, FUN_Topo)
      Storage_ls[[i]] <- do.call(rbind, Rand_ls)
    }else{
      Storage_ls[[i]] <- FUN_Topo(as.matrix(Sim_ls[[netID]][[i]]))
    }
  }
  Storage_ls
})
names(PostExt_ls) <- names(AnalysisData_ls)

PostRand_df <- do.call(rbind, lapply(PostExt_ls, "[[", "Random"))
PostStr_df <- do.call(rbind, lapply(PostExt_ls, "[[", "Strength"))
PostClim_df <- do.call(rbind, lapply(PostExt_ls, "[[", "Climate"))
PostIUCN_df <- do.call(rbind, lapply(PostExt_ls, "[[", "IUCN"))


# VISUALISATION ============================================================
Plot_df <- rbind(PreExt_df, PostStr_df, PostClim_df, PostIUCN_df, PostRand_df)
Plot_df$Study <- c(rep(rownames(PostIUCN_df), 4), 
                   rep(rownames(PostIUCN_df), each = nrow(PostRand_df)/length(AnalysisData_ls)))
Plot_df$Sim <- c(rep(c("Pre", "Centrality", "Climate", "IUCN"), each = nrow(PostIUCN_df)), 
                 rep("Random", nrow(PostRand_df)))

Plot_df <- Plot_df[Plot_df$Sim != "Random", ]

comps <- list(
  # c("Random", "Pre"),
  #             c("Random", "Centrality"),
  #             c("Random", "Climate"),
  #             c("Random", "IUCN"),
              c("Pre", "Centrality"),
              c("Pre", "Climate"),
              c("Pre", "IUCN"),
              c("Centrality", "Climate"),
              c("Centrality", "IUCN"),
              c("Climate", "IUCN")
              )

N_gplot <- ggplot(Plot_df, aes(x = factor(Sim, levels = c("Pre", "Centrality", "Climate", "IUCN", "Random")), y = n_species)) + 
  geom_boxplot() + 
  stat_compare_means(comparisons = comps, method = 't.test', label = 'p.signif') + 
  theme_bw() + labs(x = "Simulation / Proxy", y = "# of species")

NA_gplot <- ggplot(Plot_df, aes(x = factor(Sim, levels = c("Pre", "Centrality", "Climate", "IUCN", "Random")), y = n_animals)) + 
  geom_boxplot() + 
  stat_compare_means(comparisons = comps, method = 't.test', label = 'p.signif') + 
  theme_bw() + labs(x = "Simulation / Proxy", y = "# of animal species")

NP_gplot <- ggplot(Plot_df, aes(x = factor(Sim, levels = c("Pre", "Centrality", "Climate", "IUCN", "Random")), y = n_plants)) + 
  geom_boxplot() + 
  stat_compare_means(comparisons = comps, method = 't.test', label = 'p.signif') + 
  theme_bw() + labs(x = "Simulation / Proxy", y = "# of plant species")

Links_gplot <- ggplot(Plot_df, aes(x = factor(Sim, levels = c("Pre", "Centrality", "Climate", "IUCN", "Random")), y = n_links)) + 
  geom_boxplot() + 
  stat_compare_means(comparisons = comps, method = 't.test', label = 'p.signif') + 
  theme_bw() + labs(x = "Simulation / Proxy", y = "# of links")

Nest_gplot <- ggplot(Plot_df, aes(x = factor(Sim, levels = c("Pre", "Centrality", "Climate", "IUCN", "Random")), y = Nestedness)) + 
  geom_boxplot() + 
  stat_compare_means(comparisons = comps, method = 't.test', label = 'p.signif') + 
  theme_bw() + labs(x = "Simulation / Proxy", y = "Nestedness")

Mod_gplot <- ggplot(Plot_df, aes(x = factor(Sim, levels = c("Pre", "Centrality", "Climate", "IUCN", "Random")), y = Modularity)) + 
  geom_boxplot() + 
  stat_compare_means(comparisons = comps, method = 't.test', label = 'p.signif') + 
  theme_bw() + labs(x = "Simulation / Proxy", y = "Modularity")

plot_grid(N_gplot, Links_gplot, NA_gplot, NP_gplot, Nest_gplot, Mod_gplot, ncol = 2)
