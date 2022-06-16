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
print("Checking for characteristics that make networks unusable")
ErrorCheck <- rowSums(
  cbind(
    unlist(pblapply(AnalysisData_ls, function(x){sum(rowSums(x$Adjacency) > 0) > 1})),
    unlist(pblapply(AnalysisData_ls, function(x){sum(colSums(x$Adjacency) > 0) > 1}))
  )
)
AnalysisData_ls <- AnalysisData_ls[ErrorCheck == 2]

# PARALLEL EXECTIONS =======================================================
plants_sp <- rownames(plants_gowdis)
animals_sp <- rownames(animals_gowdis)

message("### REGISTERING CLUSTER")
nCores <- ifelse(parallel::detectCores()>length(AnalysisData_ls), 
                 length(AnalysisData_ls), parallel::detectCores())
cl <- parallel::makeCluster(nCores) # for parallel pbapply functions
parallel::clusterExport(cl,
                        varlist = c('FUN_Topo', "animals_gowdis", "plants_gowdis", ".ExtinctionOrderEK", "ExtinctionOrderEK", "RandomExtinctionsEK", "AnalysisData_ls", "install.load.package", "package_vec", ".DataInit", "plants_sp", "animals_sp"),
                        envir = environment()
)
clusterpacks <- clusterCall(cl, function() sapply(package_vec, install.load.package))

# PRE-EXCTINCTION ==========================================================
message("### PRE-EXTINCTION ###")
print("Network Topologies")
PreExt_df <- pblapply(lapply(AnalysisData_ls, "[[", "Adjacency"), 
                      cl = cl,
                      FUN = FUN_Topo)
PreExt_df <- do.call(rbind, PreExt_df)
PreExt_df$netID <- names(AnalysisData_ls)
PreExt_df$Proxy <- "Pre-Extinction"
PreExt_df$Simulation <- "Pre-Extinction"

# POST-EXCTINCTION =========================================================
message("### EXTINCTION SIMULATION(S) ###")

for(IS_iter in seq(0, 1, 0.05)){
  Sim_ls <- FUN_SimComp(PlantAnim = NULL, RunName = "ALL", IS = IS_iter)
  TopoComp_ls <- FUN_TopoComp(Sim_ls = Sim_ls, RunName = "ALL", IS = IS_iter)
  Sim_ls <- FUN_SimComp(PlantAnim = plants_sp, RunName = "Plants", IS = IS_iter)
  TopoComp_ls <- FUN_TopoComp(Sim_ls = Sim_ls, RunName = "Plants", IS = IS_iter)
  Sim_ls <- FUN_SimComp(PlantAnim = animals_sp, RunName = "Animals", IS = IS_iter)
  TopoComp_ls <- FUN_TopoComp(Sim_ls = Sim_ls, RunName = "Animals", IS = IS_iter)
}

# VISUALISATION ============================================================
message("### RESULT VISUALISATION ###")

RunName = "ALL"

fs <- list.files(path = Dir.Exports, pattern = paste0(RunName, "SimulationTopo"))
IS_vec <- as.numeric(unlist(regmatches(fs,
                             gregexpr("[[:digit:]]+\\.*[[:digit:]]*",fs))
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
Plot_df<- reshape(data = Eff_df,
                  idvar = "netID",
                  varying = colnames(Eff_df)[1:6],
                  v.name = "Value",
                  timevar = c("Topology"),
                  times = colnames(Eff_df)[1:6],
                  new.row.names = 1:(nrow(Eff_df)*6),
                  direction = "long"
)

Eff_gg <- ggplot(Plot_df, aes(y = Value, x = IS, col = Pry, fill = Pry)) + 
  geom_point(alpha = 0.2) + 
  stat_smooth(method = "lm") + 
  facet_wrap(~ factor(Topology, levels = c("n_species", "n_plants", "n_animals",
                                           "n_links", "Nestedness", "Modularity")), 
             scales = "free", ncol = 3) + 
  theme_bw() + 
  labs(x = "Interaction Strength % Required For Survival", y = "Effect Size")

Eff_ls <- lapply(unique(Plot_df$Pry), FUN = function(x){
  ggplot(Plot_df, aes(y = Value, x = factor(IS))) + 
      geom_boxplot() +
      facet_wrap(~factor(Topology, levels = c("n_species", "n_plants", "n_animals",
                                              "n_links", "Nestedness", "Modularity")),
                 scales = "free", ncol = 3) +
      theme_bw() +
      labs(x = "Interaction Strength % Required For Survival", y = "Effect Size", title = x)
})

Diff_df <- Topo_df[Topo_df$Simulation == "Prediction", 1:6] - PreExt_df[match(Topo_df[Topo_df$Simulation == "Prediction",]$netID, PreExt_df$netID), 1:6]
Diff_df$IS <- Topo_df[Topo_df$Simulation == "Prediction",]$IS
Diff_df$Proxy <- Topo_df[Topo_df$Simulation == "Prediction",]$Proxy

Plot_df<- reshape(data = Diff_df,
                  idvar = "netID",
                  varying = colnames(Diff_df)[1:6],
                  v.name = "Value",
                  timevar = "Topology",
                  times = colnames(Diff_df)[1:6],
                  new.row.names = 1:(nrow(Diff_df)*6),
                  direction = "long"
)

Diff_gg <- ggplot(Plot_df, aes(y = Value, x = IS, col = Proxy, fill = Proxy)) + 
  geom_point(alpha = 0.2) + 
  geom_smooth() + 
  facet_wrap(~ factor(Topology, levels = c("n_species", "n_plants", "n_animals",
                                           "n_links", "Nestedness", "Modularity")), 
             scales = "free", ncol = 3) + 
  theme_bw() + 
  labs(x = "Interaction Strength % Required For Survival", y = "Change In Network Topology Metric")


Diff_ls <- lapply(unique(Plot_df$Proxy), FUN = function(x){
  ggplot(Plot_df[Plot_df$Proxy == x, ], aes(y = Value, x = factor(IS))) +
    geom_boxplot() +
    facet_wrap(~ factor(Topology, levels = c("n_species", "n_plants", "n_animals",
                                             "n_links", "Nestedness", "Modularity")),
               scales = "free", ncol = 3) +
    theme_bw() +
    labs(x = "Interaction Strength % Required For Survival", y = "Change In Network Topology Metric", title = x)
})




pdf("Print.pdf", paper = "a4r", height = 9, width = 12)
Eff_gg
Eff_ls
Diff_gg
Diff_ls
dev.off()

# stop("work out new visualisations")
# ## By Cascade Orientation --------------------------------------------------
# FUN_PlotMod <- function(Pre_df, Post_df, RunName){
#   Plot_df <- rbind(Pre_df, Post_df)
#   Plot_df<- reshape(data = Plot_df,
#                     idvar = "netID",
#                     varying = colnames(Plot_df)[1:6],
#                     v.name = "Value",
#                     timevar = "Topology",
#                     times = colnames(Plot_df)[1:6],
#                     new.row.names = 1:(nrow(Plot_df)*6),
#                     direction = "long"
#   )
#   
#   ## PROXY COMPARISON ----
#   ProxyComp_ls <- as.list(rep(NA, length(unique(Plot_df$Topology))))
#   names(ProxyComp_ls) <- unique(Plot_df$Topology)
#   Titles <- c("Number of Species in Network following Extinction",
#               "Number of Plant Species in Network following Extinction",
#               "Number of Animal Species in Network following Extinction",
#               "Number of Realised Links in Network following Extinction",
#               "Nestedness of Network following Extinction", 
#               "Modularity of Network following Extinction")
#   Colours <- c("black", "green", "red", "purple", "blue", "orange")
#   for(i in 1:length(unique(Plot_df$Topology))){
#     ## Comparison of random vs. predicted
#     comps <- list(c("Random", "Prediction"))
#     RandPred_gg <- ggplot(Plot_df[Plot_df$Simulation != "Pre-Extinction" & 
#                                     Plot_df$Topology == unique(Plot_df$Topology)[i], ], 
#                           aes(x = factor(Simulation, levels = c("Random", "Prediction")), y = Value)) +
#       geom_boxplot(col = Colours[i]) +
#       facet_wrap(~factor(Proxy, levels = c("Strength", "Climate", "IUCN"))) + 
#       stat_compare_means(comparisons = comps, method = 't.test', label = 'p.signif') +
#       theme_bw() + labs(x = "Simulation", y = "Network Topology Metric")
#     ## Comparison of proxy effects in predictions
#     comps <- list(c("Strength", "Climate"),
#                   c("Strength", "IUCN"),
#                   c("Climate", "IUCN")
#     )
#     PredComp_gg <- ggplot(Plot_df[Plot_df$Simulation != "Pre-Extinction" & 
#                                     Plot_df$Simulation != "Random" & 
#                                     Plot_df$Topology == unique(Plot_df$Topology)[i], ], 
#                           aes(x = factor(Proxy, levels = c("Strength", "Climate", "IUCN")), y = Value)) +
#       geom_boxplot(col = Colours[i]) +
#       stat_compare_means(comparisons = comps, method = 't.test', label = 'p.signif') +
#       theme_bw() + labs(x = "Proxy", y = "Network Topology Metric")
#     ## Fusing of Plots
#     title <- ggdraw() + draw_label(paste0(Titles[i], " (", RunName, ")"), fontface='bold')
#     ProxyComp_ls[[i]] <- cowplot::plot_grid(title, RandPred_gg, PredComp_gg, ncol = 1, rel_heights = c(0.05, 0.5, 0.5)) 
#   }
#   ProxyComp_ls
#   
#   ## NETWORK TOPOLOGIES PRE AND POST EXTINCTION ----
#   comps <- list(c("Pre-Extinction", "Strength"),
#                 c("Pre-Extinction", "Climate"),
#                 c("Pre-Extinction", "IUCN")
#                 # ,
#                 # c("Strength", "Climate"),
#                 # c("Strength", "IUCN"),
#                 # c("Climate", "IUCN")
#   )
#   
#   PrePost_gg <- ggplot(Plot_df[Plot_df$Simulation != "Random", ], 
#                        aes(x = factor(Proxy, levels = c("Pre-Extinction", "Strength", "Climate", "IUCN")), y = Value)) +
#     geom_boxplot() +
#     facet_wrap(~factor(Topology, levels = c("n_species", "n_plants", "n_animals",
#                                             "n_links", "Nestedness", "Modularity")), 
#                scales = "free") + 
#     stat_compare_means(comparisons = comps, method = 't.test', label = 'p.signif') +
#     theme_bw() + labs(x = "Proxy", y = "Network Topology Metric")
#   
#   ## RETURN PLOT LIST ----
#   return(list(ProxyComp = ProxyComp_ls,
#               PrePost = PrePost_gg))
#   
# }
# 
# Plots_ALL <- FUN_PlotMod(Pre_df = PreExt_df, Post_df = TopoComp_ALL$Topo_df, RunName = "ALL")
# Plots_Plants <- FUN_PlotMod(Pre_df = PreExt_df, Post_df = TopoComp_Plants$Topo_df, RunName = "Plants")
# Plots_Animals <- FUN_PlotMod(Pre_df = PreExt_df, Post_df = TopoComp_Animals$Topo_df, RunName = "Animals")
# 
# 
# ## By Extinction Proxy -----------------------------------------------------
# 
# 
# 
# 
# 
# PreExt_df$Cascade <- "Pre-Extinction"
# TopoComp_ALL$Topo_df$Cascade <- "ALL"
# TopoComp_Plants$Topo_df$Cascade <- "Bottom-Up"
# TopoComp_Animals$Topo_df$Cascade <- "Top-Down"
# 
# 
# 
# 
# Plot_df <- rbind(PreExt_df,
#                  TopoComp_ALL$Topo_df,
#                  TopoComp_Plants$Topo_df, 
#                  TopoComp_Animals$Topo_df
# )
# Plot_df <- reshape(data = Plot_df,
#                    idvar = "netID",
#                    varying = colnames(Plot_df)[1:6],
#                    v.name = "Value",
#                    timevar = "Topology",
#                    times = colnames(Plot_df)[1:6],
#                    new.row.names = 1:(nrow(Plot_df)*6),
#                    direction = "long"
# )
# 
# comps <- list(c("Pre-Extinction", "ALL"),
#               c("Pre-Extinction", "Bottom-Up"),
#               c("Pre-Extinction", "Top-Down")
# )
# 
# Plots_ExtProxy <- as.list(rep(NA, length(unique(Plot_df$Proxy))-1))
# names(Plots_ExtProxy) <- unique(Plot_df$Proxy)[-1]
# 
# for(i in 1:length(Plots_ExtProxy)){
#   k <- unique(Plot_df$Proxy)[-1][i]
#   Plots_ExtProxy[[k]] <- ggplot(Plot_df[Plot_df$Simulation != "Random" & (Plot_df$Proxy == k | Plot_df$Proxy == "Pre-Extinction"), ], 
#                                 aes(x = factor(Cascade, levels = c("Pre-Extinction", "ALL", "Bottom-Up", "Top-Down")), 
#                                     y = Value)) + 
#     geom_violin(col = "brown") + 
#     facet_wrap(~factor(Topology, levels = c("n_species", "n_plants", "n_animals",
#                                             "n_links", "Nestedness", "Modularity")), 
#                scales = "free") + 
#     stat_compare_means(comparisons = comps, method = 't.test', label = 'p.signif') +
#     theme_bw() + labs(x = "Extinction Cascade", title = k)
# }
# Plots_ExtProxy