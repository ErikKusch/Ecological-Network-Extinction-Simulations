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
source("X - NetworkExtinctionFunsRewiring.R")
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

# POTENTIAL ASSOCIATIONS ===================================================
message("### IDENTIFYING POTENTIAL REWIRING PARTNERS ###")
plants_sp <- rownames(plants_gowdis)
animals_sp <- rownames(animals_gowdis)

## Metaweb -----------------------------------------------------------------
# metaweb_mat <- add_matrices(lapply(lapply(AnalysisData_ls, "[[", "Adjacency"), function(x){x > 0}))
# metaweb_mat <- metaweb_mat[ , colnames(metaweb_mat) %in% colnames(animals_gowdis)]
# metaweb_mat <- metaweb_mat[rownames(metaweb_mat) %in% rownames(plants_gowdis), ]
meta_df <- traits_df
meta_df$value <- meta_df$value>0
meta_df <- meta_df[meta_df$animal.phylo.id %in% animals_sp | meta_df$plant.phylo.id %in% plants_sp, ]

## Potential Partner Identification ----------------------------------------
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

# PARALLEL EXECTIONS =======================================================
CutOffs <- list(Strength = 0.75,
                Climate = 2,
                IUCN = 5)

message("### REGISTERING CLUSTER")
nCores <- ifelse(parallel::detectCores()>length(AnalysisData_ls),
                 length(AnalysisData_ls), parallel::detectCores())
cl <- parallel::makeCluster(nCores) # for parallel pbapply functions
parallel::clusterExport(cl,
                        varlist = c('FUN_Topo', "CutOffs",  "AnalysisData_ls", 
                                    "animals_gowdis", "plants_gowdis", "plants_sp", "animals_sp", 
                                    "RewClass_ls", "meta_df",
                                    "install.load.package", "package_vec", 
                                    ".DataInit", "SimulateExtinctions", "RandomExtinctions", "ExtinctionOrder"),
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
parallel::clusterExport(cl,
                        varlist = c("PreExt_df"),
                        envir = environment()
)

# POST-EXCTINCTION =========================================================
message("### EXTINCTION SIMULATION(S) ###")

if(all(unlist(
  lapply(list(file.path(Dir.Exports, "PlotTopoAll_ls.RData")
              # ,
              # file.path(Dir.Exports, "PlotTopoPlants_ls.RData"), 
              # file.path(Dir.Exports, "PlotTopoAnimals_ls.RData")
  ),
  file.exists)
))){
  print("Simulations and topologies already calculated - Loading 3 files from hard drive")
  PlotTopoAll_ls <- loadObj(file.path(Dir.Exports, "PlotTopoAll_ls.RData"))
  # PlotTopoPlants_ls <- loadObj(file.path(Dir.Exports, "PlotTopoPlants_ls.RData"))
  # PlotTopoAnimals_ls <- loadObj(file.path(Dir.Exports, "PlotTopoAnimals_ls.RData"))
}else{
  ## Extinction Simulations --------------------------------------------------
  for(Rewiring_Iter in seq(0, 1, 0.05)){
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
    }
  }
  
  # ## Sensitivity Analysis ----------------------------------------------------
  # message("Sensitivity analysis for WHICH = 'Strength' in FUN_SimComp.")
  
  ## Topology Loading and Storing as one object ------------------------------
  ## while loading in the topologies, we also compute absolute and relative change of each simulation to the pre-extinction network topologies
  PlotTopoAll_ls <- loadTopo(RunName = "ALL", CutOffs = CutOffs, Pre = PreExt_df)
  saveObj(PlotTopoAll_ls, file.name = file.path(Dir.Exports, "PlotTopoAll_ls.RData"))
  PlotTopoPlants_ls <- loadTopo(RunName = "Plants", CutOffs = CutOffs, Pre = PreExt_df)
  saveObj(PlotTopoPlants_ls, file.name = file.path(Dir.Exports, "PlotTopoPlants_ls.RData"))
  PlotTopoAnimals_ls <- loadTopo(RunName = "Animals", CutOffs = CutOffs, Pre = PreExt_df)
  saveObj(PlotTopoAnimals_ls, file.name = file.path(Dir.Exports, "PlotTopoAnimals_ls.RData"))
}

# VISUALISATION ============================================================
message("### RESULT VISUALISATION ###")
pal_lm <- c("#016392", "#E19825", "#3E8853")

## Venn-Diagram of Proxy Agreement -----------------------------------------
print("Extinction Proxy Overlap ---")
for(RunName in c("ALL", "Plants", "Animals")
){
  ExtSpecies_ls <- FUN_SimComp(PlantAnim = NULL, RunName = RunName, IS = 1, Rewiring = 1, CutOffs = CutOffs)
  
  # ## number of primary extinctions per proxy
  # Clim_ls <- lapply(lapply(ExtSpecies_ls, "[[", "Climate"), "[[", "Removed")
  # IUCN_ls <- lapply(lapply(ExtSpecies_ls, "[[", "IUCN"), "[[", "Removed")
  # Centr_ls <- lapply(lapply(ExtSpecies_ls, "[[", "Strength"), "[[", "Removed")
  # ClimIUCN_ls <- lapply(1:length(Clim_ls), function(x){sum(Clim_ls[[x]] %in% IUCN_ls[[x]])})
  # ClimCentr_ls <- lapply(1:length(Clim_ls), function(x){sum(Clim_ls[[x]] %in% Centr_ls[[x]])})
  # IUCNCentr_ls <- lapply(1:length(Clim_ls), function(x){sum(IUCN_ls[[x]] %in% Centr_ls[[x]])})
  # All_ls <- lapply(1:length(Clim_ls), function(x){
  #   sum(
  #     IUCN_ls[[x]][IUCN_ls[[x]] %in% Centr_ls[[x]]] %in% Clim_ls[[x]]
  #     )
  #   })
  # ## absolute numbers of primary extinctions per network
  # PrimaryExt_df <- data.frame(
  #   Climate = unlist(lapply(Clim_ls, length)),
  #   IUCN = unlist(lapply(IUCN_ls, length)),
  #   Centrality = unlist(lapply(Centr_ls, length)),
  #   `Climate+IUCN` = unlist(ClimIUCN_ls),
  #   `Climate+Centrality` = unlist(ClimCentr_ls),
  #   `IUCN+Centrality` = unlist(IUCNCentr_ls),
  #   All = unlist(All_ls)
  # )
  # ## relative numbers of primary extinctions per network
  # PrimaryExt_df <- apply(X = PrimaryExt_df, MARGIN = 2, function(x){
  #   x/PreExt_df$n_species
  # })
  # 
  # ## summary statistics
  # PrimExt_df <- data.frame(min = apply(PrimaryExt_df, 2, min),
  #                          mean = apply(PrimaryExt_df, 2, mean),
  #                          max = apply(PrimaryExt_df, 2, max),
  #                          sd = apply(PrimaryExt_df, 2, sd))
  # print(PrimExt_df)
  
  ## Total Venn-Diagram
  Venn_ls <- list(Climate = unlist(lapply(lapply(ExtSpecies_ls, "[[", "Climate"), "[[", "Removed")),
                  IUCN = unlist(lapply(lapply(ExtSpecies_ls, "[[", "IUCN"), "[[", "Removed")),
                  Centrality = unlist(lapply(lapply(ExtSpecies_ls, "[[", "Strength"), "[[", "Removed")))
  ggvenn(Venn_ls, fill_color = pal_lm, fill_alpha = 0.8, text_color = "white")
  ggsave(filename = file.path(Dir.Exports, paste0("PLOT_Proxy", RunName,".png")), width = 4, height = 3, units = "cm", scale = 7, dpi = 1e3)
}

## Grid-Based Visualisations -----------------------------------------------
print("Grid-based visualisations ---")
PlotTopo_ls <- list(ALL = PlotTopoAll_ls
                    # ,
                    # Animals = PlotTopoAnimals_ls,
                    # Plants = PlotTopoPlants_ls
)

### Effects relative to pre-extinction network topology ----
for(RunName in names(PlotTopo_ls)){
  ### Relative Changes ----
  Change_df <- PlotTopo_ls[[RunName]]$Change
  TopoPlots <- c("n_species", "n_links", "Modularity", "Nestedness") # "n_species", "n_animals", "n_plants", "Modularity", "n_links", "Nestedness"
  Change_df <- Change_df[Change_df$Topology %in% TopoPlots & Change_df$Proxy %in% c("Climate", "IUCN", "Strength"), ]
  
  plot_df <- aggregate(RelChange ~ Proxy+Topology+IS+RE, FUN = mean, data = Change_df[Change_df$RelChange != 0, ])
  sd_df <- aggregate(RelChange ~ Proxy+Topology+IS+RE, FUN = sd, data = Change_df[Change_df$RelChange != 0, ])
  
  # iterate my Topology, then fuse plots
  
  matplot_ls <- list(Pred = as.list(rep(NA, length(TopoPlots))),
                     SD = as.list(rep(NA, length(TopoPlots))))
  
  
  for(TopoIter in 1:length(TopoPlots)){
    plot_df2 <- plot_df[plot_df$Topology == TopoPlots[TopoIter], ]
    sd_df2 <- sd_df[sd_df$Topology == TopoPlots[TopoIter], ]
    
    matplot_ls$Pred[[TopoIter]] <- ggplot(plot_df2, aes(x = RE, y = IS)) +
      geom_tile(aes(fill = RelChange)) +
      coord_fixed() + 
      facet_wrap(~Proxy, ncol = 1) + 
      theme_bw() + 
      # xlab("Probability of Rewiring Required to Realise Novel Links") + 
      # ylab("Proportion of Initial Interaction Strength Required for Continued Existence") + 
      xlab("") + 
      ylab("") + 
      guides(fill = guide_colourbar(barwidth = 15,
                                    barheight = 1.5,
                                    title = "Proportional Network Metric Change",
                                    title.position = "bottom",
                                    # direction = "horizontal",
                                    legend.location = "bottom")) + 
      scale_fill_viridis("Mean", option = "B", direction = -1) + 
      theme(legend.position = "bottom") + ggtitle(TopoPlots[TopoIter]) + 
      theme(plot.margin = unit(c(0,0,0,0), "lines"))
    
    matplot_ls$SD[[TopoIter]] <- ggplot(sd_df2, aes(x = RE, y = IS)) +
      geom_tile(aes(fill = RelChange)) +
      coord_fixed() +
      facet_wrap(~Proxy, ncol = 1) + 
      theme_bw() + 
      # xlab("Probability of Rewiring Required to Realise Novel Links") + 
      # ylab("Proportion of Initial Interaction Strength Required for Continued Existence") + 
      xlab("") + 
      ylab("") + 
      guides(fill = guide_colourbar(barwidth = 15,
                                    barheight = 1.5,
                                    title = "Proportional Network Metric Change",
                                    title.position = "bottom",
                                    # direction = "horizontal",
                                    legend.location = "bottom")) + 
      scale_fill_viridis("SD", option = "D") + 
      theme(legend.position = "bottom") + ggtitle(TopoPlots[TopoIter]) + 
      theme(plot.margin = unit(c(0,0,0,0), "lines"))
  }
  y.grob <- textGrob("Proportion of Initial Interaction Strength Required for Continued Existence", 
                     gp=gpar(fontface="bold", col="black", fontsize=15), rot=90)
  x.grob <- textGrob("Probability of Rewiring Required to Realise Novel Links", 
                     gp=gpar(fontface="bold", col="black", fontsize=15))
  
  MatPred_plot <- plot_grid(plotlist = matplot_ls$Pred, nrow = 1)
  MatPred_plot <- grid.arrange(arrangeGrob(MatPred_plot, left = y.grob, top = x.grob))
  ggsave(MatPred_plot, filename = file.path(Dir.Exports, paste0("PLOT_MatrixChange", RunName, ".png")), width = 40/1.2, height = 34/1.2, units = "cm")
  
  MatSD_plot <- plot_grid(plotlist = matplot_ls$SD, nrow = 1)
  MatSD_plot <- grid.arrange(arrangeGrob(MatSD_plot, left = y.grob, top = x.grob))
  ggsave(MatSD_plot, filename = file.path(Dir.Exports, paste0("PLOT_MatrixChange_SD", RunName, ".png")), width = 40/1.2, height = 34/1.2, units = "cm") 
}

### Effects of bottom-up and top-down compared to ALL and each other ----
CompCasc_ls <- list(Plants = NA,
                    Animals = NA)
for(RunName in c("Plants", "Animals")){
  ### Baseline ALL Changes ----
  Change_df <- PlotTopo_ls[["ALL"]]$Change
  TopoPlots <- c("n_species", "n_links", "Modularity", "Nestedness") # "n_species", "n_animals", "n_plants", "Modularity", "n_links", "Nestedness"
  Change_df <- Change_df[Change_df$Topology %in% TopoPlots & Change_df$Proxy %in% c("Climate", "IUCN", "Strength"), ]
  base_plot_df <- aggregate(RelChange ~ Proxy+Topology+IS+RE, FUN = mean, data = Change_df[Change_df$RelChange != 0, ])
  base_sd_df <- aggregate(RelChange ~ Proxy+Topology+IS+RE, FUN = sd, data = Change_df[Change_df$RelChange != 0, ])
  
  ### Relative Changes ----
  Change_df <- PlotTopo_ls[[RunName]]$Change
  TopoPlots <- c("n_species", "n_links", "Modularity", "Nestedness") # "n_species", "n_animals", "n_plants", "Modularity", "n_links", "Nestedness"
  Change_df <- Change_df[Change_df$Topology %in% TopoPlots & Change_df$Proxy %in% c("Climate", "IUCN", "Strength"), ]
  test_plot_df <- aggregate(RelChange ~ Proxy+Topology+IS+RE, FUN = mean, data = Change_df[Change_df$RelChange != 0, ])
  test_sd_df <- aggregate(RelChange ~ Proxy+Topology+IS+RE, FUN = sd, data = Change_df[Change_df$RelChange != 0, ])
  
  CompCasc_ls[[RunName]] <- list(mean = test_plot_df,
                                 sd = test_sd_df)
  
  ### Changes compared to baseline ----
  base_combs <- with(base_plot_df, paste(Proxy, Topology, IS, RE, sep="-"))
  test_combs <- with(test_plot_df, paste(Proxy, Topology, IS, RE, sep="-"))
  
  plot_df <- base_plot_df[base_combs %in% test_combs, ]
  plot_df$RelChange <- base_plot_df[base_combs %in% test_combs, 5] - test_plot_df[test_combs %in% base_combs, 5]
  
  sd_df <- base_sd_df[base_combs %in% test_combs, ]
  sd_df$RelChange <- base_sd_df[base_combs %in% test_combs, 5] - test_sd_df[test_combs %in% base_combs, 5]
  
  # iterate my Topology, then fuse plots
  matplot_ls <- list(Pred = as.list(rep(NA, length(TopoPlots))),
                     SD = as.list(rep(NA, length(TopoPlots))))
  
  
  for(TopoIter in 1:length(TopoPlots)){
    plot_df2 <- plot_df[plot_df$Topology == TopoPlots[TopoIter], ]
    sd_df2 <- sd_df[sd_df$Topology == TopoPlots[TopoIter], ]
    
    matplot_ls$Pred[[TopoIter]] <- ggplot(plot_df2, aes(x = RE, y = IS)) +
      geom_tile(aes(fill = RelChange)) +
      coord_fixed() + 
      facet_wrap(~Proxy, ncol = 1) + 
      theme_bw() + 
      # xlab("Probability of Rewiring Required to Realise Novel Links") + 
      # ylab("Proportion of Initial Interaction Strength Required for Continued Existence") + 
      xlab("") + 
      ylab("") + 
      guides(fill = guide_colourbar(barwidth = 15,
                                    barheight = 1.5,
                                    title = "Proportional Network Metric Change",
                                    title.position = "bottom",
                                    # direction = "horizontal",
                                    legend.location = "bottom")) + 
      scale_fill_viridis("Mean", option = "B", direction = -1) + 
      theme(legend.position = "bottom") + ggtitle(TopoPlots[TopoIter]) + 
      theme(plot.margin = unit(c(0,0,0,0), "lines"))
    
    matplot_ls$SD[[TopoIter]] <- ggplot(sd_df2, aes(x = RE, y = IS)) +
      geom_tile(aes(fill = RelChange)) +
      coord_fixed() +
      facet_wrap(~Proxy, ncol = 1) + 
      theme_bw() + 
      # xlab("Probability of Rewiring Required to Realise Novel Links") + 
      # ylab("Proportion of Initial Interaction Strength Required for Continued Existence") + 
      xlab("") + 
      ylab("") + 
      guides(fill = guide_colourbar(barwidth = 15,
                                    barheight = 1.5,
                                    title = "Proportional Network Metric Change",
                                    title.position = "bottom",
                                    # direction = "horizontal",
                                    legend.location = "bottom")) + 
      scale_fill_viridis("SD", option = "D") + 
      theme(legend.position = "bottom") + ggtitle(TopoPlots[TopoIter]) + 
      theme(plot.margin = unit(c(0,0,0,0), "lines"))
  }
  y.grob <- textGrob("Proportion of Initial Interaction Strength Required for Continued Existence", 
                     gp=gpar(fontface="bold", col="black", fontsize=15), rot=90)
  x.grob <- textGrob("Probability of Rewiring Required to Realise Novel Links", 
                     gp=gpar(fontface="bold", col="black", fontsize=15))
  
  MatPred_plot <- plot_grid(plotlist = matplot_ls$Pred, nrow = 1)
  MatPred_plot <- grid.arrange(arrangeGrob(MatPred_plot, left = y.grob, top = x.grob))
  ggsave(MatPred_plot, filename = file.path(Dir.Exports, paste0("PLOT_MatrixChange", RunName, "-ComparedtoALL.png")), width = 40/1.2, height = 34/1.2, units = "cm")
  
  MatSD_plot <- plot_grid(plotlist = matplot_ls$SD, nrow = 1)
  MatSD_plot <- grid.arrange(arrangeGrob(MatSD_plot, left = y.grob, top = x.grob))
  ggsave(MatSD_plot, filename = file.path(Dir.Exports, paste0("PLOT_MatrixChange_SD", RunName, "-ComparedtoALL.png")), width = 40/1.2, height = 34/1.2, units = "cm") 
}

plants_combs <- with(CompCasc_ls$Plants$mean, paste(Proxy, Topology, IS, RE, sep="-"))
animals_combs <- with(CompCasc_ls$Animals$mean, paste(Proxy, Topology, IS, RE, sep="-"))

plot_df <- CompCasc_ls$Animals$mean[base_combs %in% test_combs, ]
plot_df$RelChange <- CompCasc_ls$Plants$mean[plants_combs %in% animals_combs, 5] - 
  CompCasc_ls$Animals$mean[animals_combs %in% plants_combs, 5]
sd_df <- CompCasc_ls$Animals$sd[base_combs %in% test_combs, ]
sd_df$RelChange <- CompCasc_ls$Plants$sd[plants_combs %in% animals_combs, 5] - 
  CompCasc_ls$Animals$sd[animals_combs %in% plants_combs, 5]

# iterate my Topology, then fuse plots
matplot_ls <- list(Pred = as.list(rep(NA, length(TopoPlots))),
                   SD = as.list(rep(NA, length(TopoPlots))))


for(TopoIter in 1:length(TopoPlots)){
  plot_df2 <- plot_df[plot_df$Topology == TopoPlots[TopoIter], ]
  sd_df2 <- sd_df[sd_df$Topology == TopoPlots[TopoIter], ]
  
  matplot_ls$Pred[[TopoIter]] <- ggplot(plot_df2, aes(x = RE, y = IS)) +
    geom_tile(aes(fill = RelChange)) +
    coord_fixed() + 
    facet_wrap(~Proxy, ncol = 1) + 
    theme_bw() + 
    # xlab("Probability of Rewiring Required to Realise Novel Links") + 
    # ylab("Proportion of Initial Interaction Strength Required for Continued Existence") + 
    xlab("") + 
    ylab("") + 
    guides(fill = guide_colourbar(barwidth = 15,
                                  barheight = 1.5,
                                  title = "Proportional Network Metric Change",
                                  title.position = "bottom",
                                  # direction = "horizontal",
                                  legend.location = "bottom")) + 
    scale_fill_viridis("Mean", option = "B", direction = -1) + 
    theme(legend.position = "bottom") + ggtitle(TopoPlots[TopoIter]) + 
    theme(plot.margin = unit(c(0,0,0,0), "lines"))
  
  matplot_ls$SD[[TopoIter]] <- ggplot(sd_df2, aes(x = RE, y = IS)) +
    geom_tile(aes(fill = RelChange)) +
    coord_fixed() +
    facet_wrap(~Proxy, ncol = 1) + 
    theme_bw() + 
    # xlab("Probability of Rewiring Required to Realise Novel Links") + 
    # ylab("Proportion of Initial Interaction Strength Required for Continued Existence") + 
    xlab("") + 
    ylab("") + 
    guides(fill = guide_colourbar(barwidth = 15,
                                  barheight = 1.5,
                                  title = "Proportional Network Metric Change",
                                  title.position = "bottom",
                                  # direction = "horizontal",
                                  legend.location = "bottom")) + 
    scale_fill_viridis("SD", option = "D") + 
    theme(legend.position = "bottom") + ggtitle(TopoPlots[TopoIter]) + 
    theme(plot.margin = unit(c(0,0,0,0), "lines"))
}
y.grob <- textGrob("Proportion of Initial Interaction Strength Required for Continued Existence", 
                   gp=gpar(fontface="bold", col="black", fontsize=15), rot=90)
x.grob <- textGrob("Probability of Rewiring Required to Realise Novel Links", 
                   gp=gpar(fontface="bold", col="black", fontsize=15))

MatPred_plot <- plot_grid(plotlist = matplot_ls$Pred, nrow = 1)
MatPred_plot <- grid.arrange(arrangeGrob(MatPred_plot, left = y.grob, top = x.grob))
ggsave(MatPred_plot, filename = file.path(Dir.Exports, "PLOT-RelChange_Plants-Animals.png"), width = 40/1.2, height = 34/1.2, units = "cm")

MatSD_plot <- plot_grid(plotlist = matplot_ls$SD, nrow = 1)
MatSD_plot <- grid.arrange(arrangeGrob(MatSD_plot, left = y.grob, top = x.grob))
ggsave(MatSD_plot, filename = file.path(Dir.Exports, "PLOT-RelChange_Plants-Animals_SD.png"), width = 40/1.2, height = 34/1.2, units = "cm") 




PlotTopo_ls2 <- lapply(names(PlotTopo_ls), function(x){
  returnme <- PlotTopo_ls[[x]]$Change
  returnme$Cascade <- x
  returnme
})

LM_df <- do.call(rbind, PlotTopo_ls2)
LM_df <- LM_df[LM_df$Proxy != "IUCN_Climate", ]

model_iter <- "n_species"

summary(lm(RelChange ~ 0 + (IS+RE)*Proxy, data = LM_df[LM_df$Topology == model_iter, ]))






















### Effect Sizes ----
EffectSize_df <- PlotTopoAll_ls$EffectSize
TopoPlots <- c("n_species", "n_links", "Modularity", "Nestedness")
EffectSize_df <- stats::reshape(data = EffectSize_df, 
                                times = colnames(EffectSize_df)[1:(ncol(EffectSize_df)-4)],
                                varying = list(colnames(EffectSize_df)[1:(ncol(EffectSize_df)-4)]),
                                timevar = "Topology",
                                direction = "long")
EffectSize_df <- EffectSize_df[EffectSize_df$Topology %in% TopoPlots, ]
colnames(EffectSize_df)[ncol(EffectSize_df)-1] <- "EffectSize"
EffectSize_df <- EffectSize_df[which(abs(EffectSize_df$EffectSize) != Inf), ]
EffectSize_df <- EffectSize_df[which(!is.na(EffectSize_df$EffectSize)), ]

plot_df <- aggregate(EffectSize ~ Pry+Topology+IS+RE, FUN = mean, data = EffectSize_df)
sd_df <- aggregate(EffectSize ~ Pry+Topology+IS+RE, FUN = sd, data = EffectSize_df)

Upper <- plot_df$EffectSize + sd_df$EffectSize*2
Lower <- plot_df$EffectSize - sd_df$EffectSize*2
plot_df$Sig <- ifelse(abs(sign(Upper) + sign(Lower)) == 2, TRUE, FALSE)

effectsizeplot <- ggplot(plot_df, aes(x = RE, y = IS)) +
  geom_tile(aes(fill = EffectSize)) +
  coord_fixed() + 
  facet_wrap(~Topology+Pry) + 
  scale_fill_gradient2(high = "darkgreen", low = "darkred") +
  geom_point(aes(shape = Sig), size = 2) +
  scale_shape_manual(values=c(32, 15), na.translate = FALSE, name = "Significance") +
  guides(shape = FALSE) + 
  theme_bw() + 
  xlab("Probability of Rewiring Required to Realise Novel Links") + 
  ylab("Proportion of Initial Interaction Strength Required for Continued Existence") + 
  guides(fill = guide_colourbar(barwidth = 2,
                                barheight = 15,
                                title = "Effect Size"))
effectsizeplot
ggsave(effectsizeplot, filename = file.path(Dir.Exports, "PLOT_MatrixChangeEffectSize_ALL.png"), width = 40/1.2, height = 34/1.2, units = "cm") 

stop("significance of effect size?")

# sd_gg <- ggplot(sd_df[sd_df$Topology %in% TopoPlots,], aes(x = RE, y = IS)) +
#   geom_tile(aes(fill = EffectSize)) +
#   coord_fixed() +
#   facet_wrap(~Topology+Pry) + 
#   scale_fill_viridis("SD", option = "D")
# 
# print(plot_grid(pred_gg, sd_gg, nrow = 2))

## Effects of each extinction proxy ----------------------------------------

stop("increase this RE selection to 1 when trying to recreate initial plots")
Plot_df <- PlotTopoAll_ls$Change[PlotTopoAll_ls$Change$RE == 1, ]


TopoPlots <- c("n_species", "n_links", "Modularity", "Nestedness")
ggplot(Plot_df[Plot_df$Topology %in% TopoPlots & Plot_df$Proxy != "IUCN_Climate",],
       aes(x = IS, y = RelChange, col = Proxy)) + 
  geom_point(alpha = 0.4) + 
  geom_smooth() + 
  facet_wrap(~Topology, scales = "free", ncol = 1) + 
  scale_color_manual(values = pal_lm) +
  theme_bw() + 
  labs(y = "Relative Change in Network Topology (Pre- vs Post-Extinction)",
       x = "Network Dependency")

ggplot(Plot_df[Plot_df$Topology %in% TopoPlots & Plot_df$Proxy != "IUCN_Climate",],
       aes(x = factor(IS), y = RelChange, col = Proxy)) + 
  geom_boxplot() + 
  facet_wrap(~Topology + Proxy, scales = "free", ncol = 3) + 
  scale_color_manual(values = pal_lm) + 
  theme_bw() + 
  labs(y = "Relative Change in Network Topology (Pre- vs Post-Extinction)",
       x = "Network Dependency")

Compare_df <- Plot_df[Plot_df$Topology %in% TopoPlots,]
Clim_df <- Compare_df[Compare_df$Proxy == "Climate", ]
IUCN_df <- Compare_df[Compare_df$Proxy == "IUCN", ]

Plot_df <- Clim_df[,1:5]
Plot_df <- cbind(Plot_df, Clim_df$RelChange - IUCN_df$RelChange)
colnames(Plot_df)[ncol(Plot_df)] <- "Difference"

ggplot(Plot_df,
       aes(x = factor(IS), y = Difference)) +
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  facet_wrap(~Topology, scales = "free", ncol = 1) +
  scale_color_manual(values = pal_lm) +
  theme_bw() +
  labs(title = "Difference between Climate- and IUCN-driven Extinctions",
       y = "Climate - IUCN Difference",
       x = "Network Dependency")


## Effects of each extinction cascade --------------------------------------

# ## Individual Cascade-Directions -------------------------------------------
# tTest.Calc <- function(Diff_df){
#   tTests <- expand.grid(unique(Diff_df$Proxy), unique(Diff_df$Topology), unique(Diff_df$IS))
#   Signif <- lapply(1:nrow(tTests), function(x){
#     vals <-Diff_df$Value[Diff_df$Proxy == tTests[x,1] & 
#                            Diff_df$Topology == tTests[x,2] & 
#                            Diff_df$IS == tTests[x,3]]
#     if(sum(!is.na(vals)) > 1){
#       ifelse(t.test( na.omit(vals)[!is.infinite(na.omit(vals))] )[["p.value"]] < 0.05, TRUE, FALSE)
#     }else{
#       FALSE
#     }
#   })
#   tTests$Signif <- unlist(Signif)
#   colnames(tTests) <- c("Proxy", "Topology", "IS", "Signif")
#   
#   Diff_df <- merge(x = Diff_df, y = tTests, by = colnames(tTests)[1:3])
#   Diff_df
# }
# 
# Diff.Calc <- function(Topo_df){
#   Diff_df <- Topo_df[Topo_df$Simulation == "Prediction", 1:6] - PreExt_df[match(Topo_df[Topo_df$Simulation == "Prediction",]$netID, PreExt_df$netID), 1:6]
#   Diff_df$IS <- Topo_df[Topo_df$Simulation == "Prediction",]$IS
#   Diff_df$Proxy <- Topo_df[Topo_df$Simulation == "Prediction",]$Proxy
#   Diff_df$netID <- unlist(lapply(strsplit(rownames(Diff_df), split = "[.]"), "[[", 2))
#   Diff_df <- reshape(data = Diff_df,
#                      idvar = "netID",
#                      varying = colnames(Diff_df)[1:6],
#                      v.name = "Value",
#                      timevar = "Topology",
#                      times = colnames(Diff_df)[1:6],
#                      new.row.names = 1:(nrow(Diff_df)*6),
#                      direction = "long"
#   )
#   
#   tTest.Calc(Diff_df)
# }
# 
# 
# Plot.Run <- function(topolist, RunName){
#   Topo_df <- topolist[[1]]
#   Eff_df <- topolist[[2]]
#   pal_signif <- c("#c32200", "#62c300")
#   
#   ## Effect size against random simulation
#   Plot_df<- reshape(data = Eff_df,
#                     idvar = "netID",
#                     varying = colnames(Eff_df)[1:6],
#                     v.name = "Value",
#                     timevar = c("Topology"),
#                     times = colnames(Eff_df)[1:6],
#                     new.row.names = 1:(nrow(Eff_df)*6),
#                     direction = "long"
#   )
#   colnames(Plot_df)[2] <- "Proxy"
#   Plot_df <- tTest.Calc(Plot_df)
#   pal_lm <- c("#016392", "#E19825", "#3E8853")
#   Eff_gg <- ggplot(Plot_df, aes(y = Value, x = IS, fill = Proxy, col = Proxy)) + 
#     geom_point(alpha = 0.2) + 
#     stat_smooth(method = "lm") + 
#     scale_fill_manual(values = pal_lm) + 
#     scale_color_manual(values = pal_lm) + 
#     facet_wrap(~ factor(Topology, levels = c("n_species", "n_plants", "n_animals",
#                                              "n_links", "Nestedness", "Modularity")), 
#                scales = "free", ncol = 3) + 
#     theme_bw() + 
#     labs(x = "Interaction Strength % Required For Survival", y = "Effect Size")
#   
#   Eff_ls <- lapply(unique(Plot_df$Proxy), FUN = function(x){
#     ggplot(Plot_df[Plot_df$Proxy == x, ], aes(y = Value, x = factor(IS), fill = Signif)) + 
#       geom_boxplot() +
#       scale_fill_manual(values = pal_signif) + 
#       facet_wrap(~factor(Topology, levels = c("n_species", "n_plants", "n_animals",
#                                               "n_links", "Nestedness", "Modularity")),
#                  scales = "free", ncol = 3) +
#       theme_bw() +
#       labs(x = "Interaction Strength % Required For Survival", y = "Effect Size", title = x)
#   })
#   
#   ## Absolute effect
#   Diff_df <- Diff.Calc(Topo_df)
#   Diff_gg <- ggplot(Diff_df, aes(y = Value, x = IS, col = Proxy, fill = Proxy)) + 
#     geom_point(alpha = 0.2) + 
#     geom_smooth() + 
#     scale_fill_manual(values = pal_lm) + 
#     scale_color_manual(values = pal_lm) + 
#     facet_wrap(~ factor(Topology, levels = c("n_species", "n_plants", "n_animals",
#                                              "n_links", "Nestedness", "Modularity")), 
#                scales = "free", ncol = 3) + 
#     theme_bw() + 
#     labs(x = "Interaction Strength % Required For Survival", y = "Change In Network Topology Metric")
#   
#   Diff_ls <- lapply(unique(Diff_df$Proxy), FUN = function(x){
#     ggplot(Diff_df[Diff_df$Proxy == x, ], aes(y = Value, x = factor(IS), col = Signif)) +
#       geom_boxplot() +
#       scale_color_manual(values = pal_signif) + 
#       facet_wrap(~ factor(Topology, levels = c("n_species", "n_plants", "n_animals",
#                                                "n_links", "Nestedness", "Modularity")),
#                  scales = "free", ncol = 3) +
#       theme_bw() +
#       labs(x = "Interaction Strength % Required For Survival", y = "Change In Network Topology Metric", title = x)
#   })
#   
#   ## Relative effect
#   Topo_df <- topolist[[1]]
#   Topo_df <- Topo_df[Topo_df$Simulation == "Prediction", ]
#   Plot_df <- merge(PreExt_df, Topo_df, by = c("netID"))
#   Rel_ls <- lapply(c("n_species", "n_plants", "n_animals", "n_links", "Nestedness", "Modularity"), function(x){
#     Change_vec <- apply(Plot_df[ , grepl(x,colnames(Plot_df))], 1, diff)
#     Plot_df <- data.frame(
#       Pre = Plot_df[,grep(x,colnames(Plot_df))[1]],
#       Value = Change_vec,
#       netID = Plot_df$netID,
#       Proxy = Plot_df$Proxy.y,
#       IS = Plot_df$IS
#     )
#     Plot_df$RelChange <- abs(Plot_df$Value)/Plot_df$Pre
#     
#     Plot_df2 <- aggregate(RelChange ~ IS + Pre + Proxy, FUN = mean, data = Plot_df)
#     ggplot(Plot_df2, aes(x = Pre , y = IS, fill = RelChange)) + 
#       geom_tile(aes(fill = RelChange)) +
#       # coord_fixed() + 
#       facet_wrap(~factor(Proxy, levels = c("Strength", "Climate", "IUCN")),
#                  scales = "free", nrow = 3) +
#       scale_fill_viridis() + 
#       theme_bw() + labs(x = paste(x, "(pre-extinction"), y = "Interaction Strength % Required For Survival")
#     
#     ggplot(Plot_df, aes(x = factor(IS), y = RelChange)) + 
#       geom_boxplot() + 
#       facet_wrap(~factor(Proxy, levels = c("Strength", "Climate", "IUCN")),
#                  scales = "free", nrow = 3) +
#       theme_bw() + labs(y = paste("Relative change in", x, "compared to pre-extinction level"), 
#                         x = "Interaction Strength % Required For Survival",
#                         title = x)
#   })
#   
#   ## Saving plots
#   pdf(file.path(Dir.Exports, paste0(RunName, "_Print.pdf")), paper = "a4r", height = 9, width = 12)
#   print(Eff_gg)
#   print(Eff_ls)
#   print(Diff_gg)
#   print(Diff_ls)
#   print(Rel_ls)
#   dev.off()
# }
# 
# Plot.Run(PlotTopoAll_ls, "ALL")
# Plot.Run(PlotTopoAnimals_ls, "Animals")
# Plot.Run(PlotTopoPlants_ls, "Plants")
# 
# ## Comparison of Cascade-Directions ----------------------------------------
# DiffAll <- Diff.Calc(PlotTopoAll_ls[[1]])
# DiffAnim <- Diff.Calc(PlotTopoAnimals_ls[[1]])
# DiffPlan <- Diff.Calc(PlotTopoPlants_ls[[1]])
# 
# df1 <- cbind(DiffAll[,-5:-6], 
#       DiffAnim[,5] - DiffAll[,5])
# colnames(df1)[5] <- "Value"
# df1 <- tTest.Calc(df1)
# df1$Comparison <- "Top-Down"
# df2 <- cbind(DiffAll[,-5:-6], 
#       DiffPlan[,5] - DiffAll[,5])
# colnames(df2)[5] <- "Value"
# df2 <- tTest.Calc(df2)
# df2$Comparison <- "Bottom-Up"
# 
# Plot_df <- na.omit(rbind(df1, df2))
# 
# # comps <- list( c("Top-Down", "Bottom-Up"))
# pal_signif <- c("#c32200", "#62c300")
# pal_cascade <- c("#880094", "#947100")
# Comps_gg <- lapply(c("n_species", "n_plants", "n_animals", "n_links", "Nestedness", "Modularity"), function(x){
#   ggplot(Plot_df[Plot_df$Topology == x, ], 
#          aes(x = factor(IS), y = Value, col = Comparison, fill = Signif)) + 
#     geom_hline(aes(yintercept = 0)) +
#     geom_boxplot() + 
#     scale_fill_manual(values = pal_signif) +
#     scale_color_manual(values = pal_cascade) +
#     # stat_compare_means(comparisons = comps, method = 't.test', label = 'p.signif') +
#     facet_wrap(~factor(Proxy, levels = c("Strength", "Climate", "IUCN")),
#                scales = "free", nrow = 3) +
#     theme_bw() +
#     labs(x = "Interaction Strength % Required For Survival", y = paste0("Change In ", x), title = x)
# })
# 
# pdf(file.path(Dir.Exports, "Comps.pdf"), paper = "a4", height = 12, width = 9)
# print(Comps_gg)
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
