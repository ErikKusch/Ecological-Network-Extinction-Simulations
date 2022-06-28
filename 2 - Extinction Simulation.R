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

# PARALLEL EXECTIONS =======================================================
plants_sp <- rownames(plants_gowdis)
animals_sp <- rownames(animals_gowdis)

CutOffs <- list(Strength = 0.75,
                Climate = 2,
                IUCN = 5)

message("### REGISTERING CLUSTER")
nCores <- ifelse(parallel::detectCores()>length(AnalysisData_ls), 
                 length(AnalysisData_ls), parallel::detectCores())
cl <- parallel::makeCluster(nCores) # for parallel pbapply functions
parallel::clusterExport(cl,
                        varlist = c('FUN_Topo', "animals_gowdis", "plants_gowdis", "AnalysisData_ls", "install.load.package", "package_vec", ".DataInit", "plants_sp", "animals_sp", "CutOffs", "ExtinctionOrder", "RandomExtinctions", "SimulateExtinctions"),
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

for(Rewiring_Iter in seq(0, 1, 0.5)){
  for(IS_iter in seq(0, 1, 0.05)){
    Sim_ls <- FUN_SimComp(PlantAnim = NULL, RunName = "ALL", 
                          IS = IS_iter, Rewiring = Rewiring_Iter,
                          CutOffs = CutOffs)
    # TopoComp_ls <- FUN_TopoComp(Sim_ls = Sim_ls, RunName = "ALL", 
    #                             IS = IS_iter, Rewiring = Rewiring_Iter,
    #                             CutOffs = CutOffs)
    # Sim_ls <- FUN_SimComp(PlantAnim = plants_sp, RunName = "Plants",
    #                       IS = IS_iter, Rewiring = Rewiring_Iter,
    #                       CutOffs = CutOffs)
    # TopoComp_ls <- FUN_TopoComp(Sim_ls = Sim_ls, RunName = "Plants",
    #                             IS = IS_iter, Rewiring = Rewiring_Iter,
    #                             CutOffs = CutOffs)
    # Sim_ls <- FUN_SimComp(PlantAnim = animals_sp, RunName = "Animals",
    #                       IS = IS_iter, Rewiring = Rewiring_Iter,
    #                       CutOffs = CutOffs)
    # TopoComp_ls <- FUN_TopoComp(Sim_ls = Sim_ls, RunName = "Animals",
    #                             IS = IS_iter, Rewiring = Rewiring_Iter,
    #                             CutOffs = CutOffs)
  }
}



message("Sensitivity analysis for WHICH = 'Strength' in FUN_SimComp.")

# VISUALISATION ============================================================
message("### RESULT VISUALISATION ###")
pal_lm <- c("#016392", "#E19825", "#3E8853")

## while loading in the topologies, we also compute absolute and relative change of each simulation to the pre-extinction network topologies
PlotTopoAll_ls <- loadTopo(RunName = "ALL", CutOffs = CutOffs, Pre = PreExt_df)
PlotTopoPlants_ls <- loadTopo(RunName = "Plants", CutOffs = CutOffs, Pre = PreExt_df)
PlotTopoAnimals_ls <- loadTopo(RunName = "Animals", CutOffs = CutOffs, Pre = PreExt_df)

## Venn-Diagram of Proxy Agreement -----------------------------------------
print("Extinction Proxy Overlap ---")
for(RunName in c("ALL", "Plants", "Animals")){
  ExtSpecies_ls <- FUN_SimComp(PlantAnim = NULL, RunName = RunName, IS = 0, CutOffs = CutOffs)
  
  ## number of primary extinctions per proxy
  Clim_ls <- lapply(lapply(ExtSpecies_ls, "[[", "Climate"), "[[", "Removed")
  IUCN_ls <- lapply(lapply(ExtSpecies_ls, "[[", "IUCN"), "[[", "Removed")
  Centr_ls <- lapply(lapply(ExtSpecies_ls, "[[", "Strength"), "[[", "Removed")
  ClimIUCN_ls <- lapply(1:length(Clim_ls), function(x){sum(Clim_ls[[x]] %in% IUCN_ls[[x]])})
  ClimCentr_ls <- lapply(1:length(Clim_ls), function(x){sum(Clim_ls[[x]] %in% Centr_ls[[x]])})
  IUCNCentr_ls <- lapply(1:length(Clim_ls), function(x){sum(IUCN_ls[[x]] %in% Centr_ls[[x]])})
  All_ls <- lapply(1:length(Clim_ls), function(x){
    sum(
      IUCN_ls[[x]][IUCN_ls[[x]] %in% Centr_ls[[x]]] %in% Clim_ls[[x]]
      )
    })
  ## absolute numbers of primary extinctions per network
  PrimaryExt_df <- data.frame(
    Climate = unlist(lapply(Clim_ls, length)),
    IUCN = unlist(lapply(IUCN_ls, length)),
    Centrality = unlist(lapply(Centr_ls, length)),
    `Climate+IUCN` = unlist(ClimIUCN_ls),
    `Climate+Centrality` = unlist(ClimCentr_ls),
    `IUCN+Centrality` = unlist(IUCNCentr_ls),
    All = unlist(All_ls)
  )
  ## relative numbers of primary extinctions per network
  PrimaryExt_df <- apply(X = PrimaryExt_df, MARGIN = 2, function(x){
    x/PreExt_df$n_species
  })
  
  ## summary statistics
  PrimExt_df <- data.frame(min = apply(PrimaryExt_df, 2, min),
                           mean = apply(PrimaryExt_df, 2, mean),
                           max = apply(PrimaryExt_df, 2, max),
                           sd = apply(PrimaryExt_df, 2, sd))
  print(PrimExt_df)
  
  ## Total Venn-Diagram
  Venn_ls <- list(Climate = unlist(lapply(lapply(ExtSpecies_ls, "[[", "Climate"), "[[", "Removed")),
                  IUCN = unlist(lapply(lapply(ExtSpecies_ls, "[[", "IUCN"), "[[", "Removed")),
                  Centrality = unlist(lapply(lapply(ExtSpecies_ls, "[[", "Strength"), "[[", "Removed")))
  ggvenn(Venn_ls, fill_color = pal_lm, fill_alpha = 0.8, text_color = "white")
  ggsave(filename = file.path(Dir.Exports, paste0("PLOT_Proxy", RunName,".png")), width = 4, height = 3, units = "cm", scale = 7, dpi = 1e3)
}

## Effects of each extinction proxy ----------------------------------------
TopoPlots <- c("n_species", "n_animals", "n_plants")
ggplot(PlotTopoAll_ls$Change[PlotTopoAll_ls$Change$Topology %in% TopoPlots,],
       aes(x = IS, y = RelChange, col = Proxy)) + 
  geom_point(alpha = 0.4) + 
  geom_smooth() + 
  facet_wrap(~Topology, scales = "free", ncol = 1) + 
  scale_color_manual(values = pal_lm) + 
  theme_bw() + 
  labs(y = "Relative Change in Network Topology (Pre- vs Post-Extinction)",
       x = "Network Dependency")

ggplot(PlotTopoAll_ls$Change[PlotTopoAll_ls$Change$Topology %in% TopoPlots,],
       aes(x = factor(IS), y = RelChange, col = Proxy)) + 
  geom_boxplot() + 
  facet_wrap(~Topology + Proxy, scales = "free", ncol = 3) + 
  scale_color_manual(values = pal_lm) + 
  theme_bw() + 
  labs(y = "Relative Change in Network Topology (Pre- vs Post-Extinction)",
       x = "Network Dependency")

Compare_df <- PlotTopoAll_ls$Change[PlotTopoAll_ls$Change$Topology %in% TopoPlots,]
Clim_df <- Compare_df[Compare_df$Proxy == "Climate", ]
IUCN_df <- Compare_df[Compare_df$Proxy == "IUCN", ]

Plot_df <- Clim_df[,1:4]
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
