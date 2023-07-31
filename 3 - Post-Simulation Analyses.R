#' ####################################################################### #
#' PROJECT: [Biodiversity Simplification & Ecological Network Topology] 
#' CONTENTS: 
#'  - Post-Simulation Analyses
#'  DEPENDENCIES:
#'  - "0 - Preamble.R"
#'  - "2 - Extinction Simulation" has to have been run
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE =================================================================
rm(list=ls())
set.seed(42)

## Sourcing ----------------------------------------------------------------
source("2 - Extinction Simulation.R")

# DATA LOADING & MANIPULATING ==============================================
message("### DATA PREPARATION ###")

## Target Metrics ----------------------------------------------------------
TopoPlots <- c("n_species", "n_links") # "n_species", "n_links, "n_animals", "n_plants", "Modularity", "Nestedness"

## Data Loading ------------------------------------------------------------
load(file.path(Dir.Data, "Networks.RData")) 
Fricke_fs <- list.files(Dir.D.Fricke, full.names = TRUE)
hush_ls <- pblapply(Fricke_fs[-grep(x = Fricke_fs, pattern = ".zip")],
                    load,
                    .GlobalEnv)
world <- ne_countries(scale = "medium", returnclass = "sf")
map_df <- networks_df[networks_df$net.id %in% names(AnalysisData_ls), ]

DT <- sf::st_as_sf(map_df[, c("net.id","longitude","latitude")], 
                   coords = c("longitude","latitude"), 
                   crs = crs(world))
sf_use_s2(FALSE)
world$samplenets <- lengths(st_intersects(world, DT))

PlotTopo_ls <- list(ALL = PlotTopoAll_ls,
                    Animals = PlotTopoAnimals_ls,
                    Plants = PlotTopoPlants_ls
)
PlotTopo_ls$ALL$Change$Casc <- "ALL"
PlotTopo_ls$Animals$Change$Casc <- "Animals"
PlotTopo_ls$Plants$Change$Casc <- "Plants"

## Data Limiting -----------------------------------------------------------
AnalysisData_ls <- AnalysisData_ls[28] # loading data for single network
RM_spec <- FUN_SimComp(PlantAnim = NULL, RunName = "FIG2", CutOffs = CutOffs, 
                       PotPartners = RewClass_ls, Traits = meta_df, 
                       IS = 1, Rewiring = 1, WHICH = "Climate")$`Arias-Arone 2016`$Climate$Removed # removed species in target network

# FIGURE 2 - One-Network Resilience Landscape ------------------------------
print("########## FIGURE 2")
### Matrix Visualisation ----
Change_df <- PlotTopoAll_ls$Change
Change_df <- Change_df[Change_df$Topology %in% "n_species" & 
                         Change_df$Proxy %in% "Climate" & 
                         Change_df$netID == "Arias-Arone 2016", ]

ggplot(Change_df, aes(x = RE, y = IS)) +
  geom_tile(aes(fill = RelChange)) +
  coord_fixed() + 
  theme_bw() + 
  xlab("Probability of Rewiring Required to Realise Novel Links") + 
  ylab("Proportion of Initial Interaction Strength Required for Continued Existence") + 
  guides(fill = guide_colourbar(barwidth = 15,
                                barheight = 1.5,
                                title = "Proportional Network Metric Change",
                                title.position = "bottom",
                                # direction = "horizontal",
                                legend.location = "bottom")) + 
  scale_fill_viridis("Mean", option = "B", direction = -1) + 
  theme(legend.position = "bottom") + ggtitle("Relative Loss of Species") + 
  theme(plot.margin = unit(c(0,0,0,0), "lines")) + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = file.path(Dir.Exports, "FIG2_Arias-Arone_2016_MAT.png"), width = 12/2, height = 12/2, units = "cm", scale = 3, dpi = 1e3)

### Pre-Extinction ----
Pre <- AnalysisData_ls["Arias-Arone 2016"][["Arias-Arone 2016"]][["Adjacency"]]
AnalysisData_ls <- AnalysisData_ls["Arias-Arone 2016"]
Pre_igraph <- igraph::graph_from_adjacency_matrix(Pre, weighted = TRUE, mode = "undirected")
g <- Pre_igraph

V(g)$type <- V(g)$name %in% rownames(plants_gowdis)
V(g)$color <- as.numeric(V(g)$name %in% RM_spec)
Pre_g <- g
V(g)$name <- paste0(substr(unlist(lapply(strsplit(V(g)$name, split = " "), "[[", 1)), start = 1, stop = 1),
                    ".",
                    substr(unlist(lapply(strsplit(V(g)$name, split = " "), "[[", 2)), start = 1, stop = 2))

png(filename = file.path(Dir.Exports, "FIG2_Arias-Arone_2016_Initial.png"), width = 16, height = 16, units = "cm", res = 1700)
plot(g, layout = layout_as_bipartite(g, hgap = 5, vgap = 0.5)[,c(2,1)], asp = 1.3, 
     vertex.color=c("grey","orange")[V(g)$color+1], 
     edge.width = E(g)$weight,
     vertex.label.color = "black",
     vertex.label.cex = 0.8,
     vertex.label.dist = 0,
     vertex.label.degree = -pi,
     vertex.size = 20
)
dev.off()

### Post-Extinction ----
plot.Concept <- function(g = Pre_g, g2 = NULL){
  layout_g <- g
  
  # colour of primary and secondary extinction vertices
  V(g)$color[V(Pre_g)$name %nin% V(g2)$name & V(Pre_g)$name %nin% RM_spec] <- 2
  
  # edge constellations
  g_df <- get.data.frame(g)
  g2_df <- get.data.frame(g2)
  
  ## lost edges
  g_el <- paste(g_df[,1], g_df[,2], sep = "_")
  g2_el <- paste(g2_df[,1], g2_df[,2], sep = "_")
  E(g)$color <- as.numeric(g_el %nin% g2_el) # 0 unaltered link, 1 deleted link
  
  ## rewired edges
  ### added to pre-existing
  match_df <- merge(g_df, g2_df, by=c("from","to")) # NA's match
  Rew_weights <- match_df$weight.y[which(apply(match_df[,3:4], 1, diff) != 0)]
  Rew_edges <- paste(match_df[,1], match_df[,2], sep = "_")[which(apply(match_df[,3:4], 1, diff) != 0)]
  E(g)$weight[na.omit(match(Rew_edges, g_el))] <- Rew_weights
  E(g)$color[na.omit(match(Rew_edges, g_el))] <- 2
  
  ### completely new
  if(length(E(g2)[g2_el %nin% g_el]) >= 1){
    test1 <- unlist(lapply(strsplit(g2_el[g2_el %nin% g_el], "_"), "[[", 1))
    test2 <- unlist(lapply(strsplit(g2_el[g2_el %nin% g_el], "_"), "[[", 2))
    g <- add_edges(g, edges = c(rbind(test1, test2)),
                   attr = list(color = 2, weight = E(g2)$weight[g2_el %nin% g_el]
                   )
    )
  }
  # node abbreviations
  V(g)$name <- paste0(substr(unlist(lapply(strsplit(V(g)$name, split = " "), "[[", 1)), start = 1, stop = 1),
                      ".",
                      substr(unlist(lapply(strsplit(V(g)$name, split = " "), "[[", 2)), start = 1, stop = 2))
  # plotting
  plot(g, layout = layout_as_bipartite(layout_g, hgap = 5, vgap = 0.5)[,c(2,1)], asp = 1.3, 
       vertex.color = c("grey","orange", "darkred")[V(g)$color+1], 
       edge.color = c("grey","darkred", "forestgreen")[E(g)$color+1], 
       edge.width = log(E(g)$weight+1)*3,
       vertex.label.color = "black",
       vertex.label.cex = 0.8,
       vertex.label.dist = 0,
       vertex.label.degree = -pi,
       vertex.size = 20
  )
}

## x axis, 0 - rewiring
Post00 <- FUN_SimComp(PlantAnim = NULL, RunName = "FIG2", CutOffs = CutOffs, 
                      PotPartners = RewClass_ls, Traits = meta_df, 
                      IS = 0, Rewiring = 0)$`Arias-Arone 2016`$Climate$Prediction
Post00_igraph <- igraph::graph_from_adjacency_matrix(Post00, weighted = TRUE, mode = "undirected")

Post00.5 <- FUN_SimComp(PlantAnim = NULL, RunName = "FIG2", CutOffs = CutOffs, 
                        PotPartners = RewClass_ls, Traits = meta_df, 
                        IS = 0, Rewiring = 0.5)$`Arias-Arone 2016`$Climate$Prediction
Post00.5_igraph <- igraph::graph_from_adjacency_matrix(Post00.5, weighted = TRUE, mode = "undirected")

Post01 <- FUN_SimComp(PlantAnim = NULL, RunName = "FIG2", CutOffs = CutOffs, 
                      PotPartners = RewClass_ls, Traits = meta_df, 
                      IS = 0, Rewiring = 1)$`Arias-Arone 2016`$Climate$Prediction
Post01_igraph <- igraph::graph_from_adjacency_matrix(Post01, weighted = TRUE, mode = "undirected")

png(filename = file.path(Dir.Exports, "FIG2_Arias-Arone_2016_0-0.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post00_igraph)
dev.off()
png(filename = file.path(Dir.Exports, "FIG2_Arias-Arone_2016_0-0.5.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post00.5_igraph)
dev.off()
png(filename = file.path(Dir.Exports, "FIG2_Arias-Arone_2016_0-1.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post01_igraph)
dev.off()

## y axis, 0 - interdependency (all the same because everything can and does rewire)
Post0.50 <- FUN_SimComp(PlantAnim = NULL, RunName = "FIG2", CutOffs = CutOffs, 
                        PotPartners = RewClass_ls, Traits = meta_df, 
                        IS = 0.5, Rewiring = 0)$`Arias-Arone 2016`$Climate$Prediction
Post0.50_igraph <- igraph::graph_from_adjacency_matrix(Post0.50, weighted = TRUE, mode = "undirected")

Post10 <- FUN_SimComp(PlantAnim = NULL, RunName = "FIG2", CutOffs = CutOffs, 
                      PotPartners = RewClass_ls, Traits = meta_df, 
                      IS = 1, Rewiring = 0)$`Arias-Arone 2016`$Climate$Prediction
Post10_igraph <- igraph::graph_from_adjacency_matrix(Post10, weighted = TRUE, mode = "undirected")

png(filename = file.path(Dir.Exports, "FIG2_Arias-Arone_2016_0.5-0.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post0.50_igraph)
dev.off()
png(filename = file.path(Dir.Exports, "FIG2_Arias-Arone_2016_1-0.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post10_igraph)
dev.off()

## combinations
Post0.50.5 <- FUN_SimComp(PlantAnim = NULL, RunName = "FIG2", CutOffs = CutOffs, 
                          PotPartners = RewClass_ls, Traits = meta_df, 
                          IS = 0.5, Rewiring = 0.5)$`Arias-Arone 2016`$Climate$Prediction
Post0.50.5_igraph <- igraph::graph_from_adjacency_matrix(Post0.50.5, weighted = TRUE, mode = "undirected")

Post10.5 <- FUN_SimComp(PlantAnim = NULL, RunName = "FIG2", CutOffs = CutOffs, 
                        PotPartners = RewClass_ls, Traits = meta_df, 
                        IS = 1, Rewiring = 0.5)$`Arias-Arone 2016`$Climate$Prediction
Post10.5_igraph <- igraph::graph_from_adjacency_matrix(Post10.5, weighted = TRUE, mode = "undirected")

Post0.51 <- FUN_SimComp(PlantAnim = NULL, RunName = "FIG2", CutOffs = CutOffs, 
                        PotPartners = RewClass_ls, Traits = meta_df,  
                        IS = 0.5, Rewiring = 1)$`Arias-Arone 2016`$Climate$Prediction
Post0.51_igraph <- igraph::graph_from_adjacency_matrix(Post0.51, weighted = TRUE, mode = "undirected")

Post11 <- FUN_SimComp(PlantAnim = NULL, RunName = "FIG2", CutOffs = CutOffs, 
                      PotPartners = RewClass_ls, Traits = meta_df,  
                      IS = 1, Rewiring = 1)$`Arias-Arone 2016`$Climate$Prediction
Post11_igraph <- igraph::graph_from_adjacency_matrix(Post11, weighted = TRUE, mode = "undirected")

png(filename = file.path(Dir.Exports, "FIG2_Arias-Arone_2016_0.5-0.5.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post0.50.5_igraph)
dev.off()
png(filename = file.path(Dir.Exports, "FIG2_Arias-Arone_2016_1-0.5.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post10.5_igraph)
dev.off()
png(filename = file.path(Dir.Exports, "FIG2_Arias-Arone_2016_0.5-1.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post0.51_igraph)
dev.off()
png(filename = file.path(Dir.Exports, "FIG2_Arias-Arone_2016_1-1.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post11_igraph)
dev.off()

# FIGURE 3, S7, S8, S9 - Multi-Network Resilience Landscape -------------
print("########## FIGURE 3, S7, S8, S9")

### Effects relative to pre-extinction network topology ----
MeanNames <- c("FIGS8", "FIG3", "FIGS7", "FIGS9")
SDNames <- c("FIGS8", "FIG3", "FIGS7", "FIGS9")
RunName = "ALL"
# for(RunName in names(PlotTopo_ls)){
Change_df <- PlotTopo_ls[[RunName]]$Change

SSP585_df <- PlotTopoClimSSP585_ls$Change
SSP585_df <- SSP585_df[SSP585_df$Proxy == "Climate", ]
SSP585_df$Proxy <- "SSP585"
SSP585_df$Casc <- "ALL"

Change_df <- rbind(Change_df, SSP585_df)

Proxy_ls <- as.list(rep(NA, length(unique(Change_df$Proxy))))
names(Proxy_ls) <- unique(Change_df$Proxy)

### Relative Changes ----
# Change_df <- PlotTopoAll_ls$Change
# Change_df <- PlotTopo_ls[[RunName]]$Change
# sink(file = file.path(Dir.Exports, paste0("MODEL_", RunName, ".txt")))
# 
# for(ModIter in TopoPlots){
#   print(paste("###################", ModIter))
#   print(
#     summary(lm(RelChange ~ 0 + (IS+RE)*Proxy, 
#                data = Change_df[model_df$Topology == ModIter, ]))
#   )
# }
# sink()

for(proxy_enum in 1:length(unique(Change_df$Proxy))){
  proxy_iter <- unique(Change_df$Proxy)[proxy_enum]
  iter_df <- Change_df[Change_df$Topology %in% TopoPlots & Change_df$Proxy %in% proxy_iter, ]
  plot_df <- aggregate(RelChange ~ Proxy+Topology+IS+RE, FUN = mean, data = iter_df)
  sd_df <- aggregate(RelChange ~ Proxy+Topology+IS+RE, FUN = sd, data = iter_df)
  
  # iterate my Topology, then fuse plots
  matplot_ls <- list(Pred = as.list(rep(NA, length(TopoPlots))), SD = as.list(rep(NA, length(TopoPlots))))
  
  for(TopoIter in 1:length(TopoPlots)){
    plot_df2 <- plot_df[plot_df$Topology == TopoPlots[TopoIter], ]
    sd_df2 <- sd_df[sd_df$Topology == TopoPlots[TopoIter], ]
    
    matplot_ls$Pred[[TopoIter]] <- ggplot(plot_df2, aes(x = RE, y = IS)) +
      geom_tile(aes(fill = RelChange)) +
      coord_fixed() + 
      facet_wrap(~Proxy+Topology, ncol = 2) + 
      theme_bw() + 
      # xlab("Probability of Rewiring Required to Realise Novel Links") + 
      # ylab("Proportion of Initial Interaction Strength Required for Continued Existence") + 
      xlab("") + 
      ylab("") + 
      guides(fill = guide_colourbar(barwidth = 2,
                                    barheight = 20,
                                    # title = "Proportional Network Metric Change",
                                    title = "",
                                    # title.position = "bottom",
                                    # direction = "horizontal",
                                    # legend.location = "bottom"
      )) + 
      scale_fill_viridis("Mean", option = "B", direction = -1) + 
      # theme(legend.position = "bottom") + 
      # ggtitle(TopoPlots[TopoIter]) + 
      theme(plot.margin = unit(c(0,0,0,0), "lines"))
    
    matplot_ls$SD[[TopoIter]] <- ggplot(sd_df2, aes(x = RE, y = IS)) +
      geom_tile(aes(fill = RelChange)) +
      coord_fixed() +
      facet_wrap(~Proxy+Topology, ncol = 2) + 
      theme_bw() + 
      # xlab("Probability of Rewiring Required to Realise Novel Links") + 
      # ylab("Proportion of Initial Interaction Strength Required for Continued Existence") + 
      xlab("") + 
      ylab("") + 
      guides(fill = guide_colourbar(barwidth = 2,
                                    barheight = 20,
                                    # title = "Proportional Network Metric Change",
                                    title = "",
                                    # title.position = "bottom",
                                    # direction = "horizontal",
                                    # legend.location = "bottom"
      )) + 
      scale_fill_viridis("SD", option = "D") + 
      # theme(legend.position = "bottom") + 
      # ggtitle(TopoPlots[TopoIter]) + 
      theme(plot.margin = unit(c(0,0,0,0), "lines"))
  }
  y.grob <- textGrob("Proportion of Initial Interaction Strength Required for Continued Existence", 
                     gp=gpar(fontface="bold", col="black", fontsize=15), rot=90)
  x.grob <- textGrob("Probability of Rewiring Required to Realise Novel Links", 
                     gp=gpar(fontface="bold", col="black", fontsize=15))
  
  MatPred_plot <- plot_grid(plotlist = matplot_ls$Pred, nrow = 1)
  # MatPred_plot <- grid.arrange(arrangeGrob(MatPred_plot, left = y.grob, bottom = x.grob))
  MatSD_plot <- plot_grid(plotlist = matplot_ls$SD, nrow = 1)
  # MatSD_plot <- grid.arrange(arrangeGrob(MatSD_plot, left = y.grob, bottom = x.grob))
  
  MatComb_plot <- plot_grid(MatPred_plot, MatSD_plot, ncol = 1, labels = "")
  MatComb_plot <- grid.arrange(arrangeGrob(MatComb_plot, left = y.grob, bottom = x.grob))
  
  # if(proxy_enum == 3){
  #   ggsave(MatPred_plot, filename = file.path(Dir.Exports, paste0(
  #     MeanNames[proxy_enum], 
  #     "_MatrixChange-", RunName, "-", proxy_iter, ".png")), width = 34/1.2, height = 32/1.2, units = "cm")
  #   ggsave(MatSD_plot, filename = file.path(Dir.Exports, paste0(
  #     SDNames[proxy_enum], 
  #     "_MatrixChange_SD-", RunName, "-", proxy_iter, ".png")), width = 34/1.2, height = 32/1.2, units = "cm") 
  #   Proxy_ls[[proxy_iter]] <- matplot_ls 
  # }else{
    ggsave(MatComb_plot, 
           filename = file.path(Dir.Exports, paste0(MeanNames[proxy_enum], 
                                                    "_MatrixChange-", RunName, "-", proxy_iter, ".png")), width = 32/1.2, height = 30/1.2, units = "cm") 
    Proxy_ls[[proxy_iter]] <- matplot_ls 
  # }
}
# RelChangePlots_ls <- Proxy_ls
# }

# FIGURE S4 - Network Data Locations -----------------------------------------
print("########## FIGURE S4")
ggplot(data = world) +
  theme_minimal_grid() + 
  geom_sf(color = "black",
          size = 0.2
          ,
          # aes(fill = factor(samplenets))
          fill = "white"
  ) + 
  # scale_fill_viridis_d(option = "B") +
  geom_point(data = map_df, aes(x = longitude, y = latitude), pch = 17, size = 3, col = "darkred") +
  xlab("Longitude") + ylab("Latitude") +
  labs(fill = "", title = "Network Sampling Sites") + 
  theme(plot.title = element_text(hjust = 0.5))
# theme(legend.position = "bottom") + 
# guides(fill = guide_legend(nrow=1))

ggsave(filename = file.path(Dir.Exports, "FIGS4_NetworksCountries.png"), width = 16/2, height = 9/2, units = "cm", scale = 3, dpi = 1e3)


# FIGURE S5 - Extinction proxy Overlap ---------------------------------------
print("########## FIGURE S5")
pal_lm <- c("#016392", "#E19825", "#3E8853")
RunName = "ALL"
# for(RunName in c("ALL", "Plants", "Animals")){
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
ggsave(filename = file.path(Dir.Exports, paste0("FIGS5_Proxy", RunName,".png")), width = 4, height = 3, units = "cm", scale = 7, dpi = 1e3)
# }

# FIGURE S10 - Effect Sizes --------------------------------------------------
print("########## FIGURE S10")
RunName = "ALL"
# for(RunName in names(PlotTopo_ls)){
## extract necessary data
# EffectSize_df <- PlotTopoAll_ls$EffectSize
# RandomSD_df <- PlotTopoAll_ls$RandomSD
EffectSize_df <- PlotTopo_ls[[RunName]]$EffectSize
RandomSD_df <- PlotTopo_ls[[RunName]]$RandomSD

## reformat to effect sizes
# Sig_df <- !(abs(EffectSize_df[,TopoPlots]) > RandomSD_df[,TopoPlots])
# Plot_df <- EffectSize_df[, TopoPlots]
# Plot_df[Sig_df] <- NA
Plot_df <- EffectSize_df[, TopoPlots] / RandomSD_df[, TopoPlots]
Plot_df$netID <- EffectSize_df$netID
Plot_df$RE <- EffectSize_df$RE
Plot_df$IS <- EffectSize_df$IS
Plot_df$Pry <- EffectSize_df$Pry

# ## make effectsizes relative to base networks
# Plot_df <- Plot_df[ , c("netID", TopoPlots, "Pry", "IS", "RE")]
# merged_df <- base::merge(Plot_df, PreExt_df[, c("netID", TopoPlots)], by = "netID")
# Rel_effSizes <- merged_df[ , 2:(length(TopoPlots)+1)] /  # predictions
#   merged_df[ , (1+2*length(TopoPlots)):ncol(merged_df)] # pre-extinctions
# colnames(Rel_effSizes) <- TopoPlots
# Rel_effSizes$netID <- Plot_df$netID
# Rel_effSizes$RE <- Plot_df$RE
# Rel_effSizes$IS <- Plot_df$IS
# Rel_effSizes$Pry <- Plot_df$Pry
# Plot_df <- Rel_effSizes

## reshape data for aggregation and formatting
Plot_df <- stats::reshape(data = Plot_df, 
                          times = colnames(Plot_df)[1:length(TopoPlots)],
                          varying = list(colnames(Plot_df)[1:length(TopoPlots)]),
                          timevar = "Topology",
                          direction = "long")
colnames(Plot_df)[ncol(Plot_df)-1] <- "EffectSize"
Plot_df <- Plot_df[which(abs(Plot_df$EffectSize) != Inf), ]
Plot_df <- Plot_df[which(!is.na(Plot_df$EffectSize)), ]

# ## removal of effectsizes where there was no deletion sequence
# ExtSpecies_ls <- FUN_SimComp(PlantAnim = NULL, RunName = RunName, IS = 1, Rewiring = 1, CutOffs = CutOffs)
# Plot_df[Plot_df$netID %in% names(which(unlist(lapply(lapply(lapply(ExtSpecies_ls, "[[", "Climate"), "[[", "Removed"), length)) == 0)) &
#           Plot_df$Pry == "Climate", "EffectSize"] <- NA
# Plot_df[Plot_df$netID %in% names(which(unlist(lapply(lapply(lapply(ExtSpecies_ls, "[[", "IUCN"), "[[", "Removed"), length)) == 0)) &
#           Plot_df$Pry == "IUCN", "EffectSize"] <- NA
# Plot_df[Plot_df$netID %in% names(which(unlist(lapply(lapply(lapply(ExtSpecies_ls, "[[", "Strength"), "[[", "Removed"), length)) == 0)) &
#           Plot_df$Pry == "Strength", "EffectSize"] <- NA
# Plot_df <- na.omit(Plot_df)         
# Plot_df <- Plot_df[, -ncol(Plot_df)]

## Plot Creation
ES_ls <- as.list(rep(NA, length(unique(Plot_df$Pry))))
names(ES_ls) <- unique(Plot_df$Pry)
for(proxy_iter in unique(Plot_df$Pry)){
  Topo_ls <- as.list(rep(NA, length(TopoPlots)))
  names(Topo_ls) <- TopoPlots
  for(topo_iter in TopoPlots){
    iter_df <- Plot_df[Plot_df$Topology == topo_iter &
                         Plot_df$Pry == proxy_iter, ]
    mean_df <- aggregate(EffectSize ~ Pry+Topology+IS+RE, FUN = mean, data = iter_df)
    Topo_ls[[topo_iter]] <- ggplot(mean_df, aes(x = RE, y = IS)) +
      geom_tile(aes(fill = EffectSize)) +
      coord_fixed() +
      facet_wrap(~Pry+Topology) +
      scale_fill_gradient2(high = "darkgreen", low = "darkred") +
      theme_bw() +  xlab("") + ylab("") +
      guides(fill = guide_colourbar(barwidth = 2, barheight = 20, title = ""))
  }
  ES_ls[[proxy_iter]] <- Topo_ls
  
  ## plotting and saving
  y.grob <- textGrob("Proportion of Initial Interaction Strength Required for Continued Existence",
                     gp=gpar(fontface="bold", col="black", fontsize=25), rot=90)
  x.grob <- textGrob("Probability of Rewiring Required to Realise Novel Links",
                     gp=gpar(fontface="bold", col="black", fontsize=25))
  MatPred_plot <- plot_grid(plotlist = ES_ls[[proxy_iter]], nrow = 1)
  ES_ls[[proxy_iter]] <- MatPred_plot
}

ggsave(ES_ls$Climate,
       width = 34/1.2, height = 30/1.2, units = "cm",
       filename = file.path(Dir.Exports, paste0("FIGS10A_EffectSizes-", RunName, ".png"))
)

ggsave(
  plot_grid(plotlist = ES_ls[c("IUCN", "Strength")], nrow = 1), 
  filename = file.path(Dir.Exports, paste0("FIGS10B+C_EffectSizes-", RunName, ".png")), 
  width = (34/1.2)*2, height = 30/1.2, units = "cm")

# EffSizePlots_ls <- Proxy_ls

## as exptected, the below produces an error. Signficiance how?!
# lm_df <- iter_df[iter_df$IS == 0.5 & iter_df$RE == 0.5, ]
# lm_df
# lme4::lmer(EffectSize ~ 1 + (1|netID), data = lm_df)
# }

# FIGURE 4 - Proxy Comparisons -----------------------------------------------
print("########## FIGURE 4")
model_df <- do.call(rbind, lapply(PlotTopo_ls, FUN = function(x){x[["Change"]]}))
model2_df <- PlotTopoClimSSP585_ls$Change
model2_df$Casc <- "ALL"
model2_df$Proxy <- "SSP585"
model_df <- rbind(model_df, model2_df)
model_df <- model_df[model_df$Casc == "ALL", ]
model_df$Proxy <- factor(model_df$Proxy)
model_df$netID <- factor(model_df$netID)

# plot_ls <- as.list(rep(NA, length(TopoPlots)))
# names(plot_ls) <- TopoPlots

Plotting_df <- data.frame(Proxy = NA, Posterior = NA, id = NA, Parameter = NA, NetworkMetric = NA)

for(model_iter in TopoPlots){
  print(model_iter)
  brms_fam <- ifelse(model_iter %in% c("Modularity", "Nestedness"), 
                     "normal", "zero_one_inflated_beta")
  iter_df <- model_df[model_df$Topology == model_iter, ]
  # iter_df <- iter_df[which(iter_df$RelChange >= 0 & iter_df$RelChange <= 1),] # necessary for nestedness models
  if(file.exists(file.path(Dir.Exports, paste0("FIG4_MODEL_", model_iter, ".RData")))){
    load(file.path(Dir.Exports, paste0("FIG4_MODEL_", model_iter, ".RData")))
  }else{
    Casc_brms <- brm(RelChange ~ 0+Proxy+IS*RE, # + (1|netID),
                     data = iter_df,
                     family = brms_fam,
                     chains = 4, cores = 4, thin = 10, iter = 1e4)
    save(Casc_brms, file = file.path(Dir.Exports, paste0("FIG4_MODEL_", model_iter, ".RData")))
  }
  
  Posterior_df <- posterior_samples(Casc_brms)
  # print(nrow(Posterior_df))
  ProxyIntercept_df <- inv_logit_scaled(Posterior_df[, 1:4]) #exp
  ProxyIntercept_df <- stats::reshape(ProxyIntercept_df,
                                      direction = "long",
                                      varying = colnames(ProxyIntercept_df),
                                      times = colnames(ProxyIntercept_df),
                                      timevar = "Proxy",
                                      v.names = "Posterior",
  )
  ProxyIntercept_df$Parameter = "Intercept"
  ProxyIntercept_df$Proxy <- gsub(ProxyIntercept_df$Proxy,
                                  pattern = "b_ProxyStrength", replacement = "Centrality")
  ProxyIntercept_df$Proxy <- gsub(ProxyIntercept_df$Proxy,
                                  pattern = "b_ProxyClimate", replacement = "SSP245")
  ProxyIntercept_df$Proxy <- gsub(ProxyIntercept_df$Proxy,
                                  pattern = "b_ProxySSP585", replacement = "SSP585")
  ProxyIntercept_df$Proxy <- gsub(ProxyIntercept_df$Proxy,
                                  pattern = "b_ProxyIUCN", replacement = "IUCN")
  ProxyIntercept_df$NetworkMetric <- model_iter
  
  ISRE_df <- inv_logit_scaled(Posterior_df[, c("b_IS", "b_RE", "b_IS:RE")])
  ISRE_df <- stats::reshape(ISRE_df,
                                      direction = "long",
                                      varying = colnames(ISRE_df),
                                      times = colnames(ISRE_df),
                                      timevar = "Proxy",
                                      v.names = "Posterior",
  )
  ISRE_df$Parameter = "Resilience"
  ISRE_df$Proxy <- gsub(ISRE_df$Proxy,
                                  pattern = "b_IS", replacement = "LLS")
  ISRE_df$Proxy <- gsub(ISRE_df$Proxy,
                                  pattern = "b_RE", replacement = "RPT")
  ISRE_df$Proxy <- gsub(ISRE_df$Proxy,
                                  pattern = "LLS:RE", replacement = "Interaction")
  ISRE_df$NetworkMetric <- model_iter
  
  Plotting_df <- rbind(Plotting_df, ProxyIntercept_df, ISRE_df)
  
  # plot_ls[[model_iter]] <- ggplot(ProxyIntercept_df, aes(x = Posterior, 
  #                                                        y = Proxy)) +
  #   stat_halfeye() +
  #   theme_bw() + 
  #   geom_vline(xintercept = 0, col = "black") +
  #   labs(title = model_iter) + 
  #   theme(text = element_text(size = 20)) + theme(plot.margin = unit(c(0,3,0,3), "lines"))
  # 
  # if(model_iter %in% c("n_links", "Nestedness")){
  #   plot_ls[[model_iter]] <- plot_ls[[model_iter]] + labs(y = NULL)
  # }
}
Plotting_df <- Plotting_df[-1,]
Plotting_df$Proxy <- factor(Plotting_df$Proxy, 
                                  levels = c("Centrality", "SSP585", "SSP245", "IUCN", "Interaction", "RPT", "LLS"))
Plot_ls <- as.list(rep(NA, length(unique(Plotting_df$Parameter))))
names(Plot_ls) <- unique(Plotting_df$Parameter)
for(i in names(Plot_ls)){
  Plot_ls[[i]] <- ggplot(Plotting_df[Plotting_df$Parameter == i, ], 
                         aes(x = Posterior, y = Proxy)
  ) +
    stat_halfeye() +
    theme_bw() +
    facet_wrap(~factor(NetworkMetric, levels = TopoPlots), scales = "free_x") + 
    geom_vline(xintercept = 0, col = "black") + 
    theme(text = element_text(size = 20))
  # labs(title = model_iter) +
  # theme(text = element_text(size = 20)) + theme(plot.margin = unit(c(0,3,0,3), "lines"))
}

ggsave(
  cowplot::plot_grid(plotlist = Plot_ls, ncol = 1),
  filename = file.path(Dir.Exports, "FIG4_ProxyComparison-ALL.png"), 
  width = 30/0.6, height = 34/1.2, units = "cm")

# FIGURE 5 - Cascade Comparisons ---------------------------------------------
print("########## FIGURE 5")
model_df <- do.call(rbind, lapply(PlotTopo_ls, FUN = function(x){x[["Change"]]}))
model_df <- model_df[model_df$Proxy == "Climate", ]
model_df$Casc <- factor(model_df$Casc)
model_df$netID <- factor(model_df$netID)

plot_ls <- as.list(rep(NA, length(TopoPlots)))
names(plot_ls) <- TopoPlots

Plotting_df <- data.frame(Cascade = NA, Posterior = NA, id = NA, Parameter = NA, NetworkMetric = NA)

for(model_iter in TopoPlots){
  print(model_iter)
  brms_fam <- ifelse(model_iter %in% c("Modularity", "Nestedness"), 
                     "normal", "zero_one_inflated_beta")
  iter_df <- model_df[model_df$Topology == model_iter, ]
  # iter_df <- iter_df[which(iter_df$RelChange >= 0 & iter_df$RelChange <= 1),] # necessary for nestedness models
  if(file.exists(file.path(Dir.Exports, paste0("FIG5_MODEL_", model_iter, ".RData")))){
    load(file.path(Dir.Exports, paste0("FIG5_MODEL_", model_iter, ".RData")))
  }else{
    Casc_brms <- brm(RelChange ~ 0+Casc+IS*RE,
                     data = iter_df,
                     family = brms_fam,
                     chains = 4, cores = 4, thin = 10, iter = 1e4)
    save(Casc_brms, file = file.path(Dir.Exports, paste0("FIG5_MODEL_", model_iter, ".RData")))
  }
  
  Posterior_df <- posterior_samples(Casc_brms)
  
  CascIntercept_df <- inv_logit_scaled(Posterior_df[, 1:3])
  CascIntercept_df <- stats::reshape(CascIntercept_df,
                                     direction = "long",
                                     varying = colnames(CascIntercept_df),
                                     times = colnames(CascIntercept_df),
                                     timevar = "Cascade",
                                     v.names = "Posterior",
  )
  CascIntercept_df$Parameter = "Intercept"
  CascIntercept_df$Cascade <- gsub(CascIntercept_df$Cascade,
                                   pattern = "b_CascALL", replacement = "Bidirectional")
  CascIntercept_df$Cascade <- gsub(CascIntercept_df$Cascade,
                                   pattern = "b_CascPlants", replacement = "Bottom-Up")
  CascIntercept_df$Cascade <- gsub(CascIntercept_df$Cascade,
                                   pattern = "b_CascAnimals", replacement = "Top-Down")
  CascIntercept_df$NetworkMetric <- model_iter
  
  ISRE_df <- inv_logit_scaled(Posterior_df[, c("b_IS", "b_RE", "b_IS:RE")])
  ISRE_df <- stats::reshape(ISRE_df,
                            direction = "long",
                            varying = colnames(ISRE_df),
                            times = colnames(ISRE_df),
                            timevar = "Cascade",
                            v.names = "Posterior",
  )
  ISRE_df$Parameter = "Resilience"
  ISRE_df$Cascade <- gsub(ISRE_df$Cascade,
                        pattern = "b_IS", replacement = "LLS")
  ISRE_df$Cascade <- gsub(ISRE_df$Cascade,
                        pattern = "b_RE", replacement = "RPT")
  ISRE_df$Cascade <- gsub(ISRE_df$Cascade,
                        pattern = "LLS:RE", replacement = "Interaction")
  ISRE_df$NetworkMetric <- model_iter
  
  Plotting_df <- rbind(Plotting_df, CascIntercept_df, ISRE_df)
  # plot_ls[[model_iter]] <- ggplot(CascIntercept_df, aes(x = Posterior, y = Cascade)) +
  #   stat_halfeye() +
  #   theme_bw() + 
  #   # geom_vline(xintercept = 0, col = "red") + 
  #   labs(title = model_iter) + 
  #   theme(text = element_text(size = 20)) + theme(plot.margin = unit(c(0,3,0,3), "lines"))
  # 
  # if(model_iter %in% c("n_links", "Nestedness")){
  #   plot_ls[[model_iter]] <- plot_ls[[model_iter]] + labs(y = NULL)
  # }
}

Plotting_df <- Plotting_df[-1,]
Plotting_df$Cascade <- factor(Plotting_df$Cascade, 
                            levels = c("Bidirectional", "Bottom-Up", "Top-Down", "Interaction", "RPT", "LLS"))
Plot_ls <- as.list(rep(NA, length(unique(Plotting_df$Parameter))))
names(Plot_ls) <- unique(Plotting_df$Parameter)
for(i in names(Plot_ls)){
  Plot_ls[[i]] <- ggplot(Plotting_df[Plotting_df$Parameter == i, ], 
                         aes(x = Posterior, y = Cascade)
  ) +
    stat_halfeye() +
    theme_bw() +
    facet_wrap(~factor(NetworkMetric, levels = TopoPlots), scales = "free_x") + 
    geom_vline(xintercept = 0, col = "black") + 
    theme(text = element_text(size = 20))
  # labs(title = model_iter) +
  # theme(text = element_text(size = 20)) + theme(plot.margin = unit(c(0,3,0,3), "lines"))
}

ggsave(
  cowplot::plot_grid(plotlist = Plot_ls, ncol = 1),
  filename = file.path(Dir.Exports, "FIG5_CascComparison-Climate.png"), 
  width = 30/0.6, height = 34/1.2, units = "cm")


# FIGURE S11 & S12 - Cascade Comparison Matrices -----------------------------
print("########## FIGURE S11 & S12")
RunNames <- c("Plants", "Animals")
FigNames <- c("FigS11", "FigS12")
for(RunIter in 1:2){
  RunName <- RunNames[RunIter]
  FigName <- FigNames[RunIter]
  ### Baseline ALL Changes ----
  Change_df <- PlotTopo_ls[["ALL"]]$Change
  Change_df <- Change_df[Change_df$Topology %in% TopoPlots & Change_df$Proxy %in% c("Climate", "IUCN", "Strength"), ]
  base_plot_df <- aggregate(RelChange ~ Proxy+Topology+IS+RE, FUN = mean, data = Change_df)
  base_sd_df <- aggregate(RelChange ~ Proxy+Topology+IS+RE, FUN = sd, data = Change_df)
  
  ### Relative Changes ----
  Change_df <- PlotTopo_ls[[RunName]]$Change
  Change_df <- Change_df[Change_df$Topology %in% TopoPlots & Change_df$Proxy %in% c("Climate", "IUCN", "Strength"), ]
  test_plot_df <- aggregate(RelChange ~ Proxy+Topology+IS+RE, FUN = mean, data = Change_df)
  test_sd_df <- aggregate(RelChange ~ Proxy+Topology+IS+RE, FUN = sd, data = Change_df)
  
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
      facet_wrap(~Proxy, ncol = 3) +
      theme_bw() +
      # xlab("Probability of Rewiring Required to Realise Novel Links") +
      # ylab("Proportion of Initial Interaction Strength Required for Continued Existence") +
      xlab("") +
      ylab("") +
      guides(fill = guide_colourbar(barwidth = 1.5,
                                    barheight = 12,
                                    title = "", #Proportional Network \n Metric Change
                                    title.position = "bottom",
                                    # direction = "horizontal",
                                    legend.location = "bottom")) +
      scale_fill_gradient2("SD", low = "#D55E00", high = "#009E73") + 
      # theme(legend.position = "bottom") + 
      ggtitle(TopoPlots[TopoIter]) 
    # +
    #   theme(plot.margin = unit(c(0,0,0,0), "lines"))
    
    matplot_ls$SD[[TopoIter]] <- ggplot(sd_df2, aes(x = RE, y = IS)) +
      geom_tile(aes(fill = RelChange)) +
      coord_fixed() +
      facet_wrap(~Proxy, ncol = 3) +
      theme_bw() +
      # xlab("Probability of Rewiring Required to Realise Novel Links") +
      # ylab("Proportion of Initial Interaction Strength Required for Continued Existence") +
      xlab("") +
      ylab("") +
      guides(fill = guide_colourbar(barwidth = 1.5,
                                    barheight = 12,
                                    title = "", #Proportional Network \n Metric Change
                                    title.position = "bottom",
                                    # direction = "horizontal",
                                    legend.location = "bottom")) +
      scale_fill_gradient2("SD", low = "#D55E00", high = "#009E73") + 
      # theme(legend.position = "bottom") + 
      ggtitle(TopoPlots[TopoIter]) 
    # +
    #   theme(plot.margin = unit(c(0,0,0,0), "lines"))
  }
  y.grob <- textGrob("Proportion of Initial Interaction Strength Required for Continued Existence",
                     gp=gpar(fontface="bold", col="black", fontsize=12), rot=90)
  x.grob <- textGrob("Probability of Rewiring Required to Realise Novel Links",
                     gp=gpar(fontface="bold", col="black", fontsize=12))
  
  MatPred_plot <- plot_grid(plotlist = matplot_ls$Pred, nrow = 2)
  MatPred_plot <- grid.arrange(arrangeGrob(MatPred_plot, left = y.grob, bottom = x.grob))
  ggsave(MatPred_plot, filename = file.path(Dir.Exports, paste0(FigName, "_MatrixChange", RunName, "-ComparedtoALL.png")), width = 28/1.2, height = 22/1.2, units = "cm")
  # MatSD_plot <- plot_grid(plotlist = matplot_ls$SD, nrow = 1)
  # MatSD_plot <- grid.arrange(arrangeGrob(MatSD_plot, left = y.grob, top = x.grob))
  # ggsave(MatSD_plot, filename = file.path(Dir.Exports, paste0("PLOT_MatrixChange_SD", RunName, "-ComparedtoALL.png")), width = 44/1.2, height = 34/1.2, units = "cm")
}


# FIGURE S13 - Cascade Species Group Vulnerability ---------------------------
## animals
Removal_df <- PlotTopo_ls[["Animals"]]$Topology
Removal_df <- Removal_df[Removal_df$Simulation == "Prediction" &
                           Removal_df$IS == 0 & Removal_df$RE == 0, c("netID", "Proxy", "Removed")]

Change_df <- PlotTopo_ls[["Animals"]]$Change
Change_df <- Change_df[Change_df$Topology == "n_plants", ]

merged_df <- merge.data.frame(x = Change_df, y = Removal_df, by = c("netID", "Proxy"))

## plants
Removal_df <- PlotTopo_ls[["Plants"]]$Topology
Removal_df <- Removal_df[Removal_df$Simulation == "Prediction" &
                           Removal_df$IS == 0 & Removal_df$RE == 0, c("netID", "Proxy", "Removed")]
Change_df <- PlotTopo_ls[["Plants"]]$Change
Change_df <- Change_df[Change_df$Topology == "n_animals", ]

merged_df2 <- merge.data.frame(x = Change_df, y = Removal_df, by = c("netID", "Proxy"))

## plotting
plot_df <- rbind(merged_df, merged_df2)
plot_df$Value <- plot_df$RelChange/plot_df$Removed
plot_df$Value[which(is.nan(plot_df$Value))] <- NA


sd_df <- aggregate(Value ~ Proxy+Casc+IS+RE, FUN = sd, data = plot_df)
plot_df <- aggregate(Value ~ Proxy+Casc+IS+RE, FUN = mean, data = plot_df)

Pred_gplot <- ggplot(plot_df, aes(x = RE, y = IS)) +
  geom_tile(aes(fill = Value)) +
  coord_fixed() +
  facet_wrap(~Proxy+Casc, ncol = 2) +
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
  theme(legend.position = "bottom")

SD_gplot <- ggplot(sd_df, aes(x = RE, y = IS)) +
  geom_tile(aes(fill = Value)) +
  coord_fixed() +
  facet_wrap(~Proxy+Casc, ncol = 2) +
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
  theme(legend.position = "bottom")

y.grob <- textGrob("Proportion of Initial Interaction Strength Required for Continued Existence",
                   gp=gpar(fontface="bold", col="black", fontsize=15), rot=90)
x.grob <- textGrob("Probability of Rewiring Required to Realise Novel Links",
                   gp=gpar(fontface="bold", col="black", fontsize=15))

Export_plot <- grid.arrange(arrangeGrob(plot_grid(Pred_gplot, SD_gplot, nrow = 1), left = y.grob, top = x.grob))
ggsave(Export_plot, filename = file.path(Dir.Exports, paste0("FIGS13", "_LossOfSpecReltoOtherGroupPrimLoss.png")), width = 40/1.2, height = 34/1.2, units = "cm")

