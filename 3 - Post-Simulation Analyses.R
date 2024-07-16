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
sf::sf_use_s2(FALSE)
world$samplenets <- lengths(st_intersects(world, DT))

PlotTopo_ls <- list(SSP245 = PlotTopoSSP245_ls,
                    SSP585 = PlotTopoSSP585_ls
)
PlotTopo_ls$SSP245$Change$SSP <- 245
PlotTopo_ls$SSP585$Change$SSP <- 585

# CONCEPT ARIAS-ARONE One-Network Resilience Landscape ---------------------
print("########## Arias-Arone (Single-Network Panels)")
AnalysisData_ls <- AnalysisData_ls[28] # loading data for single network
RM_spec <- FUN_SimComp(PlantAnim = NULL, RunName = "FIG2", CutOffs = CutOffs, 
                       PotPartners = RewClass_ls, Traits = meta_df, 
                       IS = 1, Rewiring = 1, WHICH = c("SSP245"))$`Arias-Arone 2016`$Climate$Removed # removed species in target network

### Matrix Visualisation ----
Change_df <- PlotTopo_ls$SSP245$Change
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

# RESILIENCE LANDSCAPE MATRICES --------------------------------------------
print("########## Resilience Landscape Matrices")

lapply(PlotTopo_ls, FUN = function(SSP_iter){
  Change_df <- SSP_iter$Change
  Change_df <- Change_df[Change_df$Topology %in% TopoPlots, ]
  ## also remove networks for which no species werre removed here:
  EffectSize_df <- SSP_iter$EffectSize[SSP_iter$EffectSize$netID %nin% unique(SSP_iter$Topology$netID[SSP_iter$Topology$Removed == 0]),]
  RandomSD_df <- SSP_iter$RandomSD[SSP_iter$RandomSD$netID %nin% unique(SSP_iter$Topology$netID[SSP_iter$Topology$Removed == 0]),]
  
  ## landscapes ----
  ## calc mean and sd per location and driver
  Mean_df <- aggregate(x = Change_df, RelChange ~ Topology + IS + RE + Proxy, FUN = mean)
  SD_df <- aggregate(x = Change_df, RelChange ~ Topology + IS + RE + Proxy, FUN = sd)

  Mean_gg <- ggplot(Mean_df, 
         aes(x = RE, y = IS)) +
    geom_tile(aes(fill = RelChange)) +
    coord_fixed() + 
    facet_grid(Proxy~factor(Topology, levels = TopoPlots)) + 
    theme_bw() + 
    xlab("") + 
    ylab("") + 
    guides(fill = guide_colourbar(barwidth = 23.5,
                                  barheight = 2,
                                  # title = "Proportional Network Metric Change",
                                  title = "",
                                  direction = "horizontal"
    )) + 
    scale_fill_viridis("Mean", option = "B", direction = -1) + 
    theme(plot.margin = unit(c(0,0,0,0), "lines"), legend.position = "bottom")
  
  SD_gg <- ggplot(SD_df, 
                  aes(x = RE, y = IS)) +
    geom_tile(aes(fill = RelChange)) +
    coord_fixed() + 
    facet_grid(Proxy~factor(Topology, levels = TopoPlots)) + 
    theme_bw() + 
    xlab("") + 
    ylab("") + 
    guides(fill = guide_colourbar(barwidth = 23.5,
                                  barheight = 2,
                                  # title = "Proportional Network Metric Change",
                                  title = "",
                                  direction = "horizontal"
    )) + 
    scale_fill_viridis("SD", option = "D") + 
    theme(plot.margin = unit(c(0,0,0,0), "lines"), legend.position = "bottom")
  
  y.grob <- textGrob("Proportion of Initial Interaction Strength Required for Continued Existence", 
                     gp=gpar(fontface="bold", col="black", fontsize=15), rot=90)
  x.grob <- textGrob("Probability of Rewiring Required to Realise Novel Links", 
                     gp=gpar(fontface="bold", col="black", fontsize=15))
  Mats_gg <- cowplot::plot_grid(Mean_gg, ggplot()+theme_void(), SD_gg, ncol = 3, 
                                rel_widths = c(1, 0.15, 1))
  MatComb_plot <- grid.arrange(arrangeGrob(Mats_gg, left = y.grob, top = x.grob))
  
  ggsave(MatComb_plot, 
         filename = file.path(Dir.Exports, 
                              paste0("FIG_Matrices_SSP", unique(Change_df$SSP), ".png")),
         width = 34/1.2, height = 28/1.2, units = "cm")
  
  ## comparison to contemporary ----
  ContempComparison_ls <- lapply(TopoPlots, FUN = function(Topoiter){
    Contemp_df <- Change_df[Change_df$Topology == Topoiter & 
                              Change_df$IS == 0 & Change_df$RE == 1, ]
    ContempComp <- Contemp_df[, c("RelChange", "netID", "Proxy")]
    
    Comgrid <- expand.grid(unique(Change_df$RE), unique(Change_df$IS), unique(Change_df$Proxy))
    colnames(Comgrid) <- c("RE", "IS", "Proxy")
    
    contempttest <- apply(Comgrid, 1, FUN = function(rowiter){
      # print(rowiter)
      REiter <- rowiter[1]
      ISiter <- rowiter[2]
      Proxyiter <- rowiter[3]
      
      ToComp <- Change_df[Change_df$Topology == Topoiter & 
                            Change_df$RE == as.numeric(REiter) & 
                            Change_df$IS == as.numeric(ISiter) & 
                            Change_df$Proxy == as.character(Proxyiter), 
                          c("RelChange", "netID")]
      colnames(ToComp)[1] <- "CompRel"
      testdf <- merge(ContempComp[ContempComp$Proxy == as.character(Proxyiter), ], ToComp)
      
      
      
      
      ttest <- wilcox.test(testdf$RelChange, testdf$CompRel, paired = TRUE)
      # ttest <- t.test(testdf$RelChange, testdf$CompRel, paired = TRUE)
      data.frame(
        meandiff = mean(testdf$RelChange-testdf$CompRel),
        Topo = Topoiter,
        RE = REiter,
        IS = ISiter,
        Proxy = Proxyiter,
        testval = ttest$statistic,
        testp = ttest$p.value)
    })
    do.call(rbind, contempttest)
  })
  ttestComp_df <- do.call(rbind, ContempComparison_ls)
  ttestComp_df$Sig <- ttestComp_df$testp < 0.05
  ttestComp_df$Sig[which(!ttestComp_df$Sig)] <- NA
  
  TComp_gg <- ggplot(data = ttestComp_df, 
         aes(x = as.numeric(RE), y = as.numeric(IS),
             fill = meandiff
             )
         ) + 
    geom_tile(aes(col = Sig), linewidth = 0.2) +
    scale_color_manual(values = "black", na.translate = FALSE, name = "Signifance") + 
    # geom_point(aes(shape = Sig), size = 4.5) +
    # scale_shape_manual(values=c(0), na.translate = FALSE, name = "Significance") +
    geom_point(aes(shape = factor(sign(meandiff))), size = 2) +
    scale_shape_manual(values=c(95, 61, 43), na.translate = FALSE, name = "Over-/Underprediction") + 
    coord_fixed() + 
    facet_grid(factor(Topo, levels = TopoPlots) ~ Proxy) + 
    theme_bw() + 
    xlab("") + 
    ylab("") + 
    guides(fill = guide_colourbar(barwidth = 23.5,
                                  barheight = 2,
                                  # title = "Proportional Network Metric Change",
                                  title = "",
                                  direction = "horizontal"
    )) + 
    scale_fill_gradient2(low = "darkred", high = "forestgreen") + 
    theme(plot.margin = unit(c(0,0,0,0), "lines"), legend.position = "bottom")
  TComp_gg <- grid.arrange(arrangeGrob(TComp_gg, left = y.grob, top = x.grob))
  ggsave(TComp_gg, 
         filename = file.path(Dir.Exports, 
                              paste0("FIG_ContempComp_SSP", unique(Change_df$SSP), ".png")),
         width = 34/1.2, height = 28/1.2, units = "cm")
  
  ## effect sizes ----
  ## reformat to effect sizes
  Plot_df <- EffectSize_df[, TopoPlots] / RandomSD_df[, TopoPlots]
  Plot_df$netID <- EffectSize_df$netID
  Plot_df$RE <- EffectSize_df$RE
  Plot_df$IS <- EffectSize_df$IS
  Plot_df$Proxy <- EffectSize_df$Proxy

  ## reshape data for aggregation and formatting
  Plot_df <- stats::reshape(data = Plot_df,
                            times = colnames(Plot_df)[1:length(TopoPlots)],
                            varying = list(colnames(Plot_df)[1:length(TopoPlots)]),
                            timevar = "Topology",
                            direction = "long")
  colnames(Plot_df)[ncol(Plot_df)-1] <- "EffectSize"
  Plot_df <- Plot_df[which(abs(Plot_df$EffectSize) != Inf), ]
  Plot_df <- Plot_df[which(!is.na(Plot_df$EffectSize)), ]

  ## Plot Creation
  Mean_df <- aggregate(x = Plot_df, EffectSize ~ Topology + IS + RE + Proxy, FUN = mean)
  SD_df <- aggregate(x = Plot_df, EffectSize ~ Topology + IS + RE + Proxy, FUN = sd)
  
  topo_ls <- lapply(unique(Mean_df$Topology), FUN = function(topoiter){
    prx_ls <- lapply(unique(Mean_df$Proxy), FUN = function(prox_iter){
      Mean_gg <- ggplot(Mean_df[Mean_df$Proxy == prox_iter & 
                                  Mean_df$Topology == topoiter, ], 
                        aes(x = RE, y = IS)) +
        geom_tile(aes(fill = EffectSize)) +
        coord_fixed() + 
        facet_grid(Proxy~factor(Topology, levels = TopoPlots)) + 
        theme_bw() + 
        xlab("") + 
        ylab("") + 
        guides(fill = guide_colourbar(barwidth = 23.5,
                                      barheight = 2,
                                      # title = "Proportional Network Metric Change",
                                      title = "",
                                      direction = "horizontal"
        )) + 
        scale_fill_gradient2(low = "darkred", high = "forestgreen") + 
        theme(plot.margin = unit(c(0,0,0,0), "lines"), legend.position = "bottom")
      
      SD_gg <- ggplot(SD_df[SD_df$Proxy == prox_iter & 
                              SD_df$Topology == topoiter, ], 
                      aes(x = RE, y = IS)) +
        geom_tile(aes(fill = EffectSize)) +
        coord_fixed() + 
        facet_grid(Proxy~factor(Topology, levels = TopoPlots)) + 
        theme_bw() + 
        xlab("") + 
        ylab("") + 
        guides(fill = guide_colourbar(barwidth = 23.5,
                                      barheight = 2,
                                      # title = "Proportional Network Metric Change",
                                      title = "",
                                      direction = "horizontal"
        )) + 
        scale_fill_viridis("SD", option = "D") + 
        theme(plot.margin = unit(c(0,0,0,0), "lines"), legend.position = "bottom")
      
      cowplot::plot_grid(Mean_gg, SD_gg, ncol = 2)
    })
    
    Eff_gg <- cowplot::plot_grid(plotlist = prx_ls, ncol = 1)
    MatComb_plot <- grid.arrange(arrangeGrob(Eff_gg, left = y.grob, top = x.grob))
    
    ggsave(MatComb_plot, 
           filename = file.path(Dir.Exports, 
                                paste0("FIG_EffectSizes_SSP", unique(Change_df$SSP), "_", topoiter, ".png")),
           width = 34/1.2, height = 60/1.2, units = "cm")
  })
})

# NETWORK DATA LOCATIONS ---------------------------------------------------
print("########## FIGURE S4")
ggplot(data = world) +
  theme_minimal_grid() + 
  geom_sf(color = "black",
          size = 0.2,
          fill = "white"
  ) + 
  geom_point(data = map_df, aes(x = longitude, y = latitude), pch = 17, size = 3, col = "darkred") +
  xlab("Longitude") + ylab("Latitude") +
  labs(fill = "", title = "Network Sampling Sites") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = file.path(Dir.Exports, "FIGS4_NetworksCountries.png"), width = 16/2, height = 9/2, units = "cm", scale = 3, dpi = 1e3)

# PRIMARY PROXY COMPARISONS ------------------------------------------------
print("########## PRIMARY PROXY MODELS")
model_df <- do.call(rbind, lapply(PlotTopo_ls, FUN = function(x){x[["Change"]]}))
model_df$Proxy <- factor(model_df$Proxy)
model_df$SSP <- factor(model_df$SSP)
model_df$netID <- factor(model_df$netID)

# plot_ls <- as.list(rep(NA, length(TopoPlots)))
# names(plot_ls) <- TopoPlots
Plotting_df <- data.frame(Proxy = NA, Posterior = NA, id = NA, Parameter = NA, 
                          SSP = 585, Driver = NA, NetworkMetric = NA, PlotYLabels = NA)

for(model_iter in TopoPlots){
  print(model_iter)
  brms_fam <- ifelse(model_iter %in% c("Modularity", "Nestedness"), 
                     "normal", "zero_one_inflated_beta")
  iter_df <- model_df[model_df$Topology == model_iter, ]
  iter_df$PrimProx <- paste(iter_df$SSP, iter_df$Proxy, sep = "-")
  
  if(file.exists(file.path(Dir.Exports, paste0("PROXCOMP_", model_iter, ".RData")))){
    load(file.path(Dir.Exports, paste0("PROXCOMP_", model_iter, ".RData")))
  }else{
    Prox_brms <- brm(RelChange ~ 0+PrimProx+IS*RE, # + (1|netID),
                     data = iter_df,
                     family = brms_fam,
                     chains = 4, cores = 4, thin = 10, iter = 1e4)
    save(Prox_brms, file = file.path(Dir.Exports, paste0("PROXCOMP_", model_iter, ".RData")))
  }
  
  Posterior_df <- posterior_samples(Prox_brms, pars = "^b_")
  # print(nrow(Posterior_df))
  ProxyIntercept_df <- inv_logit_scaled(Posterior_df) #exp
  ProxyIntercept_df <- stats::reshape(ProxyIntercept_df,
                                      direction = "long",
                                      varying = colnames(ProxyIntercept_df),
                                      times = colnames(ProxyIntercept_df),
                                      timevar = "Proxy",
                                      v.names = "Posterior",
  )
  ProxyIntercept_df$Parameter <- ifelse(
    startsWith(ProxyIntercept_df$Proxy, "b_PrimProx"), "Drivers", "Resilience"
  )
  ProxyIntercept_df$SSP <- factor(gsub(".*?([0-9]+).*", "\\1", ProxyIntercept_df$Proxy))
  ProxyIntercept_df$Driver <- ProxyIntercept_df$Proxy
  ProxyIntercept_df$Driver <- gsub(pattern = "b_PrimProx245M", replacement = "", ProxyIntercept_df$Driver)
  ProxyIntercept_df$Driver <- gsub(pattern = "b_PrimProx585M", replacement = "", ProxyIntercept_df$Driver)
  
  # ProxyIntercept_df$Proxy <- gsub(ProxyIntercept_df$Proxy,
  #                                 pattern = "b_ProxyStrength", replacement = "Centrality")
  # ProxyIntercept_df$Proxy <- gsub(ProxyIntercept_df$Proxy,
  #                                 pattern = "b_ProxyClimate", replacement = "SSP245")
  # ProxyIntercept_df$Proxy <- gsub(ProxyIntercept_df$Proxy,
  #                                 pattern = "b_ProxySSP585", replacement = "SSP585")
  # ProxyIntercept_df$Proxy <- gsub(ProxyIntercept_df$Proxy,
  #                                 pattern = "b_ProxyIUCN", replacement = "IUCN")
  ProxyIntercept_df$NetworkMetric <- model_iter
  
  # ISRE_df <- inv_logit_scaled(Posterior_df[, c("b_IS", "b_RE", "b_IS:RE")])
  # ISRE_df <- stats::reshape(ISRE_df,
  #                                     direction = "long",
  #                                     varying = colnames(ISRE_df),
  #                                     times = colnames(ISRE_df),
  #                                     timevar = "Proxy",
  #                                     v.names = "Posterior",
  # )
  
  ProxyIntercept_df$Proxy <- gsub(ProxyIntercept_df$Proxy,
                                  pattern = "b_IS", replacement = "LLS")
  ProxyIntercept_df$Proxy <- gsub(ProxyIntercept_df$Proxy,
                                  pattern = "b_RE", replacement = "RPT")
  ProxyIntercept_df$Proxy <- gsub(ProxyIntercept_df$Proxy,
                                  pattern = "LLS:RE", replacement = "Interaction")
  # ISRE_df$NetworkMetric <- model_iter
  
  ProxyIntercept_df$PlotYLabels <- factor(paste(ProxyIntercept_df$SSP, 
                                                ProxyIntercept_df$Driver, 
                                                sep = " - "))
  
  Plotting_df <- 
    # rbind(Plotting_df, 
    ProxyIntercept_df
  # )
  
  mod_gg <- cowplot::plot_grid(
    ggplot(Plotting_df[Plotting_df$Parameter == "Drivers",], 
           aes(x = Posterior, 
               y = factor(Driver, levels = rev(levels(factor(Driver)))))) + 
      stat_halfeye(normalize = "xy") +
      geom_vline(xintercept = 0) +
      facet_wrap(~SSP, ncol = 1, scales = "free") +
      lims(x = c(0, max(Plotting_df[Plotting_df$Parameter == "Drivers","Posterior"])+0.05)) +
      theme_bw() + labs(y = "Primary Extinction Risk Proxy") + 
      theme(strip.text = element_text(size = 20))
    ,
    ggplot() + theme_void(),
    ggplot(Plotting_df[Plotting_df$Parameter == "Resilience",], 
           aes(x = Posterior, 
               y = factor(Proxy, levels = rev(c("LLS", "RPT", "Interaction"))))) + 
      stat_halfeye(normalize = "xy") +
      geom_vline(xintercept = 0) +
      theme_bw() + labs(y = "Extinction Cascade Resilience Mechanisms"),
    ncol = 3, rel_widths = c(1, 0.05, 1.2)
  )
  ggsave(
    mod_gg,
    filename = file.path(Dir.Exports, 
                         paste0("FIG_ProxyComparison_", 
                                model_iter, "_.png")), 
    width = 13*2.4, height = 6.5*2.4, units = "cm")
  
}