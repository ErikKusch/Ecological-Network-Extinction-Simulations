#' ####################################################################### #
#' PROJECT: [Biodiversity Simplification & Ecological Network Topology] 
#' CONTENTS: 
#'  - Conceptual Visualisations
#'  DEPENDENCIES:
#'  - "0 - Preamble.R"
#'  - "0 - Fricke_Functions.R"
#'  - "1 - DataRetrieval.R" has to have been run and produced "AnalysesData.RData" in "Dir.Data" directory
#'  - "2 - Extinction Simulation.R" has to have been run
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE =================================================================
rm(list=ls())
set.seed(42)
# CutOffs <- list(Strength = 0.75,
#                 Climate = 2,
#                 IUCN = 5)

## Sourcing ----------------------------------------------------------------
source("0 - Preamble.R")
source("0 - Fricke_Functions.R")

source2 <- function(file, start, end, ...) {
  file.lines <- scan(file, what=character(), skip=start-1, nlines=end-start+1, sep='\n')
  file.lines.collapsed <- paste(file.lines, collapse='\n')
  source(textConnection(file.lines.collapsed), ...)
}

source2(file = "2 - Extinction Simulation.R", start = 1, end = 204)

# MAP OF NETWORK LOCATIONS =================================================
## DATA LOADING ------------------------------------------------------------
message("### LOADING NETWORK DATA ###")
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

## PLOTTING ----------------------------------------------------------------
ggplot(data = world) +
  theme_minimal_grid() + 
  geom_sf(color = "black",
          size = 0.2,
          aes(fill = factor(samplenets))
          # fill = "grey" 
          ) + 
  scale_fill_viridis_d(option = "B") + 
  # geom_point(data = map_df, aes(x = longitude, y = latitude), pch = 13, size = 2, col = "red") +
  xlab("Longitude") + ylab("Latitude") +
  labs(fill = "", title = "Number of Networks per Country") + 
  theme(plot.title = element_text(hjust = 0.5))
  # theme(legend.position = "bottom") + 
  # guides(fill = guide_legend(nrow=1))

ggsave(filename = file.path(Dir.Exports, "NetworksCountries.png"), width = 16/2, height = 9/2, units = "cm", scale = 3, dpi = 1e3)

# leaflet() %>%
#   # Base groups
#   addTiles() %>%
#   # Overlay groups
#   addCircles(data = map_df, ~longitude, ~latitude,
#                    group = "net.id", 
#                    label = ~study.id,
#                    color ="red",
#              radius = 1e2
#   )

# NETWORK EXTINCTION THRESHOLDS AND REWIRING ===============================
## DATA LOADING ------------------------------------------------------------
AnalysisData_ls <- AnalysisData_ls[28]
RM_spec <- FUN_SimComp(PlantAnim = NULL, RunName = "CONCEPT", CutOffs = CutOffs, 
                       PotPartners = RewClass_ls, Traits = meta_df, 
                       IS = 1, Rewiring = 1, WHICH = "Climate")$`Arias-Arone 2016`$Climate$Removed # removed species in target network
PlotTopoAll_ls <- loadObj(file.path(Dir.Exports, "PlotTopoAll_ls.RData")) # data for matrix visualisation

## PLOTTING ----------------------------------------------------------------
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
ggsave(filename = file.path(Dir.Exports, "Arias-Arone_2016_MAT.png"), width = 12/2, height = 12/2, units = "cm", scale = 3, dpi = 1e3)

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

png(filename = file.path(Dir.Exports, "Arias-Arone_2016_Initial.png"), width = 16, height = 16, units = "cm", res = 1700)
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
  plot(g, layout = layout_as_bipartite(g, hgap = 5, vgap = 0.5)[,c(2,1)], asp = 1.3, 
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
Post00 <- FUN_SimComp(PlantAnim = NULL, RunName = "CONCEPT", CutOffs = CutOffs, 
                      PotPartners = RewClass_ls, Traits = meta_df, 
                      IS = 0, Rewiring = 0)$`Arias-Arone 2016`$Climate$Prediction
Post00_igraph <- igraph::graph_from_adjacency_matrix(Post00, weighted = TRUE, mode = "undirected")

Post00.5 <- FUN_SimComp(PlantAnim = NULL, RunName = "CONCEPT", CutOffs = CutOffs, 
                      PotPartners = RewClass_ls, Traits = meta_df, 
                      IS = 0, Rewiring = 0.5)$`Arias-Arone 2016`$Climate$Prediction
Post00.5_igraph <- igraph::graph_from_adjacency_matrix(Post00.5, weighted = TRUE, mode = "undirected")

Post01 <- FUN_SimComp(PlantAnim = NULL, RunName = "CONCEPT", CutOffs = CutOffs, 
                      PotPartners = RewClass_ls, Traits = meta_df, 
                      IS = 0, Rewiring = 1)$`Arias-Arone 2016`$Climate$Prediction
Post01_igraph <- igraph::graph_from_adjacency_matrix(Post01, weighted = TRUE, mode = "undirected")

png(filename = file.path(Dir.Exports, "Arias-Arone_2016_0-0.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post00_igraph)
dev.off()
png(filename = file.path(Dir.Exports, "Arias-Arone_2016_0-0.5.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post00.5_igraph)
dev.off()
png(filename = file.path(Dir.Exports, "Arias-Arone_2016_0-1.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post01_igraph)
dev.off()


## y axis, 0 - interdependency (all the same because everything can and does rewire)
Post0.50 <- FUN_SimComp(PlantAnim = NULL, RunName = "CONCEPT", CutOffs = CutOffs, 
                        PotPartners = RewClass_ls, Traits = meta_df, 
                        IS = 0.5, Rewiring = 0)$`Arias-Arone 2016`$Climate$Prediction
Post0.50_igraph <- igraph::graph_from_adjacency_matrix(Post0.50, weighted = TRUE, mode = "undirected")

Post10 <- FUN_SimComp(PlantAnim = NULL, RunName = "CONCEPT", CutOffs = CutOffs, 
                      PotPartners = RewClass_ls, Traits = meta_df, 
                      IS = 1, Rewiring = 0)$`Arias-Arone 2016`$Climate$Prediction
Post10_igraph <- igraph::graph_from_adjacency_matrix(Post10, weighted = TRUE, mode = "undirected")


png(filename = file.path(Dir.Exports, "Arias-Arone_2016_0.5-0.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post0.50_igraph)
dev.off()
png(filename = file.path(Dir.Exports, "Arias-Arone_2016_1-0.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post10_igraph)
dev.off()


## combinations
Post0.50.5 <- FUN_SimComp(PlantAnim = NULL, RunName = "CONCEPT", CutOffs = CutOffs, 
                        PotPartners = RewClass_ls, Traits = meta_df, 
                        IS = 0.5, Rewiring = 0.5)$`Arias-Arone 2016`$Climate$Prediction
Post0.50.5_igraph <- igraph::graph_from_adjacency_matrix(Post0.50.5, weighted = TRUE, mode = "undirected")

Post10.5 <- FUN_SimComp(PlantAnim = NULL, RunName = "CONCEPT", CutOffs = CutOffs, 
                        PotPartners = RewClass_ls, Traits = meta_df, 
                        IS = 1, Rewiring = 0.5)$`Arias-Arone 2016`$Climate$Prediction
Post10.5_igraph <- igraph::graph_from_adjacency_matrix(Post10.5, weighted = TRUE, mode = "undirected")

Post0.51 <- FUN_SimComp(PlantAnim = NULL, RunName = "CONCEPT", CutOffs = CutOffs, 
                      PotPartners = RewClass_ls, Traits = meta_df,  
                      IS = 0.5, Rewiring = 1)$`Arias-Arone 2016`$Climate$Prediction
Post0.51_igraph <- igraph::graph_from_adjacency_matrix(Post0.51, weighted = TRUE, mode = "undirected")

Post11 <- FUN_SimComp(PlantAnim = NULL, RunName = "CONCEPT", CutOffs = CutOffs, 
                      PotPartners = RewClass_ls, Traits = meta_df,  
                      IS = 1, Rewiring = 1)$`Arias-Arone 2016`$Climate$Prediction
Post11_igraph <- igraph::graph_from_adjacency_matrix(Post11, weighted = TRUE, mode = "undirected")


png(filename = file.path(Dir.Exports, "Arias-Arone_2016_0.5-0.5.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post0.50.5_igraph)
dev.off()
png(filename = file.path(Dir.Exports, "Arias-Arone_2016_1-0.5.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post10.5_igraph)
dev.off()
png(filename = file.path(Dir.Exports, "Arias-Arone_2016_0.5-1.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post0.51_igraph)
dev.off()
png(filename = file.path(Dir.Exports, "Arias-Arone_2016_1-1.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post11_igraph)
dev.off()

































