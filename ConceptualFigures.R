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

ggsave(filename = file.path(Dir.Exports, "CONCEPT_NetworksCountries.png"), width = 16/2, height = 9/2, units = "cm", scale = 3, dpi = 1e3)

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

# ENVIRONMENTAL LAYERS =====================================================
world <- ne_countries(scale = "medium", returnclass = "sp")

Plot.Concept.Clim <- function(train_ERA, FName){
  Raw_df <- as.data.frame(train_ERA, xy = TRUE) # turn raster into dataframe
  world_plot <- crop(world, extent(train_ERA))
  colnames(Raw_df)[c(-1, -2)] <- c("Temperature", "Soil Moisture")
  Raw_df <- gather(data = Raw_df, key = Values, value = "value", colnames(Raw_df)[c(-1, -2)]) #  make ggplot-ready
  COL <- list(inferno(100), cividis(100))
  for(i in 1:2){
    Raw_plot <- ggplot() + # create plot
      geom_raster(data = Raw_df[Raw_df$Values == unique(Raw_df$Values)[i], ], 
                  aes(x = x, y = y, fill = value)) + # plot the covariate data
      theme_bw() +
      geom_polygon(data = world_plot, aes(x = long, y = lat, group = group), colour = "black", fill = "NA") + 
      geom_point(data = map_df, aes(x = longitude, y = latitude), pch = 17, size = 2, col = "darkred") +
      labs(x = "Longitude", y = "Latitude") + # make plot more readable
      scale_fill_gradientn(name = "", colours = COL[[i]], na.value = "transparent") + # add colour and legend
      labs(title = unique(Raw_df$Values)[i]) + 
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + # reduce margins (for fusing of plots)
      theme(legend.key.size = unit(2, 'cm')) + 
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(Raw_plot, width = 16/2, height = 7.8/2, units = "cm", scale = 3, dpi = 1e3, 
           filename = file.path(Dir.Exports, paste0("CONCEPT_", FName, "_", unique(Raw_df$Values)[i], ".png")))
  }
}

## DATA LOADING ------------------------------------------------------------
## present learning
Enviro_ras <- stack(file.path(Dir.Data, "Enviro_Pres.nc"))
nets_df <- metanet[metanet$net.id %in% names(AnalysisData_ls), ]
colnames(nets_df)[3:4] <- c("Lat", "Lon")
nets_shp <- KrigR:::buffer_Points(nets_df, Buffer = 5, ID = "net.id")
train_ERA <- Enviro_ras
train_ERA <- crop(train_ERA, extent(nets_shp))
train_ERA <- mask(train_ERA, nets_shp)

## projections raw
### SSP 
train_SSP <- lapply(c(file.path(Dir.D.Projections, "ssp245_tas_2081-2100.nc"), 
                      file.path(Dir.D.Projections, "ssp245_mrsos_2081-2100.nc")
), stack)
train_SSP <- lapply(train_SSP, FUN = function(x){
  x <- crop(x,extent(train_ERA))
  x <- mask(x, nets_shp)
  mean(x)
})
train_SSP[[1]] <- resample(x = train_SSP[[1]], y = train_SSP[[2]])
train_SSP <- stack(train_SSP)
names(train_SSP) <- c("Temperature", "Moisture")

### HISTORICAL
train_HIST <- lapply(c(file.path(Dir.D.Projections, "historical_tas_1981-2000.nc"), 
                       file.path(Dir.D.Projections, "historical_mrsos_1981-2000.nc")
), stack)
train_HIST <- lapply(train_HIST, FUN = function(x){
  x <- crop(x,extent(train_ERA))
  x <- mask(x, nets_shp)
  mean(x)
})
train_HIST[[1]] <- resample(x = train_HIST[[1]], y = train_HIST[[2]])
train_HIST <- stack(train_HIST)
names(train_HIST) <- c("Temperature", "Moisture")

## projections kriged
load(file.path(Dir.Data, "Projections.RData"))
names(Projections_stack[[1]]) <- c("Tair.Historical", "Tair.SSP245", "Tair.Diff")
names(Projections_stack[[2]]) <- c("Qsoil.Historical", "Qsoil.SSP245", "Qsoil.Diff")

## PLOTTING ----------------------------------------------------------------
Plot.Concept.Clim(train_ERA = train_ERA, FName = "Present")
Plot.Concept.Clim(train_ERA = train_SSP, FName = "SSP_Future")
Plot.Concept.Clim(train_ERA = train_HIST, FName = "SSP_Historic")
Plot.Concept.Clim(train_ERA = stack(resample(Projections_stack$Temp$Tair.Diff, Projections_stack$Water$Qsoil.Diff),
                                    Projections_stack$Water$Qsoil.Diff), 
                  FName = "SSP_KrigDiff")

# SPECIES DISTRIBUTIONS ====================================================
world <- ne_countries(scale = "medium", returnclass = "sf")
## DATA LOADING ------------------------------------------------------------
occ_fs <- list.files(Dir.D.Occurrences, full.names = TRUE, pattern = ".rds")
occ_ls <- as.list(pbsapply(occ_fs, readRDS))
occ_spec <- list.files(Dir.D.Occurrences, pattern = ".rds")
names(occ_ls) <- tools::file_path_sans_ext(occ_spec)

## PLOTTING ----------------------------------------------------------------
for(x in names(occ_ls)[1:5]){
  y <- occ_ls[[x]]
  ggplot(data = world) +
    theme_minimal_grid() + 
    geom_sf(color = "black", size = 0.2, fill = "white") + 
    geom_point(data = as.data.frame(y), aes(x = decimalLongitude, y = decimalLatitude), pch = 17, size = 2, col = "darkred") +
    xlab("Longitude") + ylab("Latitude") +
    labs(fill = "", title = x) + 
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename = file.path(Dir.Exports, paste0("CONCEPT_", x, ".jpg")), width = 16/2, height = 9/2, units = "cm", scale = 3, dpi = 1e3)
}


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
ggsave(filename = file.path(Dir.Exports, "CONCEPT_Arias-Arone_2016_MAT.png"), width = 12/2, height = 12/2, units = "cm", scale = 3, dpi = 1e3)

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

png(filename = file.path(Dir.Exports, "CONCEPT_Arias-Arone_2016_Initial.png"), width = 16, height = 16, units = "cm", res = 1700)
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

png(filename = file.path(Dir.Exports, "CONCEPT_Arias-Arone_2016_0-0.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post00_igraph)
dev.off()
png(filename = file.path(Dir.Exports, "CONCEPT_Arias-Arone_2016_0-0.5.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post00.5_igraph)
dev.off()
png(filename = file.path(Dir.Exports, "CONCEPT_Arias-Arone_2016_0-1.png"), width = 16, height = 16, units = "cm", res = 1700)
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


png(filename = file.path(Dir.Exports, "CONCEPT_Arias-Arone_2016_0.5-0.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post0.50_igraph)
dev.off()
png(filename = file.path(Dir.Exports, "CONCEPT_Arias-Arone_2016_1-0.png"), width = 16, height = 16, units = "cm", res = 1700)
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


png(filename = file.path(Dir.Exports, "CONCEPT_Arias-Arone_2016_0.5-0.5.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post0.50.5_igraph)
dev.off()
png(filename = file.path(Dir.Exports, "CONCEPT_Arias-Arone_2016_1-0.5.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post10.5_igraph)
dev.off()
png(filename = file.path(Dir.Exports, "CONCEPT_Arias-Arone_2016_0.5-1.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post0.51_igraph)
dev.off()
png(filename = file.path(Dir.Exports, "CONCEPT_Arias-Arone_2016_1-1.png"), width = 16, height = 16, units = "cm", res = 1700)
plot.Concept(g2 = Post11_igraph)
dev.off()

































