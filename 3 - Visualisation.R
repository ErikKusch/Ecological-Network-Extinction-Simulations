#' ####################################################################### #
#' PROJECT: [Biodiversity Simplification & Ecological Network Topology] 
#' CONTENTS: 
#'  - Network Visualisations
#'  DEPENDENCIES:
#'  - "1 - DataRetrieval.R"
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE =================================================================
rm(list=ls())
set.seed(42)

## Sourcing ----------------------------------------------------------------
source("0 - Preamble.R")
source("0 - Fricke_Functions.R")

message("########### STARTING VISUALISATION GENERATION ###########")

# NETWORK DATA ==========================================================
message("### LOADING NETWORK DATA ###")
load(file.path(Dir.Data, "Networks.RData")) 
Fricke_fs <- list.files(Dir.D.Fricke, full.names = TRUE)
hush_ls <- pblapply(Fricke_fs[-grep(x = Fricke_fs, pattern = ".zip")],
                    load,
                    .GlobalEnv)


### this should probably go at the end of "1 - Data Retrieval"
PreExt_df <- pblapply(names(List_ls), FUN = function(netID){
  plot_df <- List_ls[[netID]]
  data.frame(Nestedness = networklevel(web = plot_df, index = "weighted nestedness"),
             Modularity = NOS(web = plot_df)$mod,
             RobustnessAnimal = robustness(second.extinct(web = plot_df, participant = "higher")),
             RobustnessPlant = robustness(second.extinct(web = plot_df, participant = "lower"))
  )
})
leaflet_df <- cbind(networks_df, do.call(rbind, PreExt_df))
leaflet_df[, c("Nestedness", "Modularity", "RobustnessAnimal", "RobustnessPlant")] <- round(leaflet_df[, c("Nestedness", "Modularity", "RobustnessAnimal", "RobustnessPlant")], 2)

# SPECIES NAMES & ABBREVIATIONS =========================================
## Animals Abbreviations for entire data set ----
animals_spec <- unique(unlist(lapply(List_ls, FUN = function(x){colnames(x)})))
AnimalsLegend <- data.frame(Abbr = paste0(substr(unlist(lapply(strsplit(animals_spec, split = " "), "[[", 1)), start = 1, stop = 1),
                                         ".",
                                         substr(unlist(lapply(strsplit(animals_spec, split = " "), "[[", 2)), start = 1, stop = 2)),
                           Full = animals_spec
)
### manipulating duplicate abbreviations
AnimalsLegend <- AnimalsLegend %>% group_by(Abbr) %>% 
  mutate(
    Abbr = paste(Abbr, row_number()-1, sep = "_")
  )
AnimalsLegend$Abbr <- gsub(AnimalsLegend$Abbr, pattern = "_0", replacement = "")
AnimalsLegend <- Sort.DF(as.data.frame(AnimalsLegend), "Abbr")

## Plants Abbreviations for entire data set ----
plants_spec <- unique(unlist(lapply(List_ls, FUN = function(x){rownames(x)})))
PlantsLegend <- data.frame(Abbr = paste0(substr(unlist(lapply(strsplit(plants_spec, split = " "), "[[", 1)), start = 1, stop = 1),
                                         ".",
                                         substr(unlist(lapply(strsplit(plants_spec, split = " "), "[[", 2)), start = 1, stop = 2)),
                           Full = plants_spec
)
### manipulating duplicate abbreviations
PlantsLegend <- PlantsLegend %>% group_by(Abbr) %>% 
  mutate(
    Abbr = paste(Abbr, row_number()-1, sep = "_")
  )
PlantsLegend$Abbr <- gsub(PlantsLegend$Abbr, pattern = "_0", replacement = "")
PlantsLegend <- Sort.DF(as.data.frame(PlantsLegend), "Abbr")

## back-translating abbreviations into network list ----
Abbreviated_ls <- pblapply(names(List_ls), FUN = function(x){
  colnames(List_ls[[x]]) <- AnimalsLegend$Abbr[match(colnames(List_ls[[x]]), AnimalsLegend$Full)]
  rownames(List_ls[[x]]) <- PlantsLegend$Abbr[match(rownames(List_ls[[x]]), PlantsLegend$Full)]
  List_ls[[x]]
})
names(Abbreviated_ls) <- names(List_ls)

# GENERATION OF PLOT OBJECTS ============================================
Dir.E.LeafletImages <- file.path(Dir.Exports, "LeafletImages")
if(!dir.exists(Dir.E.LeafletImages)){dir.create(Dir.E.LeafletImages)}
Plots_ls <- pblapply(names(Abbreviated_ls), FUN = function(netID){
  # print(netID) # just here for debugging
  # data extraction
  plot_df <- Abbreviated_ls[[netID]]
  ## plot object generation
  plotweb(web = plot_df, y.lim = c(0.3, 1.7)) # network plot
  p <- recordPlot()
  p <- plot_grid(p, NULL, ncol = 1, rel_heights = c(1,0)) # making into ggplot-type plot object for leaflet later
  ggsave(plot = p, filename = file.path(Dir.E.LeafletImages, paste0(netID, ".png")), width = 2*max(nrow(plot_df), ncol(plot_df)), height = 18, units = "cm", limitsize = FALSE)
  p
})
names(Plots_ls) <- names(Abbreviated_ls)

# LEAFLET MAP GENERATION ================================================
# ## currently not used, but is useful for purely local visualisation
# netpaths <- file.path(Dir.E.LeafletImages, paste0(leaflet_df$net.id, ".png"))
# label_ls <- lapply(paste( "<b> Study ID: </b>" , leaflet_df$study.id, "<br>",
#                           "Network ID:" , leaflet_df$net.id, "<br>"
#                           , "Weighted Nestedness", leaflet_df$Nestedness, "<br>"
#                           , "Modularity", leaflet_df$Modularity, "<br>", 
#                           "Robustness to animal extinction", leaflet_df$RobustnessAnimal, "<br>"
#                           , "Robustness to plant extinction", leaflet_df$RobustnessPlant),
#                    htmltools::HTML)
# addPopupImages(netpaths, group = "net.id", width = 600, maxWidth = 600, height = 300)

leaflet() %>%
  # Base groups
  addTiles() %>%
  # Overlay groups
  addCircleMarkers(data = leaflet_df, ~longitude, ~latitude,
                   group = "net.id", 
                   label = ~study.id
  ) %>%
addPopupGraphs(Plots_ls, group = "net.id", width = 800, maxWidth = 800, height = 400)


save.image(file = "Shiny.RData")
