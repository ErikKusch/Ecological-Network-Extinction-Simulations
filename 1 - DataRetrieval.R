#' ####################################################################### #
#' PROJECT: [Biodiversity Simplification & Ecological Network Topology] 
#' CONTENTS: 
#'  - Data Download
#'  - Data Manipulation
#'  DEPENDENCIES:
#'  - 0 - Preamble.R
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE =================================================================
rm(list=ls())
set.seed(42)

## Sourcing ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
source("0 - Preamble.R")
source("0 - Fricke_Functions.R")

# ECOLOGICAL DATA ==========================================================
# cite this data as: https://doi.org/10.5281/zenodo.5565122
if(!file.exists(file.path(Dir.D.Fricke, "Fricke2021.zip"))){
  download.file(url = "https://zenodo.org/record/5565123/files/evancf/global-dispersal-change-v1.0.0.zip?download=1", mode = "wb", 
                destfile = file.path(Dir.D.Fricke, "Fricke2021.zip"))
  unzip(zipfile = file.path(Dir.D.Fricke, "Fricke2021.zip"), exdir = Dir.D.Fricke)
  file.copy(from = list.files(path = file.path(Dir.D.Fricke, "evancf-global-dispersal-change-52409b4/data-and-model-outputs/"), full.names = TRUE), to = Dir.D.Fricke, recursive = TRUE)
  unlink(file.path(Dir.D.Fricke, "evancf-global-dispersal-change-52409b4"), recursive = TRUE)
}
Fricke_fs <- list.files(Dir.D.Fricke, full.names = TRUE)
lapply(Fricke_fs[-grep(x = Fricke_fs, pattern = ".zip")],
       load,
       .GlobalEnv)

## Networks ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(!file.exists(file.path(Dir.Data, "Networks.RData"))){
  nets <- unique(metanet$net.id)[unique(metanet$net.id) %in% unique(int.set$net.id)]
  List_ls <- as.list(rep(NA, length(nets)))
  names(List_ls) <- nets
  
  pb <- txtProgressBar(min = 0, max = length(nets), style = 3)
  
  for(i in 1:length(nets)){
    focal.net.id <- nets[i]
    focal.set <- subset(int.set, net.id == focal.net.id)
    if(nrow(focal.set) == 1){
      next()
      setTxtProgressBar(pb, i)}
    
    focal.net <- net.spread(split.by = "net.id", 
                            split.vals = focal.net.id,
                            tax.type = "phylo.id", 
                            data.type = "quant",
                            long.df = focal.set)
    if(length(focal.net) == 0){
      next()
      setTxtProgressBar(pb, i)
    }
    focal.net <- order.net(focal.net[[1]], focal.net[[1]])
    List_ls[[i]] <- focal.net 
    setTxtProgressBar(pb, i)
  }
  List_ls <- List_ls[which(!is.na(List_ls))]
  List_ls <- List_ls[unlist(lapply(List_ls, function(x){sum(x>1)/sum(x>0)})) > 0.5]
  List_ls <- List_ls[unlist(lapply(List_ls, function(x){sum(dim(x) > 10) == 2}))]
  networks_df <- metanet[metanet$net.id %in% names(List_ls), ]
  save(List_ls, networks_df, file = file.path(Dir.Data, "Networks.RData"))
}else{
  load(file.path(Dir.Data, "Networks.RData")) 
}

# library(leaflet)
# 
# 
# Map <- leaflet(networks_df) 
# Map <- addProviderTiles(Map, providers$Esri.WorldTopoMap) 
# Map <- addMarkers(Map, ~longitude, ~latitude, label = ~net.id,
#                   labelOptions = labelOptions(textsize = "12px")
# )
# Map

## Occurrences +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
plants_spec <- unlist(lapply(List_ls, FUN = function(x){rownames(x)}))
length(unique(animals_spec))
animals_spec <- unlist(lapply(List_ls, FUN = function(x){colnames(x)}))
length(unique(plants_spec))

## Traits ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
colnames(int.set)[5:29]
traits_df <- int.set[int.set$animal.phylo.id %in% animals_spec & int.set$plant.phylo.id %in% plants_spec, ]
plant_means <- aggregate(traits_df[18:25], by=list(Species=traits_df$plant.phylo.id), FUN=mean)
animals_means <- aggregate(traits_df[18:25], by=list(Species=traits_df$plant.phylo.id), FUN=mean)

# ABIOTIC DATA =============================================================
download_ERA(
  Variable = "2m_temperature",
  DataSet = "era5-land",
  DateStart = "1982-01-01",
  DateStop = "2020-12-31",
  TResolution = "month",
  TStep = 1,
  API_User = API_User,
  API_Key = API_Key,
  Dir = Dir.Data,
  FileName = "AT_Climatology.nc",
  SingularDL = TRUE
)

download_ERA(
  Variable = "2m_temperature",
  DataSet = "era5-land",
  DateStart = "1982-01-01",
  DateStop = "2020-12-31",
  TResolution = "month",
  TStep = 1,
  API_User = API_User,
  API_Key = API_Key,
  Dir = Dir.Data,
  FileName = "AT_Climatology.nc"
)


