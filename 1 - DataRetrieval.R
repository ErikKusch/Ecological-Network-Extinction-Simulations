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

## Sourcing ----------------------------------------------------------------
source("0 - Preamble.R")
source("0 - Fricke_Functions.R")
source("0 - Data_Functions.R")

# ECOLOGICAL DATA ==========================================================
# cite this data as: https://doi.org/10.5281/zenodo.5565122
if(!file.exists(file.path(Dir.D.Fricke, "Fricke2021.zip"))){
  message("Downloading network data")
  download.file(url = "https://zenodo.org/record/5565123/files/evancf/global-dispersal-change-v1.0.0.zip?download=1", mode = "wb", 
                destfile = file.path(Dir.D.Fricke, "Fricke2021.zip"))
  unzip(zipfile = file.path(Dir.D.Fricke, "Fricke2021.zip"), exdir = Dir.D.Fricke)
  file.copy(from = list.files(path = file.path(Dir.D.Fricke, "evancf-global-dispersal-change-52409b4/data-and-model-outputs/"), full.names = TRUE), to = Dir.D.Fricke, recursive = TRUE)
  unlink(file.path(Dir.D.Fricke, "evancf-global-dispersal-change-52409b4"), recursive = TRUE)
}
message("Loading network data")
Fricke_fs <- list.files(Dir.D.Fricke, full.names = TRUE)
lapply(Fricke_fs[-grep(x = Fricke_fs, pattern = ".zip")],
       load,
       .GlobalEnv)

## Networks ----------------------------------------------------------------
if(!file.exists(file.path(Dir.Data, "Networks.RData"))){
  message("Subsetting network data")
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
  List_ls <- List_ls[unlist(lapply(List_ls, function(x){sd(x[x>0])/mean(x[x>0])})) > 0.5]
  List_ls <- List_ls[unlist(lapply(List_ls, function(x){sum(dim(x) > 7) == 2}))]
  networks_df <- metanet[metanet$net.id %in% names(List_ls), ]
  save(List_ls, networks_df, file = file.path(Dir.Data, "Networks.RData"))
}else{
  message("Network data already subsetted")
  load(file.path(Dir.Data, "Networks.RData")) 
}

message("currently only running for two networks for testing purposes")
List_ls <- List_ls[1:2]

# library(leaflet)
# 
# 
# Map <- leaflet(networks_df)
# Map <- addProviderTiles(Map, providers$Esri.WorldTopoMap)
# Map <- addMarkers(Map, ~longitude, ~latitude, label = ~net.id,
#                   labelOptions = labelOptions(textsize = "12px")
# )
# Map

## Occurrences -------------------------------------------------------------
if(file.exists(file.path(Dir.Data, "Occurrences.RData"))){
  message("Occurrences already obtained")
  load(file.path(Dir.Data, "Occurrences.RData"))
}else{
  message("Obtaining occurrence records")
  ### Plants ####
  message("PLANTS")
  plants_spec <- unlist(lapply(List_ls, FUN = function(x){rownames(x)}))
  plants_occ <- occ_data(scientificName = unique(plants_spec))
  plants_occ <- lapply(plants_occ, FUN = function(x){nrow(x[[2]])})
  Failed_plants <- names(plants_occ)[which(unlist(plants_occ) == 0)]
  if(length(Failed_plants) != 0){stop("Not all plant species are found on gbif")}
  plants_gbif <- Gbif_Species(species = plants_spec, year_vec = 2000:2020)
  
  ### Animals ####
  message("ANIMALS")
  animals_spec <- unlist(lapply(List_ls, FUN = function(x){colnames(x)}))
  animals_occ <- occ_data(scientificName = unique(animals_spec))
  animals_occ <- lapply(animals_occ, FUN = function(x){nrow(x[[2]])})
  Failed_animals <- names(animals_occ)[which(unlist(animals_occ) == 0)]
  if(length(Failed_animals) != 0){stop("Not all animal species are found on gbif")} 
  animals_gbif <- Gbif_Species(species = animals_spec, year_vec = 2000:2020)
  
  ### Saving ####
  save(animals_spec, animals_gbif, plants_spec, plants_gbif, file = file.path(Dir.Data, "Occurrences.RData"))
}

## Traits ------------------------------------------------------------------
if(file.exists(file.path(Dir.Data, "Traits.RData"))){
  message("Traits already extracted")
  load(file.path(Dir.Data, "Traits.RData"))
}else{
  message("Extracting traits")
  colnames(int.set)[5:29]
  traits_df <- int.set[int.set$animal.phylo.id %in% animals_spec & int.set$plant.phylo.id %in% plants_spec, ]
  plant_means <- aggregate(traits_df[18:25], by=list(Species=traits_df$plant.phylo.id), FUN=mean)
  animals_means <- aggregate(traits_df[18:25], by=list(Species=traits_df$plant.phylo.id), FUN=mean)
  save(traits_df, animals_means, plant_means, file = file.path(Dir.Data, "Traits.RData"))
}

# ABIOTIC DATA =============================================================
## BUFFER CREATION ---------------------------------------------------------
nets_df <- metanet[metanet$study.id %in% names(List_ls), ]
colnames(nets_df)[3:4] <- c("Lat", "Lon")
nets_shp <- KrigR:::buffer_Points(nets_df, Buffer = 10, ID = "study.id")

plot(Enviro_ras[[1]])
plot(nets_shp, add = TRUE)

## BASELINE ----------------------------------------------------------------
if(file.exists(file.path(Dir.Data, "Enviro_Pres.nc"))){
  message("Environmental data already obtained")
  Enviro_ras <- stack(file.path(Dir.Data, "Enviro_Pres.nc"))
}else{
  message("Obtaining environmental data")
  if(!file.exists(file.path(Dir.Data, "AT_Climatology.nc"))){
    download_ERA(
      Variable = "2m_temperature",
      DataSet = "era5-land",
      DateStart = "1982-01-01",
      DateStop = "1999-12-31",
      TResolution = "month",
      TStep = 1,
      API_User = API_User,
      API_Key = API_Key,
      Dir = Dir.Data,
      FileName = "AT_Climatology.nc",
      SingularDL = TRUE
    )
  }
  AT_ras <- stack(file.path(Dir.Data, "AT_Climatology.nc"))
  if(!file.exists(file.path(Dir.Data, "PP_Climatology.nc"))){
    download_ERA(
      Variable = "total_precipitation",
      PrecipFix = TRUE,
      DataSet = "era5-land",
      DateStart = "1982-01-01",
      DateStop = "1999-12-31",
      TResolution = "month",
      TStep = 1,
      API_User = API_User,
      API_Key = API_Key,
      Dir = Dir.Data,
      FileName = "PP_Climatology.nc",
      SingularDL = TRUE
    )
  }
  PP_ras <- stack(file.path(Dir.Data, "PP_Climatology.nc"))
  Enviro_ras <- stack(mean(AT_ras), mean(PP_ras))
  writeRaster(Enviro_ras, format = "CDF", filename = file.path(Dir.Data, "Enviro_Pres.nc"))
}

## PROJECTIONS -------------------------------------------------------------
message("For now, we are just using air temperature data")
train_ERA <- Enviro_ras[[1]]
train_ERA <- crop(train_ERA, extent(nets_shp))
train_ERA <- mask(train_ERA, nets_shp)

### SSP ----
train_SSP <- stack(file.path(Dir.D.Projections, "ssp585_tas_2041-2060.nc"))
train_SSP <- train_SSP[[1]]
train_SSP <- crop(train_SSP,extent(train_ERA))
train_SSP <- mask(train_SSP, nets_shp)

### HISTORICAL ----
train_HIST <- stack(file.path(Dir.D.Projections, "historical_tas_1981-2000.nc"))
train_HIST <- train_HIST[[1]]
train_HIST <- crop(train_HIST,extent(train_ERA))
train_HIST <- mask(train_HIST, nets_shp)

## KRIGING -----------------------------------------------------------------
### HISTORICAL DATA ----
GMTED <- download_DEM(
  Train_ras = train_ERA,
  Target_res = 0.008334,
  Keep_Temporary = TRUE,
  Dir = Dir.Data
)
Cov_coarse <- GMTED[[1]]
Cov_fine <- GMTED[[2]]

Output_ERA <- krigR(
  Data = train_ERA,
  Covariates_coarse = Cov_coarse, 
  Covariates_fine = Cov_fine,   
  KrigingEquation = "ERA ~ DEM",  
  Cores = 1, 
  Dir = Dir.Data,  
  FileName = "hist_nmax120", 
  Keep_Temporary = FALSE,
  nmax = 120
)

### SSP ----
GMTED <- download_DEM(
  Train_ras = train_SSP,
  Target_res = 0.008334,
  Keep_Temporary = TRUE,
  Dir = Dir.Data
)
Cov_coarse <- GMTED[[1]]
Cov_fine <- GMTED[[2]]

Output_SSP <- krigR(
  Data = train_SSP,
  Covariates_coarse = Cov_coarse, 
  Covariates_fine = Cov_fine,   
  KrigingEquation = "ERA ~ DEM",  
  Cores = 1, 
  Dir = Dir.Data,  
  FileName = "SSP585_2041-2060_nmax120", 
  Keep_Temporary = FALSE,
  nmax = 120
)

### HISTORICAL ----
GMTED <- download_DEM(
  Train_ras = train_HIST,
  Target_res = 0.008334,
  Keep_Temporary = TRUE,
  Dir = Dir.Data
)
Cov_coarse <- GMTED[[1]]
Cov_fine <- GMTED[[2]]

Output_SSP <- krigR(
  Data = train_HIST,
  Covariates_coarse = Cov_coarse, 
  Covariates_fine = Cov_fine,   
  KrigingEquation = "ERA ~ DEM",  
  Cores = 1, 
  Dir = Dir.Data,  
  FileName = "CMIP-HIST_nmax120", 
  Keep_Temporary = FALSE,
  nmax = 120
)


# EXTINCTION PROXIES =======================================================
stop("create species by site by extinction proxy df here")

## Network Centrality ------------------------------------------------------

## Safety Margins ----------------------------------------------------------

## IUCN Criteria -----------------------------------------------------------











# 
# 
# ### extinction risk trials
# animals_gbif <- Gbif_Species(species = animals_spec[3], year_vec = 2000:2020)
# test_df <- animals_gbif
# x_range <- range(test_df$decimalLongitude)
# y_range <- range(test_df$decimalLatitude)
# coordinates(test_df) <- ~decimalLongitude + decimalLatitude
# test_ras <- mean(crop(AT_ras, extent(x_range, y_range)))
# 
# 
# library(rworldmap)
# 
# newmap <- getMap(resolution = "low")
# # plot map
# plot(newmap, asp = 1, border = "darkgray", col = "black", bg = "gray95",
#      xlim = x_range, ylim = y_range
#      )
# plot(test_ras, add = TRUE)
# points(test_df, col = "red", cex = .5, pch = 20)
# message("pretty bad range for some species also means pretty bad sampling from 9x9km ERA5-Land. IDEA: could locally krig from point-buffers that extract for just the sample locations.")
# 
# test_At <- raster::extract(x = test_ras, y = test_df)
# test_At <- na.omit(test_At)
# 
# 
# set.seed(42)
# 
# at_c <- c()
# Boot <- 1e3
# pb <- txtProgressBar(min = 0, max = Boot, style = 3)
# for(i in 1:Boot){
#   at_c <- c(at_c, sample(x = test_At, size = round(length(test_At)*0.7), replace = TRUE))
#   setTxtProgressBar(pb, i)
#   }
# 
# # plot(density(at_c))
# # abline(v = median(at_c), col = "blue")
# # abline(v = median(at_c)+2*sd(at_c), col = "red")
# # abline(v = median(at_c)-2*sd(at_c), col = "red")
# 
# library(MASS)
# fit <- fitdistr(at_c, "normal")
# para <- fit$estimate
# plot(density(at_c))
# curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)
# abline(v = para[1], col = "blue")
# abline(v = para[1]+2*para[2], col = "blue", lty = 2)
# abline(v = para[1]-2*para[2], col = "blue", lty = 2)
# 






















