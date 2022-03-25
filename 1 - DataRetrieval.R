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

message("########### STARTING DATA RETRIEVAL ###########")

# NETWORK DATA ==========================================================
message("### NETWORK DATA ###")
# cite this data as: https://doi.org/10.5281/zenodo.5565122
if(!file.exists(file.path(Dir.D.Fricke, "Fricke2021.zip"))){
  print("Downloading network data")
  download.file(url = "https://zenodo.org/record/5565123/files/evancf/global-dispersal-change-v1.0.0.zip?download=1", mode = "wb", 
                destfile = file.path(Dir.D.Fricke, "Fricke2021.zip"))
  unzip(zipfile = file.path(Dir.D.Fricke, "Fricke2021.zip"), exdir = Dir.D.Fricke)
  file.copy(from = list.files(path = file.path(Dir.D.Fricke, "evancf-global-dispersal-change-52409b4/data-and-model-outputs/"), full.names = TRUE), to = Dir.D.Fricke, recursive = TRUE)
  unlink(file.path(Dir.D.Fricke, "evancf-global-dispersal-change-52409b4"), recursive = TRUE)
  unlink(file.path(Dir.D.Fricke, "gbm.mods1.RData"), recursive = TRUE)
  unlink(file.path(Dir.D.Fricke, "gbm.mods2.RData"), recursive = TRUE)
  unlink(file.path(Dir.D.Fricke, "range.matrices.RData"), recursive = TRUE)
}
print("Loading network data")
Fricke_fs <- list.files(Dir.D.Fricke, full.names = TRUE)
hush_ls <- pblapply(Fricke_fs[-grep(x = Fricke_fs, pattern = ".zip")],
                    load,
                    .GlobalEnv)

## Networks ----------------------------------------------------------------
if(!file.exists(file.path(Dir.Data, "Networks.RData"))){
  print("Subsetting network data")
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
  print("Network data already subsetted")
  load(file.path(Dir.Data, "Networks.RData")) 
}

# message("currently only running for two networks for testing purposes")
# List_ls <- List_ls[1:2]

# library(leaflet)
# 
# 
# Map <- leaflet(networks_df)
# Map <- addProviderTiles(Map, providers$Esri.WorldTopoMap)
# Map <- addMarkers(Map, ~longitude, ~latitude, label = ~net.id,
#                   labelOptions = labelOptions(textsize = "12px")
# )
# Map

# SPATIAL DATA =============================================================
## BUFFER CREATION ---------------------------------------------------------
nets_df <- metanet[metanet$study.id %in% names(List_ls), ]
colnames(nets_df)[3:4] <- c("Lat", "Lon")
nets_shp <- KrigR:::buffer_Points(nets_df, Buffer = 10, ID = "study.id")

## LANDMASK ----------------------------------------------------------------
Countries_shp <- ne_countries(scale = 10, type = "countries")
Land_shp <- crop(Countries_shp, extent(nets_shp))

## CENTROIDS ---------------------------------------------------------------
if(file.exists(file.path(Dir.D.Occurrences, "Centroids.RData"))){
  load(file.path(Dir.D.Occurrences, "Centroids.RData"))
}else{
  ### Countries and States
  Shapes_shp <- rbind(
    ne_countries(scale = 10, type = "countries")[,"name"],
    ne_countries(scale = 10, type = "sovereignty")[,"name"],
    ne_states()[,"name"]
  )
  ### Continents
  continents_vec <- c("africa", "europe", "asia", "oceania", "north america", "south america", "antarctica")
  Conts_ls <- as.list(rep(NA, length(continents_vec)))
  names(Conts_ls) <- continents_vec
  for(i in continents_vec){
    Conts_ls[[i]] <- raster::aggregate(ne_countries(scale = 10, type = "countries", continent = i))
    Conts_ls[[i]]$name <- i
  }
  Conts_shp <- do.call(rbind, Conts_ls)
  ### Final product and centroids
  Shapes_shp <- rbind(Shapes_shp, Conts_shp)
  Shapes_ct <- gCentroid(Shapes_shp, byid = TRUE)
  save(Shapes_ct, file = file.path(Dir.D.Occurrences, "Centroids.RData"))
}

# ABIOTIC DATA =============================================================
## BASELINE ----------------------------------------------------------------
message("### HISTORICAL ERA5-LAND DATA ###")
if(file.exists(file.path(Dir.Data, "Enviro_Pres.nc"))){
  print("Environmental data already obtained")
  Enviro_ras <- stack(file.path(Dir.Data, "Enviro_Pres.nc"))
}else{
  print("Obtaining environmental data")
  if(!file.exists(file.path(Dir.D.Climatologies, "AT_Climatology.nc"))){
    download_ERA(
      Variable = "2m_temperature",
      DataSet = "era5-land",
      DateStart = "1982-01-01",
      DateStop = "1999-12-31",
      TResolution = "month",
      TStep = 1,
      API_User = API_User,
      API_Key = API_Key,
      Dir = Dir.D.Climatologies,
      FileName = "AT_Climatology.nc",
      SingularDL = TRUE
    )
  }
  AT_ras <- stack(file.path(Dir.D.Climatologies, "AT_Climatology.nc"))

  if(!file.exists(file.path(Dir.D.Climatologies, "QS_Climatology.nc"))){
    download_ERA(
      Variable = "volumetric_soil_water_layer_1",
      PrecipFix = FALSE,
      DataSet = "era5-land",
      DateStart = "1982-01-01",
      DateStop = "1999-12-31",
      TResolution = "month",
      TStep = 1,
      API_User = API_User,
      API_Key = API_Key,
      Dir = Dir.D.Climatologies,
      FileName = "QS_Climatology.nc",
      SingularDL = TRUE
    )
  }
  QS_ras <- stack(file.path(Dir.D.Climatologies, "QS_Climatology.nc"))
  Enviro_ras <- stack(mean(AT_ras), mean(QS_ras))
  writeRaster(Enviro_ras, format = "CDF", filename = file.path(Dir.Data, "Enviro_Pres.nc"))
}
plot(Enviro_ras[[1]])
plot(nets_shp, add = TRUE)


## PROJECTIONS -------------------------------------------------------------
message("### PROJECTION DATA ###")
print("Loading projection data")
if(!file.exists(file.path(Dir.Data, "Projections.nc"))){
  print("Loading raw projection data")
  message("For now, we are just using air temperature data")
  train_ERA <- Enviro_ras[[1]]
  train_ERA <- crop(train_ERA, extent(nets_shp))
  train_ERA <- mask(train_ERA, nets_shp)
  
  ### SSP ----
train_SSP <- stack(c(file.path(Dir.D.Projections, "ssp245_tas_2081-2000.nc")
                       # ,
                       # file.path(Dir.D.Projections, "ssp245_qs_2081-2000.nc")
                       )
                     )
  train_SSP <- train_SSP
  train_SSP <- crop(train_SSP,extent(train_ERA))
  train_SSP <- mask(train_SSP, nets_shp)
  
  ### HISTORICAL ----
  train_HIST <- stack(c(file.path(Dir.D.Projections, "historical_tas_1981-2000.nc")
                        # ,
                        # file.path(Dir.D.Projections, "historical_qs_1981-2000.nc")
                        )
                      )
  train_HIST <- train_HIST
  train_HIST <- crop(train_HIST,extent(train_ERA))
  train_HIST <- mask(train_HIST, nets_shp)
}

## KRIGING -----------------------------------------------------------------
message("### PROJECTION KRIGING ###")
if(!file.exists(file.path(Dir.Data, "Projections.nc"))){
  print("Kriging raw projection data")
  ### TEMPERATURE ----
  #### Covariates ----
  GMTED <- download_DEM(
    Train_ras = train_SSP,
    Target_res = train_ERA,
    Shape = Land_shp,
    Keep_Temporary = FALSE,
    Dir = Dir.D.Projections
  )
  Cov_coarse <- GMTED[[1]]
  Cov_coarse <- mask(Cov_coarse, nets_shp)
  Cov_fine <- GMTED[[2]]
  Cov_fine <- mask(Cov_fine, nets_shp)
  
  #### SSP ----
  if(!file.exists(file.path(Dir.D.Projections, "K_ssp245_tas_2081-2000_nmax120.nc"))){
    Output_SSP <- krigR(
      Data = train_SSP[[1]],
      Covariates_coarse = Cov_coarse, 
      Covariates_fine = Cov_fine,   
      KrigingEquation = "ERA ~ DEM",  
      Cores = numberOfCores, 
      Dir = Dir.D.Projections,  
      FileName = "K_ssp245_tas_2081-2000_nmax120", 
      Keep_Temporary = FALSE,
      nmax = 120
    )
  }

  TAS_ras <- mean(stack(file.path(Dir.D.Projections, "K_ssp245_tas_2081-2000_nmax120.nc")))
  
  #### HISTORICAL ----
  if(!file.exists(file.path(Dir.D.Projections, "K_CMIP-HIST_tas_nmax120.nc"))){
    Output_SSP <- krigR(
      Data = train_HIST[[1]],
      Covariates_coarse = Cov_coarse, 
      Covariates_fine = Cov_fine,   
      KrigingEquation = "ERA ~ DEM",  
      Cores = numberOfCores, 
      Dir = Dir.D.Projections,  
      FileName = "K_CMIP-HIST_tas_nmax120", 
      Keep_Temporary = FALSE,
      nmax = 120
    )
  }
  CMIP_tas_ras <- mean(stack(file.path(Dir.D.Projections, "K_CMIP-HIST_tas_nmax120.nc")))
  
  ### SOIL MOISTURE ----
  #### Covariates ----
  if(!file.exists(file.path(Dir.Data, "SoilCovs_ls.RData"))){
    SoilCovs_vec <- c("tkdry", "tksat", "csol", "k_s", "lambda", "psi", "theta_s") # need these names for addressing soil covariates, documentation of these can be found here http://globalchange.bnu.edu.cn/research/soil4.jsp
    # create lists to combine soil data into one
    SoilCovs_ls <- as.list(rep(NA, length(SoilCovs_vec)))
    names(SoilCovs_ls) <- c(SoilCovs_vec)
    ## Downloading, unpacking, and stacking
    for(Soil_Iter in SoilCovs_vec){
      if(file.exists(file.path(Dir.D.Projections, paste0(Soil_Iter, ".nc")))) {
        print(paste(Soil_Iter, "already downloaded and processed."))
        SoilCovs_ls[[which(names(SoilCovs_ls) == Soil_Iter)]] <- raster(file.path(Dir.D.Projections, paste0(Soil_Iter, ".nc")))
        next()
      } # if not downloaded and processed yet
      print(paste("Handling", Soil_Iter, "data."))
      Dir.Soil <- file.path(Dir.D.Projections, Soil_Iter)
      dir.create(Dir.Soil)
      download.file(paste0("http://globalchange.bnu.edu.cn/download/data/worldptf/", Soil_Iter,".zip"),
                    destfile = file.path(Dir.Soil, paste0(Soil_Iter, ".zip"))
      ) # download data
      unzip(file.path(Dir.Soil, paste0(Soil_Iter, ".zip")), exdir = Dir.Soil) # unzip data
      File <- list.files(Dir.Soil, pattern = ".nc")[1] # only keep first soil layer
      Soil_ras <- raster(file.path(Dir.Soil, File)) # load data
      SoilCovs_ls[[which(names(SoilCovs_ls) == Soil_Iter)]] <- Soil_ras # save to list
      writeRaster(x = Soil_ras, filename = file.path(Dir.D.Projections, Soil_Iter), format = "CDF")
      plot(Soil_ras, main = Soil_Iter, colNA = "black")
      unlink(Dir.Soil, recursive = TRUE)
    }
    SoilCovs_stack <- stack(SoilCovs_ls)
    # range_m <- KrigR::mask_Shape(SoilCovs_stack, nets_shp)
    SoilCovs <- crop(SoilCovs_stack, extent(nets_shp))
    SoilCovs <- mask(SoilCovs, nets_shp)
    Cov_coarse <- resample(x = SoilCovs, y = train_SSP)
    Cov_fine <- resample(x = SoilCovs, y = train_ERA)
    SoilCovs_ls <- list(Coarse = Cov_coarse,
                        Fine = Cov_fine)
    save(SoilCovs_ls, file = file.path(Dir.D.Projections, "SoilCovs_ls.RData"))
    unlink(paste0(Dir.D.Projections, "/", SoilCovs_vec, ".nc"))
  }else{
    load(file.path(Dir.D.Projections, "SoilCovs_ls.RData"))
    Cov_coarse <- SoilCovs_ls[[1]]
    Cov_fine <- SoilCovs_ls[[2]]
  }
    
  #### SSP ----
  if(!file.exists(file.path(Dir.D.Projections, "K_ssp245_qs_2081-2000_nmax120.nc"))){
    QS_SSP <- krigR(
      Data = train_SSP[[2]],
      Covariates_coarse = Cov_coarse, 
      Covariates_fine = Cov_fine,   
      KrigingEquation = "ERA ~ DEM",  
      Cores = numberOfCores, 
      Dir = Dir.D.Projections,  
      FileName = "K_ssp245_qs_2081-2000_nmax120", 
      Keep_Temporary = FALSE,
      nmax = 120
    )
  }
  QS_ras <- mean(stack(file.path(Dir.D.Projections, "K_ssp245_qs_2081-2000_nmax120.nc")))
  
  #### HISTORICAL ----
  if(!file.exists(file.path(Dir.D.Projections, "K_CMIP-HIST_qs_nmax120.nc"))){
    Output_SSP <- krigR(
      Data = train_HIST[[2]],
      Covariates_coarse = Cov_coarse, 
      Covariates_fine = Cov_fine,   
      KrigingEquation = "ERA ~ DEM",  
      Cores = numberOfCores, 
      Dir = Dir.D.Projections,  
      FileName = "K_CMIP-HIST_qs_nmax120", 
      Keep_Temporary = FALSE,
      nmax = 120
    )
  }

  CMIP_qs_ras <- mean(stack(file.path(Dir.D.Projections, "K_CMIP-HIST_qs_nmax120.nc")))
  
  ### DIFFERENCE AND FUSING ----
  Projections_stack <- stack(CMIP_tas_ras,
                             TAS_ras,
                             TAS_ras - CMIP_tas_ras,
                             CMIP_qs_ras,
                             QS_ras,
                             QS_ras - CMIP_qs_ras)
  writeRaster(Projections_stack, filename = file.path(Dir.Data, "Projections"), format = "CDF")
}else{
  print("Raw projection data already kriged")
}
Projections_stack <- stack(file.path(Dir.Data, "Projections.nc"))
names(Projections_stack) <- c("Tair.Historical", "Tair.SSP245", "Tair.Diff",
                              "Qsoil.Historical", "Qsoil.SSP245", "Qsoil.Diff")
# plot(Projections_stack)

# OCCURRENCE DATA ==========================================================
message("### OCCURRENCE DATA ###")

if(file.exists(file.path(Dir.D.Occurrences, "NonOcc_spec.RData"))){
  load(file.path(Dir.D.Occurrences, "NonOcc_spec.RData"))
}else{
  NonOcc_spec <- c()
}

occ_fs <- list.files(Dir.D.Occurrences)
occ_spec <- c(gsub(occ_fs, pattern = ".rds", replacement = ""), NonOcc_spec)
### Plants ####
print("PLANTS")
plants_spec <- unlist(lapply(List_ls, FUN = function(x){rownames(x)}))
if(sum(unique(plants_spec) %nin% occ_spec) > 0){
  plants_occ <- occ_data(scientificName = unique(plants_spec))
  plants_occ <- lapply(plants_occ, FUN = function(x){nrow(x[[2]])})
  Failed_plants <- names(plants_occ)[which(unlist(plants_occ) == 0)]
  if(length(Failed_plants) != 0){stop("Not all plant species are found on gbif")}
  plants_gbif <- Gbif_Species(species = plants_spec, year_vec = 1982:1999)
  print("Identifying outliers & saving occurrence data")
  NonOcc_spec <- c(NonOcc_spec, names(plants_gbif[lapply(plants_gbif, nrow) == 0]))
  save(NonOcc_spec, file = file.path(Dir.D.Occurrences, "NonOcc_spec.RData"))
  plants_gbif <- plants_gbif[lapply(plants_gbif, nrow) != 0] # remove species for which no records are present
  ## remove species for which 20 records or less are present
  n_occ <- unlist(lapply(plants_gbif, nrow))
  NonOcc_spec <- c(NonOcc_spec, names(plants_gbif)[n_occ <= 20])
  save(NonOcc_spec, file = file.path(Dir.D.Occurrences, "NonOcc_spec.RData"))
  plants_gbif <- plants_gbif[n_occ > 20]
  ## saving data
  hush_ls <- pblapply(names(plants_gbif), function(df){
    x <- na.omit(plants_gbif[[df]][, c("key", "decimalLatitude", "decimalLongitude")])
    x <- Gbif_Outliers(x = x, Enviro_ras = Enviro_ras, Centroids = Shapes_ct)
    saveRDS(x, file = file.path(Dir.D.Occurrences, paste0(df, ".rds")))
  }) 
}else{
  print("No new plant species")
}

### Animals ####
print("ANIMALS")
animals_spec <- unlist(lapply(List_ls, FUN = function(x){colnames(x)}))
if(sum(unique(animals_spec) %nin% occ_spec) > 0){
  animals_occ <- occ_data(scientificName = unique(animals_spec))
  animals_occ <- lapply(animals_occ, FUN = function(x){nrow(x[[2]])})
  Failed_animals <- names(animals_occ)[which(unlist(animals_occ) == 0)]
  if(length(Failed_animals) != 0){stop("Not all animal species are found on gbif")} 
  animals_gbif <- Gbif_Species(species = animals_spec, year_vec = 1982:1999)
  print("Saving occurrence data")
  NonOcc_spec <- c(NonOcc_spec, names(animals_gbif[lapply(animals_gbif, nrow) == 0]))
  save(NonOcc_spec, file = file.path(Dir.D.Occurrences, "NonOcc_spec.RData"))
  animals_gbif <- animals_gbif[lapply(animals_gbif, nrow) != 0] # remove species for which no records are present
  ## remove species for which 20 records or less are present
  n_occ <- unlist(lapply(animals_gbif, nrow))
  NonOcc_spec <- c(NonOcc_spec, names(animals_gbif)[n_occ <= 20])
  save(NonOcc_spec, file = file.path(Dir.D.Occurrences, "NonOcc_spec.RData"))
  animals_gbif <- animals_gbif[n_occ > 20]
  ## saving data
  hush_ls <- pblapply(names(animals_gbif), function(df){
    x <- na.omit(animals_gbif[[df]][, c("key", "decimalLatitude", "decimalLongitude")])
    x <- Gbif_Outliers(x = x, Enviro_ras = Enviro_ras, Centroids = Shapes_ct)
    saveRDS(x, file = file.path(Dir.D.Occurrences, paste0(df, ".rds")))
  }) 
}else{
  print("No new animal species")
}

print("Loading occurrence data")
occ_fs <- list.files(Dir.D.Occurrences, full.names = TRUE, pattern = ".rds")
occ_ls <- as.list(pbsapply(occ_fs, readRDS))
occ_spec <- list.files(Dir.D.Occurrences, pattern = ".rds")
names(occ_ls) <- gsub(occ_spec, pattern = ".rds", replacement = "")

## Removing species for which no occurrences were found --------------------
message("## CLEANING NETWORKS OFF SPECIES WITHOUT OCCURRENCES")
List_ls <- pblapply(List_ls, FUN = function(x, spec){x[rownames(x) %nin% NonOcc_spec , colnames(x) %nin% NonOcc_spec]}, spec = NonOcc_spec)
plants_spec <- unlist(lapply(List_ls, FUN = function(x){rownames(x)}))
animals_spec <- unlist(lapply(List_ls, FUN = function(x){colnames(x)}))


# TRAIT DATA ===============================================================
message("### TRAIT DATA ###")
if(file.exists(file.path(Dir.Data, "Traits.RData"))){
  print("Traits already extracted")
  load(file.path(Dir.Data, "Traits.RData"))
}else{
  print("Extracting traits")
  colnames(int.set)[5:29]
  traits_df <- int.set[int.set$animal.phylo.id %in% animals_spec & int.set$plant.phylo.id %in% plants_spec, ]
  plant_means <- aggregate(traits_df[18:25], by=list(Species=traits_df$plant.phylo.id), FUN=mean)
  animals_means <- aggregate(traits_df[18:25], by=list(Species=traits_df$plant.phylo.id), FUN=mean)
  save(traits_df, animals_means, plant_means, file = file.path(Dir.Data, "Traits.RData"))
}

# CLIMATE PREFERENCES ======================================================
message("### CLIMATE PREFERENCES")
if(file.exists(file.path(Dir.Data, "ClimPrefs.RData"))){
  load(file.path(Dir.Data, "ClimPrefs.RData"))
  NewPrefs_pos <- which(names(occ_ls) %nin% Preferences_df$spec)
  if(length(NewPrefs_pos)>0){
    occ_ls2 <- occ_ls[NewPrefs_pos]
    Preferences_df <- Clim_Preferences(data = occ_ls, Enviro_ras = Enviro_ras, Outliers = TRUE, Boot = 1e3)
    save(Preferences_df, file = file.path(Dir.Data, "ClimPrefs.RData"))
  }else{
    print("Already calculated climatic preferences for all species")
  }
}else{
  Preferences_df <- Clim_Preferences(data = occ_ls, Enviro_ras = Enviro_ras, Outliers = TRUE, Boot = 1e3)
  save(Preferences_df, file = file.path(Dir.Data, "ClimPrefs.RData"))
}

# EXTINCTION PROXIES =======================================================
message("### EXTINCTION PROXIES")
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






















