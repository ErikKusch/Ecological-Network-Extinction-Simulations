#' ####################################################################### #
#' PROJECT: [Biodiversity Simplification & Ecological Network Topology] 
#' CONTENTS: 
#'  - Data Download
#'  - Data Manipulation
#'  - Extinction Risk Calculation
#'  DEPENDENCIES:
#'  - "0 - Preamble.R"
#'  - "0 - Fricke_Functions.R"
#'  - "0 - Data_Functions.R"
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
message("### SPATIAL DATA ###")
## BUFFER CREATION ---------------------------------------------------------
nets_df <- metanet[metanet$net.id %in% names(List_ls), ]
colnames(nets_df)[3:4] <- c("Lat", "Lon")
nets_shp <- KrigR:::buffer_Points(nets_df, Buffer = 5, ID = "net.id")

## LANDMASK ----------------------------------------------------------------
Countries_shp <- ne_countries(scale = 10, type = "countries") # may need to run: devtools::install_github("ropensci/rnaturalearthhires")
CountriesCoarse_shp <- ne_countries(scale = 50, type = "countries")
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

## PROTECTED AREAS ---------------------------------------------------------
# CITATION: UNEP-WCMC and IUCN (2022), Protected Planet: The World Database on Protected Areas (WDPA) and World Database on Other Effective Area-based Conservation Measures (WD-OECM) [Online], March 2022, Cambridge, UK: UNEP-WCMC and IUCN. Available at: www.protectedplanet.net.
Dir.D.Protected <- file.path(Dir.Data, "Protected Areas")
if(!dir.exists(Dir.D.Protected)){dir.create(Dir.D.Protected)}
if(!file.exists(file.path(Dir.D.Protected, "ProtectedAreas.RData"))){
  download.file(url = "https://d1gam3xoknrgr2.cloudfront.net/current/WDPA_WDOECM_Mar2022_Public_all_shp.zip", 
                destfile = file.path(Dir.D.Protected, "WDPA_WDOECM_Mar2022_Public_all_shp.zip"))
  unzip(zipfile = file.path(Dir.D.Protected, "WDPA_WDOECM_Mar2022_Public_all_shp.zip"), 
        exdir = Dir.D.Protected,
        files = c("WDPA_WDOECM_Mar2022_Public_all_shp_0.zip",
                  "WDPA_WDOECM_Mar2022_Public_all_shp_1.zip",
                  "WDPA_WDOECM_Mar2022_Public_all_shp_2.zip"))
  ## repo1
  unzip(zipfile = file.path(Dir.D.Protected, "WDPA_WDOECM_Mar2022_Public_all_shp_0.zip"),
        exdir = Dir.D.Protected)
  ProtectedAreas_shp <- readOGR(file.path(Dir.D.Protected, "WDPA_WDOECM_Mar2022_Public_all_shp-polygons.shp"))[,1]
  save(ProtectedAreas_shp, file = file.path(Dir.D.Protected, "ProtectedAreas.RData"))
  ## repo2
  unzip(zipfile = file.path(Dir.D.Protected, "WDPA_WDOECM_Mar2022_Public_all_shp_1.zip"),
        exdir = Dir.D.Protected)
  ProtectedAreas_shp <- rbind(ProtectedAreas_shp,
                              readOGR(file.path(Dir.D.Protected, "WDPA_WDOECM_Mar2022_Public_all_shp-polygons.shp"))[,1]
  )
  save(ProtectedAreas_shp, file = file.path(Dir.D.Protected, "ProtectedAreas.RData"))
  ## repo3
  unzip(zipfile = file.path(Dir.D.Protected, "WDPA_WDOECM_Mar2022_Public_all_shp_0.zip"),
        exdir = Dir.D.Protected)
  ProtectedAreas_shp <- rbind(ProtectedAreas_shp,
                              readOGR(file.path(Dir.D.Protected, "WDPA_WDOECM_Mar2022_Public_all_shp-polygons.shp"))[,1]
  )
  save(ProtectedAreas_shp, file = file.path(Dir.D.Protected, "ProtectedAreas.RData"))
  ## combine shapes
  # ProtectedAreas_shp <- raster::aggregate(ProtectedAreas_shp)
  # save(ProtectedAreas_shp, file = file.path(Dir.D.Protected, "ProtectedAreas.RData"))
  ## clean directory
  unlink(list.files(Dir.D.Protected, full.names = TRUE, 
                    pattern = "WDPA_WDOECM_Mar2022_Public_all_shp"))
}
load(file.path(Dir.D.Protected, "ProtectedAreas.RData"))

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
# plot(Enviro_ras[[1]])
# plot(nets_shp, add = TRUE)

## PROJECTIONS -------------------------------------------------------------
message("### PROJECTION DATA ###")
print("Loading projection data")
if(!file.exists(file.path(Dir.Data, "Projections.nc"))){
  print("Loading raw projection data")
  train_ERA <- Enviro_ras[[1]]
  train_ERA <- crop(train_ERA, extent(nets_shp))
  train_ERA <- mask(train_ERA, nets_shp)
  
  ### SSP ----
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
  
  ### HISTORICAL ----
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
}

## KRIGING -----------------------------------------------------------------
message("### PROJECTION KRIGING ###")
if(!file.exists(file.path(Dir.Data, "Projections.RData"))){
  print("Kriging raw projection data")
  ### TEMPERATURE ----
  #### Covariates ----
  GMTED <- download_DEM(
    Train_ras = train_SSP,
    Target_res = train_ERA,
    Shape = nets_shp,
    Keep_Temporary = FALSE,
    Dir = Dir.D.Projections
  )
  Cov_coarse <- GMTED[[1]]
  # Cov_coarse <- mask(Cov_coarse, nets_shp)
  Cov_fine <- GMTED[[2]]
  # Cov_fine <- mask(Cov_fine, nets_shp)
  
  #### SSP ----
  if(!file.exists(file.path(Dir.D.Projections, "K_ssp245_tas_2081-2100_nmax120.nc"))){
    Output_SSP <- krigR(
      Data = train_SSP$Temperature,
      Covariates_coarse = Cov_coarse, 
      Covariates_fine = Cov_fine,   
      KrigingEquation = "ERA ~ DEM",  
      Cores = numberOfCores, 
      Dir = Dir.D.Projections,  
      FileName = "K_ssp245_tas_2081-2100_nmax120", 
      Keep_Temporary = FALSE,
      nmax = 120
    )
  }
  TAS_ras <- raster(file.path(Dir.D.Projections, "K_ssp245_tas_2081-2100_nmax120.nc"))
  
  #### HISTORICAL ----
  if(!file.exists(file.path(Dir.D.Projections, "K_CMIP-HIST_tas_nmax120.nc"))){
    Output_SSP <- krigR(
      Data = train_HIST$Temperature,
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
  CMIP_tas_ras <- raster(file.path(Dir.D.Projections, "K_CMIP-HIST_tas_nmax120.nc"))
  
  ### SOIL MOISTURE ----
  #### Covariates ----
  if(!file.exists(file.path(Dir.D.Projections, "SoilCovs_ls.RData"))){
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
    Cov_coarse2 <- crop(SoilCovs_stack, extent(Cov_coarse))
    Cov_coarse2 <- resample(x = Cov_coarse2, y = Cov_coarse)
    range_m <- KrigR::mask_Shape(Cov_coarse2[[1]], nets_shp)
    Cov_coarse2 <- mask(Cov_coarse2, range_m)
    Cov_fine2 <- crop(SoilCovs_stack, extent(Cov_fine))
    Cov_fine2 <- resample(x = Cov_fine2, y = Cov_fine)
    range_m <- KrigR::mask_Shape(Cov_fine2[[1]], nets_shp)
    Cov_fine2 <- mask(Cov_fine2, range_m)
    SoilCovs_ls <- list(Coarse = Cov_coarse2,
                        Fine = Cov_fine2)
    save(SoilCovs_ls, file = file.path(Dir.D.Projections, "SoilCovs_ls.RData"))
    unlink(paste0(Dir.D.Projections, "/", SoilCovs_vec, ".nc"))
  }else{
    load(file.path(Dir.D.Projections, "SoilCovs_ls.RData"))
    Cov_coarse2 <- SoilCovs_ls[[1]]
    Cov_fine2 <- SoilCovs_ls[[2]]
  }
  
  #### SSP ----
  if(!file.exists(file.path(Dir.D.Projections, "K_ssp245_mrsos_2081-2100_nmax120.nc"))){
    QS_SSP <- krigR(
      Data = train_SSP$Moisture,
      Covariates_coarse = Cov_coarse2, 
      Covariates_fine = Cov_fine2,   
      KrigingEquation = "ERA ~ tkdry+tksat+csol+k_s+lambda+psi+theta_s",  
      Cores = numberOfCores, 
      Dir = Dir.D.Projections,  
      FileName = "K_ssp245_mrsos_2081-2100_nmax120", 
      Keep_Temporary = FALSE,
      nmax = 120
    )
  }
  MRSOS_ras <- mean(stack(file.path(Dir.D.Projections, "K_ssp245_mrsos_2081-2100_nmax120.nc"))) * 0.013
  MRSOS_ras <- reclassify(MRSOS_ras, cbind(-Inf, 0, 0))
  
  #### HISTORICAL ----
  if(!file.exists(file.path(Dir.D.Projections, "K_CMIP-HIST_mrsos_nmax120.nc"))){
    Output_SSP <- krigR(
      Data = train_HIST$Moisture,
      Covariates_coarse = Cov_coarse2, 
      Covariates_fine = Cov_fine2,   
      KrigingEquation = "ERA ~ tkdry+tksat+csol+k_s+lambda+psi+theta_s",  
      Cores = numberOfCores, 
      Dir = Dir.D.Projections,  
      FileName = "K_CMIP-HIST_mrsos_nmax120", 
      Keep_Temporary = FALSE,
      nmax = 120
    )
  }
  CMIP_mrsos_ras <- mean(stack(file.path(Dir.D.Projections, "K_CMIP-HIST_mrsos_nmax120.nc"))) * 0.013
  CMIP_mrsos_ras <- reclassify(CMIP_mrsos_ras, cbind(-Inf, 0, 0))
  
  ## Resolve slightly larger TAS rasters
  # stop("resolve extent issues")
  # TAS_ras <- crop(TAS_ras, extent(MRSOS_ras))
  # CMIP_tas_ras <- crop(CMIP_tas_ras, extent(MRSOS_ras))
  
  ### DIFFERENCE AND FUSING ----
  Projections_stack <- list(Temp = stack(CMIP_tas_ras,
                                         TAS_ras,
                                         TAS_ras - CMIP_tas_ras),
                            Water = stack(
                              CMIP_mrsos_ras,
                              MRSOS_ras,
                              MRSOS_ras - CMIP_mrsos_ras)
  )
  save(Projections_stack, file = file.path(Dir.Data, "Projections.RData"))
}else{
  print("Raw projection data already kriged")
}
load(file.path(Dir.Data, "Projections.RData"))
names(Projections_stack[[1]]) <- c("Tair.Historical", "Tair.SSP245", "Tair.Diff")
names(Projections_stack[[2]]) <- c("Qsoil.Historical", "Qsoil.SSP245", "Qsoil.Diff")
# plot(Projections_stack)

# OCCURRENCE DATA ==========================================================
message("### OCCURRENCE DATA ###")
## species for which no data is presented and no data could be obtained previously
if(file.exists(file.path(Dir.D.Occurrences, "NonOcc_spec.RData"))){
  load(file.path(Dir.D.Occurrences, "NonOcc_spec.RData"))
}else{
  NonOcc_spec <- c()
}
## species for which data is present
occ_fs <- list.files(Dir.D.Occurrences)
occ_spec <- c(tools::file_path_sans_ext(occ_fs), NonOcc_spec)
### Plants ####
print("PLANTS")
plants_spec <- unique(unlist(lapply(List_ls, FUN = function(x){rownames(x)})))
if(sum(plants_spec %nin% occ_spec) > 0){
  plants_occ <- occ_data(scientificName = plants_spec[plants_spec %nin% occ_spec])
  if(sum(plants_spec %nin% occ_spec) == 1){
    plants_occ <- list(nrow(plants_occ$data)) 
    names(plants_occ) <- plants_spec[plants_spec %nin% occ_spec]
  }else{
    plants_occ <- lapply(plants_occ, FUN = function(x){nrow(x[[2]])})
  }
  Failed_plants <- names(plants_occ)[which(unlist(plants_occ) == 0)]
  if(length(Failed_plants) != 0){stop("Not all plant species are found on gbif")}
  plants_gbif <- Gbif_Species(species = plants_spec[plants_occ %nin% occ_spec], year_vec = 1982:1999)
  if(sum(plants_spec %nin% occ_spec) == 1){
    plants_gbif <- list(plants_gbif)
    names(plants_gbif) <- plants_spec[plants_spec %nin% occ_spec]
  }
  print("Identifying outliers & saving occurrence data")
  NonOcc_spec <- c(NonOcc_spec, names(plants_gbif[lapply(plants_gbif, nrow) == 0]))
  save(NonOcc_spec, file = file.path(Dir.D.Occurrences, "NonOcc_spec.RData"))
  plants_gbif <- plants_gbif[lapply(plants_gbif, nrow) != 0] # remove species for which no records are present
  ## remove species for which 20 records or less are present
  n_occ <- unlist(lapply(plants_gbif, nrow))
  NonOcc_spec <- c(NonOcc_spec, names(plants_gbif)[n_occ <= 20])
  save(NonOcc_spec, file = file.path(Dir.D.Occurrences, "NonOcc_spec.RData"))
  plants_gbif <- plants_gbif[n_occ > 20]
  ## outliers
  # print("Identifying Outliers")
  # hush_ls <- lapply(names(plants_gbif), function(df){
  #   x <- na.omit(as.data.frame(plants_gbif[[df]])[, c("key", "decimalLatitude", "decimalLongitude")])
  #   x$species <- df
  #   x
  # })
  # plants_df <- na.omit(do.call(rbind, hush_ls))
  # plants_df <- clean_coordinates(plants_df,
  #                                lat = "decimalLatitude",
  #                                lon = "decimalLongitude",
  #                                verbose = FALSE)
  # plants_gbif <- split(plants_df, f = plants_df$species)
  ## saving data
  hush_ls <- pblapply(names(plants_gbif), function(df){
    x <- plants_gbif[[df]][,-4] # remove species column
    x <- Gbif_Outliers(x = x, Enviro_ras = Enviro_ras, Centroids = Shapes_ct)
    saveRDS(x, file = file.path(Dir.D.Occurrences, paste0(df, ".rds")))
  }) 
}else{
  print("No new plant species")
}

### Animals ####
print("ANIMALS")
animals_spec <- unique(unlist(lapply(List_ls, FUN = function(x){colnames(x)})))
if(sum(animals_spec %nin% occ_spec) > 0){
  animals_occ <- occ_data(scientificName = animals_spec[animals_spec %nin% occ_spec])
  if(sum(animals_spec %nin% occ_spec) == 1){
    animals_occ <- list(nrow(animals_occ$data)) 
    names(animals_occ) <- animals_spec[animals_spec %nin% occ_spec]
  }else{
    animals_occ <- lapply(animals_occ, FUN = function(x){nrow(x[[2]])})
  }
  Failed_animals <- names(animals_occ)[which(unlist(animals_occ) == 0)]
  if(length(Failed_animals) != 0){stop("Not all animal species are found on gbif")} 
  animals_gbif <- Gbif_Species(species = animals_spec[animals_spec %nin% occ_spec], year_vec = 1982:1999)
  if(sum(animals_spec %nin% occ_spec) == 1){
    animals_gbif <- list(animals_gbif)
    names(animals_gbif) <- animals_spec[animals_spec %nin% occ_spec]
  }
  NonOcc_spec <- c(NonOcc_spec, names(animals_gbif[lapply(animals_gbif, nrow) == 0]))
  save(NonOcc_spec, file = file.path(Dir.D.Occurrences, "NonOcc_spec.RData"))
  animals_gbif <- animals_gbif[lapply(animals_gbif, nrow) != 0] # remove species for which no records are present
  ## remove species for which 20 records or less are present
  n_occ <- unlist(lapply(animals_gbif, nrow))
  NonOcc_spec <- c(NonOcc_spec, names(animals_gbif)[n_occ <= 20])
  save(NonOcc_spec, file = file.path(Dir.D.Occurrences, "NonOcc_spec.RData"))
  animals_gbif <- animals_gbif[n_occ > 20]
  ## outliers
  # print("Identifying Outliers")
  # hush_ls <- lapply(names(animals_gbif), function(df){
  #   x <- na.omit(as.data.frame(animals_gbif[[df]])[, c("key", "decimalLatitude", "decimalLongitude")])
  #   x$species <- df
  #   x
  # })
  # animals_df <- na.omit(do.call(rbind, hush_ls))
  # animals_df <- clean_coordinates(animals_df,
  #                                 lat = "decimalLatitude",
  #                                 lon = "decimalLongitude",
  #                                 verbose = FALSE)
  # animals_gbif <- split(animals_df, f = animals_df$species)
  ## saving data
  print("Saving occurrence data")
  hush_ls <- pblapply(names(animals_gbif), function(df){
    x <- animals_gbif[[df]][,-4] # remove species column
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
names(occ_ls) <- tools::file_path_sans_ext(occ_spec)

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
  traits_df <- int.set[int.set$animal.phylo.id %in% animals_spec & int.set$plant.phylo.id %in% plants_spec, ]
  colnames(traits_df)
  plant_means <- aggregate(traits_df[22:29], by=list(Species=traits_df$plant.phylo.id), FUN=mean)
  plants_gowdis <- gowdis(plant_means)
  plants_gowdis <- as.matrix(plants_gowdis)
  dimnames(plants_gowdis) <- list(plant_means$Species, plant_means$Species)
  animal_means <- aggregate(traits_df[7:15], by=list(Species=traits_df$animal.phylo.id), FUN=mean)
  animals_gowdis <- gowdis(animal_means)
  animals_gowdis <- as.matrix(animals_gowdis)
  dimnames(animals_gowdis) <- list(animal_means$Species, animal_means$Species)
  save(traits_df, animal_means, animals_gowdis, plant_means, plants_gowdis, file = file.path(Dir.Data, "Traits.RData"))
}

# CLIMATE PREFERENCES ======================================================
message("### CLIMATE PREFERENCES ###")
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
Preferences_df$Temp_median <- as.numeric(Preferences_df$Temp_median)
Preferences_df$Temp_sd <- as.numeric(Preferences_df$Temp_sd)
Preferences_df$Water_median <- as.numeric(Preferences_df$Water_median)
Preferences_df$Water_sd <- as.numeric(Preferences_df$Water_sd)

# EXTINCTION PROXIES =======================================================
message("### EXTINCTION PROXIES ###")

## Network Centrality ------------------------------------------------------
message("## Network Centrality")
if(!file.exists(file.path(Dir.Data, "Prox_NetworkCentrality.RData"))){
  Prox.Centrality_ls <- pblapply(names(List_ls), FUN = function(x){
    # print(x)
    graph <- graph_from_incidence_matrix(List_ls[[x]], weighted = TRUE)
    sort(igraph::strength(graph), decreasing = TRUE)
  })
  names(Prox.Centrality_ls) <- names(List_ls)
  save(Prox.Centrality_ls, file = file.path(Dir.Data, "Prox_NetworkCentrality.RData"))
}else{
  print("Already computed")
  load(file.path(Dir.Data, "Prox_NetworkCentrality.RData"))
}

## Safety Margins ----------------------------------------------------------
message("## Climate Criteria")

if(!file.exists(file.path(Dir.Data, "Prox_Climate.RData"))){
  Prox.Climate_ls <- pblapply(names(List_ls), function(netID){
    # print(netID)
    ## Species Identities
    Plants_spec <- rownames(List_ls[[netID]])
    Animals_spec <- colnames(List_ls[[netID]])
    
    ## Network position
    extract_df <- networks_df[networks_df$net.id == netID, ]
    coordinates(extract_df) <- ~ longitude + latitude
    
    ## Environmental differences at network location
    Present <- extract(Enviro_ras$X1, extract_df, method = "bilinear")
    if(is.na(Present)){
      Present <- mean(unlist(extract(Enviro_ras$X1, extract_df, buffer = 1e4)), na.rm = TRUE)
    }
    TairDiff <- Present +
      extract(Projections_stack[[1]]$Tair.Diff, extract_df, method = "bilinear")
    Present <- extract(Enviro_ras$X2, extract_df, method = "bilinear")
    if(is.na(Present)){
      Present <- mean(unlist(extract(Enviro_ras$X2, extract_df, buffer = 1e4)), na.rm = TRUE)
    }
    QsoilDiff <- Present +
      extract(Projections_stack[[2]]$Qsoil.Diff, extract_df, method = "bilinear")
    
    ## calculation of climate stress for each species
    Prox_df <- data.frame(species = c(Plants_spec, Animals_spec),
                          Tair = NA,
                          Qsoil = NA)
    for(speciesIter in Prox_df$species){
      PrefIter_df <- Preferences_df[Preferences_df$spec == speciesIter, ]
      Prox_df[Prox_df$species == speciesIter, 2:3] <- c((PrefIter_df$Temp_median - TairDiff) / PrefIter_df$Temp_sd,
                                                        (PrefIter_df$Water_median - QsoilDiff) / PrefIter_df$Water_sd)
    }
    
    ## creating order of extinction risk /climate stress
    Order_df <- Prox_df
    Order_df[, 2:3] <- abs(Order_df[, 2:3])
    WhichMax <- apply(Order_df[, 2:3], MARGIN = 1, FUN = which.max)
    Order_vec <- sapply(1:nrow(Order_df), FUN = function(x){
      Order_df[x, WhichMax[x]+1]
    })
    names(Order_vec) <- Order_df$species
    
    ## saving climate proxies
    list(Order = sort(Order_vec, decreasing = TRUE),
         ClimRisk = Prox_df)
  })
  names(Prox.Climate_ls) <- names(List_ls)
  save(Prox.Climate_ls, file = file.path(Dir.Data, "Prox_Climate.RData"))
}
load(file.path(Dir.Data, "Prox_Climate.RData"))


## IUCN Criteria -----------------------------------------------------------
message("## IUCN Criteria")
if(!file.exists(file.path(Dir.Data, "Prox_IUCNCriteria.rds"))){
  Prox.IUCN_df <- data.frame(Species = names(occ_ls),
                             Category = NA,
                             Source = NA,
                             N = NA)
}else{
  Prox.IUCN_df <- readRDS(file.path(Dir.Data, "Prox_IUCNCriteria.rds"))
}
IUCN_spec <- names(occ_ls)[names(occ_ls) %nin% Prox.IUCN_df$Species[!is.na(Prox.IUCN_df$Category)]] # identify species for which IUCN records need to be retrieved
if(length(IUCN_spec) > 0){
  ### Direct IUCN retrieval ----
  print("Direct retrieval of IUCN criteria")
  IUCN_ls <- pblapply(IUCN_spec, FUN = function(x){
    rl_search(x, key = IUCN_Key)$result
  }
  )
  names(IUCN_ls) <- IUCN_spec
  Ident_df <- do.call(rbind, IUCN_ls[which(unlist(lapply(IUCN_ls, class)) == "data.frame")])
  if(!is.null(Ident_df)){
    Prox.IUCN_df$Category[match(Ident_df$scientific_name, Prox.IUCN_df$Species)] <- Ident_df$category
    Prox.IUCN_df$Source[match(Ident_df$scientific_name, Prox.IUCN_df$Species)] <- "IUCN"
    saveRDS(Prox.IUCN_df, file = file.path(Dir.Data, "Prox_IUCNCriteria.rds"))
  }
  
  ## Calculating IUCN criteria ----
  print("Computation of IUCN criteria where needed")
  Calc_spec <- Prox.IUCN_df$Species[is.na(Prox.IUCN_df$Category)]
  Calc_ls <- occ_ls[IUCN_spec %in% Calc_spec]
  Calc_ls <- lapply(names(Calc_ls), FUN = function(x){
    Calc_ls[[x]]$tax <- x
    Calc_ls[[x]]
  })
  Calc_df <- do.call(rbind, Calc_ls)
  Calc_df <- Calc_df[!Calc_df$Out_Enviro & !Calc_df$Out_Centroid, ]
  Calc_df <- as.data.frame(Calc_df)
  Calc_df <- Calc_df[,-c(1, 4, 5)]
  colnames(Calc_df) <- c("ddlat", "ddlon", "tax")
  ## check for which species longitude range does not exceed 180° (a limitation of the Conr package)
  print("Figuring out for which species entire data range can be used and subsequent IUCN criteria computation")
  Longi_check <- pbsapply(unique(Calc_df$tax), FUN = function(i){
    diff(range(Calc_df$ddlon[Calc_df$tax == i])) >= 180}
  )
  
  # save.image(file = "PreIUCNCalculation.RData") ## can be removed when all IUCN criteria have been computeds
  
  IUCN_calc <- IUCN.eval(DATA = Calc_df[Calc_df$tax %in% names(Longi_check)[!Longi_check], ], 
                         parallel = TRUE, NbeCores = ifelse(numberOfCores>10, 10, numberOfCores),
                         protec.areas = ProtectedAreas_shp, ID_shape_PA = "WDPAID"
  )
  Prox.IUCN_df$Source[match(IUCN_calc$tax, Prox.IUCN_df$Species)] <- "ConR"
  Prox.IUCN_df$Category[match(IUCN_calc$tax, Prox.IUCN_df$Species)] <- IUCN_calc$Category_CriteriaB
  Prox.IUCN_df$N[match(IUCN_calc$tax, Prox.IUCN_df$Species)] <- table(Calc_df$tax)
  saveRDS(Prox.IUCN_df, file = file.path(Dir.Data, "Prox_IUCNCriteria.rds"))
  ## limitation of species occurrences to the 180° longitude range for which there is the most data
  print("Figuring out for which species parts of the data range need to be used, selecting the most data-rich range, and subsequent IUCN criteria computation")
  Calc_df2 <- Calc_df[Calc_df$tax %in% names(Longi_check)[Longi_check], ]
  for(i in names(Longi_check)[Longi_check]){
    num_vec <- rep(NA, length(-180:1))
    names(num_vec) <- -180:1
    for(k in -180:1){
      num_vec[which(k == names(num_vec))] <- 
        sum(Calc_df2$ddlon[Calc_df2$tax == i] >= k & Calc_df2$ddlon[Calc_df2$tax == i] < k+180)
    }
    x <- as.numeric(names(which.max(num_vec)))
    Calc_df2$ddlon[Calc_df2$tax == i][
      Calc_df2$ddlon[Calc_df2$tax == i] < x |
        Calc_df2$ddlon[Calc_df2$tax == i] >= x+180] <- NA
  }
  Calc_df2 <- na.omit(Calc_df2)
  IUCN_calc2 <- IUCN.eval(DATA = Calc_df2, 
                          parallel = TRUE, NbeCores = ifelse(numberOfCores>10, 10, numberOfCores),
                          protec.areas = ProtectedAreas_shp, ID_shape_PA = "WDPAID"
  )
  Prox.IUCN_df$Source[match(IUCN_calc2$tax, Prox.IUCN_df$Species)] <- "ConR - Reduced"
  Prox.IUCN_df$Category[match(IUCN_calc2$tax, Prox.IUCN_df$Species)] <- IUCN_calc2$Category_CriteriaB
  Prox.IUCN_df$N[match(IUCN_calc2$tax, Prox.IUCN_df$Species)] <- table(Calc_df$tax)
  saveRDS(Prox.IUCN_df, file = file.path(Dir.Data, "Prox_IUCNCriteria.rds"))
}else{
  print("All IUCN criteria already obtained")
}

# SAVING ALL DATA AS ONE OBJECT ============================================
save(Prox.Climate_ls, Prox.IUCN_df, Prox.Centrality_ls, networks_df, List_ls, traits_df, animals_gowdis, plants_gowdis, animal_means, plant_means,
     file = file.path(Dir.Data, "AnalysesData.RData"))
