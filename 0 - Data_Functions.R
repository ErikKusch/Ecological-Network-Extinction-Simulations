#' ####################################################################### #
#' PROJECT: [Biodiversity Simplification & Ecological Network Topology] 
#' CONTENTS: 
#'  - Gbif Data Retrieval
#'  DEPENDENCIES:
#'  - 
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# GBIF Occurrences =========================================================
Gbif_Species <- function(species = NULL, year_vec = 2000:2020){
  n_species <- length(unique(species))
  for(k in year_vec){
    message(paste("Year", k, "/" , max(year_vec)))
    Gbif <- NA
    while(is.na(Gbif)){
      try(Gbif <- occ_data(scientificName = unique(species),
                           hasGeospatialIssue = FALSE,
                           hasCoordinate = TRUE,
                           year = k,
                           limit = 1e5
      ))
    }
    if(n_species > 1){
      Data_ls <- lapply(Gbif, '[[', 2) # extract downloaded gbif data into list object, each element corresponds to a genuskey
      if(k == year_vec[[1]]){
        Data <- Data_ls
      }else{
        Data <- pblapply(X = 1:n_species, FUN = function(x){
          data.frame(
            rbindlist(list(Data[[x]],
                           Data_ls[[x]]),
                      fill = TRUE)
          )
        }
        )
      }
      names(Data) <- names(Data_ls)
    }else{
      if(k == year_vec[[1]]){
        Data <- as.data.frame(Gbif$data)
      }else{
        Data <- data.frame(rbindlist(list(Data, as.data.frame(Gbif$data)), fill = TRUE))
      }
    }
  }
  return(Data)
}

# GBIF Outliers ============================================================
Gbif_Outliers <- function(x = NULL, Enviro_ras = Enviro_ras, Centroids = Shapes_ct){
  coordinates(x) <- ~ decimalLongitude + decimalLatitude
  ## Environmental outliers, classified as beyond mean +/- 4*sd
  enviro_vals <- data.frame(raster::extract(Enviro_ras, x))
  enviro_vals$key <- x$key
  enviro_vals <- na.omit(enviro_vals)
  x <- x[x$key %in% enviro_vals$key, ]
  sd_vec <- apply(enviro_vals[,-3], MARGIN = 2, FUN = sd, na.rm = TRUE)
  mean_vec <- apply(enviro_vals[,-3], MARGIN = 2, FUN = mean, na.rm = TRUE)
  x$Out_Enviro <- FALSE
  Out_pos <- which(enviro_vals[,1] > mean_vec[1]+4*sd_vec[1] | enviro_vals[,1] < mean_vec[1]-4*sd_vec[1] |
                     enviro_vals[,2] > mean_vec[2]+4*sd_vec[2] | enviro_vals[,2] < mean_vec[2]-4*sd_vec[2])
  if(length(Out_pos)>0){x$Out_Enviro[Out_pos] <- TRUE}
  ## Centroids, matched to centroids defined above and coordinates rounded to fourth digit
  x$Out_Centroid <- FALSE
  Out_pos <- which(round(coordinates(x), 4)[,1] %in% round(coordinates(Centroids), 4)[,1] &
                     round(coordinates(x), 4)[,2] %in% round(coordinates(Centroids), 4)[,2])
  if(length(Out_pos)>0){x$Out_Centroid[Out_pos] <- TRUE}
  return(x)
}

# Climate Preferences ======================================================
Clim_Preferences <- function(data = occ_ls, Enviro_ras = Enviro_ras, Outliers = TRUE, Boot = 1e3){
  occ_spd <- do.call(rbind, data)
  occ_spd$spec <- rep(names(data), unlist(lapply(data, nrow)))
  spec_vec <- unique(occ_spd$spec)
  preferences_df <- data.frame(spec = NA,
                               Temp_median = NA, 
                               Temp_sd = NA,
                               Water_median = NA,
                               Water_sd = NA) 
  set.seed(42)
  counter <- 1
  pb <- txtProgressBar(min = 0, max = length(spec_vec), style = 3)
  
  for(spec in spec_vec){
    occ_iter <- occ_spd[occ_spd$spec == spec, ]
    if(Outliers){
      NonOut_pos <- rowSums(as.data.frame(occ_iter)[, 4:5]) == 0 # 4:5 are the outlier colums here
      occ_iter <- occ_iter[NonOut_pos,]
    }
    extract_df <- na.omit(raster::extract(Enviro_ras, occ_iter))
    values_df <- data.frame(
      X1 = NA,
      X2 = NA
    )
    for(i in 1:Boot){
      rows <- sample(x = 1:nrow(extract_df), size = round(nrow(extract_df)*0.7), replace = TRUE)
      values_df <- rbind(values_df, extract_df[rows, ])
    }
    values_df <- na.omit(values_df)
    
    preferences_df <- rbind(preferences_df, 
                            c(spec,
                              apply(values_df, 2, median)[1], apply(values_df, 2, sd)[1],
                              apply(values_df, 2, median)[2], apply(values_df, 2, sd)[2])
    )
    setTxtProgressBar(pb, counter)
    counter <- counter + 1
  }
  preferences_df <- na.omit(preferences_df)
  return(preferences_df)
}
