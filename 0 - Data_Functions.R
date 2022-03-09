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
