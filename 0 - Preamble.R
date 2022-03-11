#' ####################################################################### #
#' PROJECT: [Biodiversity Simplification & Ecological Network Topology] 
#' CONTENTS: 
#'  - Directory Establishment
#'  - Package Loading
#'  - Helper Functions
#'  DEPENDENCIES:
#'  - 
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# DIRECTORIES ===============================================================
Dir.Base <- getwd() # read out the project directory
## DATA ---------------------------------------------------------------------
Dir.Data <- file.path(Dir.Base, "Data")
Dir.D.Fricke <- file.path(Dir.Data, "Fricke2021")
Dir.D.Climatologies <- file.path(Dir.Data, "Climatologies")
Dir.D.Projections <- file.path(Dir.Data, "Projections")
Dir.D.Occurrences <- file.path(Dir.Data, "Occurrences")
DataDirs <- c(Dir.Data, Dir.D.Climatologies, Dir.D.Fricke, Dir.D.Occurrences)
CreateDir <- sapply(DataDirs, function(x) if(!dir.exists(x)) dir.create(x))
## EXPORTS ------------------------------------------------------------------
Dir.Exports <- file.path(Dir.Base, "Exports")
# DirEx.Observations <- file.path(Dir.Exports, "Observations")
ExportDirs <- c(Dir.Exports)
CreateDir <- sapply(ExportDirs, function(x) if(!dir.exists(x)) dir.create(x))
rm(list = c("CreateDir", "ExportDirs", "DataDirs"))

# PACKAGES ================================================================
## CRAN -------------------------------------------------------------------
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}

package_vec <- c(
  "devtools", # needed for non-cran packages further down
  "rgeos", # for loading shapefiles
  "tidyverse", # for data handling
  "rgbif", # for occurrence retrieval
  "pbapply", # for apply with progress bar
  "data.table", # for data handling
  "rnaturalearth", # for landmask in projection kriging
  "rnaturalearthdata" # for landmask in projection kriging
  
)
sapply(package_vec, install.load.package)

## KrigR ------------------------------------------------------------------
if("KrigR" %in% rownames(installed.packages()) == FALSE){ # KrigR check
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
  devtools::install_github("https://github.com/ErikKusch/KrigR", force = TRUE)
}
library(KrigR) 
try(source("X - PersonalSettings.R")) # I do this here to specify number of cores and API credentials and am thus not sharing this file
# CDS API (needed for ERA5-Land downloads)
if(!exists("API_Key") | !exists("API_User")){ # CS API check: if CDS API credentials have not been specified elsewhere
  API_User <- readline(prompt = "Please enter your Climate Data Store API user number and hit ENTER.")
  API_Key <- readline(prompt = "Please enter your Climate Data Store API key number and hit ENTER.")
} # end of CDS API check
# NUMBER OF CORES
if(!exists("numberOfCores")){ # Core check: if number of cores for parallel processing has not been set yet
  numberOfCores <- as.numeric(readline(prompt = paste("How many cores do you want to allocate to these processes? Your machine has", parallel::detectCores())))
} # end of Core check

# FUNCTIONALITY =============================================================
`%nin%` <- Negate(`%in%`)

hush <- function(code){
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}

Sort.DF <- function(Data = NULL, Column = NULL){
  Data[order(Data[ , Column] ), ]
}