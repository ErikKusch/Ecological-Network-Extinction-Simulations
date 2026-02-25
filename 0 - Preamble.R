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
CreateDir <- sapply(DataDirs, function(x) if (!dir.exists(x)) dir.create(x))
## EXPORTS ------------------------------------------------------------------
Dir.Exports <- file.path(Dir.Base, "Exports")
# DirEx.Observations <- file.path(Dir.Exports, "Observations")
ExportDirs <- c(Dir.Exports)
CreateDir <- sapply(ExportDirs, function(x) if (!dir.exists(x)) dir.create(x))
rm(list = c("CreateDir", "ExportDirs", "DataDirs"))

# PACKAGES ================================================================
try(source("X - PersonalSettings.R"))

## CRAN -------------------------------------------------------------------
# devtools::install_github("ErikKusch/NetworkExtinction", ref = "ErikDevel")
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}

package_vec <- c(
  "devtools", # needed for non-cran packages further down
  # "rgeos", # for loading shapefiles
  "tidyverse", # for data handling
  "rgbif", # for occurrence retrieval
  "pbapply", # for apply with progress bar
  "parallel", # for clusterCall
  "data.table", # for data handling
  "rnaturalearth", # for landmask in projection kriging
  "rnaturalearthdata", # for landmask in projection kriging
  "rinat", # for downloading INaturalist data (needed for revision)
  # "rredlist", # for IUCN risk retrieval
  # "ConR", # for computation of IUCN risks
  # "CoordinateCleaner", # for additional occurrence cleaning; NOT NEEDED ANYMORE
  "igraph", # for graph operations
  "FD", # for gower distance of trait data
  "reshape2", # for making network matrices into plottable data frames via melt()
  "bipartite", # for bipartite network analyses
  "leaflet", # for html map products to investigate networks separately
  "leafpop", # for graph popups in leaflet output
  "cowplot", # for arranging of plots
  "gridExtra", # for table grobs as legends in plots
  "dplyr", # for data cleaning
  "ggpubr", # for t-test comparisons in ggplots
  "viridis", # for extra colours in ggplot
  "ggvenn", # for venn diagrams
  "randomForest", # for classification of associations potential
  "brms", # for post-simulation bayesian models with zero-inflated beta
  "tidybayes", # for brms output visualisations
  "NetworkExtinction", # for network extinction simulations
  "boot", # for bootstrapping of niches
  "iterators",
  "sf",
  "terra",
  "sp"
)
sapply(package_vec, install.load.package)

## KrigR ------------------------------------------------------------------
if ("KrigR" %in% rownames(installed.packages()) == FALSE) { # KrigR check
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
  devtools::install_github("https://github.com/ErikKusch/KrigR", force = TRUE)
}
library(KrigR)

## API Credentials --------------------------------------------------------
# CDS API (needed for ERA5-Land downloads)
if (!exists("API_Key") | !exists("API_User")) { # CS API check: if CDS API credentials have not been specified elsewhere
  API_User <- readline(prompt = "Please enter your Climate Data Store API user number and hit ENTER.")
  API_Key <- readline(prompt = "Please enter your Climate Data Store API key number and hit ENTER.")
} # end of CDS API check
# NUMBER OF CORES
if (!exists("numberOfCores")) { # Core check: if number of cores for parallel processing has not been set yet
  numberOfCores <- as.numeric(readline(prompt = paste("How many cores do you want to allocate to these processes? Your machine has", parallel::detectCores())))
} # end of Core check

if (as.character(options("gbif_user")) == "NULL") {
  Register <- readline(prompt = "Please enter your GBIF user name")
  options(gbif_user = Register)
}

if (as.character(options("gbif_email")) == "NULL") {
  Register <- readline(prompt = "Please enter your GBIF email")
  options(gbif_email = Register)
}

if (as.character(options("gbif_user")) == "NULL") {
  Register <- readline(prompt = "Please enter your GBIF password")
  options(gbif_pwd = "Register")
}

# if(!exists("IUCN_Key")){
#   IUCN_Key <- readline(prompt = "Please enter your IUCN API key.")
# }

# FUNCTIONALITY =============================================================
`%nin%` <- Negate(`%in%`)

hush <- function(code) {
  sink("NUL") # use /dev/null in UNIX
  tmp <- code
  sink()
  return(tmp)
}

Sort.DF <- function(Data = NULL, Column = NULL, decreasing = FALSE) {
  Data[order(Data[, Column], decreasing = decreasing), ]
}

substrRight <- function(x, n) {
  substr(x, nchar(x) - n + 1, nchar(x))
}

# loads .RData objects to a specified named object
loadRData <- function(fileName) {
  # loads an RData file, and returns it
  get(load(fileName))
}

# return name of an object
objName <- function(z) {
  deparse(substitute(z))
}

# add a list of matrices together
add_matrices <- function(a) {
  cols <- sort(unique(unlist(lapply(a, colnames))))
  rows <- sort(unique(unlist(lapply(a, rownames))))
  out <- array(0, dim = c(length(rows), length(cols)), dimnames = list(rows, cols))
  for (m in a) out[rownames(m), colnames(m)] <- out[rownames(m), colnames(m)] + m
  out
}

# source only part of a file
source2 <- function(file, start, end, ...) {
  file.lines <- scan(file, what = character(), skip = start - 1, nlines = end - start + 1, sep = "\n")
  file.lines.collapsed <- paste(file.lines, collapse = "\n")
  source(textConnection(file.lines.collapsed), ...)
}
