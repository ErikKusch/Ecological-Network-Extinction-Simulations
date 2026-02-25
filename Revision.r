#' ####################################################################### #
#' PROJECT: [Biodiversity Simplification & Ecological Network Topology]
#' CONTENTS:
#'  - INaturalist taxon matching
#'  - INaturalist distribution download
#'  - Extinction Risk Calculation
#'  - Comparison of satefy margins with different approaches
#'  DEPENDENCIES:
#'  - "0 - Preamble.R"
#'  - "0 - Fricke_Functions.R"
#'  - "0 - Data_Functions.R"
#'  - "1- DataRetrieval.R" run
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE =================================================================
rm(list = ls())
set.seed(42)

## Sourcing ----------------------------------------------------------------
source("0 - Preamble.R")
source("0 - Fricke_Functions.R")
source("0 - Data_Functions.R")

## Data Loading ------------------------------------------------------------
load(file.path(Dir.Data, "AnalysesData.RData")) # load the data needed for the analysis
## data now loaded:
#' - List_ls: List object of which each element is a frugivory adjacency matrix with plants as rows and animals as columns
#' - networks_df: Metadata for all networks expressed as adjacency matrices in List_ls (column "net.id" in networks_df can be matched to List_ls names)
#' - Prox.Centrality_ls: List object of which each element is a named vector that sorts the species for each network in our study from most to least central (as measured by connection strength)
rm(Prox.Centrality)
#' - Prox.Climate_ls: List object of which each element a list containing a named vector that sorts the species for each network in our study from most to least at risk from a climate-standpoint and a dataframe showing the climate risk for each species for temperature and soil moisture
#' - Prox.IUCN_df: data frame containing IUCN criteria for all species in our analyses. NOT ANYMORE
rm(Prox.IUCN_df)
#' - traits_df: data frame of trait expressions per species in our analysis
#' - animals_gowdis: square matrix of animal species dissimilarity in trait space
#' - plants_gowdis: square matrix of plant species dissimilarity in trait space


# CLIMATE SAFETY MARGIN REVISION ===========================================
## Paper-Contained Climate Preferences -------------------------------------
load(file.path(Dir.Data, "ClimPrefs.RData")) # loads Preferences_df
load(file.path(Dir.Data, "Prox_Climate.RData")) # loads Prox.Climate_ls; a list of two elements (ssp245 and ssp585) that are lists where each element corresponds to a network and containg `Order` (order of primary extinctions) and `ClimRisk` (a data frame reporting climate safety margins for each species)

## INaturalist Climate Safety Margins --------------------------------------
### INaturalist ranges ++++++++++
if (file.exists(file.path(Dir.Data, "Revision_INAT.RData"))) {
    load(file.path(Dir.Data, "Revision_INAT.RData"))
} else {
    INatIDs_ls <- pblapply(
        Preferences_df$spec,
        # cl = cl,
        FUN = function(SpecIter) {
            # message(SpecIter)
            FNAME <- paste0("INAT_", SpecIter, ".RData")
            if (file.exists(FNAME)) {
                load(FNAME)
            } else {
                inat_id <- tryCatch(
                    {
                        get_inat_obs(taxon_name = SpecIter, maxresults = 1e2)
                    }, # which(inat_df$species_guess == SpecIter)[1]
                    error = function(e) {
                        e
                    }
                )
                Sys.sleep(10) # try not to overload INaturalist servers
                save(inat_id, file = FNAME)
                # print(inat_id)
            }
            inat_id
        }
    )
    names(INatIDs_ls) <- Preferences_df$spec
    ObjectClasses <- unlist(lapply(INatIDs_ls, FUN = function(iter) {
        class(iter)[1]
    }))

    Success_spec <- Preferences_df$spec[which(ObjectClasses == "data.frame")]

    TaxonIDs_ls <- pblapply(Success_spec, FUN = function(Iter_spec) {
        # Iter_spec <- Success_spec[1]
        Inat_df <- INatIDs_ls[[Iter_spec]]

        ExactMatch <- which(Inat_df$scientific_name == Iter_spec)[1]

        if (is.na(ExactMatch)) {
            data.frame(
                taxonID = names(table(Inat_df$taxon_id))[which.max(table(Inat_df$taxon_id))],
                matchtype = "INatAssigned"
            )
        } else {
            data.frame(
                taxonID = Inat_df$taxon_id[ExactMatch],
                matchtype = "Exact"
            )
        }
    })
    names(TaxonIDs_ls) <- Success_spec
    INatTaxon_df <- do.call(rbind, TaxonIDs_ls)

    ## INaturalist Maps --------------------------------------------------------
    INatMaps_ls <- pblapply(INatTaxon_df$taxonID[INatTaxon_df$matchtype == "Exact"], FUN = function(taxonIter) {
        print(taxonIter)
        tryCatch(
            {
                geojson_data <- st_read(paste0("https://inaturalist-open-data.s3.us-east-1.amazonaws.com/geomodel/geojsons/latest/", taxonIter, ".geojson"), quiet = TRUE)
            }, # which(inat_df$species_guess == SpecIter)[1]
            error = function(e) {
                e
            }
        )
    })
    names(INatMaps_ls) <- names(TaxonIDs_ls)[INatTaxon_df$matchtype == "Exact"]
    ObjectClasses <- unlist(lapply(INatMaps_ls, FUN = function(iter) {
        class(iter)[1]
    }))

    Success_spec <- names(INatMaps_ls)[which(ObjectClasses == "sf")]

    INatRanges_ls <- INatMaps_ls[which(ObjectClasses == "sf")]

    save(INatRanges_ls, INatTaxon_df, file = file.path(Dir.Data, "Revision_INAT.RData"))
    unlink(paste0("INAT_", Preferences_df$spec, ".RData"))
}

### Safety Margin Computation ++++++++++
ssps <- c("ssp245", "ssp585")

FNAME <- file.path(Dir.Data, "Revision_INAT_ClimateSafetyMargins.RData")
if (file.exists(FNAME)) {
    load(FNAME)
} else {
    ## Loading Environmental Data
    Enviro_ras <- raster::stack(file.path(Dir.Data, "Enviro_Pres.nc"))
    krigs_ls <- as.list(c(NA, NA))
    names(krigs_ls) <- ssps
    for (ssp in ssps) {
        load(file.path(Dir.Data, paste0("Projections", ssp, ".RData")))
        names(Projections_stack[[1]]) <- c("Tair.Historical", "Tair.SSP", "Tair.Diff")
        names(Projections_stack[[2]]) <- c("Qsoil.Historical", "Qsoil.SSP", "Qsoil.Diff")
        krigs_ls[[ssp]] <- Projections_stack
    }

    ## calculate preferences
    preferences_df <- pblapply(names(INatRanges_ls), FUN = function(Spec_Iter) {
        # Spec_Iter <- names(INatRanges_ls)[1]

        extracted_df <- raster::extract(Enviro_ras, INatRanges_ls[[Spec_Iter]])
        extracted_df <- do.call(rbind, (extracted_df))
        preferences_df <- data.frame(
            spec = Spec_Iter,
            Temp_mean = mean(extracted_df[, 1], na.rm = TRUE),
            Temp_sd = sd(extracted_df[, 1], na.rm = TRUE),
            Water_mean = mean(extracted_df[, 2], na.rm = TRUE),
            Water_sd = sd(extracted_df[, 2], na.rm = TRUE)
        )
    })
    preferences_df <- do.call(rbind, preferences_df)

    ProxClim_ls <- as.list(c(NA, NA))
    names(ProxClim_ls) <- ssps
    for (ssp in ssps) {
        if (!file.exists(file.path(Dir.Data, paste0("Revision_INAT_Prox_Climate", ssp, ".RData")))) {
            Prox.Climate_ls <- pblapply(names(List_ls), function(netID) {
                # netID <- names(List_ls)[1]
                print(netID)
                ## Species Identities
                Plants_spec <- rownames(List_ls[[netID]])
                Animals_spec <- colnames(List_ls[[netID]])

                ## Network position
                extract_df <- networks_df[networks_df$net.id == netID, ]
                coordinates(extract_df) <- ~ longitude + latitude

                ## Environmental differences at network location
                Present <- raster::extract(Enviro_ras$X1, extract_df, method = "bilinear")
                if (is.na(Present)) {
                    Present <- mean(unlist(raster::extract(Enviro_ras$X1, extract_df, buffer = 1e4)), na.rm = TRUE)
                }
                TairDiff <- Present +
                    raster::extract(krigs_ls[[ssp]][[1]]$Tair.Diff, extract_df, method = "bilinear")
                Present <- raster::extract(Enviro_ras$X2, extract_df, method = "bilinear")
                if (is.na(Present)) {
                    Present <- mean(unlist(raster::extract(Enviro_ras$X2, extract_df, buffer = 1e4)), na.rm = TRUE)
                }
                QsoilDiff <- Present +
                    raster::extract(krigs_ls[[ssp]][[2]]$Qsoil.Diff, extract_df, method = "bilinear")

                ## calculation of climate stress for each species
                Prox_df <- data.frame(
                    species = c(Plants_spec, Animals_spec),
                    Tair = NA,
                    Qsoil = NA
                )
                for (speciesIter in Prox_df$species) {
                    PrefIter_df <- preferences_df[preferences_df$spec == speciesIter, ]
                    if (nrow(PrefIter_df) != 0) {
                        Prox_df[Prox_df$species == speciesIter, 2:3] <- c(
                            (PrefIter_df$Temp_mean - TairDiff) / PrefIter_df$Temp_sd,
                            (PrefIter_df$Water_mean - QsoilDiff) / PrefIter_df$Water_sd
                        )
                    }
                }

                ## creating order of extinction risk /climate stress
                Order_df <- na.omit(Prox_df)
                Order_df$Qsoil[Order_df$species %in% Animals_spec] <- 0 # no consideration for Qsoil effects on animals
                Order_df[, 2:3] <- abs(Order_df[, 2:3])
                WhichMax <- apply(Order_df[, 2:3], MARGIN = 1, FUN = which.max)
                Order_vec <- sapply(1:nrow(Order_df), FUN = function(x) {
                    Order_df[x, WhichMax[x] + 1]
                })
                names(Order_vec) <- Order_df$species

                ## saving climate proxies
                list(
                    Order = sort(Order_vec, decreasing = TRUE),
                    ClimRisk = Prox_df
                )
            })
            names(Prox.Climate_ls) <- names(List_ls)
            save(Prox.Climate_ls, file = file.path(Dir.Data, "Prox_Climate.RData"))
        }
        load(file.path(Dir.Data, "Prox_Climate.RData"))
        ProxClim_ls[[ssp]] <- Prox.Climate_ls
    }
    Prox.INat.Climate_ls <- ProxClim_ls
    save(Prox.INat.Climate_ls, preferences_df, file = FNAME)
}

## Safety Margin Comparison ------------------------------------------------



# PHYLOGENETIC SIGNAL ON TRAITS ============================================
# compare gowdis (which is already calculated) to phylogenetic distance of animals and plants respectively
