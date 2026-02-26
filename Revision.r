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
rm(Prox.Climate_ls)
#' - ProxClim_ls: List object containing Prox.Climate_ls for both ssps
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
### INaturalist ranges ----------
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

### Safety Margin Computation ----------
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
### Data Combination ----------
approaches_ls <- lapply(list(ProxClim_ls, Prox.INat.Climate_ls), FUN = function(Approach) {
    # Approach <- ProxClim_ls

    data_ls <- lapply(names(Approach), FUN = function(ssp) {
        # ssp <- names(Approach)[1]
        Nets_ls <- Approach[[ssp]]
        df_ls <- lapply(names(Nets_ls), FUN = function(index) {
            df <- Nets_ls[[index]]$ClimRisk
            # df$Approach <- "GBIF"
            df$ssp <- ssp
            df$net.id <- index
            df
        })
        do.call(rbind, df_ls)
    })
    do.call(rbind, data_ls)
})
# revis_climprox_df <- do.call(rbind, ssp_ls)
revis_climprox_df <- full_join(
    approaches_ls[[1]], approaches_ls[[2]],
    by = c("species", "net.id", "ssp")
)
colnames(revis_climprox_df) <- c("species", "Tair_GBIF", "Qsoil_GBIF", "ssp", "net.id", "Tair_INat", "Qsoil_INat")
revis_climprox_df <- na.omit(revis_climprox_df)
revis_climprox_df[, c(2:3, 6:7)] <- abs(revis_climprox_df[, c(2:3, 6:7)])
head(revis_climprox_df)

revis_climprox_df$group <- "Animals"
revis_climprox_df$group[revis_climprox_df$species %in% rownames(plants_gowdis)] <- "Plants"

Exts_df <- do.call(
    rbind,
    apply(revis_climprox_df, 1, FUN = function(x) {
        # x <- revis_climprox_df[1,]
        # print(x)

        if (tail(x, 1) == "Plants") {
            GBIF_ext <- any(abs(as.numeric(x[2:3])) > 2)
            INat_ext <- any(abs(as.numeric(x[6:7])) > 2)
        } else {
            GBIF_ext <- abs(as.numeric(x[2])) > 2
            INat_ext <- abs(as.numeric(x[6])) > 2
        }
        data.frame(
            GBIF_ext = GBIF_ext,
            INat_ext = INat_ext
        )
    })
)
revis_climprox_df <- cbind(revis_climprox_df, Exts_df)

### Correlation of Margins ----------
TairSsp245 <- cor.test(
    revis_climprox_df$Tair_GBIF[revis_climprox_df$ssp == "ssp245"],
    revis_climprox_df$Tair_INat[revis_climprox_df$ssp == "ssp245"]
)
TairSsp585 <- cor.test(
    revis_climprox_df$Tair_GBIF[revis_climprox_df$ssp == "ssp585"],
    revis_climprox_df$Tair_INat[revis_climprox_df$ssp == "ssp585"]
)
QsoilSsp245 <- cor.test(
    revis_climprox_df$Qsoil_GBIF[revis_climprox_df$ssp == "ssp245" & revis_climprox_df$species %in% rownames(plants_gowdis)],
    revis_climprox_df$Qsoil_INat[revis_climprox_df$ssp == "ssp245" & revis_climprox_df$species %in% rownames(plants_gowdis)]
)
QsoilSsp585 <- cor.test(
    revis_climprox_df$Qsoil_GBIF[revis_climprox_df$ssp == "ssp585" & revis_climprox_df$species %in% rownames(plants_gowdis)],
    revis_climprox_df$Qsoil_INat[revis_climprox_df$ssp == "ssp585" & revis_climprox_df$species %in% rownames(plants_gowdis)]
)

labels_df <- data.frame(
    label =
        c(
            TairSsp245$estimate, TairSsp585$estimate,
            QsoilSsp245$estimate, QsoilSsp585$estimate
        ),
    X = c(30, 30, 7.5, 7.5),
    Y = c(2.5, 2.5, 1, 1),
    Margin = c("Tair", "Tair", "Qsoil", "Qsoil"),
    ssp = c("ssp245", "ssp585", "ssp245", "ssp585")
)

Tair_gg <- ggplot(revis_climprox_df, aes(x = Tair_GBIF, y = Tair_INat)) +
    geom_point() +
    stat_smooth(method = "lm") +
    geom_vline(xintercept = 2, linetype = "dashed") +
    geom_hline(yintercept = 2, linetype = "dashed") +
    geom_label(
        data = labels_df[labels_df$Margin == "Tair", ],
        aes(label = paste("Correlation =", round(label, 3)), x = X, y = Y),
        size = 3, fill = "white"
    ) +
    geom_text(
        data = subset(
            revis_climprox_df,
            Tair_GBIF > 30
        ),
        aes(label = species),
        hjust = 1.1,
        size = 3
    ) +
    theme_bw() +
    facet_wrap(~ssp) +
    labs(y = "INaturalist", x = "GBIF")

Qsoil_gg <- ggplot(
    revis_climprox_df[revis_climprox_df$species %in% rownames(plants_gowdis), ],
    aes(x = Qsoil_GBIF, y = Qsoil_INat)
) +
    geom_point() +
    stat_smooth(method = "lm") +
    geom_vline(xintercept = 2, linetype = "dashed") +
    geom_hline(yintercept = 2, linetype = "dashed") +
    geom_label(
        data = labels_df[labels_df$Margin == "Qsoil", ],
        aes(label = paste("Correlation =", round(label, 3)), x = X, y = Y),
        size = 3, fill = "white"
    ) +
    ggrepel::geom_text_repel(
        data = subset(
            revis_climprox_df,
            species %in% rownames(plants_gowdis) & Qsoil_GBIF > 10
        ),
        aes(label = species),
        hjust = 1,
        nudge_x = -0.1,
        size = 3,
        box.padding = 0.5, # space around text
        point.padding = 0.3, # space between text and point
        max.overlaps = Inf # ensures all labels attempt to appear
    ) +
    theme_bw() +
    facet_wrap(~ssp) +
    labs(y = "INaturalist", x = "GBIF")

### Congruency of Primary Extinctions ----------
Congruency_df <- do.call(
    rbind,
    lapply(ssps, FUN = function(ssp) {
        Iter_df <- revis_climprox_df[revis_climprox_df$ssp == ssp, ]

        plants_total <- sum(Iter_df$species %in% rownames(plants_gowdis))
        Ext_plants_g <- sum(Iter_df$GBIF_ext[Iter_df$species %in% rownames(plants_gowdis)])
        Ext_plants_i <- sum(Iter_df$INat_ext[Iter_df$species %in% rownames(plants_gowdis)])
        Ext_plants_both <- sum(Iter_df$GBIF_ext[Iter_df$species %in% rownames(plants_gowdis)] + Iter_df$INat_ext[Iter_df$species %in% rownames(plants_gowdis)] == 2)

        animals_total <- sum(Iter_df$species %in% rownames(animals_gowdis))
        Ext_animals_g <- sum(Iter_df$GBIF_ext[Iter_df$species %in% rownames(animals_gowdis)])
        Ext_animals_i <- sum(Iter_df$INat_ext[Iter_df$species %in% rownames(animals_gowdis)])
        Ext_animals_both <- sum(Iter_df$GBIF_ext[Iter_df$species %in% rownames(animals_gowdis)] + Iter_df$INat_ext[Iter_df$species %in% rownames(animals_gowdis)] == 2)

        data.frame(
            Values = c(
                plants_total, Ext_plants_g, Ext_plants_i, Ext_plants_both,
                animals_total, Ext_animals_g, Ext_animals_i, Ext_animals_both
            ),
            Groups = rep(c("Plants", "Animals"), each = 4),
            Approach = rep(c("Number of Species", "GBIF", "INaturalist", "Shared"), 2),
            ssp = ssp
        )
    })
)

venn_counts <- Congruency_df %>%
    filter(Approach %in% c("GBIF", "INaturalist", "Shared")) %>%
    tidyr::pivot_wider(
        names_from = Approach,
        values_from = Values
    )

groups <- unique(venn_counts$Groups)
ssps <- unique(venn_counts$ssp)

venn_ls <- lapply(seq_along(ssps), FUN = function(i) {
    ret_ls <- lapply(seq_along(groups), FUN = function(j) {
        row <- venn_counts %>%
            filter(Groups == groups[j], ssp == ssps[i])

        venn_gg <- draw.pairwise.venn(
            area1 = row$GBIF,
            area2 = row$INaturalist,
            cross.area = row$Shared,
            category = c("GBIF", "iNat"),
            fill = c("#1f78b4", "#33a02c"),
            alpha = 0.75,
            ext.text = FALSE,
            scaled = TRUE,
            cex = 0.9, # shrink counts
            cat.cex = 1,
            ind = FALSE
        )
        ggdraw(venn_gg)
    })
    names(ret_ls) <- groups
    ret_ls
})
names(venn_ls) <- ssps
Cong_gg <- cowplot::plot_grid(
    cowplot::plot_grid(
        label_row("     I - Animals"),
        venn_ls[["ssp245"]]$Animals,
        label_row("     II - Plants"),
        venn_ls[["ssp245"]]$Plants,
        rel_heights = c(0.1, 1, 0.1, 1),
        ncol = 1
    ),
    cowplot::plot_grid(
        label_row("     I - Animals"),
        venn_ls[["ssp585"]]$Animals,
        label_row("     II - Plants"),
        venn_ls[["ssp585"]]$Plants,
        rel_heights = c(0.1, 1, 0.1, 1),
        ncol = 1
    ),
    ncol = 2
)

### Saving Plot ----------
p <- cowplot::plot_grid(
    label_row("(A) Air Temperature"),
    Tair_gg,
    label_row("(B) Soil Moisture"),
    Qsoil_gg,
    label_row("(C) Primary Extinction Congruency"),
    Cong_gg,
    ncol = 1, rel_heights = c(0.1, 1, 0.1, 1, 0.1, 1.5)
)
ggsave(p, file = file.path(Dir.Exports, "Revision_ClimateSafetyMargins.png"), width = 20 / 1.8, height = 22 / 1.8)

# PHYLOGENETIC SIGNAL ON TRAITS ============================================
# compare gowdis (which is already calculated) to phylogenetic distance of animals and plants respectively
# install.packages("rotl")
library(ape)
library(rotl)

phylo_trait_corr <- lapply(list(
    animals_gowdis,
    plants_gowdis
), FUN = function(dist_obj) {
    # print(head(dist_obj))
    # dist_obj <- plants_gowdis

    species <- rownames(dist_obj)
    # Match species in Open Tree of Life
    matched <- tnrs_match_names(names = species)
    # Keep only species with valid OTT IDs
    valid_ids <- matched$ott_id[!is.na(matched$ott_id)]
    tree <- tol_induced_subtree(
        ott_ids = valid_ids[valid_ids %nin% c(5524967, 3907331)] # these give errors
    )
    # plot(tree)
    # tree

    # Extract tip labels from tree
    tree_tips <- tree$tip.label

    # Remove OTT IDs for easier matching
    tip_names <- gsub(pattern = "_", replacement = " ", sub("_ott.*$", "", tree_tips))

    # Match your species
    my_species <- rownames(dist_obj)
    present <- my_species[my_species %in% tip_names]
    absent <- my_species[!my_species %in% tip_names]

    # cat("Present species:\n")
    # print(length(present))
    # cat("Absent species:\n")
    # print(length(absent))

    # Find tip labels in the tree corresponding to present species
    present_tips <- tree_tips[tip_names %in% present]

    # Prune tree to only these tips
    pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, present_tips))

    # calculate distance
    phylo_dist <- cophenetic(pruned_tree) # tree has no branch length, so only topological distances
    # phylo_dist

    colnames(phylo_dist) <- gsub(pattern = "_", replacement = " ", sub("_ott.*$", "", colnames(phylo_dist)))
    rownames(phylo_dist) <- gsub(pattern = "_", replacement = " ", sub("_ott.*$", "", rownames(phylo_dist)))

    # dim(phylo_dist)
    # dim(dist_obj)

    # Get species names from each matrix
    phylo_species <- rownames(phylo_dist)
    trait_species <- rownames(dist_obj)

    # Identify shared species
    shared_species <- intersect(phylo_species, trait_species)
    # length(shared_species) # optional: see how many overlap

    # Subset phylogenetic distance matrix
    phylo_dist_sub <- phylo_dist[shared_species, shared_species]

    # Subset trait distance matrix
    trait_dist_sub <- dist_obj[shared_species, shared_species]

    # Extract lower triangle values
    phylo_vec <- phylo_dist_sub[lower.tri(phylo_dist_sub)]
    trait_vec <- trait_dist_sub[lower.tri(trait_dist_sub)]

    # plot(phylo_vec, trait_vec)

    # Quick correlation
    cortest <- cor.test(phylo_vec, trait_vec, method = "pearson")
    print(cortest)
    return(cortest)
})
names(phylo_trait_corr) <- c("Animals", "Plants")
phylo_trait_corr
