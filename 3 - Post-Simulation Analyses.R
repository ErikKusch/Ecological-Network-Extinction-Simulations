#' ####################################################################### #
#' PROJECT: [Biodiversity Simplification & Ecological Network Topology] 
#' CONTENTS: 
#'  - Post-Simulation Analyses
#'  DEPENDENCIES:
#'  - "0 - Preamble.R"
#'  - "2 - Extinction Simulation" has to have been run
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE =================================================================
rm(list=ls())
set.seed(42)

## Sourcing ----------------------------------------------------------------
source("0 - Preamble.R")

# DATA LOADING & MANIPULATING ==============================================
message("### DATA PREPARATION ###")

## Data Loading ------------------------------------------------------------
source("2 - Extinction Simulations.R")
PlotTopo_ls <- list(ALL = PlotTopoAll_ls,
                    Animals = PlotTopoAnimals_ls,
                    Plants = PlotTopoPlants_ls
)
TopoPlots <- c("n_species", "n_links", "Modularity", "Nestedness") # "n_species", "n_animals", "n_plants", "Modularity", "n_links", "Nestedness"

# LINEAR MODELS ============================================================
for(model_iter in TopoPlots){
  model_df <- PlotTopoAll_ls$Change
  sink(file = file.path(Dir.Exports, paste0("MODEL_", model_iter, ".txt")))
  print(
    summary(lm(RelChange ~ 0 + (IS+RE)*Proxy, data = model_df[model_df$Topology == model_iter, ]))
  )
  sink()
}

# VISUALISATION ============================================================
message("### RESULT VISUALISATION ###")
pal_lm <- c("#016392", "#E19825", "#3E8853")

## Venn-Diagram of Proxy Agreement -----------------------------------------
print("Extinction Proxy Overlap ---")
for(RunName in c("ALL", "Plants", "Animals")
){
  ExtSpecies_ls <- FUN_SimComp(PlantAnim = NULL, RunName = RunName, IS = 1, Rewiring = 1, CutOffs = CutOffs)
  
  # ## number of primary extinctions per proxy
  # Clim_ls <- lapply(lapply(ExtSpecies_ls, "[[", "Climate"), "[[", "Removed")
  # IUCN_ls <- lapply(lapply(ExtSpecies_ls, "[[", "IUCN"), "[[", "Removed")
  # Centr_ls <- lapply(lapply(ExtSpecies_ls, "[[", "Strength"), "[[", "Removed")
  # ClimIUCN_ls <- lapply(1:length(Clim_ls), function(x){sum(Clim_ls[[x]] %in% IUCN_ls[[x]])})
  # ClimCentr_ls <- lapply(1:length(Clim_ls), function(x){sum(Clim_ls[[x]] %in% Centr_ls[[x]])})
  # IUCNCentr_ls <- lapply(1:length(Clim_ls), function(x){sum(IUCN_ls[[x]] %in% Centr_ls[[x]])})
  # All_ls <- lapply(1:length(Clim_ls), function(x){
  #   sum(
  #     IUCN_ls[[x]][IUCN_ls[[x]] %in% Centr_ls[[x]]] %in% Clim_ls[[x]]
  #     )
  #   })
  # ## absolute numbers of primary extinctions per network
  # PrimaryExt_df <- data.frame(
  #   Climate = unlist(lapply(Clim_ls, length)),
  #   IUCN = unlist(lapply(IUCN_ls, length)),
  #   Centrality = unlist(lapply(Centr_ls, length)),
  #   `Climate+IUCN` = unlist(ClimIUCN_ls),
  #   `Climate+Centrality` = unlist(ClimCentr_ls),
  #   `IUCN+Centrality` = unlist(IUCNCentr_ls),
  #   All = unlist(All_ls)
  # )
  # ## relative numbers of primary extinctions per network
  # PrimaryExt_df <- apply(X = PrimaryExt_df, MARGIN = 2, function(x){
  #   x/PreExt_df$n_species
  # })
  # 
  # ## summary statistics
  # PrimExt_df <- data.frame(min = apply(PrimaryExt_df, 2, min),
  #                          mean = apply(PrimaryExt_df, 2, mean),
  #                          max = apply(PrimaryExt_df, 2, max),
  #                          sd = apply(PrimaryExt_df, 2, sd))
  # print(PrimExt_df)
  
  ## Total Venn-Diagram
  Venn_ls <- list(Climate = unlist(lapply(lapply(ExtSpecies_ls, "[[", "Climate"), "[[", "Removed")),
                  IUCN = unlist(lapply(lapply(ExtSpecies_ls, "[[", "IUCN"), "[[", "Removed")),
                  Centrality = unlist(lapply(lapply(ExtSpecies_ls, "[[", "Strength"), "[[", "Removed")))
  ggvenn(Venn_ls, fill_color = pal_lm, fill_alpha = 0.8, text_color = "white")
  ggsave(filename = file.path(Dir.Exports, paste0("PLOT_Proxy", RunName,".png")), width = 4, height = 3, units = "cm", scale = 7, dpi = 1e3)
}

## Grid-Based Visualisations -----------------------------------------------
print("Grid-based visualisations ---")

### Effects relative to pre-extinction network topology ----
RunName = "ALL"
# for(RunName in names(PlotTopo_ls)){

  Proxy_ls <- as.list(rep(NA, length(unique(Change_df$Proxy))))
  names(Proxy_ls) <- unique(Change_df$Proxy)

  ### Relative Changes ----
  Change_df <- PlotTopoAll_ls$Change
# Change_df <- PlotTopo_ls[[RunName]]$Change
  
  for(proxy_iter in unique(Change_df$Proxy)){
    iter_df <- Change_df[Change_df$Topology %in% TopoPlots & Change_df$Proxy %in% proxy_iter, ]
    plot_df <- aggregate(RelChange ~ Proxy+Topology+IS+RE, FUN = mean, data = iter_df[iter_df$RelChange != 0, ])
    sd_df <- aggregate(RelChange ~ Proxy+Topology+IS+RE, FUN = sd, data = iter_df[iter_df$RelChange != 0, ])
    
    # iterate my Topology, then fuse plots
    matplot_ls <- list(Pred = as.list(rep(NA, length(TopoPlots))), SD = as.list(rep(NA, length(TopoPlots))))
    
    for(TopoIter in 1:length(TopoPlots)){
      plot_df2 <- plot_df[plot_df$Topology == TopoPlots[TopoIter], ]
      sd_df2 <- sd_df[sd_df$Topology == TopoPlots[TopoIter], ]
      
      matplot_ls$Pred[[TopoIter]] <- ggplot(plot_df2, aes(x = RE, y = IS)) +
        geom_tile(aes(fill = RelChange)) +
        coord_fixed() + 
        facet_wrap(~Proxy+Topology, ncol = 1) + 
        theme_bw() + 
        # xlab("Probability of Rewiring Required to Realise Novel Links") + 
        # ylab("Proportion of Initial Interaction Strength Required for Continued Existence") + 
        xlab("") + 
        ylab("") + 
        guides(fill = guide_colourbar(barwidth = 2,
                                      barheight = 20,
                                      # title = "Proportional Network Metric Change",
                                      title = "",
                                      # title.position = "bottom",
                                      # direction = "horizontal",
                                      # legend.location = "bottom"
                                      )) + 
        scale_fill_viridis("Mean", option = "B", direction = -1) + 
        # theme(legend.position = "bottom") + 
        # ggtitle(TopoPlots[TopoIter]) + 
        theme(plot.margin = unit(c(0,0,0,0), "lines"))
      
      matplot_ls$SD[[TopoIter]] <- ggplot(sd_df2, aes(x = RE, y = IS)) +
        geom_tile(aes(fill = RelChange)) +
        coord_fixed() +
        facet_wrap(~Proxy+Topology, ncol = 1) + 
        theme_bw() + 
        # xlab("Probability of Rewiring Required to Realise Novel Links") + 
        # ylab("Proportion of Initial Interaction Strength Required for Continued Existence") + 
        xlab("") + 
        ylab("") + 
        guides(fill = guide_colourbar(barwidth = 2,
                                      barheight = 20,
                                      # title = "Proportional Network Metric Change",
                                      title = "",
                                      # title.position = "bottom",
                                      # direction = "horizontal",
                                      # legend.location = "bottom"
                                      )) + 
        scale_fill_viridis("SD", option = "D") + 
        # theme(legend.position = "bottom") + 
        # ggtitle(TopoPlots[TopoIter]) + 
        theme(plot.margin = unit(c(0,0,0,0), "lines"))
    }
    y.grob <- textGrob("Proportion of Initial Interaction Strength Required for Continued Existence", 
                       gp=gpar(fontface="bold", col="black", fontsize=15), rot=90)
    x.grob <- textGrob("Probability of Rewiring Required to Realise Novel Links", 
                       gp=gpar(fontface="bold", col="black", fontsize=15))
    
    MatPred_plot <- plot_grid(plotlist = matplot_ls$Pred, nrow = 2)
    MatPred_plot <- grid.arrange(arrangeGrob(MatPred_plot, left = y.grob, top = x.grob))
    ggsave(MatPred_plot, filename = file.path(Dir.Exports, paste0("PLOT_MatrixChange-", RunName, "-", proxy_iter, ".png")), width = 34/1.2, height = 32/1.2, units = "cm")
    
    MatSD_plot <- plot_grid(plotlist = matplot_ls$SD, nrow = 2)
    MatSD_plot <- grid.arrange(arrangeGrob(MatSD_plot, left = y.grob, top = x.grob))
    ggsave(MatSD_plot, filename = file.path(Dir.Exports, paste0("PLOT_MatrixChange_SD-", RunName, "-", proxy_iter, ".png")), width = 34/1.2, height = 32/1.2, units = "cm") 
    
    Proxy_ls[[proxy_iter]] <- matplot_ls
  }
  
  RelChangePlots_ls <- Proxy_ls
  
  
  
  
# }

### Effect Sizes ----
RunName = "ALL"

## extract necessary data
EffectSize_df <- PlotTopoAll_ls$EffectSize
RandomSD_df <- PlotTopoAll_ls$RandomSD

## reformat to effect sizes
# Sig_df <- !(abs(EffectSize_df[,TopoPlots]) > RandomSD_df[,TopoPlots])
# Plot_df <- EffectSize_df[, TopoPlots]
# Plot_df[Sig_df] <- NA
Plot_df <- EffectSize_df[, TopoPlots] / RandomSD_df[, TopoPlots]
Plot_df$netID <- EffectSize_df$netID
Plot_df$RE <- EffectSize_df$RE
Plot_df$IS <- EffectSize_df$IS
Plot_df$Pry <- EffectSize_df$Pry

# ## make effectsizes relative to base networks
# Plot_df <- Plot_df[ , c("netID", TopoPlots, "Pry", "IS", "RE")]
# merged_df <- base::merge(Plot_df, PreExt_df[, c("netID", TopoPlots)], by = "netID")
# Rel_effSizes <- merged_df[ , 2:(length(TopoPlots)+1)] /  # predictions
#   merged_df[ , (1+2*length(TopoPlots)):ncol(merged_df)] # pre-extinctions
# colnames(Rel_effSizes) <- TopoPlots
# Rel_effSizes$netID <- Plot_df$netID
# Rel_effSizes$RE <- Plot_df$RE
# Rel_effSizes$IS <- Plot_df$IS
# Rel_effSizes$Pry <- Plot_df$Pry
# Plot_df <- Rel_effSizes

## reshape data for aggregation and formatting
Plot_df <- stats::reshape(data = Plot_df, 
                          times = colnames(Plot_df)[1:length(TopoPlots)],
                          varying = list(colnames(Plot_df)[1:length(TopoPlots)]),
                          timevar = "Topology",
                          direction = "long")
colnames(Plot_df)[ncol(Plot_df)-1] <- "EffectSize"
Plot_df <- Plot_df[which(abs(Plot_df$EffectSize) != Inf), ]
Plot_df <- Plot_df[which(!is.na(Plot_df$EffectSize)), ]

# ## removal of effectsizes where there was no deletion sequence
# ExtSpecies_ls <- FUN_SimComp(PlantAnim = NULL, RunName = RunName, IS = 1, Rewiring = 1, CutOffs = CutOffs)
# Plot_df[Plot_df$netID %in% names(which(unlist(lapply(lapply(lapply(ExtSpecies_ls, "[[", "Climate"), "[[", "Removed"), length)) == 0)) &
#           Plot_df$Pry == "Climate", "EffectSize"] <- NA
# Plot_df[Plot_df$netID %in% names(which(unlist(lapply(lapply(lapply(ExtSpecies_ls, "[[", "IUCN"), "[[", "Removed"), length)) == 0)) &
#           Plot_df$Pry == "IUCN", "EffectSize"] <- NA
# Plot_df[Plot_df$netID %in% names(which(unlist(lapply(lapply(lapply(ExtSpecies_ls, "[[", "Strength"), "[[", "Removed"), length)) == 0)) &
#           Plot_df$Pry == "Strength", "EffectSize"] <- NA
# Plot_df <- na.omit(Plot_df)         
# Plot_df <- Plot_df[, -ncol(Plot_df)]

## Plot Creation
Proxy_ls <- as.list(rep(NA, length(unique(Plot_df$Pry))))
names(Proxy_ls) <- unique(Plot_df$Pry)
for(proxy_iter in unique(Plot_df$Pry)){
  Topo_ls <- as.list(rep(NA, length(TopoPlots)))
  names(Topo_ls) <- TopoPlots
  for(topo_iter in TopoPlots){
    iter_df <- Plot_df[Plot_df$Topology == topo_iter &
                         Plot_df$Pry == proxy_iter, ]
    mean_df <- aggregate(EffectSize ~ Pry+Topology+IS+RE, FUN = mean, data = iter_df)
    Topo_ls[[topo_iter]] <- ggplot(mean_df, aes(x = RE, y = IS)) +
      geom_tile(aes(fill = EffectSize)) +
      coord_fixed() +
      facet_wrap(~Pry+Topology) +
      scale_fill_gradient2(high = "darkgreen", low = "darkred") +
      theme_bw() +  xlab("") + ylab("") +
      guides(fill = guide_colourbar(barwidth = 2, barheight = 20, title = ""))
  }
  Proxy_ls[[proxy_iter]] <- Topo_ls
  
  ## plotting and saving
  y.grob <- textGrob("Proportion of Initial Interaction Strength Required for Continued Existence",
                     gp=gpar(fontface="bold", col="black", fontsize=15), rot=90)
  x.grob <- textGrob("Probability of Rewiring Required to Realise Novel Links",
                     gp=gpar(fontface="bold", col="black", fontsize=15))
  MatPred_plot <- plot_grid(plotlist = Proxy_ls[[proxy_iter]], nrow = 2)
  MatPred_plot <- grid.arrange(arrangeGrob(MatPred_plot, left = y.grob, top = x.grob))
  ggsave(MatPred_plot, filename = file.path(Dir.Exports, paste0("PLOT_EffectSizes-", RunName, "-", proxy_iter, ".png")), width = 34/1.2, height = 30/1.2, units = "cm")
  
}
EffSizePlots_ls <- Proxy_ls



## as exptected, the below produces an error. Signficiance how?!
# lm_df <- iter_df[iter_df$IS == 0.5 & iter_df$RE == 0.5, ]
# lm_df
# lme4::lmer(EffectSize ~ 1 + (1|netID), data = lm_df)





# 
# 
# 
# # for(Pry_iter in unique(Plot_df$Pry)){
# #   iter_df <- mean_df[mean_df$Pry == Pry_iter,]
# #   itersd_df <- sd_df[sd_df$Pry == Pry_iter,]
# #   iter_num_df <- numsig_df[numsig_df$Pry == Pry_iter,]
# #   
# #   effectsizeplot <- ggplot(iter_df, aes(x = RE, y = IS)) +
# #     geom_tile(aes(fill = EffectSize)) +
# #     coord_fixed() + 
# #     facet_wrap(~factor(Topology, levels = unique(Plot_df$Topology))) + 
# #     scale_fill_gradient2(high = "darkgreen", low = "darkred") +
# #     theme_bw() +  xlab("") + ylab("") + 
# #     guides(fill = guide_colourbar(barwidth = 2, barheight = 15, title = "Mean"))
# #   
# #   effsdplot <- ggplot(itersd_df, aes(x = RE, y = IS)) +
# #     geom_tile(aes(fill = EffectSize)) +
# #     coord_fixed() + 
# #     facet_wrap(~factor(Topology, levels = unique(Plot_df$Topology))) + 
# #     scale_fill_viridis_c(option = "B", direction = -1) + 
# #     theme_bw() + xlab("") + ylab("") + 
# #     guides(fill = guide_colourbar(barwidth = 2, barheight = 15, title = "SD"))
# #   
# #   numsigplot <- ggplot(iter_num_df, aes(x = RE, y = IS)) +
# #     geom_tile(aes(fill = EffectSize)) +
# #     coord_fixed() + 
# #     facet_wrap(~factor(Topology, levels = unique(Plot_df$Topology))) + 
# #     scale_fill_viridis_c(option = "E", direction = -1) + 
# #     theme_bw() + xlab("") + ylab("") + 
# #     guides(fill = guide_colourbar(barwidth = 2, barheight = 15, title = "#"))
# #   
# #   y.grob <- textGrob("Proportion of Initial Interaction Strength Required for Continued Existence", 
# #                      gp=gpar(fontface="bold", col="black", fontsize=15), rot=90)
# #   x.grob <- textGrob("Probability of Rewiring Required to Realise Novel Links", 
# #                      gp=gpar(fontface="bold", col="black", fontsize=15))
# #   
# #   MatPred_plot <- plot_grid(plotlist = list(effectsizeplot, effsdplot, numsigplot), nrow = 3)
# #   MatPred_plot <- grid.arrange(arrangeGrob(MatPred_plot, left = y.grob, top = x.grob))
# #   ggsave(MatPred_plot, filename = file.path(Dir.Exports, paste0("PLOT_EffectSizes-", RunName, "-", Pry_iter, ".png")), width = 34/1.2, height = 36/1.2, units = "cm")
# # }
# 
# ### Effects of bottom-up and top-down compared to ALL and each other ----
# CompCasc_ls <- list(Plants = NA,
#                     Animals = NA)
# for(RunName in c("Plants", "Animals")){
#   ### Baseline ALL Changes ----
#   Change_df <- PlotTopo_ls[["ALL"]]$Change
#   TopoPlots <- c("n_species", "n_links", "Modularity", "Nestedness") # "n_species", "n_animals", "n_plants", "Modularity", "n_links", "Nestedness"
#   Change_df <- Change_df[Change_df$Topology %in% TopoPlots & Change_df$Proxy %in% c("Climate", "IUCN", "Strength"), ]
#   base_plot_df <- aggregate(RelChange ~ Proxy+Topology+IS+RE, FUN = mean, data = Change_df[Change_df$RelChange != 0, ])
#   base_sd_df <- aggregate(RelChange ~ Proxy+Topology+IS+RE, FUN = sd, data = Change_df[Change_df$RelChange != 0, ])
#   
#   ### Relative Changes ----
#   Change_df <- PlotTopo_ls[[RunName]]$Change
#   TopoPlots <- c("n_species", "n_links", "Modularity", "Nestedness") # "n_species", "n_animals", "n_plants", "Modularity", "n_links", "Nestedness"
#   Change_df <- Change_df[Change_df$Topology %in% TopoPlots & Change_df$Proxy %in% c("Climate", "IUCN", "Strength"), ]
#   test_plot_df <- aggregate(RelChange ~ Proxy+Topology+IS+RE, FUN = mean, data = Change_df[Change_df$RelChange != 0, ])
#   test_sd_df <- aggregate(RelChange ~ Proxy+Topology+IS+RE, FUN = sd, data = Change_df[Change_df$RelChange != 0, ])
#   
#   CompCasc_ls[[RunName]] <- list(mean = test_plot_df,
#                                  sd = test_sd_df)
#   
#   ### Changes compared to baseline ----
#   base_combs <- with(base_plot_df, paste(Proxy, Topology, IS, RE, sep="-"))
#   test_combs <- with(test_plot_df, paste(Proxy, Topology, IS, RE, sep="-"))
#   
#   plot_df <- base_plot_df[base_combs %in% test_combs, ]
#   plot_df$RelChange <- base_plot_df[base_combs %in% test_combs, 5] - test_plot_df[test_combs %in% base_combs, 5]
#   
#   sd_df <- base_sd_df[base_combs %in% test_combs, ]
#   sd_df$RelChange <- base_sd_df[base_combs %in% test_combs, 5] - test_sd_df[test_combs %in% base_combs, 5]
#   
#   # iterate my Topology, then fuse plots
#   matplot_ls <- list(Pred = as.list(rep(NA, length(TopoPlots))),
#                      SD = as.list(rep(NA, length(TopoPlots))))
#   
#   
#   for(TopoIter in 1:length(TopoPlots)){
#     plot_df2 <- plot_df[plot_df$Topology == TopoPlots[TopoIter], ]
#     sd_df2 <- sd_df[sd_df$Topology == TopoPlots[TopoIter], ]
#     
#     matplot_ls$Pred[[TopoIter]] <- ggplot(plot_df2, aes(x = RE, y = IS)) +
#       geom_tile(aes(fill = RelChange)) +
#       coord_fixed() + 
#       facet_wrap(~Proxy, ncol = 1) + 
#       theme_bw() + 
#       # xlab("Probability of Rewiring Required to Realise Novel Links") + 
#       # ylab("Proportion of Initial Interaction Strength Required for Continued Existence") + 
#       xlab("") + 
#       ylab("") + 
#       guides(fill = guide_colourbar(barwidth = 15,
#                                     barheight = 1.5,
#                                     title = "Proportional Network Metric Change",
#                                     title.position = "bottom",
#                                     # direction = "horizontal",
#                                     legend.location = "bottom")) + 
#       scale_fill_viridis("Mean", option = "B", direction = -1) + 
#       theme(legend.position = "bottom") + ggtitle(TopoPlots[TopoIter]) + 
#       theme(plot.margin = unit(c(0,0,0,0), "lines"))
#     
#     matplot_ls$SD[[TopoIter]] <- ggplot(sd_df2, aes(x = RE, y = IS)) +
#       geom_tile(aes(fill = RelChange)) +
#       coord_fixed() +
#       facet_wrap(~Proxy, ncol = 1) + 
#       theme_bw() + 
#       # xlab("Probability of Rewiring Required to Realise Novel Links") + 
#       # ylab("Proportion of Initial Interaction Strength Required for Continued Existence") + 
#       xlab("") + 
#       ylab("") + 
#       guides(fill = guide_colourbar(barwidth = 15,
#                                     barheight = 1.5,
#                                     title = "Proportional Network Metric Change",
#                                     title.position = "bottom",
#                                     # direction = "horizontal",
#                                     legend.location = "bottom")) + 
#       scale_fill_viridis("SD", option = "D") + 
#       theme(legend.position = "bottom") + ggtitle(TopoPlots[TopoIter]) + 
#       theme(plot.margin = unit(c(0,0,0,0), "lines"))
#   }
#   y.grob <- textGrob("Proportion of Initial Interaction Strength Required for Continued Existence", 
#                      gp=gpar(fontface="bold", col="black", fontsize=15), rot=90)
#   x.grob <- textGrob("Probability of Rewiring Required to Realise Novel Links", 
#                      gp=gpar(fontface="bold", col="black", fontsize=15))
#   
#   MatPred_plot <- plot_grid(plotlist = matplot_ls$Pred, nrow = 1)
#   MatPred_plot <- grid.arrange(arrangeGrob(MatPred_plot, left = y.grob, top = x.grob))
#   ggsave(MatPred_plot, filename = file.path(Dir.Exports, paste0("PLOT_MatrixChange", RunName, "-ComparedtoALL.png")), width = 40/1.2, height = 34/1.2, units = "cm")
#   
#   MatSD_plot <- plot_grid(plotlist = matplot_ls$SD, nrow = 1)
#   MatSD_plot <- grid.arrange(arrangeGrob(MatSD_plot, left = y.grob, top = x.grob))
#   ggsave(MatSD_plot, filename = file.path(Dir.Exports, paste0("PLOT_MatrixChange_SD", RunName, "-ComparedtoALL.png")), width = 40/1.2, height = 34/1.2, units = "cm") 
# }
# 
# plants_combs <- with(CompCasc_ls$Plants$mean, paste(Proxy, Topology, IS, RE, sep="-"))
# animals_combs <- with(CompCasc_ls$Animals$mean, paste(Proxy, Topology, IS, RE, sep="-"))
# 
# plot_df <- CompCasc_ls$Animals$mean[base_combs %in% test_combs, ]
# plot_df$RelChange <- CompCasc_ls$Plants$mean[plants_combs %in% animals_combs, 5] - 
#   CompCasc_ls$Animals$mean[animals_combs %in% plants_combs, 5]
# sd_df <- CompCasc_ls$Animals$sd[base_combs %in% test_combs, ]
# sd_df$RelChange <- CompCasc_ls$Plants$sd[plants_combs %in% animals_combs, 5] - 
#   CompCasc_ls$Animals$sd[animals_combs %in% plants_combs, 5]
# 
# # iterate my Topology, then fuse plots
# matplot_ls <- list(Pred = as.list(rep(NA, length(TopoPlots))),
#                    SD = as.list(rep(NA, length(TopoPlots))))
# 
# 
# for(TopoIter in 1:length(TopoPlots)){
#   plot_df2 <- plot_df[plot_df$Topology == TopoPlots[TopoIter], ]
#   sd_df2 <- sd_df[sd_df$Topology == TopoPlots[TopoIter], ]
#   
#   matplot_ls$Pred[[TopoIter]] <- ggplot(plot_df2, aes(x = RE, y = IS)) +
#     geom_tile(aes(fill = RelChange)) +
#     coord_fixed() + 
#     facet_wrap(~Proxy, ncol = 1) + 
#     theme_bw() + 
#     # xlab("Probability of Rewiring Required to Realise Novel Links") + 
#     # ylab("Proportion of Initial Interaction Strength Required for Continued Existence") + 
#     xlab("") + 
#     ylab("") + 
#     guides(fill = guide_colourbar(barwidth = 15,
#                                   barheight = 1.5,
#                                   title = "Proportional Network Metric Change",
#                                   title.position = "bottom",
#                                   # direction = "horizontal",
#                                   legend.location = "bottom")) + 
#     scale_fill_viridis("Mean", option = "B", direction = -1) + 
#     theme(legend.position = "bottom") + ggtitle(TopoPlots[TopoIter]) + 
#     theme(plot.margin = unit(c(0,0,0,0), "lines"))
#   
#   matplot_ls$SD[[TopoIter]] <- ggplot(sd_df2, aes(x = RE, y = IS)) +
#     geom_tile(aes(fill = RelChange)) +
#     coord_fixed() +
#     facet_wrap(~Proxy, ncol = 1) + 
#     theme_bw() + 
#     # xlab("Probability of Rewiring Required to Realise Novel Links") + 
#     # ylab("Proportion of Initial Interaction Strength Required for Continued Existence") + 
#     xlab("") + 
#     ylab("") + 
#     guides(fill = guide_colourbar(barwidth = 15,
#                                   barheight = 1.5,
#                                   title = "Proportional Network Metric Change",
#                                   title.position = "bottom",
#                                   # direction = "horizontal",
#                                   legend.location = "bottom")) + 
#     scale_fill_viridis("SD", option = "D") + 
#     theme(legend.position = "bottom") + ggtitle(TopoPlots[TopoIter]) + 
#     theme(plot.margin = unit(c(0,0,0,0), "lines"))
# }
# y.grob <- textGrob("Proportion of Initial Interaction Strength Required for Continued Existence", 
#                    gp=gpar(fontface="bold", col="black", fontsize=15), rot=90)
# x.grob <- textGrob("Probability of Rewiring Required to Realise Novel Links", 
#                    gp=gpar(fontface="bold", col="black", fontsize=15))
# 
# MatPred_plot <- plot_grid(plotlist = matplot_ls$Pred, nrow = 1)
# MatPred_plot <- grid.arrange(arrangeGrob(MatPred_plot, left = y.grob, top = x.grob))
# ggsave(MatPred_plot, filename = file.path(Dir.Exports, "PLOT-RelChange_Plants-Animals.png"), width = 40/1.2, height = 34/1.2, units = "cm")
# 
# MatSD_plot <- plot_grid(plotlist = matplot_ls$SD, nrow = 1)
# MatSD_plot <- grid.arrange(arrangeGrob(MatSD_plot, left = y.grob, top = x.grob))
# ggsave(MatSD_plot, filename = file.path(Dir.Exports, "PLOT-RelChange_Plants-Animals_SD.png"), width = 40/1.2, height = 34/1.2, units = "cm") 
# 
# 
# 
# 
# PlotTopo_ls2 <- lapply(names(PlotTopo_ls), function(x){
#   returnme <- PlotTopo_ls[[x]]$Change
#   returnme$Cascade <- x
#   returnme
# })
# 
# LM_df <- do.call(rbind, PlotTopo_ls2)
# LM_df <- LM_df[LM_df$Proxy != "IUCN_Climate", ]
# 
# model_iter <- "n_species"
# summary(lm(RelChange ~ 0 + (IS+RE)*Proxy, data = LM_df[LM_df$Topology == model_iter, ]))
# 
# 
# 
# 
# # 
# # 
# # stop("significance of effect size?")
# # 
# # # sd_gg <- ggplot(sd_df[sd_df$Topology %in% TopoPlots,], aes(x = RE, y = IS)) +
# # #   geom_tile(aes(fill = EffectSize)) +
# # #   coord_fixed() +
# # #   facet_wrap(~Topology+Pry) + 
# # #   scale_fill_viridis("SD", option = "D")
# # # 
# # # print(plot_grid(pred_gg, sd_gg, nrow = 2))
# # 
# # ## Effects of each extinction proxy ----------------------------------------
# # 
# # stop("increase this RE selection to 1 when trying to recreate initial plots")
# # Plot_df <- PlotTopoAll_ls$Change[PlotTopoAll_ls$Change$RE == 1, ]
# # 
# # 
# # TopoPlots <- c("n_species", "n_links", "Modularity", "Nestedness")
# # ggplot(Plot_df[Plot_df$Topology %in% TopoPlots & Plot_df$Proxy != "IUCN_Climate",],
# #        aes(x = IS, y = RelChange, col = Proxy)) + 
# #   geom_point(alpha = 0.4) + 
# #   geom_smooth() + 
# #   facet_wrap(~Topology, scales = "free", ncol = 1) + 
# #   scale_color_manual(values = pal_lm) +
# #   theme_bw() + 
# #   labs(y = "Relative Change in Network Topology (Pre- vs Post-Extinction)",
# #        x = "Network Dependency")
# # 
# # ggplot(Plot_df[Plot_df$Topology %in% TopoPlots & Plot_df$Proxy != "IUCN_Climate",],
# #        aes(x = factor(IS), y = RelChange, col = Proxy)) + 
# #   geom_boxplot() + 
# #   facet_wrap(~Topology + Proxy, scales = "free", ncol = 3) + 
# #   scale_color_manual(values = pal_lm) + 
# #   theme_bw() + 
# #   labs(y = "Relative Change in Network Topology (Pre- vs Post-Extinction)",
# #        x = "Network Dependency")
# # 
# # Compare_df <- Plot_df[Plot_df$Topology %in% TopoPlots,]
# # Clim_df <- Compare_df[Compare_df$Proxy == "Climate", ]
# # IUCN_df <- Compare_df[Compare_df$Proxy == "IUCN", ]
# # 
# # Plot_df <- Clim_df[,1:5]
# # Plot_df <- cbind(Plot_df, Clim_df$RelChange - IUCN_df$RelChange)
# # colnames(Plot_df)[ncol(Plot_df)] <- "Difference"
# # 
# # ggplot(Plot_df,
# #        aes(x = factor(IS), y = Difference)) +
# #   geom_hline(yintercept = 0) +
# #   geom_boxplot() +
# #   facet_wrap(~Topology, scales = "free", ncol = 1) +
# #   scale_color_manual(values = pal_lm) +
# #   theme_bw() +
# #   labs(title = "Difference between Climate- and IUCN-driven Extinctions",
# #        y = "Climate - IUCN Difference",
# #        x = "Network Dependency")
# # 
# # 
# ## Effects of each extinction cascade --------------------------------------
# # 
# # ## Individual Cascade-Directions -------------------------------------------
# # tTest.Calc <- function(Diff_df){
# #   tTests <- expand.grid(unique(Diff_df$Proxy), unique(Diff_df$Topology), unique(Diff_df$IS))
# #   Signif <- lapply(1:nrow(tTests), function(x){
# #     vals <-Diff_df$Value[Diff_df$Proxy == tTests[x,1] & 
# #                            Diff_df$Topology == tTests[x,2] & 
# #                            Diff_df$IS == tTests[x,3]]
# #     if(sum(!is.na(vals)) > 1){
# #       ifelse(t.test( na.omit(vals)[!is.infinite(na.omit(vals))] )[["p.value"]] < 0.05, TRUE, FALSE)
# #     }else{
# #       FALSE
# #     }
# #   })
# #   tTests$Signif <- unlist(Signif)
# #   colnames(tTests) <- c("Proxy", "Topology", "IS", "Signif")
# #   
# #   Diff_df <- merge(x = Diff_df, y = tTests, by = colnames(tTests)[1:3])
# #   Diff_df
# # }
# # 
# # Diff.Calc <- function(Topo_df){
# #   Diff_df <- Topo_df[Topo_df$Simulation == "Prediction", 1:6] - PreExt_df[match(Topo_df[Topo_df$Simulation == "Prediction",]$netID, PreExt_df$netID), 1:6]
# #   Diff_df$IS <- Topo_df[Topo_df$Simulation == "Prediction",]$IS
# #   Diff_df$Proxy <- Topo_df[Topo_df$Simulation == "Prediction",]$Proxy
# #   Diff_df$netID <- unlist(lapply(strsplit(rownames(Diff_df), split = "[.]"), "[[", 2))
# #   Diff_df <- reshape(data = Diff_df,
# #                      idvar = "netID",
# #                      varying = colnames(Diff_df)[1:6],
# #                      v.name = "Value",
# #                      timevar = "Topology",
# #                      times = colnames(Diff_df)[1:6],
# #                      new.row.names = 1:(nrow(Diff_df)*6),
# #                      direction = "long"
# #   )
# #   
# #   tTest.Calc(Diff_df)
# # }
# # 
# # 
# # Plot.Run <- function(topolist, RunName){
# #   Topo_df <- topolist[[1]]
# #   Eff_df <- topolist[[2]]
# #   pal_signif <- c("#c32200", "#62c300")
# #   
# #   ## Effect size against random simulation
# #   Plot_df<- reshape(data = Eff_df,
# #                     idvar = "netID",
# #                     varying = colnames(Eff_df)[1:6],
# #                     v.name = "Value",
# #                     timevar = c("Topology"),
# #                     times = colnames(Eff_df)[1:6],
# #                     new.row.names = 1:(nrow(Eff_df)*6),
# #                     direction = "long"
# #   )
# #   colnames(Plot_df)[2] <- "Proxy"
# #   Plot_df <- tTest.Calc(Plot_df)
# #   pal_lm <- c("#016392", "#E19825", "#3E8853")
# #   Eff_gg <- ggplot(Plot_df, aes(y = Value, x = IS, fill = Proxy, col = Proxy)) + 
# #     geom_point(alpha = 0.2) + 
# #     stat_smooth(method = "lm") + 
# #     scale_fill_manual(values = pal_lm) + 
# #     scale_color_manual(values = pal_lm) + 
# #     facet_wrap(~ factor(Topology, levels = c("n_species", "n_plants", "n_animals",
# #                                              "n_links", "Nestedness", "Modularity")), 
# #                scales = "free", ncol = 3) + 
# #     theme_bw() + 
# #     labs(x = "Interaction Strength % Required For Survival", y = "Effect Size")
# #   
# #   Eff_ls <- lapply(unique(Plot_df$Proxy), FUN = function(x){
# #     ggplot(Plot_df[Plot_df$Proxy == x, ], aes(y = Value, x = factor(IS), fill = Signif)) + 
# #       geom_boxplot() +
# #       scale_fill_manual(values = pal_signif) + 
# #       facet_wrap(~factor(Topology, levels = c("n_species", "n_plants", "n_animals",
# #                                               "n_links", "Nestedness", "Modularity")),
# #                  scales = "free", ncol = 3) +
# #       theme_bw() +
# #       labs(x = "Interaction Strength % Required For Survival", y = "Effect Size", title = x)
# #   })
# #   
# #   ## Absolute effect
# #   Diff_df <- Diff.Calc(Topo_df)
# #   Diff_gg <- ggplot(Diff_df, aes(y = Value, x = IS, col = Proxy, fill = Proxy)) + 
# #     geom_point(alpha = 0.2) + 
# #     geom_smooth() + 
# #     scale_fill_manual(values = pal_lm) + 
# #     scale_color_manual(values = pal_lm) + 
# #     facet_wrap(~ factor(Topology, levels = c("n_species", "n_plants", "n_animals",
# #                                              "n_links", "Nestedness", "Modularity")), 
# #                scales = "free", ncol = 3) + 
# #     theme_bw() + 
# #     labs(x = "Interaction Strength % Required For Survival", y = "Change In Network Topology Metric")
# #   
# #   Diff_ls <- lapply(unique(Diff_df$Proxy), FUN = function(x){
# #     ggplot(Diff_df[Diff_df$Proxy == x, ], aes(y = Value, x = factor(IS), col = Signif)) +
# #       geom_boxplot() +
# #       scale_color_manual(values = pal_signif) + 
# #       facet_wrap(~ factor(Topology, levels = c("n_species", "n_plants", "n_animals",
# #                                                "n_links", "Nestedness", "Modularity")),
# #                  scales = "free", ncol = 3) +
# #       theme_bw() +
# #       labs(x = "Interaction Strength % Required For Survival", y = "Change In Network Topology Metric", title = x)
# #   })
# #   
# #   ## Relative effect
# #   Topo_df <- topolist[[1]]
# #   Topo_df <- Topo_df[Topo_df$Simulation == "Prediction", ]
# #   Plot_df <- merge(PreExt_df, Topo_df, by = c("netID"))
# #   Rel_ls <- lapply(c("n_species", "n_plants", "n_animals", "n_links", "Nestedness", "Modularity"), function(x){
# #     Change_vec <- apply(Plot_df[ , grepl(x,colnames(Plot_df))], 1, diff)
# #     Plot_df <- data.frame(
# #       Pre = Plot_df[,grep(x,colnames(Plot_df))[1]],
# #       Value = Change_vec,
# #       netID = Plot_df$netID,
# #       Proxy = Plot_df$Proxy.y,
# #       IS = Plot_df$IS
# #     )
# #     Plot_df$RelChange <- abs(Plot_df$Value)/Plot_df$Pre
# #     
# #     Plot_df2 <- aggregate(RelChange ~ IS + Pre + Proxy, FUN = mean, data = Plot_df)
# #     ggplot(Plot_df2, aes(x = Pre , y = IS, fill = RelChange)) + 
# #       geom_tile(aes(fill = RelChange)) +
# #       # coord_fixed() + 
# #       facet_wrap(~factor(Proxy, levels = c("Strength", "Climate", "IUCN")),
# #                  scales = "free", nrow = 3) +
# #       scale_fill_viridis() + 
# #       theme_bw() + labs(x = paste(x, "(pre-extinction"), y = "Interaction Strength % Required For Survival")
# #     
# #     ggplot(Plot_df, aes(x = factor(IS), y = RelChange)) + 
# #       geom_boxplot() + 
# #       facet_wrap(~factor(Proxy, levels = c("Strength", "Climate", "IUCN")),
# #                  scales = "free", nrow = 3) +
# #       theme_bw() + labs(y = paste("Relative change in", x, "compared to pre-extinction level"), 
# #                         x = "Interaction Strength % Required For Survival",
# #                         title = x)
# #   })
# #   
# #   ## Saving plots
# #   pdf(file.path(Dir.Exports, paste0(RunName, "_Print.pdf")), paper = "a4r", height = 9, width = 12)
# #   print(Eff_gg)
# #   print(Eff_ls)
# #   print(Diff_gg)
# #   print(Diff_ls)
# #   print(Rel_ls)
# #   dev.off()
# # }
# # 
# # Plot.Run(PlotTopoAll_ls, "ALL")
# # Plot.Run(PlotTopoAnimals_ls, "Animals")
# # Plot.Run(PlotTopoPlants_ls, "Plants")
# # 
# # ## Comparison of Cascade-Directions ----------------------------------------
# # DiffAll <- Diff.Calc(PlotTopoAll_ls[[1]])
# # DiffAnim <- Diff.Calc(PlotTopoAnimals_ls[[1]])
# # DiffPlan <- Diff.Calc(PlotTopoPlants_ls[[1]])
# # 
# # df1 <- cbind(DiffAll[,-5:-6], 
# #       DiffAnim[,5] - DiffAll[,5])
# # colnames(df1)[5] <- "Value"
# # df1 <- tTest.Calc(df1)
# # df1$Comparison <- "Top-Down"
# # df2 <- cbind(DiffAll[,-5:-6], 
# #       DiffPlan[,5] - DiffAll[,5])
# # colnames(df2)[5] <- "Value"
# # df2 <- tTest.Calc(df2)
# # df2$Comparison <- "Bottom-Up"
# # 
# # Plot_df <- na.omit(rbind(df1, df2))
# # 
# # # comps <- list( c("Top-Down", "Bottom-Up"))
# # pal_signif <- c("#c32200", "#62c300")
# # pal_cascade <- c("#880094", "#947100")
# # Comps_gg <- lapply(c("n_species", "n_plants", "n_animals", "n_links", "Nestedness", "Modularity"), function(x){
# #   ggplot(Plot_df[Plot_df$Topology == x, ], 
# #          aes(x = factor(IS), y = Value, col = Comparison, fill = Signif)) + 
# #     geom_hline(aes(yintercept = 0)) +
# #     geom_boxplot() + 
# #     scale_fill_manual(values = pal_signif) +
# #     scale_color_manual(values = pal_cascade) +
# #     # stat_compare_means(comparisons = comps, method = 't.test', label = 'p.signif') +
# #     facet_wrap(~factor(Proxy, levels = c("Strength", "Climate", "IUCN")),
# #                scales = "free", nrow = 3) +
# #     theme_bw() +
# #     labs(x = "Interaction Strength % Required For Survival", y = paste0("Change In ", x), title = x)
# # })
# # 
# # pdf(file.path(Dir.Exports, "Comps.pdf"), paper = "a4", height = 12, width = 9)
# # print(Comps_gg)
# # dev.off()
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
