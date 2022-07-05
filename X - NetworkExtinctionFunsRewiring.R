.DataInit <- function(x){
  if(class(x) == "network"){
    x
  }
  if(class(x) == "matrix"){
    x <- as.network(x, loops = TRUE)
  }
  return(x)
}

ExtinctionOrder <- function(Network, Order, NetworkType = "Trophic", clust.method = "cluster_infomap",
                            IS = 0, 
                            Rewiring = FALSE, RewiringDist, RewiringProb = 0.5,
                            verbose = TRUE
                            ){
  # Setting up Objects for function run ++++++++++ ++++++++++ ++++++++++ ++++++++++
  Link_density <- Modularity <- Grado <- NULL
  Network <- .DataInit(x = Network)
  edgelist <- as.matrix.network.edgelist(Network,matrix.type="edgelist") #Prey - Predator
  Conected <- data.frame(ID = 1:network.size(Network), Grado = degree(edgelist, c("total")))
  Conected1 <- Order
  
  ## Interaction Strength Loss Preparations ++++++++++ ++++++++++
  if(length(IS )== 1){ # when the same dependency is to be applied for all species
    IS <- rep(IS, network.size(Network)) # repeat IS argument for all species
    names(IS) <- get.vertex.attribute(Network, "vertex.names") # assign species names to IS arguments
  }else{
    if(sum(!(get.vertex.attribute(Network, "vertex.names") %in% names(IS))) != 0){stop("Please ensure that the names of the nodes in your Network are matched by names of the elements in your IS argument vector.")}
  }
  ## Rewiring Preparations ++++++++++ ++++++++++
  if(!isFALSE(Rewiring)){ # if rewiring has been specified
    if(!exists("RewiringDist")){stop("To execute rewiring simulations, you need to specify the RewiringDist argument as well as the Rewiring argument.")}
    diag(RewiringDist)<- NA # set diagonal to NA as one can't rewire to oneself
    if(length(Rewiring) == 1){ # when the same Rewiring happen for all species
      fun <- deparse1(Rewiring) # turn function into string
      Rewiring <- rep(fun, network.size(Network)) # repeat function for each species
      names(Rewiring) <- get.vertex.attribute(Network, "vertex.names") # assign names of species in network to rewiring functions
    }
  }else{
    if(sum(!(get.vertex.attribute(Network, "vertex.names") %in% names(Rewiring))) != 0){stop("Please ensure that the names of the nodes in your Network are matched by names of the elements in your Rewiring argument vector.")}
  }
  
  # Base net calculations ++++++++++ ++++++++++ ++++++++++ ++++++++++
  ## base interaction strengths per node ++++++++++ ++++++++++
  options(warn=-1) # turn warning off
  Weight_mat <- net <- as.matrix.network.adjacency(Network, attrname = "weight")
  options(warn=0) # turn warnings on
  if(sum(IS) != 0){
    # if(sum(get.edge.attribute(Network, "weight"), na.rm = TRUE) == 0){
    #   stop("Either your network does not contain any edges with weights or your network does not have the edge attribute `weight` required for calculation of extinctions based on relative interaction strength loss.")
    # }
    netgraph <- suppressMessages(graph_from_adjacency_matrix(net, weighted = TRUE))
    strengthbasenet <- igraph::strength(netgraph)
  }
  
  ## identification of producers and top predators ++++++++++ ++++++++++
  indegreebasenet <- degree(Network, cmode = "indegree")
  indegreebasenetzeros <- sum(degree(Network, cmode = "indegree") == 0)
  indegreetopnetzeros <- sum(degree(Network, cmode = "outdegree") == 0)
  Producers <- get.vertex.attribute(Network, "vertex.names")[degree(Network, cmode = "indegree") == 0]
  TopPredators <- get.vertex.attribute(Network, "vertex.names")[degree(Network, cmode = "outdegree") == 0]
  
  ## output object ++++++++++ ++++++++++
  DF <- data.frame(Spp = rep(NA, length(Order)), 
                   S = rep(NA, length(Order)), 
                   L = rep(NA, length(Order)), 
                   C = rep(NA, length(Order)), 
                   Link_density = rep(NA, length(Order)),
                   SecExt = rep(NA,length(Order)), 
                   Pred_release = rep(NA,length(Order)), 
                   Iso_nodes = rep (NA,length(Order)))
  Secundaryext <- c()
  Predationrel <- c()
  accExt <- c()
  totalExt <- c()
  FinalExt <- list()
  Conected3 <- c()

  # Sequential extinction simulation ++++++++++ ++++++++++ ++++++++++ ++++++++++
  if(verbose){ProgBar <- txtProgressBar(max = length(Order), style = 3)}
  for (i in 1:length(Order)){
    # print(i)
    ### creating temporary network + deleting vertices if they have been set to go extinct ++++++++++ ++++++++++
    if(length(accExt)==0){ # on first iteration
      Temp <- Network
      DF$Spp[i] <- Conected1[i]
      delete.vertices(Temp, c(DF$Spp[1:i]))
    }
    if (length(accExt)>0){ # on any subsequent iteration
      Temp <- Network
      Temp <- delete.vertices(Temp, c(accExt))
      edgelist <- as.matrix.network.edgelist(Temp,matrix.type="edgelist")
      DF$Spp[i] <- Conected1[i]
      Temp <- Network
      delete.vertices(Temp, unique(c(c(DF$Spp[1:i]),accExt)))
    }
    
    ## network metrics to output object  ++++++++++ ++++++++++
    DF$S[i] <- network.size(Temp)
    DF$L[i] <- network.edgecount(Temp)
    DF$C[i] <- network.density(Temp)
    DF$Link_density[i] <- DF$L[i]/DF$S[i]
    
    ## premature complete annihilation message ++++++++++ ++++++++++
    if(i > 1){
      if(DF$L[i-1] == 0){
        if(verbose){setTxtProgressBar(ProgBar, length(Order))}
        warning(paste("All species in network went extinct through secondary extinction before all primary extinctions were simulated. This happened at extinction step", i-1, "out of", length(Order)))
        break
      }
    }
  
    ## calculating modularity ++++++++++ ++++++++++
    Networkclass = class(Temp)
    if (Networkclass[1] == "matrix"){
      netgraph = graph_from_adjacency_matrix(Temp, mode = "directed", weighted = TRUE)
    }
    if (Networkclass[1] == "network"){
      net = as.matrix.network.adjacency(Temp)
      netgraph = suppressMessages(graph_from_adjacency_matrix(net, mode = "directed", weighted = TRUE))
    }
    if (clust.method == "cluster_edge_betweenness"){
      Membership = suppressWarnings(cluster_edge_betweenness(netgraph, weights = TRUE, directed = TRUE, edge.betweenness = TRUE,
                                                             merges = TRUE, bridges = TRUE, modularity = TRUE, membership = TRUE))
    } else if (clust.method == "cluster_spinglass"){
      spins = 107#network.size(Temp)
      Membership = suppressWarnings(cluster_spinglass(netgraph, spins=spins)) #spins could be the Richness
    }else if (clust.method == "cluster_label_prop"){
      Membership = suppressWarnings(cluster_label_prop(netgraph, weights = TRUE, initial = NULL,
                                                       fixed = NULL))
    }else if (clust.method == "cluster_infomap"){
      nb.trials = 107#network.size(Temp)
      Membership = suppressWarnings(cluster_infomap(as.undirected(netgraph),
                                                    e.weights = E(netgraph)$weight,
                                                    v.weights = NULL,
                                                    nb.trials = nb.trials,
                                                    modularity = TRUE))
      
    } else if (clust.method == "none"){
      Membership = NA
    }else stop('Select a valid method for clustering. ?SimulateExtinction')
    DF$Modularity[i] <- Membership$modularity
    
    ## rewiring ++++++++++ ++++++++++
    accExt <- unique(append(accExt, DF$Spp[1:i]))
    if(!isFALSE(Rewiring)){
      ### identify rewiring potential ++++++++++
      Rewiring_df <- data.frame(Direction = NA,
                                Species = NA,
                                NewPartner = NA,
                                LostPartner = NA,
                                IS = NA)
      Rewiring_df <- na.omit(Rewiring_df)
      #### loop over all deleted vertices and the connections lost because of their exclusion
      for(Iter_PrimaryExt in 1:length(accExt)){
        # Iter_PrimaryExt = 1
        LostPartner <- get.vertex.attribute(Network, "vertex.names")[accExt[Iter_PrimaryExt]] # name of primary extinction species
        LostISCol <- Weight_mat[, LostPartner] # lost interaction strength with nodes now slated for secondary extinction
        LostISRow <- Weight_mat[LostPartner, ]
        Lost_df <- data.frame(LostIS = c(LostISCol, LostISRow),
                              Direction = rep(c(1,2), c(length(LostISCol), length(LostISRow))),
                              names = c(names(LostISCol), names(LostISRow))
        )
        Lost_df <- Lost_df[Lost_df$LostIS != 0, ]
        if(nrow(Lost_df)!=0){
          for(Iter_LostIS in 1:nrow(Lost_df)){ ## looping over all species that were linked to the current primary extinction
            # Iter_LostIS = 1
            LostPartnerSim <- eval(str2lang(Rewiring[which(names(Rewiring) == Lost_df$names[Iter_LostIS])]))(RewiringDist[,LostPartner]) # probability of rewiring too each node in network given rewiring function and species similraity
            RewiringCandidates <- LostPartnerSim[LostPartnerSim > RewiringProb & names(LostPartnerSim) %in% get.vertex.attribute(Temp, "vertex.names")] # rewiring probability for nodes still in temporary network and having a higher rewiring probability than 0.5
            RewiredPartner <- names(which.max(RewiringCandidates)) # most likely rewiring partner
            if(!is.null(RewiredPartner)){ # if a rewired partner has been found
              Rewiring_df <- rbind(Rewiring_df, 
                                   data.frame(Direction = Lost_df$Direction[Iter_LostIS],
                                              Species = Lost_df$names[Iter_LostIS],
                                              NewPartner = RewiredPartner,
                                              LostPartner = LostPartner,
                                              IS = Lost_df$LostIS[Iter_LostIS])
              )
            }
          }
        }
      }
      
      ### shift rewired interaction strengths ++++++++++
      if(nrow(Rewiring_df) != 0){
        #### shift interaction weights in Weight_mat
        for(Iter_Rewiring in 1:nrow(Rewiring_df)){
          # Iter_Rewiring = 1
          ## assigning shifted interaction strength
          ColSpec <- Rewiring_df[Iter_Rewiring,4-Rewiring_df[Iter_Rewiring,"Direction"]]
          RowSpec <- Rewiring_df[Iter_Rewiring,1+Rewiring_df[Iter_Rewiring,"Direction"]]
          Weight_mat[RowSpec, ColSpec] <- Weight_mat[RowSpec, ColSpec] + Rewiring_df[Iter_Rewiring,"IS"]
          ## deleting shiften interaction strength
          
          ColLost <- ifelse(Rewiring_df[Iter_Rewiring, "Direction"] == 1, 
                            Rewiring_df[Iter_Rewiring, "LostPartner"],
                            Rewiring_df[Iter_Rewiring, "Species"])
          RowLost <- ifelse(Rewiring_df[Iter_Rewiring, "Direction"] == 1, 
                            Rewiring_df[Iter_Rewiring, "Species"],
                            Rewiring_df[Iter_Rewiring, "LostPartner"])
          Weight_mat[RowLost, ColLost] <- 0
        }
        #### establishing rewired network and deleting primary extinction nodes
        Network <- as.network(Weight_mat, matrix.type = "adjacency", ignore.eval=FALSE, names.eval='weight')
        Temp <- Network
        delete.vertices(Temp, unique(c(c(DF$Spp[1:i]),accExt)))
      }
    }
    
    ## identify secondary extinctions ++++++++++ ++++++++++
    ### Relative Interaction Strength loss ++++++++++
      if(sum(IS) == 0){
        SecundaryextNames <- get.vertex.attribute(Temp, "vertex.names")[which(degree(Temp) == 0)]
        Secundaryext <- match(SecundaryextNames, get.vertex.attribute(Network, "vertex.names"))
      }else{
        AbsIS <- igraph::strength(suppressMessages(graph_from_adjacency_matrix(
          as.matrix.network.adjacency(Temp, attrname = "weight"),
          weighted = TRUE)
        ))
        RelISloss <-  AbsIS / strengthbasenet[names(strengthbasenet) %in% get.vertex.attribute(Temp, "vertex.names")]
        SecundaryextNames <- which(AbsIS == 0 | (1-RelISloss) > IS[match(names(RelISloss), names(IS))])
        Secundaryext <- match(names(SecundaryextNames), get.vertex.attribute(Network, "vertex.names"))
      }
      
    ### for trophic networks ++++++++++
    if(NetworkType == "Trophic"){
      DF$SecExt[i] <- length(SecundaryextNames[!(names(SecundaryextNames) %in% Producers)])
      DF$Pred_release[i] <- length(SecundaryextNames[!(names(SecundaryextNames) %in% TopPredators)])
    }
    ### for mutualistic networks ++++++++++
    if(NetworkType == "Mutualistic"){
      DF$SecExt[i] <- length(Secundaryext)
    }
    DF$Iso_nodes[i] <- sum(degree(Temp) == 0)
  
    ## Return of objects ++++++++++ ++++++++++
    FinalExt[[i]] <- Secundaryext
    accExt <- append(accExt, DF$Spp[1:i])
    accExt <- unique(append(accExt,Secundaryext))
    if(verbose){setTxtProgressBar(ProgBar, i)}
  }
  
  # return of final data objects ++++++++++ ++++++++++
  DF <- DF[!is.na(DF$Spp),]
  DF$AccSecExt <- cumsum(DF$SecExt)
  DF$NumExt <- 1:nrow(DF)
  DF$TotalExt <- DF$AccSecExt + DF$NumExt
  DF <- relocate(DF, Modularity, .after = Link_density)
  class(DF) <- c("data.frame", "SimulateExt")
  return(list(sims = DF,
              Network = Temp))
}

RandomExtinctions <- function(Network, nsim = 10, 
                              Record = FALSE, plot = FALSE,
                              SimNum = NULL,
                              NetworkType = "Trophic", clust.method = "cluster_infomap",
                              parallel = FALSE, ncores,
                              IS = 0, 
                              Rewiring = FALSE, RewiringDist, RewiringProb = 0.5,
                              verbose = TRUE){
  ## setting up objects
  NumExt <- sd <- AccSecExt <- AccSecExt_95CI <- AccSecExt_mean <- Lower <- NULL
  network <- .DataInit(x = Network)
  if(is.null(SimNum)){
    SimNum <- network.size(network)
  }
  
  ## simulations
  if(verbose){ProgBar <- txtProgressBar(max = nsim, style = 3)}
  if(parallel){
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    sims <- foreach(i=1:nsim, .packages = "NetworkExtinction")%dopar%{
      sims <- try(ExtinctionOrder(Network = network, Order = sample(1:network.size(network), size = SimNum), 
                                  IS = IS, 
                                  Rewiring = Rewiring, RewiringDist = RewiringDist,
                                  verbose = FALSE, RewiringProb = RewiringProb), silent = TRUE)
      try({sims$simulation <- i}, silent = TRUE)
      sims
    }
    stopCluster(cl)
  }else{
    sims <- list()
    for(i in 1:nsim){
      sims[[i]] <- try(ExtinctionOrder(Network = network, Order = sample(1:network.size(network), size = SimNum), 
                                       IS = IS, 
                                       Rewiring = Rewiring, RewiringDist = RewiringDist,
                                       verbose = FALSE, RewiringProb = RewiringProb), silent = TRUE)
      try({sims[[i]]$simulation <- i}, silent = TRUE)
      if(verbose){setTxtProgressBar(ProgBar, i)}
    }
  }
  
  ## extract objects
  temps <- lapply(sims, "[[", 2)
  sims <- lapply(sims, "[[", 1)
  cond <- sapply(sims, function(x) "data.frame" %in% class(x))
  cond <- c(1:length(cond))[cond]
  sims <- sims[cond]
  sims <- do.call(rbind, sims)
  if(Record == TRUE){
    FullSims <- sims
  }
  
  sims <- sims %>% group_by(NumExt) %>% summarise(AccSecExt_95CI = 1.96*sd(AccSecExt), AccSecExt_mean = mean(AccSecExt)) %>% mutate(Upper = AccSecExt_mean + AccSecExt_95CI, Lower = AccSecExt_mean - AccSecExt_95CI, Lower = ifelse(Lower < 0, 0, Lower))
  
  ## plot output
  if(plot == TRUE){
    g <- ggplot(sims, aes_string(x = "NumExt", y = "AccSecExt_mean")) + geom_ribbon(aes_string(ymin = "Lower", ymax = "Upper"), fill = muted("red")) + geom_line() + ylab("Acc. Secondary extinctions") + xlab("Primary extinctions") + theme_bw()
    g
  }
  
  ## object output
  if(Record == T & plot == T){
    return(list(sims = sims, graph = g, FullSims = FullSims, nets = temps))
  }else if(Record == F & plot == T){
    return(list(sims = sims, graph = g, nets = temps))
  }else if(Record == F & plot == F){
    return(list(sims = sims, nets = temps))
  }else if(Record == T & plot == F){
    return(list(sims = sims, FullSims = FullSims, nets= temps))
  }
}

SimulateExtinctions <- function(Network, Method, Order = NULL, 
                                NetworkType = "Trophic", clust.method = "cluster_infomap",
                                IS = 0, 
                                Rewiring = FALSE, RewiringDist, RewiringProb = 0.5,
                                verbose = TRUE){
  Network <- .DataInit(x = Network)
  
  if(!is.null(Order)){Method <- "Ordered"}
  
  '%ni%'<- Negate('%in%')
  if(Method %ni% c("Mostconnected", "Ordered")) stop('Choose the right method. See ?SimulateExtinction.')
  
  if(Method == "Mostconnected"){
    edgelist <- as.matrix.network.edgelist(Network,matrix.type="edgelist") #Prey - Predator
    Conected <- data.frame(ID = 1:network.size(Network), Grado = degree(edgelist, c("total")))
    Conected <- arrange(Conected, desc(Grado))
    DF <- ExtinctionOrder(Network = Network, Order = Conected$ID, clust.method = clust.method,
                          IS = IS, Rewiring = Rewiring, RewiringDist = RewiringDist, 
                          verbose = verbose, RewiringProb = RewiringProb, NetworkType = NetworkType)
  }
  if(Method == "Ordered"){
    DF <- ExtinctionOrder(Network = Network, Order = Order, clust.method = clust.method,
                          IS = IS, Rewiring = Rewiring, RewiringDist = RewiringDist, 
                          verbose = verbose, RewiringProb = RewiringProb, NetworkType = NetworkType)
  }
  
  return(DF)
}

