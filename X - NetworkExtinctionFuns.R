.DataInit <- function(x){
  if(class(x) == "network"){
    x
  }
  if(class(x) == "matrix"){
    x <- as.network(x, loops = TRUE)
  }
  return(x)
}

.ExtinctionOrderEK <- function(Network, Order, clust.method = "cluster_infomap", IS = 0.7){
  Link_density <- Modularity <- Grado <- NULL
  edgelist <- as.matrix.network.edgelist(Network,matrix.type="edgelist") #Prey - Predator
  Conected <- data.frame(ID = 1:network.size(Network), Grado = degree(edgelist, c("total")))
  
  Conected1<-  Order
  
  net <- as.matrix.network.adjacency(Network, attrname = "weight")
  netgraph <- suppressMessages(graph_from_adjacency_matrix(net, weighted = TRUE)) #
  strengthbasenet <- igraph::strength(netgraph)
  
  indegreebasenet <- degree(Network, cmode = "indegree")
  
  indegreebasenetzeros <- sum(degree(Network, cmode = "indegree") == 0)
  indegreetopnetzeros <- sum(degree(Network, cmode = "outdegree") == 0)
  
  Producers <- (1:length(degree(Network, cmode = "indegree")))[degree(Network, cmode = "indegree") == 0]
  TopPredators <- (1:length(degree(Network, cmode = "outdegree")))[degree(Network, cmode = "outdegree") == 0]
  
  DF <- data.frame(Spp = rep(NA, network.size(Network)), S = rep(NA, network.size(Network)), L = rep(NA, network.size(Network)), C = rep(NA, network.size(Network)), Link_density = rep(NA, network.size(Network)),SecExt = rep(NA,network.size(Network)), Pred_release = rep(NA,network.size(Network)), Iso_nodes =rep (NA,network.size(Network)))
  
  Secundaryext <- c()
  Predationrel <- c()
  accExt <- c()
  totalExt <- c()
  FinalExt <- list()
  Conected3 <- c()
  
  for (i in 1:length(Order)){
    
    if (length(accExt)==0){
      Temp <- Network
      DF$Spp[i] <- Conected1[i]
      delete.vertices(Temp, c(DF$Spp[1:i]))
    }
    if (length(accExt)>0){
      Temp <- Network
      Temp <- delete.vertices(Temp, c(accExt))
      edgelist <- as.matrix.network.edgelist(Temp,matrix.type="edgelist")
      Conected2 <- data.frame(ID =1:network.size(Temp), Grado = degree(edgelist, c("total")))
      for(j in sort(accExt)){
        Conected2$ID <- ifelse(Conected2$ID < j, Conected2$ID, Conected2$ID + 1)
      }
      
      DF$Spp[i] <- Conected1[i]
      Temp <- Network
      
      delete.vertices(Temp, unique(c(c(DF$Spp[1:i]),accExt)))
    }
    
    DF$S[i] <- network.size(Temp)
    DF$L[i] <- network.edgecount(Temp)
    DF$C[i] <- network.density(Temp)
    DF$Link_density [i] <- DF$L[i]/DF$S[i]
    
    Networkclass = class(Temp)
    
    if (Networkclass[1] == "matrix"){
      netgraph = graph_from_adjacency_matrix(Temp, weighted = TRUE) # !!! somehow bring weighted = NULL back for non-weighted networks
    }
    
    if (Networkclass[1] == "network"){
      net = as.matrix.network.adjacency(Temp, attrname = "weight")
      netgraph = suppressMessages(graph_from_adjacency_matrix(net, weighted = TRUE)) # !!! somehow bring weighted = NULL back for non-weighted networks
    }
    
    if (clust.method == "cluster_edge_betweenness"){
      Membership = suppressWarnings(cluster_edge_betweenness(netgraph, weights=NULL, directed = TRUE, edge.betweenness = TRUE,
                                                             merges = TRUE, bridges = TRUE, modularity = TRUE, membership = TRUE))
    } else if (clust.method == "cluster_spinglass"){
      spins = 107#network.size(Temp)
      Membership = suppressWarnings(cluster_spinglass(netgraph, spins=spins)) #spins could be the Richness
    }else if (clust.method == "cluster_label_prop"){
      Membership = suppressWarnings(cluster_label_prop(netgraph, weights = NULL, initial = NULL,
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
    if(is.na(Membership)){
      DF$Modularity[i] <- NA
    }else{
      DF$Modularity[i] <- suppressWarnings(modularity(Membership))
    }
    
    ## Producer Extinction Secondary
    SecundaryextTemp <- (1:length(degree(Temp, cmode = "indegree")))[degree(Temp, cmode = "indegree") == 0]
    for(j in sort(unique(c(c(DF$Spp[1:i]),accExt)))){ # !!! what does this step do?!?!?!
      SecundaryextTemp <- ifelse(SecundaryextTemp < j, 
                                 SecundaryextTemp, 
                                 SecundaryextTemp + 1)
    }
    Secundaryext <- SecundaryextTemp
    Secundaryext <- Secundaryext[!(Secundaryext %in% Producers)]
    
    
    ## Predator Extinction Secondary
    PredationrelTemp <- (1:length(degree(Temp, cmode = "outdegree")))[degree(Temp, cmode = "outdegree") == 0]
    for(j in sort(unique(c(c(DF$Spp[1:i]),accExt)))){
      PredationrelTemp <- ifelse(PredationrelTemp < j, PredationrelTemp, PredationrelTemp + 1)
    }
    Predationrel <- PredationrelTemp
    Predationrel <- Predationrel[!(Predationrel %in% TopPredators)]
    DF$Pred_release[i] <- length(Predationrel)
    
    
    RelISloss <- igraph::strength(suppressMessages(graph_from_adjacency_matrix(
      as.matrix.network.adjacency(Temp, attrname = "weight"), 
      weighted = TRUE)
      )) / strengthbasenet[names(strengthbasenet) %in% get.vertex.attribute(Temp, "vertex.names")]
    RelISloss < IS
    
    DF$Iso_nodes[i] <- sum(degree(Temp) == 0)
    
    Secundaryext <- which(RelISloss < IS)
    DF$SecExt[i]<- length(Secundaryext)
    
    
    # message(i)
    FinalExt[[i]] <- (Secundaryext)
    accExt <- append(accExt, DF$Spp[1:i])
    accExt <- unique(append(accExt,Secundaryext))
    
    if (DF$L[i] == 0) break
  }
  DF <- DF[complete.cases(DF),]
  DF$AccSecExt <- cumsum(DF$SecExt)
  DF$NumExt <- 1:nrow(DF)
  DF$TotalExt <- DF$AccSecExt + DF$NumExt
  DF <- relocate(DF, Modularity, .after = Link_density)
  class(DF) <- c("data.frame", "ExtinctionOrder")
  return(list(DF, Temp))
}

ExtinctionOrderEK <- function(Network, Order, clust.method = "cluster_infomap", IS = 0.3){
  Grado <- NULL
  Network <- .DataInit(x = Network)
  edgelist <- as.matrix.network.edgelist(Network,matrix.type="edgelist") #Prey - Predator
  Conected <- data.frame(ID = 1:network.size(Network), Grado = degree(edgelist, c("total")))
  
  Conected1<-  Order
  
  indegreebasenet <- degree(Network, cmode = "indegree")
  
  
  
  
  
  net <- as.matrix.network.adjacency(Network, attrname = "weight")
  netgraph <- suppressMessages(graph_from_adjacency_matrix(net, weighted = TRUE)) #
  strengthbasenet <- igraph::strength(netgraph)
  
  
  
  
  
  indegreebasenetzeros <- sum(degree(Network, cmode = "indegree") == 0)
  indegreetopnetzeros <- sum(degree(Network, cmode = "outdegree") == 0)
  
  Producers <- (1:length(degree(Network, cmode = "indegree")))[degree(Network, cmode = "indegree") == 0]
  TopPredators <- (1:length(degree(Network, cmode = "outdegree")))[degree(Network, cmode = "outdegree") == 0]
  
  DF <- data.frame(Spp = rep(NA, network.size(Network)), S = rep(NA, network.size(Network)), L = rep(NA, network.size(Network)), C = rep(NA, network.size(Network)),  SecExt = rep(NA,network.size(Network)), Pred_release = rep(NA,network.size(Network)))
  
  Secundaryext <- c()
  Predationrel <- c()
  accExt <- c()
  totalExt <- c()
  FinalExt <- list()
  Conected3 <- c()
  
  for (i in 1:length(Order)){
    
    if (length(accExt)==0){
      Temp <- Network
      DF$Spp[i] <- Conected1[i]
      delete.vertices(Temp, c(DF$Spp[1:i]))
    }
    if (length(accExt)>0){
      Temp <- Network
      Temp <- delete.vertices(Temp, c(accExt))
      edgelist <- as.matrix.network.edgelist(Temp,matrix.type="edgelist")
      Conected2 <- data.frame(ID =1:network.size(Temp), Grado = degree(edgelist, c("total")))
      for(j in sort(accExt)){
        Conected2$ID <- ifelse(Conected2$ID < j, Conected2$ID, Conected2$ID + 1)
      }
      
      DF$Spp[i] <- Conected1[i]
      Temp <- Network
      
      delete.vertices(Temp, unique(c(c(DF$Spp[1:i]),accExt)))
    }
    
    DF$S[i] <- network.size(Temp)
    DF$L[i] <- network.edgecount(Temp)
    DF$C[i] <- network.density(Temp)
    
    SecundaryextTemp <- (1:length(degree(Temp, cmode = "indegree")))[degree(Temp, cmode = "indegree") == 0]
    for(j in sort(unique(c(c(DF$Spp[1:i]),accExt)))){
      SecundaryextTemp <- ifelse(SecundaryextTemp < j, SecundaryextTemp, SecundaryextTemp + 1)
    }
    Secundaryext <- SecundaryextTemp
    Secundaryext <- Secundaryext[!(Secundaryext %in% Producers)]
    DF$SecExt[i]<- length(Secundaryext)
    
    PredationrelTemp <- (1:length(degree(Temp, cmode = "outdegree")))[degree(Temp, cmode = "outdegree") == 0]
    for(j in sort(unique(c(c(DF$Spp[1:i]),accExt)))){
      PredationrelTemp <- ifelse(PredationrelTemp < j, PredationrelTemp, PredationrelTemp + 1)
    }
    Predationrel <- PredationrelTemp
    Predationrel <- Predationrel[!(Predationrel %in% TopPredators)]
    DF$Pred_release[i] <- length(Predationrel)
    
    
    
    
    
    RelISloss <- igraph::strength(suppressMessages(graph_from_adjacency_matrix(
      as.matrix.network.adjacency(Temp, attrname = "weight"), 
      weighted = TRUE)
    )) / strengthbasenet[names(strengthbasenet) %in% get.vertex.attribute(Temp, "vertex.names")]
    Secundaryext <- which(RelISloss < IS)
    DF$SecExt[i]<- length(Secundaryext)
    
    
    
    
    
    
    # message(i)
    FinalExt[[i]] <-(Secundaryext)
    accExt <- append(accExt, DF$Spp[1:i])
    accExt <- unique(append(accExt,Secundaryext))
    
    if (DF$L[i] == 0) break
  }
  DF <- DF[complete.cases(DF),]
  DF$AccSecExt <- cumsum(DF$SecExt)
  DF$NumExt <- 1:nrow(DF)
  DF$TotalExt <- DF$AccSecExt + DF$NumExt
  class(DF) <- c("data.frame", "ExtinctionOrder")
  return(list(DF, Temp))
}

RandomExtinctionsEK <- function(Network, nsim = 10, parallel = FALSE, ncores, Record = F, plot = F, SimExt = NULL, IS = 0.3){
  NumExt <- sd <- AccSecExt <- AccSecExt_95CI <- AccSecExt_mean <- Lower <- NULL
  network <- .DataInit(x = Network)
  if(is.null(SimExt)){
    SimExt <- network.size(network)
  }
  
  if(parallel){
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    sims <- foreach(i=1:nsim, .packages = c("broom", "doParallel", "dplyr", "foreach", "ggplot2", "igraph", "magrittr", "network", "scales", "sna", "stats", "tidyr", "MASS", "parallel", "purrr"), .export = ".ExtinctionOrderEK")%dopar%{
      sims <- try(.ExtinctionOrderEK(Network = network, Order = sample(1:network.size(network), size = SimExt), IS = IS), silent = T)
      # try({sims$simulation <- i}, silent = T)
      sims
    }
    stopCluster(cl)
  }else{
    sims <- list()
    for(i in 1:nsim){
      sims[[i]] <- try(.ExtinctionOrderEK(Network = network, Order = sample(1:network.size(network), size = SimExt), IS = IS), silent = TRUE)
      try({sims[[i]]$simulation <- i}, silent = TRUE)
      message(paste("Simulation", i, "of", nsim, "ready"))
    }
  }
  
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
  if(plot == T){
    g <- ggplot(sims, aes_string(x = "NumExt", y = "AccSecExt_mean")) + geom_ribbon(aes_string(ymin = "Lower", ymax = "Upper"), fill = muted("red")) + geom_line() + ylab("Acc. Secondary extinctions") + xlab("Primary extinctions") + theme_bw()
    g
  }
  
  if(Record == T & plot == T){
    return(list(sims = sims, graph = g, FullSims = FullSims, nets = temps))
  }else if(Record == F & plot == T){
    return(list(sims = sims, graph = g, nets = temps))
  }else if(Record == F & plot == F){
    return(list(sims, nets = temps))
  }else if(Record == T & plot == F){
    return(list(sims = sims, FullSims = FullSims, nets= temps))
  }
}