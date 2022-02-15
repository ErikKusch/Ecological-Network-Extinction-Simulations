#' ####################################################################### #
#' PROJECT: [Biodiversity Simplification & Ecological Network Topology] 
#' CONTENTS: 
#'  - Network Data Loading and Formatting from Fricke 2021
#'  DEPENDENCIES:
#'  - 
#' AUTHOR: [Evan Fricke]
#' ####################################################################### #

# BIPARTITE HELPERS ========================================================
l.unique <- function(x) length(unique(x))

as.num.char <- function(x) as.numeric(as.character(x))

# NETWORK SPREADING ========================================================
net.spread <- function(split.by, split.vals, 
                       tax.type, data.type, 
                       long.df = long.df,
                       min.taxa = 2){
  if(min.taxa <= 1) stop("cannot make network with only one plant or animal taxon") # This may be obvious
  zz <- apply(cbind(split.by, as.character(split.vals), tax.type, data.type), 1, 
              function(xx) net.spread.inside(split.by = xx[1], split.vals = xx[2], tax.type = xx[3], data.type = xx[4], long.df = long.df))
  names(zz) <- split.vals
  zz <- zz[!unlist(lapply(zz, function(zzz) any(dim(zzz) <= min.taxa) | is.null(dim(zzz))))]
  zz
}

# NETWORK SPREADING WTHIN NETWORK ==========================================
net.spread.inside <- function(split.by = split.by, split.vals = split.vals, 
                              tax.type = tax.type, data.type = data.type,
                              long.df = long.df){
  
  x <- long.df[which(long.df[,split.by] == split.vals),]
  a <- x[ ,paste("animal", tax.type, sep = ".")]
  p <- x[ ,paste("plant", tax.type, sep = ".")]
  identifier <- paste(a, p, sep = " ~ ")
  y <- tapply(x$value, identifier, sum)
  
  if(data.type == "bin"){
    y <- ifelse(y > 0, 1, 0)
  }
  
  z <- spread(as.data.frame(cbind(str_split_fixed(names(y), " ~ ", 2), y)),
              value = 3,
              key = 1,
              fill = 0)
  #print(split.vals) # In case you want to track down problems
  
  rownames(z) <- z[,1]
  
  z <- z[which(!rownames(z) == "NA"), ]
  z <- z[, which(!colnames(z) == "NA")]
  
  
  # This chunk was added to make it so cases where there is 1 or 0 of the 
  # focal plant/frug taxa (at whatever focal taxa level, e.g., accepted species)
  # the first column doesn't have to get removed (because it's sort of irrelevant because these will not be added to the list).
  if(!is.null(dim(z))){
    
    z <- z[,-1] # Can do z[,-1, drop = F] to retain a 1 column dataframe rather than vector
    
    if(!is.null(dim(z))){ # The next chunk wont work if there is only 1 
      z[] <- mutate_all(z, list(as.num.char)) # Make numeric while retaining row names...
      
      z <- z[which(!rowSums(z) == 0), ] # Get rid of plant or animal taxa with no interactions
      z <- z[, which(!colSums(z) == 0)]
    }
  }
  
  z
  
}

# NETWORK ORDERING =========================================================
order.net <- function(net.to.order, ref.net){
  
  rs <- rowSums(ref.net > 0)
  cs <- colSums(ref.net > 0)
  
  row.ordered <- rownames(ref.net)[order(rs, decreasing = T)]
  col.ordered <- colnames(ref.net)[order(cs, decreasing = T)]
  
  return(net.to.order[row.ordered, col.ordered])
  
}



