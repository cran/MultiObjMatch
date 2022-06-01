#' Internal helper to build infinity sparse matrix
#'
#' @description Formats the data and make a call to
#'   \link[optmatch]{InfinitySparseMatrix-class}
#'
#' @param data the input numeric vector of cost
#' @param cols the input numeric vector corresponding to control units
#' @param rows the input numeric vector corresponding to treated units
#' @param colnames (optional) vector containing names for all control units
#' @param rownames (optional) vectpr containing names for all treated units
#' @param dimension (optional) vector of number of treated and control units
#' @param call (optional) funtion call used to create the InfinitySparseMatrix
#'
#' @return an InfinitySparseMatrix object
makeInfinitySparseMatrix <-
  function(data,
           cols,
           rows,
           colnames = NULL,
           rownames = NULL,
           dimension = NULL,
           call = NULL) {
    if (!all.equal(length(data), length(cols), length(rows))) {
      stop("Data and column/row ids must be vectors of the same length")
    }
    
    if (is.null(dimension)) {
      dimension <-
        c(max(c(rows, length(rownames))), max(c(cols, length(colnames))))
    }
    
    # if(is.null(rownames)) {
    #   rownames <- paste("T", 1:dimension[1], sep = "")
    # }
    
    # if(is.null(colnames)) {
    #   colnames <- paste("C", 1:dimension[2], sep = "")
    # }
    
    
    return(
      new(
        "InfinitySparseMatrix",
        data,
        cols = as.integer(cols),
        rows = as.integer(rows),
        colnames = colnames,
        rownames = rownames,
        dimension = as.integer(dimension),
        call = call
      )
    )
  }




#' Create network flow structure
#'
#' @param z a vector of treatment vectors
#' @param IDs (optional) the name of the units
#'
#' @return a networks structure
netFlowMatch <- function(z, IDs = NULL) {
  if (sum((z != 0) & (z != 1)) > 0)
    stop('Vector z must be binary')
  net <- structure(list(), class = 'netFlowMatch')
  net$treatedNodes <- c(1:sum(z))
  net$controlNodes <- sum(z) + 1:sum(z == 0)
  if (is.null(IDs)) {
    names(net$treatedNodes) <- which(z == 1)
    names(net$controlNodes) <- which(z == 0)
  } else{
    if (length(z) != length(IDs))
      stop("Vectors z and IDs must have same length")
    names(net$treatedNodes) <- IDs[z == 1]
    names(net$controlNodes) <- IDs[z == 0]
  }
  net$dense = TRUE
  net$sink <- length(z) + 1
  net$totalNodes <- length(z) + 1
  net
}


#' Helper function to mask edges
#'
#' @description Remove some of the treatment-control edges from a network flow
#'   representation of a match (forbidding those pairings)
#' @param net the network object
#' @param mask a matrix indicating whether to exclude the corresponding edge
#' @param replaceMask (optional) whether to mask
#'
#' @return the masked network structure object
makeSparse <- function(net, mask, replaceMask = TRUE) {
  matrix.mask <- FALSE
  sparse.mask <- TRUE
  
  #validate network flow problem
  if (!inherits(net, "netFlowMatch"))
    stop("Argument net is not of class netFlowMatch")
  
  
  #rcbalance-style edgelist case
  if (!inherits(mask, "list"))
    stop(
      'Argument mask should be a matrix or the output of a call
      to distance-building functions in optmatch or rcbalance packages'
    )
  
  if (length(mask) != length(net$treatedNodes)) {
    stop('Dimensions of mask argument do not match dimensions of problem net')
  }
  
  if (all(laply(mask, length) == length(net$controlNodes)))
    sparse.mask <- FALSE
  #transform to structure-only edgelist used in this package
  mask <-
    lapply(mask, function(x)
      as.numeric(names(x)) + length(net$treatedNodes))
  order(c(0.1, 0.1, 0.1, 0.3, 0))
  #if mask is not sparse update accordingly and return
  #	if(!sparse.mask){
  #		if(replaceMask)	net$dense <- TRUE
  #		return(net)
  #	}
  
  if (replaceMask || net$dense) {
    #convert mask to internal netFlowMatch format and return
    if (matrix.mask) {
      mask <- matrix2edgelist(mask)
    }
    net$edgeStructure <- mask
    net$dense <- FALSE
    return(net)
  } else {
    existing <- net$edgeStructure
    net$edgeStructure <- meldMask(existing, mask)
    return(net)
  }
}




#' Helper function to convert matrix to list
#'
#' @description Convert between a matrix representation of distances between
#'   treated and control units and a list of vectors (default format for
#'   build.dist.struct function in rcbalance package)
#' @param mat matrix representation of distances between treated and control
#'   units
#'
#' @return list of vector representation of distances
matrix2edgelist <- function(mat) {
  edgeList <- list()
  col.select <- rep(TRUE, dim(mat)[2])
  row.stub <- rep(FALSE, dim(mat)[1])
  for (i in 1:dim(mat)[1]) {
    step.i <- row.stub
    step.i[i] <- TRUE
    list.entry <- as.matrix(subset(mat, step.i, col.select))
    save.names <- dimnames(list.entry)
    edgeList[[i]] <- which(is.finite(list.entry))
    if (length(edgeList[[i]]) > 0) {
      names(edgeList[[i]]) <-
        save.names[[2]][which(is.finite(list.entry))]
    }
    names(edgeList)[i] <- save.names[[1]][[1]]
    
  }
  return(edgeList)
}


#' Helper function to combine two sparse distances
#'
#' @description Combine two sparse distances, 
#' allowing only pairings allowed by both
#'
#' @param mask1 matrix of the first mask
#' @param mask2 matrix of the second mask
#'
#' @return combined mask structure
meldMask <- function(mask1, mask2) {
  if (length(mask1) != length(mask2))
    stop('New distance and old disagree on number of treated')
  index.mat <- as.data.frame(rbind(cbind(rep(
    1:length(mask1), laply(mask1, length)
  ), unlist(mask1)),
  cbind(rep(
    1:length(mask2), laply(mask2, length)
  ), unlist(mask2))))
  colnames(index.mat) <- c('rows', 'cols')
  index.mat <- unique(index.mat[duplicated(index.mat),])
  #add back treated with no matches into matrix so they don't get dropped.
  dropped.treated <-
    setdiff(1:length(mask1), unique(index.mat$rows))
  index.mat <-
    rbind(index.mat, data.frame('rows' = dropped.treated, 
                                'cols' = rep(NA, length(dropped.treated))))
  #avoid compilation error
  rows <- NULL
  out.list <- dlply(index.mat, .(rows), function(x)
    x$cols)
  #replace NAs with 0-length lists
  out.list <-
    llply(out.list, function(x)
      if (is.na(x[1])) {
        c()
      } else{
        x
      })
  out.list
}



#' Add fine balance edges
#'
#' @param net the network object created for the network flow problem
#' @param treatedVals the balance value for treated nodes
#' @param controlVals the balance value for control nodes
#' @param replaceExisting (optional) whether or not to replace the existing net;
#'   TRUE by default
#'
#' @return the network structure with balance edges added
addBalance <-
  function(net,
           treatedVals,
           controlVals,
           replaceExisting = TRUE) {
    if (length(net$treatedNodes) != length(treatedVals) ||
        length(net$controlNodes) != length(controlVals)) {
      stop('Wrong numbers of control and treated values in arguments provided')
    }
    new.layer <- c(treatedVals, controlVals)
    #define a treatment assignment vector to track locations in vector
    z <- c(rep(1, length(treatedVals)), rep(0, length(controlVals)))
    if ('balanceStructure' %in% names(net) && !replaceExisting) {
      n.levels <- ncol(net$balanceStructure$variables)
      parent.layer <- NA
      for (i in c(n.levels:1)) {
        nest.tab <- table(net$balanceStructure$variables[, i], new.layer)
        ncoarse.by.fine <- apply(nest.tab, 2, function(x)
          sum(x > 0))
        #check if new.layer nests inside this layer
        if (all(ncoarse.by.fine <= 1)) {
          if (i == n.levels) {
            parent.layer <- i
            net$balanceStructure$variables <-
              cbind(net$balanceStructure$variables, new.layer)
            colnames(net$balanceStructure$variables)[i + 1] <-
              paste('f', i + 1, sep = '.')
            n.levels <- n.levels + 1
            break
          }
          #for i < n.levels, check nesting in lower layer
          nest.tab2 <-
            table(new.layer, net$balanceStructure$variables[, i + 1])
          ncoarse.by.fine2 <- apply(nest.tab2, 2, function(x)
            sum(x > 0))
          if (all(ncoarse.by.fine2 <= 1)) {
            parent.layer <- i
            net$balanceStructure$variables <-
              cbind(
                net$balanceStructure$variables[, c(1:i)],
                new.layer,
                net$balanceStructure$variables[, c((i + 1):n.levels)]
              )
            n.levels <- n.levels + 1
            colnames(net$balanceStructure$variables) <-
              paste('f', c(1:n.levels), sep = '.')
            break
          }
        }
      }
      if (is.na(parent.layer)) {
        stop(
          'New balance variable does not nest in existing balance structure; 
          try running with replaceExisting = TRUE'
        )
      }
      #adjust index to point to newly-added layer
      i <- parent.layer + 1
    } else {
      #if we're replacing an existing structure, account for nodes we just
      #scrapped
      if ('balanceStructure' %in% names(net)) {
        net$totalNodes <-
          net$totalNodes - length(unlist(net$balanceStructure$nodes))
      }
      #initialize a variables matrix
      net$balanceStructure$variables <-
        cbind(rep(1, length(new.layer)), new.layer)
      i = 2
      net$balanceStructure$nodes <- list()
      net$balanceStructure$edges <- list()
    }
    
    
    #find parent nodes for nodes in current layer
    ztab <- table(net$balanceStructure$variables[, i], z)
    zerobins <- which(apply(ztab, 1, function(x)
      all(x == 0)))
    if (length(zerobins) > 0) {
      ztab <- ztab[-zerobins,]
    }
    
    nest.tab <-
      table(net$balanceStructure$variables[, i - 1],
            net$balanceStructure$variables[, i])
    parent.categories <-
      apply(nest.tab, 2, function(x)
        rownames(nest.tab)[which (x > 0)])
    if (length(zerobins) > 0) {
      parent.categories <- parent.categories[-zerobins]
    }
    
    #create nodes for the new layer and store them in a nodes object (part of a
    #list)
    nodes <-
      matrix(net$totalNodes + 1:(3 * length(unique(new.layer))), ncol = 3)
    colnames(nodes) <-
      c('lambda', 'lambda.prime', 'lambda.double.prime')
    rownames(nodes) <- rownames(ztab)
    nodes <- as.data.frame(nodes)
    net$balanceStructure$nodes <-
      append(net$balanceStructure$nodes, list(nodes), i - 2)
    oldNodeCount <- net$totalNodes
    net$totalNodes <- max(nodes)
    
    
    #create edges internal to triangles and store them in an edges
    #object (part of a list)
    internalEdges <-
      data.frame('startn' = as.vector(as.matrix(nodes[, c(1, 1, 2)])),
                 'endn' = as.vector(as.matrix((nodes[, c(2, 3, 3)]))))
    #add edge capacities
    internalEdges$ucap <-
      c(rep(Inf, nrow(nodes)), ztab[, 2], rep(Inf, nrow(nodes)))
    
    #create new edges to layer above, overwriting old ones
    if (i - 1 > 1) {
      #ignore this step if layer above is the sink
      node.lookup <-
        match(parent.categories,
              rownames(net$balanceStructure$nodes[[i - 2]]))
      parent.nodes <-
        net$balanceStructure$nodes[[i - 2]]$lambda[node.lookup]
      nodes[[i - 2]]$externalEdges <-
        cbind(nodes$lambda.double.prime, parent.nodes, rep(Inf, nrow(nodes)))
    }
    
    #create edges to layer below
    if (i - 1 == length(net$balanceStructure$nodes)) {
      #connect new lowest layer to controls
      parent.categories <-
        net$balanceStructure$variables[net$controlNodes, i]
      #look up indices of parent nodes
      parent.nodes <-
        net$balanceStructure$nodes[[i - 1]]$lambda[match(parent.categories,
                              rownames(net$balanceStructure$nodes[[i - 1]]))]
      externalEdges <-
        cbind(net$controlNodes, parent.nodes, rep(1, length(parent.nodes)))
    } else {
      tab <- table(net$balanceStructure$variables[, i], z)
      zerobins <- which(apply(ztab, 1, function(x)
        all(x == 0)))
      if (length(zerobins) > 0) {
        ztab <- ztab[-zerobins,]
      }
      nest.tab <-
        table(net$balanceStructure$variables[, i - 1],
              net$balanceStructure$variables[, i])
      parent.categories <-
        apply(nest.tab, 2, function(x)
          rownames(nest.tab)[which (x > 0)])
      pc.found <- sapply(parent.categories, length)
      zerobins <- which(pc.found == 0)
      if (length(zerobins) > 0) {
        parent.categories <- parent.categories[-zerobins]
      }
      node.lookup <-
        match(parent.categories,
              rownames(net$balanceStructure$nodes[[i - 1]]))
      parent.nodes <-
        net$balanceStructure$nodes[[i - 1]]$lambda[node.lookup]
      externalEdges <-
        cbind(
          net$balanceStructure$nodes[[i]]$lambda.double.prime,
          parent.nodes,
          rep(Inf, length(parent.nodes))
        )
    }
    net$balanceStructure$edges <-
      append(net$balanceStructure$edges, list(
        list('internalEdges' = internalEdges, 'externalEdges' = externalEdges)
      ), i - 2)
    
    if ('exclusionEdges' %in% names(net))
      net <- addExclusion(net)
    
    net
  }




#' Add exclusion edges
#'
#' @param net the input network structure
#' @param remove (optional) whether to exclude edges; FALSE by default
#'
#' @return the network structure with exclusion edges added to allow for trdeoff
#'   for the exclusion cost
addExclusion <- function(net, remove = FALSE) {
  if (!inherits(net, "netFlowMatch"))
    stop("Argument net is not of class netFlowMatch")
  if (remove) {
    if ('exclusionEdges' %in% names(net))
      net <- net[-which(names(net) == 'exclusionEdges')]
    return(net)
  }
  if ('balanceStructure' %in% names(net)) {
    nlayer <- length(net$balanceStructure$nodes)
    ntreat <- length(net$treatedNodes)
    endn <- net$balanceStructure$nodes[[nlayer]]$lambda[match(
      as.character(net$balanceStructure$variables[1:ntreat, nlayer + 1]),
      rownames(net$balanceStructure$nodes[[nlayer]])
    )]
  } else{
    endn <- rep(net$sink, length(net$treatedNodes))
  }
  net$exclusionEdges <-
    data.frame('startn' = net$treatedNodes, 'endn' = endn)
  net
}

#' Change the edgelist to the infinity sparse matrix
#'
#' @param elist the vector of the edges
#'
#' @return the infinity sparse matrix object
edgelist2ISM <- function(elist) {
  row.idx.list <- list()
  for (i in 1:length(elist)) {
    row.idx.list[[i]] <- rep(i, length(elist[[i]]))
  }
  
  makeInfinitySparseMatrix(
    rep(1, length(unlist(elist))),
    cols = unlist(elist),
    rows = unlist(row.idx.list),
    rownames = as.character(1:length(elist)),
    colnames = as.character(1:max(unlist(row.idx.list)))
  )
}


#' change the distance matrix to cost
#'
#' @param net the network structure
#' @param distance distance matrix
#'
#' @return the vector of cost
matrix2cost <- function(net, distance) {
  if (dim(distance)[1] != length(net$treatedNodes) ||
      dim(distance)[2] != length(net$controlNodes)) {
    stop('Dimensions of distance object and problem do not agree')
  }
  #if internal network structure is sparse, apply sparsity mask to distance
  if (!net$dense) {
    mask <- edgelist2ISM(net$edgeStructure)
    distance <- distance * mask
  }
  #vectorize distance
  if (inherits(distance, 'InfinitySparseMatrix')) {
    #ISMs index row-first
    distance.v <- as.vector(distance)
  } else {
    #regular matrices index column-first and don't automatically drop infinities
    distance.v <- as.vector(t(distance))
    distance.v <- distance.v[is.finite(distance.v)]
  }
  distance.v
  #TODO -- add zeroes if other edges are present in network.
}



#' Create cost skeleton
#'
#' @description Create a more user-friendly data structure to represent the edge
#'   costs in a network.  Internally the network object used by the optmiization
#'   routine represents all the edge costs in a single vector. The "skeleton"
#'   structure decomposes this vector into a list of components, each
#'   corresponding to a different role in the network: "pairings" are edges
#'   between treated and control, exclusion" are direct links between treated
#'   units and a sink that allows them to be excluded, "balance" refers to edges
#'   that count marginal balance between groups, and "sink" indicates edges that
#'   connect control nodes to the sink. Skeletons are created so these various
#'   features can be combined (or switched on and off) easily into objective
#'   functions, and the interface to the main tradeoff function expects to see
#'   each function represented in skeleton format.
#' @param net the network structure
#'
#' @return the skeleton
costSkeleton <- function(net) {
  if (!inherits(net, "netFlowMatch"))
    stop("Argument net is not of class netFlowMatch")
  skeleton <- list()
  if (net$dense) {
    skeleton$pairings <-
      rep(0, length(net$treatedNodes) * length(net$controlNodes))
  } else {
    skeleton$pairings <- rep(0, length(unlist(net$edgeStructure)))
  }
  if ('exclusionEdges' %in% names(net)) {
    skeleton$exclusion <- rep(0, length(net$treatedNodes))
  }
  if ('balanceStructure' %in% names(net)) {
    skeleton$balance <-
      rep(0, sum(laply(net$balanceStructure$edges, function(x)
        nrow(x$internalEdges) + nrow(x$externalEdges))))
    skeleton$sink <- rep(0, nrow(net$balanceStructure$nodes[[1]]))
  } else {
    skeleton$sink <- rep(0, length(net$controlNodes))
  }
  return(skeleton)
}



#' Create a skeleton representation of the edge costs
#'
#' @param dist.struct the distance structure
#' @param net the net structure
#'
#' @return the skeleton representation of the given distance pairs and the net
pairCosts <- function(dist.struct, net) {
  if (!inherits(net, "netFlowMatch"))
    stop("Argument net is not of class netFlowMatch")
  if (length(dist.struct) != length(net$treatedNodes))
    stop('Number of treated units does not agree 
         between distance object and network')
  if (max(unlist(llply(dist.struct, function(x)
    as.numeric(names(x))))) > length(net$controlNodes))
    stop('Number of control units does not 
         agree between distance object and network')
  
  #thin down dist.struct to sparsity of network
  if (!net$dense) {
    discrep.v <-
      laply(dist.struct, length) != laply(net$edgeStructure, length)
    for (i in which(discrep.v)) {
      dist.struct[[i]] <-
        dist.struct[[i]][(as.numeric(names(dist.struct[[i]])) +
                            length(net$treatedNodes)) %in% net$edgeStructure[[i]]]
    }
  }
  
  skeleton <- costSkeleton(net)
  pair.v <- unlist(dist.struct)
  if (length(skeleton$pairings) != length(pair.v))
    stop('Distance structure is sparser than network,
         please run makeSparse so they agree')
  skeleton$pairings <- pair.v
  return(skeleton)
}


#' Turns a skeleton representation of edge costs in a network
#'
#' @description Turns a skeleton representation of edge costs in a network back
#'   into the vector representation expected by the optimization routine. See
#'   comment on the costSkeleton function for more details.
#' @param skeleton the skeleton structure
#'
#' @return vector representation expected by the optimization routine.
flattenSkeleton <- function(skeleton) {
  outv <-
    c(skeleton$pairings,
      skeleton$exclusion,
      skeleton$balance,
      skeleton$sink)
  if (is.null(names(outv)))
    names(outv) <- rep('', length(outv))
  new.name.v <- which(names(outv) == '')
  names(outv)[new.name.v] <-
    paste('edge', 1:length(new.name.v), sep = '.')
  outv
}


#' Create a skeleton representation of the exclusion edge costs
#'
#' @description Create a skeleton representation of the exclusion edge costs
#'   associated with pairings for a given distance and network
#' @param net the network structure
#' @param exclude.penalty (optional) numeric penalty for excluding a treated
#'   unit; 1 by default
#'
#' @return the skeleton with exclusion edge cost
excludeCosts <- function(net, exclude.penalty = 1) {
  if (!inherits(net, "netFlowMatch"))
    stop("Argument net is not of class netFlowMatch")
  if (!('exclusionEdges' %in% names(net)))
    stop("No exclusion edges to penalize")
  
  skeleton <- costSkeleton(net)
  skeleton$exclusion <- skeleton$exclusion + exclude.penalty
  return(skeleton)
}


#' Create a skeleton representation of the balance edge costs
#' @description Create a skeleton representation of the balance edge costs
#'   associated with pairings for a given distance and network
#' @param net the network structure
#' @param balance.penalty (optional) the numeric value for balance; 1 by default
#'
#' @return the skeleton with balance edge cost
balanceCosts <- function(net, balance.penalty = 1) {
  if (!inherits(net, "netFlowMatch"))
    stop("Argument net is not of class netFlowMatch")
  if (!('balanceStructure' %in% names(net)))
    stop("No balance structure to penalize")
  if (length(balance.penalty) != length(net$balanceStructure$nodes))
    stop('Incorrect number of penalties given for balance structure in network')
  
  costv <- c()
  for (i in 1:length(net$balanceStructure$nodes)) {
    lambda.prime <- net$balanceStructure$nodes[[i]][, 2]
    internal.subv <-
      rep(0, nrow(net$balanceStructure$edges[[i]]$internalEdges))
    internal.subv[which(net$balanceStructure$edges[[i]]$internalEdges$endn
                        %in% lambda.prime)] <-
      balance.penalty[i]
    costv <-
      c(costv, internal.subv, rep(0, nrow(
        net$balanceStructure$edges[[i]]$externalEdges
      )))
  }
  
  skeleton <- costSkeleton(net)
  skeleton$balance <- costv
  return(skeleton)
}






#add VR edges

#functions to construct distance masks for these extra structures

#think through how to handle IDs

#function to run iterative solutions and package up output.

#EASY WRAPPERS
#e.g. distanceVsDistance <- function(dist1, dist2, nsols)
# distanceVsBalance <- function(dist1, treat.vals, ctrl.vals, nsols)
# distanceVsExclusion <- function(dist1, nsols)
# balanceVsExclusion <- function(treat.vals, ctrl.vals, nsols)
# distanceVsVR <- function(dist, ubound.ratio, lbound.ratio = NULL, nsols)

#' Extract edges from the network
#'
#' @param net the network representation
#'
#' @return the list of edges
extractEdges <- function(net) {
  if (net$dense) {
    startn <- rep(net$treatedNodes, length(net$controlNodes))
    endn <-
      rep(net$controlNodes, rep(length(net$treatedNodes),
                                length(net$controlNodes)))
    ucap <-
      rep(1, length(net$treatedNodes) * length(net$controlNodes))
  } else {
    startn <- rep(net$treatedNodes, laply(net$edgeStructure, length))
    endn <- unlist(net$edgeStructure)
    ucap <- rep(1, length(endn))
  }
  edge.frame <-
    data.frame('startn' = startn,
               'endn' = endn,
               'ucap' = ucap)
  
  if ('exclusionEdges' %in% names(net)) {
    new.frame <- net$exclusionEdges
    new.frame$ucap <- 1
    edge.frame <- rbind(edge.frame, new.frame)
  }
  if ('balanceStructure' %in% names(net)) {
    new.frame <- NULL
    for (i in 1:length(net$balanceStructure$edges)) {
      x <- net$balanceStructure$edges[[i]]
      external <- x$externalEdges
      colnames(external) <- colnames(x$internalEdges)
      stub.frame <- rbind(x$internalEdges, external)
      if (is.null(new.frame)) {
        new.frame <- stub.frame
      } else {
        new.frame <- rbind(new.frame, stub.frame)
      }
    }
    edge.frame <- rbind(edge.frame, new.frame)
    
    new.frame <-
      data.frame('startn' = net$balanceStructure$nodes[[1]]$lambda.double.prime)
    new.frame$endn <- net$sink
    new.frame$ucap <- Inf
    edge.frame <- rbind(edge.frame, new.frame)
  } else {
    new.frame <- data.frame('startn' = net$controlNodes)
    new.frame$endn <- net$sink
    new.frame$ucap <- 1
    edge.frame <- rbind(edge.frame, new.frame)
  }
  return(edge.frame)
}

#' Extract the supply nodes from the net
#'
#' @param net the network representation
#'
#' @return the vector of the supply nodes
extractSupply <- function(net) {
  b <- rep(0, net$totalNodes)
  b[1:length(net$treatedNodes)] <- 1
  b[net$sink] <- -sum(b)
  return(b)
}


#' Solve the network flow problem - basic version
#'
#' @param net the network representation
#' @param f1.list the list of the first objective functions values for each node
#' @param f2.list the list of the second objective functions values for each
#'   node
#' @param rho the penalty coefficient
#' @param tol the tolerance value for precision
#'
#' @return the solution represented in a named list
solveP <- function(net, f1.list, f2.list, rho, tol = 1e-5) {
  if (!inherits(net, "netFlowMatch"))
    stop("Argument net is not of class netFlowMatch")
  
  #turn network into the format taken by RELAX code
  edge.info <- extractEdges(net)
  old.format <-
    list(
      'startn' = edge.info$startn,
      'endn' = edge.info$endn,
      'ucap' = edge.info$ucap
    )
  old.format$b <- extractSupply(net)
  
  #convert infinite capacities to large finite ones
  old.format$ucap[!is.finite(old.format$ucap)] <-
    sum(pmax(old.format$b, 0))
  
  #set network costs via f1, f2, and rho
  old.format$cost <-
    rowSums(sapply(f1.list, flattenSkeleton)) +
    rowSums(sapply(f2.list, flattenSkeleton)) * rho
  
  
  #convert costs to integers of appropriate precision, if necessary
  #optmatch developers for ideas about how to do this nicely
  cost <- old.format$cost
  if (any(cost != round(cost))) {
    intcost <- round(cost / tol)
    #use smaller scaling factor if possible
    searchtol <- 10 ^ (-c(1:floor(log10(
      .Machine$integer.max
    ))))
    searchtol <- searchtol[searchtol > tol]
    for (newtol in searchtol) {
      new.intcost <- round(intcost * tol / newtol)
      if (any(new.intcost != intcost * tol / newtol))
        break
      tol <- newtol
      intcost <- new.intcost
    }
    cost <- intcost
  }
  if (any(is.na(as.integer(cost)))) {
    stop('Integer overflow in edge costs!  Use smaller costs or penalties.')
  }
  
  old.format$cost <- cost
  o <- callrelax(old.format)
  if (o$feasible == 0) {
    stub.message <-
      'Match is infeasible or edge costs are too large for RELAX to process! Consider checking network, reducing edge costs, or raising tolerance.'
    stop(paste(stub.message, '.', sep = ''))
  }
  
  #return optimal flow solution
  return(list('x' = o$x, 'net' = old.format))
}




#' Call relax on the network
#'
#' @description this function is copied from the rcbalance package
#' @param net the network structure
#' @param solver (optional) the solver; by default, "rlemon"
#'
#' @return list of the result from the call to relax solver
callrelax <- function (net, solver = 'rlemon') {
  startn <- net$startn
  endn <- net$endn
  ucap <- net$ucap
  b <- net$b
  cost <- net$cost
  stopifnot(length(startn) == length(endn))
  stopifnot(length(startn) == length(ucap))
  stopifnot(length(startn) == length(cost))
  stopifnot(min(c(startn, endn)) >= 1)
  stopifnot(max(c(startn, endn)) <= length(b))
  stopifnot(all(startn != endn))
  
  nnodes <- length(b)
  if (solver == 'rrelaxiv') {
    if (requireNamespace('rrelaxiv', quietly = TRUE)) {
      rout <- rrelaxiv::RELAX_IV(
        startnodes = as.integer(startn),
        endnodes = as.integer(endn),
        arccosts = as.integer(cost),
        arccapacity = as.integer(ucap),
        supply = as.integer(b)
      )
      return.obj <- list(crash = 0,
                         feasible = !all(rout == 0),
                         x = rout)
      return(return.obj)
    } else {
      solver = 'rlemon'
      warning('Package rrelaxiv not available, using rlemon instead.')
    }
  }
  if (solver == 'rlemon') {
    lout <- rlemon::MinCostFlow(
      arcSources = as.integer(startn),
      arcTargets = as.integer(endn),
      arcCapacities = as.integer(ucap),
      arcCosts = as.integer(cost),
      nodeSupplies = as.integer(b),
      numNodes = max(c(startn, endn)),
      algorithm = 'CycleCancelling'
    )
    return.obj <- list(crash = 0,
                       feasible = !all(lout[[1]] == 0),
                       x = lout[[1]])
    return(return.obj)
  } else{
    stop('Argument to solver not recognized:
         please use one of rlemon and rrelaxiv')
  }
}

#' Solve the network flow problem - twoDistMatch
#'
#' @param net the network representation
#' @param f1.list the list of the first objective functions values for
#'   each node
#' @param f2.list the list of the second objective functions values for
#'   each node
#' @param f3.list the list of the third objective functions values for
#'   each node
#' @param rho1 the penalty coefficient for the second objective
#' @param rho2 the penalty coefficient for the third objective
#' @param tol the tolerance value for precision
#'
#' @return the solution represented in a named list
solveP1 <-
  function(net,
           f1.list,
           f2.list,
           f3.list,
           rho1,
           rho2 = 0,
           tol = 1e-5) {
    if (!inherits(net, "netFlowMatch"))
      stop("Argument net is not of class netFlowMatch")
    
    #turn network into the format taken by RELAX code
    edge.info <- extractEdges(net)
    old.format <-
      list(
        'startn' = edge.info$startn,
        'endn' = edge.info$endn,
        'ucap' = edge.info$ucap
      )
    old.format$b <- extractSupply(net)
    
    #convert infinite capacities to large finite ones
    old.format$ucap[!is.finite(old.format$ucap)] <-
      sum(pmax(old.format$b, 0))
    
    #set network costs via f1, f2, and rho
    old.format$cost <-
      rowSums(sapply(f1.list, flattenSkeleton)) +
      rowSums(sapply(f2.list, flattenSkeleton)) * rho1 +
      rowSums(sapply(f3.list, flattenSkeleton)) * rho2
    
    #convert costs to integers of appropriate precision, if necessary
    #optmatch developers for ideas about how to do this nicely
    cost <- old.format$cost
    if (any(cost != round(cost))) {
      intcost <- round(cost / tol)
      #use smaller scaling factor if possible
      searchtol <- 10 ^ (-c(1:floor(log10(
        .Machine$integer.max
      ))))
      searchtol <- searchtol[searchtol > tol]
      for (newtol in searchtol) {
        new.intcost <- round(intcost * tol / newtol)
        if (any(new.intcost != intcost * tol / newtol))
          break
        tol <- newtol
        intcost <- new.intcost
      }
      cost <- intcost
    }
    if (any(is.na(as.integer(cost)))) {
      stop('Integer overflow in edge costs!  Use smaller costs or penalties.')
    }
    
    old.format$cost <- cost
    o <- callrelax(old.format)
    if (o$feasible == 0) {
      stub.message <-
        'Match is infeasible or edge costs are too large for RELAX to process! Consider checking network, reducing edge costs, or raising tolerance.'
      stop(paste(stub.message, '.', sep = ''))
    }
    
    #return optimal flow solution
    return(list('x' = o$x, 'net' = old.format))
  }
