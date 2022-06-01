#' Generate penalty coefficient pairs
#' @description An internal helper function used for automatically generating
#' the set of rho values used for grid search in exploring the Pareto optimal
#' set of solutions.
#' @param paircosts.list a vector of pair-wise distance.
#' @param rho.max.factor a numeric value indicating the maximal rho values.
#' @param rho1old a vector of numeric values of rho1 used before.
#' @param rho2old a vector of numeric values of rho2 used before.
#' @param rho.min smallest rho value to consider.
#'
#' @return a vector of pairs of rho values for future search.
rho_proposition <-
  function(paircosts.list,
           rho.max.factor = 10,
           rho1old,
           rho2old,
           rho.min = 1e-2) {
    max.dist <- max(paircosts.list)
    min.dist <- min(paircosts.list)
    
    rho1 <- c(
      rho.min,
      min.dist / 2,
      min.dist,
      quantile(paircosts.list)[c(2, 4)],
      max.dist,
      rho.max.factor * max.dist
    )
    rho2 <- seq(0, max(rho1), max(rho1) / 5)
    rho1 <- rho1[rho1 >= 0]
    rho2 <- rho2[rho2 >= 0]
    rho1old <- rho1old[rho1old >= 0]
    rho2old <- rho2old[rho2old >= 0]
    
    if (is.null(rho1old) & is.null(rho2old)) {
      result <- generate_rhos(rho1, rho2)
    } else if (is.null(rho1old)) {
      result <- generate_rhos(rho1, rho2old)
    } else if (is.null(rho2old)) {
      result <- generate_rhos(rho1old, rho2)
    } else{
      result <- generate_rhos(rho1old, rho2old)
    }
    return(result)
    
  }

#' Generate rho pairs
#' @description An internal helper function that generates the set of rho value
#' pairs used for matching. This function is used when exploring the Pareto
#' optimality of the solutions to the multi-objective optimization in matching.
#'
#' @param rho1.list a vector of rho 1 values
#' @param rho2.list a vector of rho 2 values
#'
#' @return a vector of (rho1, rho2) pairs
generate_rhos <- function(rho1.list, rho2.list) {
  result <- list()
  ind = 1
  for (i in 1:length(rho1.list)) {
    for (j in 1:length(rho2.list)) {
      result[[ind]] <- c(rho1.list[i], rho2.list[j])
      ind = ind + 1
    }
  }
  return(result)
  
}


#' An internal helper function that transforms the output from the RELAX
#' algorithm to a data structure that is more interpretable for the output of
#' the main matching function
#'
#' @param out.elem a named list whose elements are: (1) the net structure (2)
#'   the edge weights of pair-wise distance (3) the edge weights of marginal
#'   balance (4) the list of rho value pairs (5) the named list of solutions
#'   from the RELAX algorithm
#' @param already.done a factor indicating the index of matches already been
#'   transformed
#' @param prev.obj an object of previously transformed matches
#'
#' @return a named list with elements containing matching information useful for
#'   the main matching function
obj.to.match <-
  function(out.elem,
           already.done = NULL,
           prev.obj = NULL) {
    tcarcs <- length(unlist(out.elem$net$edgeStructure))
    edge.info <- extractEdges(out.elem$net)
    one.sol <- function(sol) {
      x <- sol[1:tcarcs]
      match.df <-
        data.frame(
          treat = as.factor(edge.info$startn[1:tcarcs]),
          x = x,
          control = edge.info$endn[1:tcarcs]
        )
      matched.or.not <-
        daply(match.df, .(match.df$treat), function(treat.edges)
          c(as.numeric(as.character(
            treat.edges$treat[1]
          )),
          sum(treat.edges$x)), .drop_o = FALSE)
      if (any(matched.or.not[, 2] == 0)) {
        match.df <-
          match.df[-which(match.df$treat %in% 
                            matched.or.not[which(matched.or.not[,2] == 0), 1]),]
      }
      match.df$treat <- as.factor(as.character(match.df$treat))
      matches <-
        as.matrix(daply(match.df, .(match.df$treat), function(treat.edges)
          treat.edges$control[treat.edges$x ==
                                1], .drop_o = FALSE))
      matches - length(out.elem$net$treatedNodes)
    }
    if (is.null(already.done))
      return(llply(out.elem$solutions, one.sol))
    new.ones <- setdiff(1:length(out.elem$solutions), already.done)
    out.list <- list()
    out.list[already.done] <- prev.obj
    out.list[new.ones] <-
      llply(out.elem$solutions[new.ones], one.sol)
    return(out.list)
  }



#' Fit propensity scores using logistic regression.
#'
#' @param df dataframe that contains a column named "treat", the treatment
#'   vector, and columns of covariates specified.
#' @param covs factor of column names of covariates used for fitting a
#'   propensity score model.
#'
#' @return vector of estimated propensity scores (on the probability scale).
getPropensityScore <- function(df, covs) {
  pscoreModel <- glm(treat ~ .,
                     data = df[c('treat', covs)], family = binomial("logit"))
  return(predict(pscoreModel, type = "response"))
}


#' Generate a factor for exact matching.
#'
#' @param dat dataframe containing all the variables in exactList
#' @param exactList factor of names of the variables on which we want exact or
#'   close matching.
#'
#' @return factor on which to match exactly, with labels given by concatenating
#'   labels for input variables.
getExactOn <- function(dat, exactList) {
  if (length(exactList) == 1) {
    exactOn = dat[, exactList[1]]
  } else {
    exactOn = paste(dat[, exactList[1]], dat[, exactList[2]])
    if (length(exactList) >= 3) {
      for (i in 3:length(exactList)) {
        exactOn = paste(exactOn, dat[, exactList[i]])
      }
    }
  }
  return(exactOn)
}

#' Optimal tradeoffs among distance, exclusion and marginal imbalance
#' @description Explores tradeoffs among three important objective functions in
#'   an optimal matching problem:the sum of covariate distances within matched
#'   pairs, the number of treated units included in the match, and the marginal
#'   imbalance on pre-specified covariates (in total variation distance).
#'
#' @family main matching function
#'
#' @param df data frame that contain columns indicating treatment, outcome and
#'   covariates.
#' @param treatCol character of name of the column indicating treatment
#'   assignment.
#' @param myBalCol character of column name of the variable on which to evaluate
#'   marginal balance.
#' @param rhoExclude (optional) numeric vector of values of exclusion penalty.
#'   Default value is c(1).
#' @param rhoBalance (optional) factor of values of marginal balance penalty.
#'   Default value is c(1,2,3).
#' @param distMatrix (optional) a matrix that specifies the pair-wise distances
#'   between any two objects.
#' @param distList (optional) character vector of variable names used for
#'   calculating within-pair distance.
#' @param exactlist (optional) character vector, variable names that we want
#'   exact matching on; NULL by default.
#' @param propensityCols (optional) character vector, variable names on which to
#'   fit a propensity score (to supply a caliper).
#' @param pScores (optional) character, giving the variable name for the fitted
#'   propensity score.
#' @param ignore (optional) character vector of variable names that should be
#'   ignored when constructing the internal matching. NULL by default.
#' @param maxUnMatched (optional) numeric, the maximum proportion of unmatched
#'   units that can be accepted; default is 0.25.
#' @param caliperOption (optional) numeric, the propensity score caliper value
#'   in standard deviations of the estimated propensity scores; default is NULL,
#'   which is no caliper.
#' @param toleranceOption (optional) numeric, tolerance of close match distance;
#'   default is 1e-2.
#' @param maxIter (optional) integer,  maximum number of iterations to use in
#'   searching for penalty combintions that improve the matching; default is 0.
#' @param rho.max.f (optional) numeric, the scaling factor used in proposal for
#'   rhos; default is 10.
#' @importFrom plyr laply
#' @return a named list whose elements are: * "rhoList": list of penalty
#'   combinations for each match * "matchList": list of matches indexed by
#'   number 
#'   * "treatmentCol": character of treatment variable 
#'   * "covs":
#'   character vector of names of the variables used for calculating within-pair
#'   distance 
#'   * "exactCovs": character vector of names of variables that we want
#'   exact or close match on * "idMapping": numeric vector of row indices for
#'   each observation in the sorted data frame for internal use 
#'   * "stats": data
#'   frame of important statistics (total variation distance) for variable on
#'   which marginal balance is measured 
#'   * "b.var": character, name of variable
#'   on which marginal balance is measured * "dataTable": data frame sorted by
#'   treatment value 
#'   * "t": a treatment vector 
#'   * "df": the original dataframe
#'   input by the user 
#'   * "pair_cost1": list of pair-wise distance sum using the
#'   first distance measure 
#'   * "pair_cost2": list of pair-wise distance sum using
#'   the second distance measure (left NULL since only one distance measure is
#'   used here). 
#'   * "version": (for internal use) the version of the matching
#'   function called; "Basic" indicates the matching comes from distBalMatch and
#'   "Advanced" from twoDistMatch. 
#'   * "fPair": a vector of values for the first
#'   objective function; it corresponds to the pair-wise distance sum according
#'   to the first distance measure. 
#'   * "fExclude": a vector of values for the
#'   second objective function; it corresponds to the number of treated units
#'   being unmatched. 
#'   * "fMarginal": a vector of values for the third objective
#'   function; it corresponds to the marginal balanced distance for the
#'   specified variable(s).
#'
#' @details Matched designs generated by this function are Pareto optimal for
#'   the three objective functions.  The degree of relative emphasis among the
#'   three objectives in any specific solution is controlled by the penalties,
#'   denoted by Greek letter rho. Larger values of `rhoExclude` corresponds to
#'   increased emphasis on retaining treated units (all else being equal), while
#'   larger values of `rhoBalance` corresponds to increased emphasis on marginal
#'   balance. Additional details: 
#'   * Users may either specify their own distance
#'   matrix via the `distMatrix` argument or ask the function to create a
#'   Mahalanobis distance matrix internally on a set of covariates specified by
#'   the `distList` argument; if neither argument is specified an error will
#'   result.  User-specified distance matrices should have row count equal to
#'   the number of treated units and column count equal to the number of
#'   controls. 
#'   * If the `caliperOption` argument is specified, a propensity
#'   score caliper will be imposed, forbidding matches between units more than a
#'   fixed distance apart on the propensity score.  The caliper will be based
#'   either on a user-fit propensity score, identified in the input dataframe by
#'   argument `pScores`, or by an internally-fit propensity score based on
#'   logistic regression against the variables named in `psoreCols`.  If
#'   `caliperOption` is non-NULL and neither of the other arguments is specified
#'   an error will result. 
#'   * `toleranceOption` controls the precision at which
#'   the objective functions is evaluated. When matching problems are especially
#'   large or complex it may be necessary to increase toleranceOption in order
#'   to prevent integer overflows in the underlying network flow solver;
#'   generally this will be suggested in appropariate warning messages. 
#'   * While
#'   by default tradeoffs are only assessed at penalty combinations provided by
#'   the user, the user may ask for the algorithm to search over additional
#'   penalty values in order to identify additional Pareto optimal solutions.
#'   `rho.max.f` is a multiplier applied to initial penalty values to discover
#'   new solutions, and setting it larger leads to wider exploration; similarly,
#'   `maxIter` controls how long the exploration routine runs, with larger
#'   values leading to more exploration.
#'
#' @export
#' @examples
#' data("lalonde", package="cobalt")
#' psCols <- c("age", "educ", "married", "nodegree")
#' treatVal <- "treat"
#' responseVal <- "re78"
#' pairDistVal <- c("age", "married", "educ", "nodegree")
#' exactVal <- c("educ")
#' myBalVal <- c("race")
#' r1s <- c( 0.1, 0.3, 0.5, 0.7, 0.9,1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7)
#' r2s <- c(0.01)
#' matchResult <- distBalMatch(df=lalonde, treatCol=treatVal, myBalCol=myBalVal,
#' rhoExclude =r1s, rhoBalance=r2s,
#' distList=pairDistVal, exactlist=exactVal,
#' propensityCols = psCols,ignore = c(responseVal), maxUnMatched = 0.1,
#' caliperOption=NULL, toleranceOption=1e-1, maxIter=0, rho.max.f = 10)
distBalMatch <-
  function(df,
           treatCol,
           myBalCol,
           rhoExclude = c(1),
           rhoBalance = c(1, 2, 3),
           distMatrix = NULL,
           distList = NULL,
           exactlist = NULL,
           propensityCols = NULL,
           pScores = NULL,
           ignore = NULL,
           maxUnMatched = 0.25,
           caliperOption = NULL,
           toleranceOption = 1e-2,
           maxIter = 0,
           rho.max.f = 10) {
    if (is.null(treatCol) | is.null(df) | is.null(myBalCol)) {
      stop("You must input some data. See df, treatCol and myBalCol arguments.")
    }
    rho1 <- rhoExclude
    rho2 <- rhoBalance
    if (!is.null(distMatrix)) {
      distMatrix <- distanceFunctionHelper(df[[treatCol]], distMatrix)
    }
    
    
    if (is.null(exactlist)) {
      df["pseudo-exact"] = rep(1, nrow(df))
      exactlist = c("pseudo-exact")
    }
    
    
    ## 0. Data preprocessing
    result <- structure(list(), class = 'multiObjMatch')
    dat = df
    result$dataTable = dat
    dat$originalID = 1:nrow(df)
    dat$tempSorting = df[, treatCol]
    dat = dat[order(-dat$tempSorting), ]
    rownames(dat) <- as.character(1:nrow(dat))
    columnNames <- colnames(df)
    xNames <-
      columnNames[!columnNames %in% 
                    c(treatCol, ignore, "originalID", "tempSorting")]
    
    #dat['response'] = df[responseCol]
    dat['treat'] = df[treatCol]
    
    ## 1. Fit a propensity score model if the propensity score is not passed in
    if (is.null(pScores)) {
      if (is.null(propensityCols)) {
        pscore = getPropensityScore(dat, xNames)
      } else {
        pscore = getPropensityScore(dat, propensityCols)
      }
    } else{
      pscore = as.numeric(list(dat[pScores]))
    }
    
    
    ## 2. Construct nets
    base.net <- netFlowMatch(dat$treat)
    
    ## 3. Define distance and cost
    caliperType = 'none'
    if (!is.null(caliperOption)) {
      caliperType = "user"
    } else {
      calipterOption = 0
    }
    
    
    
    
    if (is.null(distMatrix)) {
      dist.i <- build.dist.struct(
        z = dat$treat,
        X = dat[distList],
        dist.type = "Mahalanobis",
        exact = getExactOn(dat, exactlist),
        calip.option = caliperType,
        calip.cov = pscore,
        caliper = caliperOption
      )
      
      distEuc.i <- build.dist.struct(
        z = dat$treat,
        X = dat[distList],
        dist.type = "Euclidean",
        exact = getExactOn(dat, exactlist),
        calip.option = caliperType,
        calip.cov = pscore,
        caliper = caliperOption
      )
    } else {
      dist.i <-
        build.dist.struct(
          z = dat$treat,
          X = dat[distList] ,
          distMat = distMatrix,
          exact = getExactOn(dat, exactlist),
          dist.type = "User",
          calip.option = caliperType,
          calip.cov = pscore,
          caliper = caliperOption
        )
      distEuc.i <- NULL
    }
    
    
    
    
    
    
    net.i <- addExclusion(makeSparse(base.net, dist.i))
    if (!is.null(distEuc.i)) {
      netEuc.i <- addExclusion(makeSparse(base.net, distEuc.i))
    }
    my.bal <-
      apply(dat[, match(myBalCol, colnames(dat)), drop = FALSE], 1, function(x)
        paste(x, collapse = "."))
    net.i <- addBalance(net.i, treatedVals =
                          my.bal[dat$treat == 1], 
                        controlVals = my.bal[dat$treat == 0])
    if (!is.null(distEuc.i)) {
      netEuc.i <- addBalance(netEuc.i,
                             treatedVals =
                               my.bal[dat$treat == 1],
                             controlVals = my.bal[dat$treat == 0])
    }
    paircost.v <- unlist(dist.i)
    if (!is.null(distEuc.i)) {
      paircostEuc.v <- unlist(distEuc.i)
    }
    rho.max.factor = rho.max.f
    ### Create edge costs for each of the objective functions we will trade off
    pair.objective.edge.costs <- pairCosts(dist.i, net.i)
    excl.objective.edge.costs <-  excludeCosts(net.i, 1)
    bal.objective.edge.costs <-
      balanceCosts(net.i, max(paircost.v) * rho.max.factor)
    if (!is.null(distEuc.i)) {
      pairEuc.objective.edge.costs <- pairCosts(distEuc.i, netEuc.i)
    }
    f1list = list(pair.objective.edge.costs)
    f2list = list(excl.objective.edge.costs)
    f3list = list(bal.objective.edge.costs)
    if (!is.null(distEuc.i)) {
      f4list = list(pairEuc.objective.edge.costs)
    }
    
    
    ## Propose possible rho values for multi-objective optimization problem
    rho_list <- list()
    rho_list <- rho_proposition(paircost.v, rho.max.f, rho1, rho2,
                                rho.min = toleranceOption)
    
    # #### EXPENSIVE PART ####
    #
    solutions <- list()
    solution.nets <- list()
    rho_counts = 1
    for (rho in rho_list) {
      temp.sol <- solveP1(
        net.i,
        f1list,
        f2list,
        f3list,
        rho1 = rho[1],
        rho2 = rho[2],
        tol = toleranceOption
      )
      
      solutions[[as.character(rho_counts)]] <- temp.sol$x
      solution.nets[[as.character(rho_counts)]] <- temp.sol$net
      print(paste('Matches finished for rho1 = ', rho[1], 
                  " and rho2 = ", rho[2]))
      rho_counts = rho_counts + 1
    }
    
    match.list <-
      obj.to.match(
        list(
          'net' = net.i,
          'costs' = pair.objective.edge.costs,
          'balance' = bal.objective.edge.costs,
          'rho.v' = 1:length(rho_list),
          'solutions' = solutions
        )
      )
    match.list <- match.list[order(as.numeric(names(match.list)))]
    
    ## remove the matching results that result in zero matches
    temp.length <- laply(match.list, nrow)
    temp.names <- names(match.list)
    nonzero.names <- temp.names[temp.length != 0]
    zero.names <- temp.names[temp.length == 0]
    
    if (length(zero.names) == length(temp.names)) {
      stop(
        'Error: Initial rho values result in non-matching. 
        Please try a new set of initial rho values.'
      )
    }
    match.list <-
      match.list[names(match.list) %in% nonzero.names == TRUE]
    
    
    ## Iterative Search for Best Rhos (if necessary)
    numIter = 1
    
    
    tempRhos <- as.numeric(names(match.list))
    percentageUnmatched <-
      1 - as.vector(laply(match.list, nrow) / sum(dat[, treatCol]))
    minInd <- which.min(percentageUnmatched)
    bestRhoInd <- tempRhos[minInd]
    bestPercentageSoFar <- as.vector(percentageUnmatched)[minInd]
    bestRho1 <- rho_list[[bestRhoInd]][1]
    bestRho2 <- rho_list[[bestRhoInd]][2]
    ind = length(rho_list) + 1
    while ((bestPercentageSoFar > maxUnMatched) &
           (numIter <= maxIter)) {
      if (is.finite(log10(bestRho1)) & is.finite(log10(bestRho2))) {
        proposedNewRho1s <-
          c(runif(2, 0, 2), 0.5 * 
              (-c((abs(rnorm(1)) * 0.01 + bestRho1 / 10):floor(log10(bestRho1))
          )), 0.5 * (-c(
            abs(rnorm(1)) * 0.01 + bestRho1:log10(bestRho1 + 
                                                    0.01 * abs(rnorm(1))) *
              10
          )))
        
        proposedNewRho2s <-
          c(runif(2, 0, 2), 0.5 * 
              (-c((abs(rnorm(1)) * 0.01 + bestRho2 / 10):floor(log10(bestRho2))
          )), 0.5 * (-c(
            abs(rnorm(1)) * 0.01 + bestRho2:log10(bestRho2 + 
                                                    0.01 * abs(rnorm(1))) *
              10
          )))
      } else {
        proposedNewRho1s <-
          c(rnorm(1) * 1e-1, rnorm(1) * 1e2, rnorm(1) * 1e3, runif(2, 0, 2))
        proposedNewRho2s <-
          c(rnorm(1) * 1e-1, rnorm(1) * 1e2, rnorm(1) * 1e3, runif(2, 0, 2))
      }
      
      proposedNewRho1s <- proposedNewRho1s[proposedNewRho1s >= 0]
      proposedNewRho2s <- proposedNewRho2s[proposedNewRho2s >= 0]
      
      for (r1 in proposedNewRho1s) {
        rho_list[[ind]] = c(r1, bestRho2)
        temp.sol <- solveP1(
          net.i,
          f1list,
          f2list,
          f3list,
          rho1 = r1,
          rho2 = bestRho2,
          tol = toleranceOption
        )
        
        solutions[[as.character(ind)]] <- temp.sol$x
        solution.nets[[as.character(ind)]] <- temp.sol$net
        print(paste('Matches finished for rho1 = ', 
                    r1, " and rho2 = ", bestRho2))
        if (is.null(temp.sol$x)) {
          print(
            paste(
              'However, rho1 = ',
              r1,
              " and rho2 = ",
              bestRho2,
              " results in empty matching."
            )
          )
        }
        ind = ind + 1
      }
      
      for (r2 in proposedNewRho2s) {
        rho_list[[ind]] = c(bestRho1, r2)
        temp.sol <- solveP1(
          net.i,
          f1list,
          f2list,
          f3list,
          rho1 = bestRho1,
          rho2 = r2,
          tol = toleranceOption
        )
        
        solutions[[as.character(ind)]] <- temp.sol$x
        solution.nets[[as.character(ind)]] <- temp.sol$net
        print(paste('Matches finished for rho1 = ', 
                    bestRho1, " and rho2 = ", r2))
        if (is.null(temp.sol$x)) {
          print(
            paste(
              'However, rho1 = ',
              bestRho1,
              " and rho2 = ",
              r2,
              " results in empty matching."
            )
          )
        }
        ind = ind + 1
      }
      
      match.list <-
        obj.to.match(
          list(
            'net' = net.i,
            'costs' = pair.objective.edge.costs,
            'balance' = bal.objective.edge.costs,
            'rho.v' = 1:length(solutions),
            'solutions' = solutions
          )
        )
      match.list <- match.list[order(as.numeric(names(match.list)))]
      
      ## remove the matching results that result in zero matches
      temp.length <- laply(match.list, nrow)
      temp.names <- names(match.list)
      nonzero.names <- temp.names[temp.length != 0]
      match.list <-
        match.list[names(match.list) %in% nonzero.names == TRUE]
      
      
      tempRhos <- as.numeric(names(match.list))
      percentageUnmatched <-
        1 - as.vector(laply(match.list, nrow) / sum(dat[, treatCol]))
      oldMinInd <- minInd
      minInd <- which.min(percentageUnmatched[-oldMinInd])
      if (minInd >= oldMinInd) {
        minInd = minInd + 1
      }
      bestRhoInd <- tempRhos[minInd]
      bestPercentageSoFar <- as.vector(percentageUnmatched)[minInd]
      if (sum(percentageUnmatched == bestPercentageSoFar) > 1) {
        bestRhoInd <-
          sample(which(percentageUnmatched == bestPercentageSoFar),
                 size = 1)
      }
      bestRho1 <- rho_list[[bestRhoInd]][1]
      bestRho2 <- rho_list[[bestRhoInd]][2]
      
      numIter = numIter + 1
    }
    my.stats1 <-
      t(
        laply(
          match.list,
          descr.stats_general,
          df = dat,
          treatCol = treatCol,
          b.vars = myBalCol,
          pair.vars = distList,
          extra = TRUE
        )
      )
    
    
    solutions.old <- solutions
    solution.nets.old <- solution.nets
    
    pair_cost_sum <- c()
    pair_cost_euc_sum <- c()
    f1_sum <- c()
    f2_sum <- c()
    f3_sum <- c()
    total_treated = sum(df[[treatCol]])
    i = 1
    for (ind in names(match.list)) {
      x <- solutions[[ind]]
      pair_cost_sum[[ind]] <- sum(paircost.v * x[1:length(paircost.v)])
      if (!is.null(distEuc.i)) {
        pair_cost_euc_sum[[ind]] <-
          sum(paircostEuc.v * x[1:length(paircostEuc.v)])
      }
      unmatch_ind = total_treated - length(match.list[[ind]])
      
      
      f1_sum[[ind]] <- sum(paircost.v * x[1:length(paircost.v)])
      f2_sum[[ind]] <- unmatch_ind
      f3_sum[[ind]] <- my.stats1[1, i]
      i = i + 1
      
      
    }
    
    rho_list_final = c()
    for (ind in names(match.list)) {
      rho_list_final[[ind]] <- rho_list[[as.numeric(ind)]]
    }
    
    result$rhoList <- rho_list_final
    result$matchList <-  match.list
    ## Store some information about treatment column, balance covariates and
    ## exact match column
    result$treatmentCol = treatCol
    result$covs = distList
    if (is.null(distList)) {
      result$covs = xNames
    }
    result$exactCovs = exactlist
    result$idMapping = dat$originalID
    result$stats = my.stats1
    result$b.var = myBalCol
    result$df = df
    result$pair_cost1 = pair_cost_sum
    result$pair_cost2 = pair_cost_euc_sum
    result$version = 'Basic'
    result$fPair <- f1_sum
    result$fExclude <- f2_sum
    result$fMarginal <- f3_sum
    return(result)
  }




#' Optimal tradeoffs among two distances and exclusion
#' @description Explores tradeoffs among three objective functions in
#' multivariate matching: sums of two different user-specified  covariate
#' distances within matched pairs, and the number of treated units included in
#' the match.
#'
#' @family main matching function
#' @param dType1 One of ("Euclidean", "Mahalanobis", "user") indicating the type
#'   of distance that are used for the first distance objective functions. NULL
#'   by default.
#' @param dType2 One of ("Euclidean", "Mahalanobis", "user")  charactor
#'   indicating the type of distance that are used for the second distance
#'   objective functions. NULL by default.
#' @param dMat1 (optional) matrix object that represents the distance matrix
#'   using the first distance measure; `dType` must be passed in as "user" if
#'   dMat is non-empty
#' @param dMat2 (optional) matrix object that represents the distance matrix
#'   using the second distance measure; `dType1` must be passed in as "user" if
#'   dMat is non-empty
#' @param df (optional) data frame that contain columns indicating treatment,
#'   outcome and covariates
#' @param treatCol (optional) character, name of the column indicating treatment
#'   assignment.
#' @param distList1 (optional) character vector names of the variables used for
#'   calculating covariate distance using first distance measure specified by
#'   dType
#' @param distList2 (optional) character vector, names of the variables used for
#'   calculating covariate distance using second distance measure specified by
#'   dType1
#' @param rhoExclude (optional) numeric vector, penalty values associated with
#'   the distance specified by `dMat` or `dType`. Default value is c(1).
#' @param rhoDistance (optional) numeric vector, penalty values associated with
#'   the distance specified by `dMat1` or `dType1`. Default value is c(1,2,3).
#' @param myBalCol (optional) character, column name of the variable on which to
#'   evaluate balance.
#' @param exactlist (optional) character vector, names of the variables on which
#'   to match exactly; NULL by default.
#' @param propensityCols character vector, names of columns on which to fit a
#'   propensity score model.
#' @param pScores (optional) character, name of the column containing fitted
#'   propensity scores; default is NULL.
#' @param ignore (optional) character vector of variable names that should be
#'   ignored when constructing the internal matching. NULL by default.
#' @param maxUnMatched (optional) numeric, maximum proportion of unmatched units
#'   that can be accepted; default is 0.25.
#' @param caliperOption (optional) numeric, the propensity score caliper value
#'   in standard deviations of the estimated propensity scores; default is NULL,
#'   which is no caliper.
#' @param toleranceOption (optional) numeric, tolerance of close match distance;
#'   default is 1e-2.
#' @param maxIter (optional) integer,  maximum number of iterations to use in
#'   searching for penalty combintions that improve the matching; default is 0.
#' @param rho.max.f (optional) numeric, the scaling factor used in proposal for
#'   rhos; default is 10.
#'
#' @return a named list whose elements are: 
#'   * "rhoList": list of penalty
#'   combinations for each match 
#'   * "matchList": list of matches indexed by
#'   number 
#'   * "treatmentCol": character of treatment variable 
#'   * "covs":character vector of names of the variables used for calculating within-pair
#'   distance 
#'   * "exactCovs": character vector of names of variables that we want
#'   exact or close match on 
#'   * "idMapping": numeric vector of row indices for
#'   each observation in the sorted data frame for internal use 
#'   * "stats": data
#'   frame of important statistics (total variation distance) for variable on
#'   which marginal balance is measured 
#'   * "b.var": character, name of variable
#'   on which marginal balance is measured (left NULL since no balance
#'   constraint is imposed here). 
#'   * "dataTable": data frame sorted by treatment
#'   value 
#'   * "t": a treatment vector 
#'   * "df": the original dataframe input by the
#'   user 
#'   * "pair_cost1": list of pair-wise distance sum using the first
#'   distance measure 
#'   * "pair_cost2": list of pair-wise distance sum using the
#'   second distance measure 
#'   * "version": (for internal use) the version of the
#'   matching function called; "Basic" indicates the matching comes from
#'   distBalMatch and "Advanced" from twoDistMatch. 
#'   * "fDist1": a vector of
#'   values for the first objective function; it corresponds to the pair-wise
#'   distance sum according to the first distance measure. 
#'   * "fExclude": a
#'   vector of values for the second objective function; it corresponds to the
#'   number of treated units being unmatched. 
#'   * "fDist2": a vector of values for
#'   the third objective function; it corresponds to the pair-wise distance sum
#'   corresponds to the
#'
#' @details Matched designs generated by this function are Pareto optimal for
#' the three objective functions.  The degree of relative emphasis among the
#' three objectives in any specific solution is controlled by the penalties,
#' denoted by Greek letter rho. Larger values for the penalties associated with
#' the two distances correspond to increased emphasis close matching on these
#' distances, at the possible cost of excluding more treated units. Additional
#' details: 
#' * Users may either specify their own distance matrices (specifying
#' the `user` option in `dType1` and/or `dType2` and supplying arguments to
#' `dMat1` and/or `dMat2` respectively) or ask the function to create
#' Mahalanobis or Euclidean distances on sets of covariates specified by the
#' `distList1` and `distList2` arguments. If `dType1` or `dType2` is not
#' specified, if one of these is set to `user` and the corresponding `dMat1`
#' argument is not provided, or if one is NOT set to `user` and the
#' corresponding `distList1` argument is not provided, an error will result. 
#' * User-specified distance matrices passed to `dMat1` or `dMat2` should have row
#' count equal to the number of treated units and column count equal to the
#' number of controls. 
#' * If the `caliperOption` argument is specified, a
#' propensity score caliper will be imposed, forbidding matches between units
#' more than a fixed distance apart on the propensity score.  The caliper will
#' be based either on a user-fit propensity score, identified in the input
#' dataframe by argument `pScores`, or by an internally-fit propensity score
#' based on logistic regression against the variables named in `psoreCols`.  If
#' `caliperOption` is non-NULL and neither of the other arguments is specified
#' an error will result. 
#' * `toleranceOption` controls the precision at which the
#' objective functions is evaluated. When matching problems are especially large
#' or complex it may be necessary to increase toleranceOption in order to
#' prevent integer overflows in the underlying network flow solver; generally
#' this will be suggested in appropariate warning messages. 
#' * While by default
#' tradeoffs are only assessed at penalty combinations provided by the user, the
#' user may ask for the algorithm to search over additional penalty values in
#' order to identify additional Pareto optimal solutions. `rho.max.f` is a
#' multiplier applied to initial penalty values to discover new solutions, and
#' setting it larger leads to wider exploration; similarly, `maxIter` controls
#' how long the exploration routine runs, with larger values leading to more
#' exploration.
#'
#' @export
#'
#' @examples
#' x1 = rnorm(100, 0, 0.5)
#' x2 = rnorm(100, 0, 0.1)
#' x3 = rnorm(100, 0, 1)
#' x4 = rnorm(100, x1, 0.1)
#' r1ss <- seq(0.1,50, 10)
#' r2ss <- seq(0.1,50, 10)
#' x = cbind(x1, x2, x3,x4)
#' z = sample(c(rep(1, 50), rep(0, 50)))
#' e1 = rnorm(100, 0, 1.5)
#' e0 = rnorm(100, 0, 1.5)
#' y1impute = x1^2 + 0.6*x2^2 + 1 + e1
#' y0impute = x1^2 + 0.6*x2^2 + e0
#' treat = (z==1)
#' y = ifelse(treat, y1impute, y0impute)
#' names(x) <- c("x1", "x2", "x3", "x4")
#' df <- data.frame(cbind(z, y, x))
#' df$x5 <- 1
#' names(x) <- c("x1", "x2", "x3", "x4")
#' df <- data.frame(cbind(z, y, x))
#' df$x5 <- 1
#' d1 <- as.matrix(dist(df["x1"]))
#' d2 <- as.matrix(dist(df["x2"]))
#' idx <- 1:length(z)
#' treatedUnits <- idx[z==1]
#' controlUnits <- idx[z==0]
#' d1 <- as.matrix(d1[treatedUnits, controlUnits])
#' d2 <- as.matrix(d2[treatedUnits, controlUnits])
#' matchResult1 <- twoDistMatch(df, "z", "y", dMat1=d1, dType1= "User", dMat2=d2,
#' dType2="User", myBalCol=c("x5"), rhoExclude=r1ss, rhoDistance=r2ss,
#' propensityCols = c("x1")) 
twoDistMatch <-
  function(dType1 = "user",
           dType2 = "user",
           dMat1 = NULL,
           df = NULL,
           dMat2 = NULL,
           treatCol = NULL,
           distList1 = NULL,
           distList2 = NULL,
           rhoExclude = c(1),
           rhoDistance = c(1, 2, 3),
           myBalCol = NULL,
           exactlist = NULL,
           propensityCols = NULL,
           pScores = NULL,
           ignore = NULL,
           maxUnMatched = 0.25,
           caliperOption = NULL,
           toleranceOption = 1e-2,
           maxIter = 0,
           rho.max.f = 10) {
    rho1 <- rhoExclude
    rho2 <- rhoDistance
    
    
    
    if (is.null(df) && is.null(dMat1)) {
      stop("You must input some data for matching to be formed.")
    }
    
    if (is.null(df) && (!is.null(dMat1))) {
      df <- data.frame(c(rep(1, dim(dMat1)[1]), rep(0, dim(dMat1)[2])))
      colnames(df) <- c("pseudo-treat")
      treatCol <- "pseudo-treat"
    }
    
    
    if (!is.null(dMat1)) {
      dMat1 <- distanceFunctionHelper(df[[treatCol]], dMat1)
      dType1 <- "User"
    }
    if (!is.null(dMat2)) {
      dMat2 <- distanceFunctionHelper(df[[treatCol]], dMat2)
      dType2 <- "User"
    }
    
    if (is.null(exactlist)) {
      df["pseudo-exact"] = rep(1, nrow(df))
      exactlist = c("pseudo-exact")
      
    }
    
    if (!is.null(caliperOption)) {
      calip.option <- "user"
    } else {
      calip.option <- "none"
      caliperOption <- 0
    }
    
    
    ## 0. Data preprocessing
    dat = df
    dat$originalID = 1:nrow(df)
    dat$tempSorting = df[, treatCol]
    dat = dat[order(-dat$tempSorting), ]
    rownames(dat) <- as.character(1:nrow(dat))
    columnNames <- colnames(df)
    xNames <-
      columnNames[!columnNames %in% 
                    c(treatCol, ignore, "originalID", "tempSorting")]
    
    #dat['response'] = df[responseCol]
    dat['treat'] = df[treatCol]
    
    result <- structure(list(), class = 'multiObjMatch')
    
    ## 1. Fit a propensity score model if the propensity score is not passed in
    if (is.null(pScores)) {
      if (is.null(propensityCols)) {
        pscore = getPropensityScore(dat, xNames)
      } else {
        pscore = getPropensityScore(dat, propensityCols)
      }
    } else{
      pscore = as.numeric(list(dat[pScores]))
    }
    
    if (is.null(dMat1)) {
      distanceMatrix1 = NULL
    } else {
      distanceMatrix1 = dMat1[dat$originalID, dat$originalID]
      rownames(distanceMatrix1) <- 1:nrow(distanceMatrix1)
      colnames(distanceMatrix1) <- 1:nrow(distanceMatrix1)
    }
    
    if (is.null(dMat2)) {
      distanceMatrix2 = NULL
    } else {
      distanceMatrix2 = dMat2[dat$originalID, dat$originalID]
      rownames(distanceMatrix2) <- 1:nrow(distanceMatrix2)
      colnames(distanceMatrix2) <- 1:nrow(distanceMatrix2)
    }
    
    
    
    ## 2. Construct nets
    base.net <- netFlowMatch(dat$treat)
    
    ## 3. Define distance and cost
    
    dist.i <- build.dist.struct(
      z = dat$treat,
      X = dat[distList1],
      distMat = distanceMatrix1,
      dist.type = dType1,
      exact = getExactOn(dat, exactlist),
      calip.option = calip.option,
      calip.cov = pscore,
      caliper = caliperOption
    )
    
    distEuc.i <- build.dist.struct(
      z = dat$treat,
      X = dat[distList2],
      distMat = distanceMatrix2,
      dist.type = dType2,
      exact = getExactOn(dat, exactlist),
      calip.option = calip.option,
      calip.cov = pscore,
      caliper = caliperOption
    )
    
    
    #net.i <- addExclusion(makeSparse(base.net, mega_dist))
    net.i <- addExclusion(makeSparse(base.net, dist.i))
    #netEuc.i <- addExclusion(makeSparse(base.net, distEuc.i))
    my.bal <-
      apply(dat[, match(myBalCol, colnames(dat)), drop = FALSE], 1, function(x)
        paste(x, collapse = "."))
    net.i <- addBalance(net.i, treatedVals =
                          my.bal[dat$treat == 1], 
                        controlVals = my.bal[dat$treat == 0])
    #netEuc.i <- addBalance(netEuc.i, treatedVals = my.bal[dat$treat==1],
    #controlVals = my.bal[dat$treat==0]) return(net.i)
    paircost.v <- unlist(dist.i)
    paircostEuc.v <- unlist(distEuc.i)
    rho.max.factor = rho.max.f
    ### Create edge costs for each of the objective functions we will trade off
    pair.objective.edge.costs <- pairCosts(dist.i, net.i)
    excl.objective.edge.costs <-  excludeCosts(net.i, max(paircost.v))
    bal.objective.edge.costs <-
      balanceCosts(net.i, max(paircost.v) * rho.max.factor)
    pairEuc.objective.edge.costs <- pairCosts(distEuc.i, net.i)
    f1list = list(pair.objective.edge.costs)
    f2list = list(excl.objective.edge.costs)
    f3list = list(bal.objective.edge.costs)
    f4list = list(pairEuc.objective.edge.costs)
    
    ## Propose possible rho values for multi-objective optimization problem
    rho_list <- list()
    rho_list <- rho_proposition(paircost.v, rho.max.f, rho1, rho2)
    
    # #### EXPENSIVE PART ####
    #
    solutions <- list()
    solution.nets <- list()
    rho_counts = 1
    for (rho in rho_list) {
      temp.sol <- solveP1(
        net.i,
        f1list,
        f2list,
        f4list,
        rho1 = rho[1],
        rho2 = rho[2],
        tol = toleranceOption
      )
      
      solutions[[as.character(rho_counts)]] <- temp.sol$x
      solution.nets[[as.character(rho_counts)]] <- temp.sol$net
      print(paste('Matches finished for rho1 = ', rho[1], 
                  " and rho2 = ", rho[2]))
      rho_counts = rho_counts + 1
    }
    
    match.list <-
      obj.to.match(
        list(
          'net' = net.i,
          'costs' = pair.objective.edge.costs,
          'balance' = bal.objective.edge.costs,
          'rho.v' = 1:length(rho_list),
          'solutions' = solutions
        )
      )
    match.list <- match.list[order(as.numeric(names(match.list)))]
    
    ## remove the matching results that result in zero matches
    temp.length <- laply(match.list, nrow)
    temp.names <- names(match.list)
    nonzero.names <- temp.names[temp.length != 0]
    zero.names <- temp.names[temp.length == 0]
    
    if (length(zero.names) == length(temp.names)) {
      stop(
        'Error: Initial rho values result in non-matching. 
        Please try a new set of initial rho values.'
      )
    }
    match.list <-
      match.list[names(match.list) %in% nonzero.names == TRUE]
    
    
    ## Iterative Search for Best Rhos (if necessary)
    numIter = 1
    
    
    tempRhos <- as.numeric(names(match.list))
    percentageUnmatched <-
      1 - as.vector(laply(match.list, nrow) / sum(dat[, treatCol]))
    minInd <- which.min(percentageUnmatched)
    bestRhoInd <- tempRhos[minInd]
    bestPercentageSoFar <- as.vector(percentageUnmatched)[minInd]
    bestRho1 <- rho_list[[bestRhoInd]][1]
    bestRho2 <- rho_list[[bestRhoInd]][2]
    ind = length(rho_list) + 1
    while ((bestPercentageSoFar > maxUnMatched) &
           (numIter <= maxIter)) {
      if (is.finite(log10(bestRho1)) & is.finite(log10(bestRho2))) {
        proposedNewRho1s <-
          c(runif(2, 0, 2), 0.5 * 
              (-c((abs(rnorm(1)) * 0.01 + bestRho1 / 10):floor(log10(bestRho1))
          )), 0.5 * 
            (-c(abs(rnorm(1)) * 0.01 + bestRho1:log10(bestRho1 + 
                                                        0.01 * abs(rnorm(1))) *
              10
          )))
        
        proposedNewRho2s <-
          c(runif(2, 0, 2), 0.5 * 
              (-c((abs(rnorm(1)) * 0.01 + bestRho2 / 10):floor(log10(bestRho2))
          )), 0.5 * (-c(
            abs(rnorm(1)) * 0.01 + bestRho2:log10(bestRho2 + 
                                                    0.01 * abs(rnorm(1))) *
              10
          )))
      } else {
        proposedNewRho1s <-
          c(rnorm(1) * 1e-1, rnorm(1) * 1e2, rnorm(1) * 1e3, runif(2, 0, 2))
        proposedNewRho2s <-
          c(rnorm(1) * 1e-1, rnorm(1) * 1e2, rnorm(1) * 1e3, runif(2, 0, 2))
      }
      
      proposedNewRho1s <- proposedNewRho1s[proposedNewRho1s >= 0]
      proposedNewRho2s <- proposedNewRho2s[proposedNewRho2s >= 0]
      
      for (r1 in proposedNewRho1s) {
        rho_list[[ind]] = c(r1, bestRho2)
        temp.sol <- solveP1(
          net.i,
          f1list,
          f2list,
          f3list,
          rho1 = r1,
          rho2 = bestRho2,
          tol = toleranceOption
        )
        
        solutions[[as.character(ind)]] <- temp.sol$x
        solution.nets[[as.character(ind)]] <- temp.sol$net
        print(paste('Matches finished for rho1 = ', r1, 
                    " and rho2 = ", bestRho2))
        if (is.null(temp.sol$x)) {
          print(
            paste(
              'However, rho1 = ',
              r1,
              " and rho2 = ",
              bestRho2,
              " results in empty matching."
            )
          )
        }
        ind = ind + 1
      }
      
      for (r2 in proposedNewRho2s) {
        rho_list[[ind]] = c(bestRho1, r2)
        temp.sol <- solveP1(
          net.i,
          f1list,
          f2list,
          f3list,
          rho1 = bestRho1,
          rho2 = r2,
          tol = toleranceOption
        )
        
        solutions[[as.character(ind)]] <- temp.sol$x
        solution.nets[[as.character(ind)]] <- temp.sol$net
        print(paste('Matches finished for rho1 = ', 
                    bestRho1, " and rho2 = ", r2))
        if (is.null(temp.sol$x)) {
          print(
            paste(
              'However, rho1 = ',
              bestRho1,
              " and rho2 = ",
              r2,
              " results in empty matching."
            )
          )
        }
        ind = ind + 1
      }
      
      match.list <-
        obj.to.match(
          list(
            'net' = net.i,
            'costs' = pair.objective.edge.costs,
            'balance' = bal.objective.edge.costs,
            'rho.v' = 1:length(solutions),
            'solutions' = solutions
          )
        )
      match.list <- match.list[order(as.numeric(names(match.list)))]
      
      ## remove the matching results that result in zero matches
      temp.length <- laply(match.list, nrow)
      temp.names <- names(match.list)
      nonzero.names <- temp.names[temp.length != 0]
      match.list <-
        match.list[names(match.list) %in% nonzero.names == TRUE]
      
      
      tempRhos <- as.numeric(names(match.list))
      percentageUnmatched <-
        1 - as.vector(laply(match.list, nrow) / sum(dat[, treatCol]))
      oldMinInd <- minInd
      minInd <- which.min(percentageUnmatched[-oldMinInd])
      if (minInd >= oldMinInd) {
        minInd = minInd + 1
      }
      bestRhoInd <- tempRhos[minInd]
      bestPercentageSoFar <- as.vector(percentageUnmatched)[minInd]
      if (sum(percentageUnmatched == bestPercentageSoFar) > 1) {
        bestRhoInd <-
          sample(which(percentageUnmatched == bestPercentageSoFar),
                 size = 1)
      }
      bestRho1 <- rho_list[[bestRhoInd]][1]
      bestRho2 <- rho_list[[bestRhoInd]][2]
      
      numIter = numIter + 1
    }
    
    solutions.old <- solutions
    solution.nets.old <- solution.nets
    
    pair_cost_sum <- c()
    pair_cost_euc_sum <- c()
    total_treated = sum(df[[treatCol]])
    f1_sum <- c()
    f2_sum <- c()
    f3_sum <- c()
    
    for (ind in names(match.list)) {
      x <- solutions[[ind]]
      pair_cost_sum[[ind]] <- sum(paircost.v * x[1:length(paircost.v)])
      pair_cost_euc_sum[[ind]] <-
        sum(paircostEuc.v * x[1:length(paircostEuc.v)])
      unmatch_ind = total_treated - length(match.list[[ind]])
      f1_sum[[ind]] <- sum(paircost.v * x[1:length(paircost.v)])
      f2_sum[[ind]] <- unmatch_ind
      f3_sum[[ind]] <- sum(paircostEuc.v * x[1:length(paircostEuc.v)])
      
      
      
    }
    
    rho_list_final = c()
    for (ind in names(match.list)) {
      rho_list_final[[ind]] <- rho_list[[as.numeric(ind)]]
    }
    
    result$rhoList <- rho_list_final
    result$matchList <-  match.list
    ## Store some information about treatment column, balance covariates and
    ## exact match column
    result$treatmentCol = treatCol
    result$covs = distList1
    if (is.null(distList1)) {
      result$covs = xNames
    }
    result$exactCovs = exactlist
    result$idMapping = dat$originalID
    result$stats = NULL
    result$b.var = myBalCol
    result$dataTable = dat
    result$t = treatCol
    result$df = df
    result$pair_cost1 = pair_cost_sum
    result$pair_cost2 = pair_cost_euc_sum
    result$version = 'Advanced'
    result$fDist1 <- f1_sum
    result$fExclude <- f2_sum
    result$fDist2 <- f3_sum
    return(result)
  }








#' Generate balance table
#' @description The helper function can generate tabular analytics that quantify
#' covariate imbalance after matching. It only works for the 'Basic' version of
#' matching (produced by `distBalMatch`).
#'
#' @family numerical analysis helper functions
#' @param matchingResult an object returned by the main matching function
#'   distBalMatch
#' @param covList (optional) a vector of names of covariates used for evaluating
#'   covariate imbalance; NULL by default.
#' @param display.all (optional) a boolean value indicating whether or not to
#'   show the data for all possible matches; TRUE by default
#' @param statList (optional) a vector of statistics that are calculated for
#'   evaluating the covariate imbalance between treated and control group. The
#'   types that are supported can be found here: \link[cobalt]{bal.tab}.
#'
#' @return a named list object containing covariate balance table and statistics
#'   for numer of units being matched for each match; the names are the
#'   character of index for each match in the `matchResult`.
#'
#' @details The result can be either directly used by indexing into the list, or
#'   post-processing by the function `compareTables` that summarizes the
#'   covariate balance information in a tidier table. Users can specify the
#'   arguments as follows: * `covList`: if it is set of NULL, all the covariates
#'   are included for the covariate balance table; otherwise, only the specified
#'   covariates will be included in the tabular result. * `display.all`: by
#'   default, the summary statistics for each match are included when the
#'   argument is set to TRUE. If the user only wants to see the summary
#'   statistics for matches with varying performance on three different
#'   objective values, the function would only display the matches with number
#'   of treated units being excluded at different quantiles. User can switch to
#'   the brief version by setting the parameter to FALSE. * `statList` is the
#'   list of statistics used for measuring balance. The argument is the same as
#'   `stats` argument in \link[cobalt]{bal.tab}, which is the function that is
#'   used for calculating the statistics. By default, only standardized
#'   difference in means is calculated.
#' @export
#'
#' @examples
#' ## Generate matches
#' data("lalonde", package="cobalt")
#' psCols <- c("age", "educ", "married", "nodegree")
#' treatVal <- "treat"
#' responseVal <- "re78"
#' pairDistVal <- c("age", "married", "educ", "nodegree")
#' exactVal <- c("educ")
#' myBalVal <- c("race")
#' r1s <- c( 0.1, 0.3, 0.5, 0.7, 0.9,1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7)
#' r2s <- c(0.01)
#' matchResult <- distBalMatch(df=lalonde, treatCol=treatVal, myBalCol=myBalVal,
#' rhoExclude =r1s, rhoBalance=r2s,
#' distList=pairDistVal, exactlist=exactVal,
#' propensityCols = psCols,ignore = c(responseVal), maxUnMatched = 0.1,
#' caliperOption=NULL, toleranceOption=1e-1, maxIter=0, rho.max.f = 10)
#' 
#' ## Generate summary table for balance 
#' balanceTables <- generateBalanceTable(matchResult)
#' balanceTableMatch10 <- balanceTables$'10'
generateBalanceTable <-
  function(matchingResult,
           covList = NULL,
           display.all = TRUE,
           statList = c("mean.diffs")) {
    if (matchingResult$version != "Basic") {
      stop("This graphing function only 
           works for result return by distBalMatch")
    }
    originalDF <- matchingResult$dataTable
    
    numTreated = sum(originalDF$treat)
    percentageUnmatched <-
      1 - as.vector(laply(matchingResult$matchList, nrow) / numTreated)
    sortedMatchList <-
      matchingResult$matchList[order(percentageUnmatched)]
    
    matchIndex <- names(sortedMatchList)
    if (display.all == TRUE & length(matchIndex) >= 5) {
      matchIndex <- matchIndex[as.integer(quantile(1:length(matchIndex)))]
    }
    balanceTable <- list()
    
    for (ind in matchIndex) {
      matchTab <- matchingResult$matchList[[ind]]
      treatedUnits <- as.numeric(rownames(matchTab))
      controlUnits <- as.vector(matchTab[, 1]) + numTreated
      if (!is.null(covList)) {
        covariates <- originalDF[c(treatedUnits, controlUnits), covList]
      } else{
        covariates <-
          originalDF[c(treatedUnits, controlUnits), matchingResult$covs]
      }
      balanceTable[[ind]] <-
        bal.tab(
          covariates,
          s.d.denom = "pooled",
          treat = as.vector(originalDF[c(treatedUnits, controlUnits), 'treat']),
          stats = statList
        )
    }
    return(balanceTable)
  }


#' Summarize covariate balance table
#'
#' @description This function would take the result of `generateBalanceTable`
#' function and combine the results in a single table. It only works for 'Basic'
#' version of the matching.
#'
#' @param balanceTable a named list, which is the result from the function
#'   `generateBalanceTable`
#'
#' @return a dataframe with combined information
#' @export
#'
#' @examples
#' ## Generate matches 
#' data("lalonde", package="cobalt")
#' psCols <- c("age", "educ", "married", "nodegree")
#' treatVal <- "treat"
#' responseVal <- "re78"
#' pairDistVal <- c("age", "married", "educ", "nodegree")
#' exactVal <- c("educ")
#' myBalVal <- c("race")
#' r1s <- c( 0.1, 0.3, 0.5, 0.7, 0.9,1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7)
#' r2s <- c(0.01)
#' matchResult <- distBalMatch(df=lalonde, treatCol=treatVal, myBalCol=myBalVal,
#' rhoExclude =r1s, rhoBalance=r2s,
#' distList=pairDistVal, exactlist=exactVal,
#' propensityCols = psCols,ignore = c(responseVal), maxUnMatched = 0.1,
#' caliperOption=NULL, toleranceOption=1e-1, maxIter=0, rho.max.f = 10)
#' 
#' ## Generate summary table for comparing matches 
#' compareTables(generateBalanceTable(matchResult))
compareTables <- function(balanceTable) {
  inds <- names(balanceTable)
  rnames <- rownames(balanceTable[[inds[1]]]$Balance)
  result <- data.frame(balanceTable[[inds[1]]]$Balance[, 1])
  colnames(result) <- c("type")
  count = 1
  for (ind in inds) {
    result[ind] <- balanceTable[[ind]]$Balance[rnames, 2]
    count = count + 1
  }
  rownames(result) <- rnames
  return(result)
}

#' Generate covariate balance in different matches
#'
#' @description This is a wrapper function for use in evaluating covariate
#' balance across different matches. The function calls `compareTables` on the
#' output from the function `generateBalanceTable`. It only works for 'Basic'
#' version of matching (using `distBalMatch`).
#'
#'
#' @param matchingResult an object returned by the main matching function
#'   distBalMatch
#' @param covList (optional) factor of names of covariates that we want to
#'   evaluate covariate balance on; default is NULL. When set to NULL, the
#'   program will compare the covariates that have been used to construct a
#'   propensity model.
#' @param display.all (optional) boolean value of whether to display all the
#'   matches; default is TRUE, where matches at each quantile is displayed
#' @param stat (optional) character of the name of the statistic used for
#'   measuring covariate balance; default is "mean.diff". This argument is the
#'   same as used in "cobalt" package, see: \link[cobalt]{bal.tab}
#'
#' @return a dataframe that shows covariate balance in different matches
#' @export
#' @examples
#' ## Generate matches 
#' data("lalonde", package="cobalt")
#' psCols <- c("age", "educ", "married", "nodegree")
#' treatVal <- "treat"
#' responseVal <- "re78"
#' pairDistVal <- c("age", "married", "educ", "nodegree")
#' exactVal <- c("educ")
#' myBalVal <- c("race")
#' r1s <- c( 0.1, 0.3, 0.5, 0.7, 0.9,1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7)
#' r2s <- c(0.01)
#' matchResult <- distBalMatch(df=lalonde, treatCol=treatVal, myBalCol=myBalVal,
#' rhoExclude =r1s, rhoBalance=r2s,
#' distList=pairDistVal, exactlist=exactVal,
#' propensityCols = psCols,ignore = c(responseVal), maxUnMatched = 0.1,
#' caliperOption=NULL, toleranceOption=1e-1, maxIter=0, rho.max.f = 10)
#' 
#' ## Generate table for comparing matches
#' compareMatching(matchResult, display.all = TRUE)
compareMatching <-
  function(matchingResult,
           covList = NULL,
           display.all = TRUE,
           stat = "mean.diff") {
    if (matchingResult$version != "Basic") {
      stop("This graphing function only works 
           for result return by distBalMatch")
    }
    return(compareTables(
      generateBalanceTable(
        matchingResult,
        covList,
        display.all = display.all,
        statList = c(stat)
      )
    ))
  }


#' An internal helper function that gives the index of matching with a wide
#' range of number of treated units left unmatched
#'
#' @param matchingResult an object returned by the main matching function
#'   distBalMatch
#'
#' @return a vector of five matching indices with the number of treated units
#'   excluded at 0th, 25th, 50th, 75th and 100th percentiles respectively.
get_five_index <- function(matchingResult) {
  df <- compareMatching(matchingResult)
  matches <- colnames(df)[2:length(colnames(df))]
  return(matches)
}



#' Marginal imbalance vs. exclusion
#' @description Plotting function that visualizes the tradeoff between the total
#' variation imbalance on a specified variable and the number of unmatched
#' treated units. This function only works for the 'Basic' version of matching
#' (conducted using `distBalMatch`).
#'
#' @family Graphical helper functions for analysis
#' @param matchingResult an object returned by the main matching function
#'   distBalMatch
#'
#' @return No return value, called for visualization of match result
#' @export
#' @examples
#' ## Generate matches 
#' data("lalonde", package="cobalt")
#' psCols <- c("age", "educ", "married", "nodegree")
#' treatVal <- "treat"
#' responseVal <- "re78"
#' pairDistVal <- c("age", "married", "educ", "nodegree")
#' exactVal <- c("educ")
#' myBalVal <- c("race")
#' r1s <- c( 0.1, 0.3, 0.5, 0.7, 0.9,1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7)
#' r2s <- c(0.01)
#' matchResult <- distBalMatch(df=lalonde, treatCol=treatVal, myBalCol=myBalVal,
#' rhoExclude =r1s, rhoBalance=r2s,
#' distList=pairDistVal, exactlist=exactVal,
#' propensityCols = psCols,ignore = c(responseVal), maxUnMatched = 0.1,
#' caliperOption=NULL, toleranceOption=1e-1, maxIter=0, rho.max.f = 10)
#' 
#' ## Generate visualization of tradeoff between total variation distance and 
#' ## number of treated units left unmatched
#' generateTVGraph(matchResult)
generateTVGraph <- function(matchingResult) {
  if (matchingResult$version != "Basic") {
    stop("This graphing function only works for result return by distBalMatch")
  }
  inds <- get_five_index(matchingResult)
  graph_labels <- names(matchingResult$matchList)
  for (i in 1:length(graph_labels)) {
    if (sum(as.numeric(inds) == as.numeric(graph_labels[i])) == 0) {
      graph_labels[i] = " "
    }
  }
  treatedSize = sum(matchingResult$dataTable$treat)
  samp.size <- laply(matchingResult$matchList, nrow)
  f1 <- matchingResult$stats[1, ] * samp.size
  f2 <- treatedSize  - samp.size
  plot(
    f2,
    f1,
    pch = 20,
    xlab = 'Treated Units Unmatched',
    ylab = 'Total Variation Imbalance',
    cex = 1.2,
    xlim = c(0, treatedSize),
    ylim = c(0, max(f1)),
    cex.lab = 1.2,
    cex.axis = 1.2,
    col = rgb(0, 0, 0, 0.3)
  )
  abline(h = 0, lty = 'dotted')
  points(f2, f1, col = rgb(0, 0, 1, 1), pch = 20)
  text(f2, f1, labels = graph_labels, pos = 2)
}

#' Distance vs. exclusion
#' @description Plotting function that generate sum of pair-wise distance vs.
#' number of unmatched treated units
#'
#' @family Graphical helper functions for analysis
#' @param matchingResult an object returned by the main matching function
#'   distBalMatch
#'
#' @return No return value, called for visualization of match result
#' @export
#' @examples
#' ## Generate matches 
#' data("lalonde", package="cobalt")
#' psCols <- c("age", "educ", "married", "nodegree")
#' treatVal <- "treat"
#' responseVal <- "re78"
#' pairDistVal <- c("age", "married", "educ", "nodegree")
#' exactVal <- c("educ")
#' myBalVal <- c("race")
#' r1s <- c( 0.1, 0.3, 0.5, 0.7, 0.9,1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7)
#' r2s <- c(0.01)
#' matchResult <- distBalMatch(df=lalonde, treatCol=treatVal, myBalCol=myBalVal,
#' rhoExclude =r1s, rhoBalance=r2s,
#' distList=pairDistVal, exactlist=exactVal,
#' propensityCols = psCols,ignore = c(responseVal), maxUnMatched = 0.1,
#' caliperOption=NULL, toleranceOption=1e-1, maxIter=0, rho.max.f = 10)
#' 
#' ## Generate visualization of tradeoff between pari-wise distance sum and 
#' ## number of treated units left unmatched
#' generatePairdistanceGraph(matchResult)
generatePairdistanceGraph <- function(matchingResult) {
  if (matchingResult$version != "Basic") {
    stop("This graphing function only works for result return by distBalMatch")
  }
  inds <- get_five_index(matchingResult)
  graph_labels <- names(matchingResult$matchList)
  for (i in 1:length(graph_labels)) {
    if (sum(as.numeric(inds) == as.numeric(graph_labels[i])) == 0) {
      graph_labels[i] = " "
    }
  }
  treatedSize = sum(matchingResult$dataTable$treat)
  samp.size <- laply(matchingResult$matchList, nrow)
  f1 <- as.numeric(matchingResult$pair_cost1)
  f2 <- treatedSize  - samp.size
  plot(
    f2,
    f1,
    pch = 20,
    xlab = 'Treated Units Unmatched',
    ylab = 'Pair-wise Distance Sum',
    cex = 1.2,
    xlim = c(0, treatedSize),
    ylim = c(0, max(f1)),
    cex.lab = 1.2,
    cex.axis = 1.2,
    col = rgb(0, 0, 0, 0.3)
  )
  abline(h = 0, lty = 'dotted')
  points(f2, f1, col = rgb(0, 0, 1, 1), pch = 20)
  text(f2, f1, labels = graph_labels, pos = 2)
}

#' Total variation imbalance vs. marginal imbalance
#' @description Plotting function that generate sum of pairwise distance vs.
#' total variation imbalance on specified balance variable. This function only
#' works for 'Basic' version of matching (conducted using `distBalMatch`).
#'
#' @family Graphical helper functions for analysis
#' @param matchingResult an object returned by the main matching function
#'   distBalMatch
#'
#' @return No return value, called for visualization of match result
#' @export
#' @examples
#' ## Generate matches 
#' data("lalonde", package="cobalt")
#' psCols <- c("age", "educ", "married", "nodegree")
#' treatVal <- "treat"
#' responseVal <- "re78"
#' pairDistVal <- c("age", "married", "educ", "nodegree")
#' exactVal <- c("educ")
#' myBalVal <- c("race")
#' r1s <- c( 0.1, 0.3, 0.5, 0.7, 0.9,1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7)
#' r2s <- c(0.01)
#' matchResult <- distBalMatch(df=lalonde, treatCol=treatVal, myBalCol=myBalVal,
#' rhoExclude =r1s, rhoBalance=r2s,
#' distList=pairDistVal, exactlist=exactVal,
#' propensityCols = psCols,ignore = c(responseVal), maxUnMatched = 0.1,
#' caliperOption=NULL, toleranceOption=1e-1, maxIter=0, rho.max.f = 10)
#' 
#' ## Visualize the tradeoff between the pair-wise distance sum and 
#' ## total variation distance 
#' generatePairdistanceBalanceGraph(matchResult)
generatePairdistanceBalanceGraph <- function(matchingResult) {
  if (matchingResult$version != "Basic") {
    stop("This graphing function only works for result return by distBalMatch")
  }
  inds <- get_five_index(matchingResult)
  graph_labels <- names(matchingResult$matchList)
  for (i in 1:length(graph_labels)) {
    if (sum(as.numeric(inds) == as.numeric(graph_labels[i])) == 0) {
      graph_labels[i] = " "
    }
  }
  treatedSize = sum(matchingResult$dataTable$treat)
  samp.size <- laply(matchingResult$matchList, nrow)
  f1 <- as.numeric(matchingResult$pair_cost1)
  f2 <- matchingResult$stats[1, ] * samp.size
  plot(
    f2,
    f1,
    pch = 20,
    xlab = 'Total Variation Imbalance on Balance Variable',
    ylab = 'Pair-wise Distance Sum',
    cex = 1.2,
    xlim = c(0, max(f2)),
    ylim = c(0, max(f1)),
    cex.lab = 1.2,
    cex.axis = 1.2,
    col = rgb(0, 0, 0, 0.3)
  )
  abline(h = 0, lty = 'dotted')
  points(f2, f1, col = rgb(0, 0, 1, 1), pch = 20)
  text(f2, f1, labels = graph_labels, pos = 2)
}

#' An internal helper function that translates the matching index in the sorted
#' data frame to the original dataframe's row index
#'
#' @param matchingResult an object returned by the main matching function
#'   distBalMatch
#'
#' @return NULL
convert_index <- function(matchingResult) {
  idMap <- matchingResult$idMapping
  numTreated = sum(matchingResult$dataTable$treat)
  matchIndex <- names(matchingResult$matchList)
  result <- list()
  
  for (ind in matchIndex) {
    matchTab <- matchingResult$matchList[[ind]]
    treatedUnits <- as.numeric(rownames(matchTab))
    controlUnits <- as.vector(matchTab[, 1]) + numTreated
    treatedOriginalID <- idMap[treatedUnits]
    controlOriginalID <- idMap[controlUnits]
    tempMatch <- list()
    tempMatch$treated = treatedOriginalID
    tempMatch$control = controlOriginalID
    result[[ind]] <- tempMatch
  }
  return(result)
}

#' An internal helper function that translate the matching index in the sorted
#' data frame to the original dataframe's row index
#'
#' @param matchingResult an object returned by the main matching function
#'   distBalMatch
#'
#' @return NULL
matched_index <- function(matchingResult) {
  return(convert_index(matchingResult))
}


#' Get matched dataframe
#' @description A function that returns the dataframe that contains only matched
#' pairs from the original data frame with specified match index
#'
#' @param matchingResult an object returned by the main matching function
#'   distBalMatch
#' @param match_num Integer index of match that the user want to extract paired
#'   observations from
#'
#' @return dataframe that contains only matched pair data
#' @export
#' @examples
#' ## Generate Matches
#' x1 = rnorm(100, 0, 0.5)
#' x2 = rnorm(100, 0, 0.1)
#' x3 = rnorm(100, 0, 1)
#' x4 = rnorm(100, x1, 0.1)
#' r1ss <- seq(0.1,50, 10)
#' r2ss <- seq(0.1,50, 10)
#' x = cbind(x1, x2, x3,x4)
#' z = sample(c(rep(1, 50), rep(0, 50)))
#' e1 = rnorm(100, 0, 1.5)
#' e0 = rnorm(100, 0, 1.5)
#' y1impute = x1^2 + 0.6*x2^2 + 1 + e1
#' y0impute = x1^2 + 0.6*x2^2 + e0
#' treat = (z==1)
#' y = ifelse(treat, y1impute, y0impute)
#' names(x) <- c("x1", "x2", "x3", "x4")
#' df <- data.frame(cbind(z, y, x))
#' df$x5 <- 1
#' names(x) <- c("x1", "x2", "x3", "x4")
#' df <- data.frame(cbind(z, y, x))
#' df$x5 <- 1
#' d1 <- as.matrix(dist(df["x1"]))
#' d2 <- as.matrix(dist(df["x2"]))
#' idx <- 1:length(z)
#' treatedUnits <- idx[z==1]
#' controlUnits <- idx[z==0]
#' d1 <- as.matrix(d1[treatedUnits, controlUnits])
#' d2 <- as.matrix(d2[treatedUnits, controlUnits])
#' matchResult1 <- twoDistMatch(df, "z", "y", dMat1=d1, dType1= "User", dMat2=d2,
#' dType2="User", myBalCol=c("x5"), rhoExclude=r1ss, rhoDistance=r2ss,
#' propensityCols = c("x1")) 
#' matchedData(matchResult1, 1)
matchedData <- function(matchingResult, match_num) {
  matched_idx <- matched_index(matchingResult)
  treated_idx <- matched_idx[[as.character(match_num)]]$treated
  control_idx <- matched_idx[[as.character(match_num)]]$control
  return(matchingResult$df[c(treated_idx, control_idx), ])
}


#' Get unmatched percentage
#' @description A function that generate the percentage of unmatched units for
#' each match.
#'
#' @family numerical analysis helper functions
#' @param matchingResult matchingResult object that contains information for all
#'   matches
#'
#' @return data frame with three columns, one containing the matching index, one
#'   containing the number of matched units, and one conatining the percentage
#'   of matched units (out of original treated group size).
#' @export
#' @examples
#' \dontrun{
#' getUnmatched(matchResult)
#' }
getUnmatched <- function(matchingResult) {
  matchingIndex <- names(matchingResult$matchList)
  matchedUnits <- laply(matchingResult$matchList, nrow)
  matchedPercentage <-
    matchedUnits / sum(matchingResult$dataTable$treat)
  result <-
    data.frame(matchingIndex, matchedUnits, matchedPercentage)
  colnames(result) <-
    c("Matching Index",
      "Number of Matched Units",
      "Percentage of Matched Units")
  return(result)
  
}


#' Penalty and objective values summary
#' @description Helper function to generate a dataframe with matching number,
#' penalty (rho) values, and objective function values.
#'
#' @family numerical analysis helper functions
#' @param matchingResult matchingResult object that contains information for all
#'   matches.
#'
#' @return a dataframe that contains objective function values and rho values
#'   corresponding coefficients before each objective function.
#' @export
#' @examples
#' ## Generate matches
#' x1 = rnorm(100, 0, 0.5)
#' x2 = rnorm(100, 0, 0.1)
#' x3 = rnorm(100, 0, 1)
#' x4 = rnorm(100, x1, 0.1)
#' r1ss <- seq(0.1,50, 10)
#' r2ss <- seq(0.1,50, 10)
#' x = cbind(x1, x2, x3,x4)
#' z = sample(c(rep(1, 50), rep(0, 50)))
#' e1 = rnorm(100, 0, 1.5)
#' e0 = rnorm(100, 0, 1.5)
#' y1impute = x1^2 + 0.6*x2^2 + 1 + e1
#' y0impute = x1^2 + 0.6*x2^2 + e0
#' treat = (z==1)
#' y = ifelse(treat, y1impute, y0impute)
#' names(x) <- c("x1", "x2", "x3", "x4")
#' df <- data.frame(cbind(z, y, x))
#' df$x5 <- 1
#' names(x) <- c("x1", "x2", "x3", "x4")
#' df <- data.frame(cbind(z, y, x))
#' df$x5 <- 1
#' d1 <- as.matrix(dist(df["x1"]))
#' d2 <- as.matrix(dist(df["x2"]))
#' idx <- 1:length(z)
#' treatedUnits <- idx[z==1]
#' controlUnits <- idx[z==0]
#' d1 <- as.matrix(d1[treatedUnits, controlUnits])
#' d2 <- as.matrix(d2[treatedUnits, controlUnits])
#' matchResult1 <- twoDistMatch(df, "z", "y", dMat1=d1, dType1= "User", dMat2=d2,
#' dType2="User", myBalCol=c("x5"), rhoExclude=r1ss, rhoDistance=r2ss,
#' propensityCols = c("x1")) 
#' 
#' ## Generate tabular summary 
#' generateRhoObj(matchResult1)
generateRhoObj <- function(matchingResult) {
  matching_index <- names(matchingResult$matchList)
  if (matchingResult$version != "Basic") {
    dat_tab <-
      cbind(
        as.character(matching_index),
        as.numeric(matchingResult$fDist1),
        as.numeric(matchingResult$fExclude),
        as.numeric(matchingResult$fDist2),
        laply(matchingResult$rhoList, function(x)
          x[1]),
        laply(matchingResult$rhoList, function(x)
          x[2])
      )
    dat_tab <- data.frame(dat_tab)
    colnames(dat_tab) <-
      c("match_index",
        "fDist1",
        "fExclude",
        "fDist2",
        "rhoExclude",
        "rhoDist2")
    dat_tab["rhoDist1"] <- rep(1, nrow(dat_tab))
    dat_tab = dat_tab[, c(
      "match_index",
      "rhoDist1",
      "rhoExclude",
      "rhoDist2",
      "fDist1",
      "fExclude",
      "fDist2"
    )]
  } else {
    dat_tab <-
      cbind(
        as.character(matching_index),
        as.numeric(matchingResult$fPair),
        as.numeric(matchingResult$fExclude),
        as.numeric(matchingResult$fMarginal),
        laply(matchingResult$rhoList, function(x)
          x[1]),
        laply(matchingResult$rhoList, function(x)
          x[2])
      )
    dat_tab <- data.frame(dat_tab)
    colnames(dat_tab) <-
      c("match_index",
        "fPair",
        "fExclude",
        "fMarginal",
        "rhoExclude",
        "rhoMarginal")
    dat_tab["rhoPair"] <- rep(1, nrow(dat_tab))
    dat_tab = dat_tab[, c(
      "match_index",
      "rhoPair",
      "rhoExclude",
      "rhoMarginal",
      "fPair",
      "fExclude",
      "fMarginal"
    )]
  }
  if (matchingResult$version == "Basic") {
    print(
      paste0(
        "rhoPair: penalty for distance objective function;",
        "rhoExclude: penalty for exclusion cost;",
        "rhoMarginal: penalty for marginal balance."
      )
    )
    print(
      paste0(
        "fPair: pair-wise distance sum objective function;",
        "fExclude: number of treated units left unmatched;",
        "fMarginal: marginal imbalance measured as the total
        variation distance on the marginal distribution of speicified variable."
      )
    )
  }
  
  if (matchingResult$version == "Advanced") {
    print(
      paste0(
        "rhoDist1: penalty for the first distance objective function;",
        "rhoExclude: penalty for exclusion cost;",
        "rhoDist2: penalty for the second distance objective function."
      )
    )
    print(
      paste0(
        "fDist1: pair-wise distance sum 
        objective function of first distance measure;",
        "fExclude: exclusion cost;",
        "fDist2: pair-wise distance sum objective 
        function of second distance measure."
      )
    )
  }
  
  
  return(dat_tab)
}



#' Internal helper function that converts axis name to internal variable name
#'
#' @param x the user input character for x-axis value
#' @param y the user input character for y-axis value
#'
#' @return a named list with variable names for visualization for internal use
convert_names <- function(x, y) {
  result <- c()
  if (x == "dist1") {
    result$x_label <- "fDist1"
    result$x_name <- "Distance 1"
  } else if (x == "exclude") {
    result$x_label <- "fExclude"
    result$x_name <- "Number of Treated Units Unmatched"
  } else if (x == "dist2") {
    result$x_label <- "fDist2"
    result$x_name <- "Distance 2"
  } else if (x == "pair") {
    result$x_label <- "fPair"
    result$x_name <- "Pairwise Distance"
  } else if (x == "marginal") {
    result$x_label <- "fMarginal"
    result$x_name <- "Marginal Balance"
  } else {
    stop(
      "Wrong argument for objective function name. 
      It must be one of the following:
         dist1, dist2, exclude, pair, marginal."
    )
  }
  
  if (y == "dist1") {
    result$y_label <- "fDist1"
    result$y_name <- "Distance 1"
  } else if (y == "exclude") {
    result$y_label <- "fExclude"
    result$y_name <- "Number of Treated Units Unmatched"
  } else if (y == "dist2") {
    result$y_label <- "fDist2"
    result$y_name <- "Distance 2"
  } else if (y == "pair") {
    result$y_label <- "fPair"
    result$y_name <- "Pairwise Distance"
  } else if (y == "marginal") {
    result$y_label <- "fMarginal"
    result$y_name <- "Marginal Balance"
  } else {
    stop(
      "Wrong argument for objective function name.
      It must be one of the following:
         dist1, dist2, exclude, distPair, distMarginal."
    )
  }
  return(result)
  
  
  
}



#' Visualize tradeoffs
#' @description Main visualization functions for showing the tradeoffs between
#' two of the three objective functions.
#'
#' @param matchingResult the matching result returned by either distBalMatch or
#'   twoDistMatch.
#' @param x_axis character, naming the objective function shown on x-axis; one
#'   of ("pair", "marginal", "dist1", "dist2", "exclude"), "dist1" by default.
#' @param y_axis character, naming the objective function shown on y-axis; one
#'   of ("pair", "marginal", "dist1", "dist2", "exclude"), "dist2" by default.
#' @param xlab (optional) the axis label for x-axis; NULL by default.
#' @param ylab (optional) the axis label for y-axis; NULL by default.
#' @param main (optional) the title of the graph; NULL by default.
#' @param display_all (optional) whether to show all the labels for match index;
#'   FALSE by default, which indicates the visualization function only labels
#'   matches at quantiles of number of treated units being excluded.
#' @param cond (optional) NULL by default, which denotes all the matches are
#'   shown; otherwise, takes a list of boolean values indicating whether to
#'   include each match
#'
#' @return No return value, called for visualization of match result
#' @details   By default, the plotting function will show the tradeoff between
#'   the first distance objective function and the marginal balance (if
#'   distBalMatch) is used; or simply the second distance objective function, if
#'   twoDistMatch is used.
#' @export
#'
#' @examples
#' ## Generate matches
#' data("lalonde", package="cobalt")
#' psCols <- c("age", "educ", "married", "nodegree")
#' treatVal <- "treat"
#' responseVal <- "re78"
#' pairDistVal <- c("age", "married", "educ", "nodegree")
#' exactVal <- c("educ")
#' myBalVal <- c("race")
#' r1s <- c( 0.1, 0.3, 0.5, 0.7, 0.9,1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7)
#' r2s <- c(0.01)
#' matchResult <- distBalMatch(df=lalonde, treatCol=treatVal, myBalCol=myBalVal,
#' rhoExclude =r1s, rhoBalance=r2s,
#' distList=pairDistVal, exactlist=exactVal,
#' propensityCols = psCols,ignore = c(responseVal), maxUnMatched = 0.1,
#' caliperOption=NULL, toleranceOption=1e-1, maxIter=0, rho.max.f = 10)
#' 
#' ## Visualization
#' visualize(matchResult, "marginal", "exclude")
#' visualize(matchResult, "pair", "exclude")
visualize <-
  function(matchingResult,
           x_axis = "dist1",
           y_axis = "dist2",
           xlab = NULL,
           ylab = NULL,
           main = NULL,
           display_all = FALSE,
           cond = NULL) {
    if (x_axis == "dist1" &&
        y_axis == "dist2" && matchingResult$version == "Basic") {
      x_axis = "pair"
      y_axis = "marginal"
    }
    if (!x_axis %in% c("pair", "marginal", "dist1", "dist2", "exclude")) {
      stop("Wrong name for argument x_axis")
    }
    if (!y_axis %in% c("pair", "marginal", "dist1", "dist2", "exclude")) {
      stop("Wrong name for argument y_axis")
    }
    
    if (((!x_axis %in% c("dist1", "dist2", "exclude")) ||
         (!y_axis %in% c("dist1", "dist2", "exclude"))) &&
        matchingResult$version == "Advanced") {
      stop("Please choose among 'dist1', 'dist2', 'exclude' for axis.")
    }
    
    if (((!x_axis %in% c("pair", "marginal", "exclude")) ||
         (!y_axis %in% c("pair", "marginal", "exclude"))) &&
        matchingResult$version == "Basic") {
      stop("Please choose among 'pair', 'marginal', 'exclude'for axis.")
    }
    
    
    naming <- convert_names(x_axis, y_axis)
    x_axis <- naming$x_label
    y_axis <- naming$y_label
    
    if (is.null(xlab)) {
      xlab <- naming$x_name
    }
    if (is.null(ylab)) {
      ylab <- naming$y_name
    }
    
    rho_obj_table <- generateRhoObj(matchingResult)
    if (!is.null(cond)) {
      rho_obj_table <- rho_obj_table[cond, ]
    }
    sorted_table <- rho_obj_table %>% arrange(fExclude)
    inds <- as.vector(sorted_table[["match_index"]])
    if (nrow(sorted_table) > 5 && display_all == FALSE) {
      inds <- inds[as.integer(quantile(1:length(inds)))]
    }
    
    graph_labels <- as.character(sorted_table[["match_index"]])
    
    for (i in 1:length(graph_labels)) {
      if (sum(as.numeric(inds) == as.numeric(graph_labels[i])) == 0) {
        graph_labels[i] = " "
      }
    }
    
    f_x <- as.numeric(as.vector(sorted_table[[x_axis]]))
    f_y <- as.numeric(as.vector(sorted_table[[y_axis]]))
    plot(
      f_x,
      f_y,
      pch = 20,
      xlab = xlab,
      ylab = ylab,
      main = main,
      cex = 1.2,
      xlim = c(0, max(f_x)),
      ylim = c(0, max(f_y)),
      cex.lab = 1.2,
      cex.axis = 1.2,
      col = rgb(0, 0, 0, 0.3)
    )
    abline(h = 0, lty = 'dotted')
    points(f_x, f_y, col = rgb(0, 0, 1, 1), pch = 20)
    text(f_x, f_y, labels = graph_labels, pos = 2)
  }




#' Helper function that change input distance matrix
#'
#' @param z the treatment vector
#' @param distMat the user input distance matrix
#'
#' @return a distance matrix where (i,j) element is the distance between unit i
#'   and j in the same order as z
distanceFunctionHelper <- function(z, distMat) {
  tmpInd <- 1:length(z)
  treatedInd <- tmpInd[z == 1]
  controlInd <- tmpInd[z == 0]
  matrixInd <- c(treatedInd, controlInd)
  finalMatrix <-
    matrix(rep(Inf, length(z) * length(z)),
           nrow = length(z),
           ncol = length(z))
  n1 = length(treatedInd)
  n0 = length(controlInd)
  for (i in 1:n1) {
    tInd <- treatedInd[i]
    for (j in 1:n0) {
      cInd <- controlInd[j]
      finalMatrix[tInd, cInd] = distMat[i, j]
      finalMatrix[cInd, tInd] = distMat[i, j]
    }
  }
  return(finalMatrix)
}
