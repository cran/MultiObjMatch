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
    mean.dist <- mean(paircosts.list)
    
    rho1 <- c(
      rho.min,
      min.dist / 2,
      min.dist,
      quantile(paircosts.list)[c(2, 4)],
      max.dist,
      rho.max.factor * max.dist,
      mean.dist
    )
    rho2 <- c(1)
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
  ind_tracker = 1
  for (i in 1:length(rho1.list)) {
    for (j in 1:length(rho2.list)) {
      result[[ind_tracker]] <- c(rho1.list[i], rho2.list[j])
      ind_tracker = ind_tracker + 1
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
#' @param data dataframe that contains a column named "treat", the treatment
#'   vector, and columns of covariates specified.
#' @param covs factor of column names of covariates used for fitting a
#'   propensity score model.
#'
#' @return vector of estimated propensity scores (on the probability scale).
getPropensityScore <- function(data, covs) {
  df = data
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
#' @param data data frame that contain columns indicating treatment, outcome and
#'   covariates.
#' @param treat_col character of name of the column indicating treatment
#'   assignment.
#' @param marg_bal_col character of column name of the variable on which to 
#' evaluate marginal balance.
#' @param exclusion_penalty (optional) numeric vector of values of exclusion 
#' penalty. Default is c(), which would trigger the auto grid search. 
#' @param balance_penalty (optional) factor of values of marginal balance
#'  penalty. Default value is c(), which would trigger the auto grid search.
#' @param dist_matrix (optional) a matrix that specifies the pair-wise distances
#'   between any two objects.
#' @param dist_col (optional) character vector of variable names used for
#'   calculating within-pair distance.
#' @param exact_col (optional) character vector, variable names that we want
#'   exact matching on; NULL by default.
#' @param propensity_col (optional) character vector, variable names on which to
#'   fit a propensity score (to supply a caliper).
#' @param pscore_name (optional) character, giving the variable name for the 
#' fitted propensity score.
#' @param ignore_col (optional) character vector of variable names that should be
#'   ignored when constructing the internal matching. NULL by default.
#' @param max_unmatched (optional) numeric, the maximum proportion of unmatched
#'   units that can be accepted; default is 0.25.
#' @param caliper_option (optional) numeric, the propensity score caliper value
#'   in standard deviations of the estimated propensity scores; default is NULL,
#'   which is no caliper.
#' @param tol (optional) numeric, tolerance of close match distance;
#'   default is 1e-2.
#' @param max_iter (optional) integer,  maximum number of iterations to use in
#'   searching for penalty combintions that improve the matching; default is 1,
#'   where the algorithm searches for one round.
#' @param rho_max_factor (optional) numeric, the scaling factor used in proposal 
#' for penalties; default is 10.
#' @param max_pareto_search_iter (optional) numeric, the number of tries to 
#' search for the tol that yield pareto optimal solutions; default is 5.
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
#'   function called; "Basic" indicates the matching comes from dist_bal_match and
#'   "Advanced" from two_dist_match. 
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
#'   denoted by Greek letter rho. Larger values of `exclusion_penalty` 
#'   corresponds to
#'   increased emphasis on retaining treated units (all else being equal), while
#'   larger values of `balance_penalty` corresponds to increased emphasis 
#'   on marginal
#'   balance. Additional details: 
#'   * Users may either specify their own distance
#'   matrix via the `dist_matrix` argument or ask the function to create a
#'   robust Mahalanobis distance matrix internally on a set of covariates specified by
#'   the `dist_col` argument; if neither argument is specified an error will
#'   result.  User-specified distance matrices should have row count equal to
#'   the number of treated units and column count equal to the number of
#'   controls. 
#'   * If the `caliper_option` argument is specified, a propensity
#'   score caliper will be imposed, forbidding matches between units more than a
#'   fixed distance apart on the propensity score.  The caliper will be based
#'   either on a user-fit propensity score, identified in the input dataframe by
#'   argument `pscore_name`, or by an internally-fit propensity score based on
#'   logistic regression against the variables named in `propensity_col`.  If
#'   `caliper_option` is non-NULL and neither of the other arguments is specified
#'   an error will result. 
#'   * `tol` controls the precision at which
#'   the objective functions is evaluated. When matching problems are especially
#'   large or complex it may be necessary to increase toleranceOption in order
#'   to prevent integer overflows in the underlying network flow solver;
#'   generally this will be suggested in appropariate warning messages. 
#'   * While
#'   by default tradeoffs are only assessed at penalty combinations provided by
#'   the user, the user may ask for the algorithm to search over additional
#'   penalty values in order to identify additional Pareto optimal solutions.
#'   `rho_max_factor` is a multiplier applied to initial penalties to discover
#'   new solutions, and setting it larger leads to wider exploration; similarly,
#'   `max_iter` controls how long the exploration routine runs, with larger
#'   values leading to more exploration.
#'
#' @export
#' @examples
#' data("lalonde", package="cobalt")
#' ps_cols <- c("age", "educ", "married", "nodegree", "race")
#' treat_val <- "treat"
#' response_val <- "re78"  
#' pair_dist_val <- c("age", "married", "educ", "nodegree", "race")
#' my_bal_val <- c("race")
#' r1s <- c(0.01,1,2,4,4.4,5.2,5.4,5.6,5.8,6)
#' r2s <- c(0.001)
#' match_result <- dist_bal_match(data=lalonde, treat_col= treat_val, 
#' marg_bal_col = my_bal_val, exclusion_penalty=r1s, balance_penalty=r2s, 
#' dist_col = pair_dist_val, 
#' propensity_col = ps_cols, max_iter=0)
dist_bal_match <-
  function(data,
           treat_col,
           marg_bal_col,
           exclusion_penalty = c(),
           balance_penalty = c(),
           dist_matrix = NULL,
           dist_col = NULL,
           exact_col = NULL,
           propensity_col = NULL,
           pscore_name = NULL,
           ignore_col = NULL,
           max_unmatched = 0.25,
           caliper_option = NULL,
           tol = 1e-2,
           max_iter = 1,
           rho_max_factor = 10,
           max_pareto_search_iter=5) {
    

    initial_match_result <- dist_bal_match_skeleton(data,
                                                    treat_col,
                                                    marg_bal_col,
                                                    exclusion_penalty,
                                                    balance_penalty,
                                                    dist_matrix,
                                                    dist_col,
                                                    exact_col,
                                                    propensity_col,
                                                    pscore_name,
                                                    ignore_col,
                                                    max_unmatched,
                                                    caliper_option,
                                                    tol,
                                                    max_iter,
                                                    rho_max_factor
                                                    )
    ## Base case - when no iteration needed 
    tol_iter_cnt = 1
    while((!check_pareto_optimality(initial_match_result)) && 
          (max_pareto_search_iter >= 1)) {
      max_pareto_search_iter <- max_pareto_search_iter - 1
      initial_match_result <- dist_bal_match_skeleton(data,
                                                      treat_col,
                                                      marg_bal_col,
                                                      exclusion_penalty,
                                                      balance_penalty,
                                                      dist_matrix,
                                                      dist_col,
                                                      exact_col,
                                                      propensity_col,
                                                      pscore_name,
                                                      ignore_col,
                                                      max_unmatched,
                                                      caliper_option,
                                                      tol/(10^tol_iter_cnt),
                                                      max_iter,
                                                      rho_max_factor
      )
      tol_iter_cnt <- tol_iter_cnt + 1
    }
    if(check_pareto_optimality(initial_match_result)){
      return(initial_match_result)
    } else {
      warning("Ordering of Pareto optimal solutions by 
              objective function value may be incorrect due to numerical 
              approximation errors.  Using a smaller value for tol may help.")
      return(initial_match_result)
    }

  }

dist_bal_match_skeleton <-
  function(data,
           treat_col,
           marg_bal_col,
           exclusion_penalty = c(),
           balance_penalty = c(),
           dist_matrix = NULL,
           dist_col = NULL,
           exact_col = NULL,
           propensity_col = NULL,
           pscore_name = NULL,
           ignore_col = NULL,
           max_unmatched = 0.25,
           caliper_option = NULL,
           tol = 1e-2,
           max_iter = 1 ,
           rho_max_factor = 10) {
    
    
    ## fix the naming; new in version 1.0.0
    df = data
    treatCol = treat_col
    myBalCol = marg_bal_col
    
    rhoExclude = exclusion_penalty
    rhoBalance = balance_penalty
    exactlist = exact_col
    propensityCols = propensity_col
    toleranceOption = tol
    maxIter = max_iter
    rho.max.f = rho_max_factor
    pScores = pscore_name
    distMatrix = dist_matrix
    distList = dist_col
    maxUnMatched = max_unmatched
    caliperOption = caliper_option
    ignore = ignore_col
    
    
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
    result$dataTable = dat
    result$numTreat <- sum(dat$treat)
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
    
    # 
    # if(is.null(distList)){
    #   distList <- xNames
    # }
    # 
    ## 2. Construct nets
    base.net <- netFlowMatch(dat$treat)
    
    ## 3. Define distance and cost
    caliperType = 'none'
    if (!is.null(caliperOption)) {
      caliperType = "user"
    } else {
      caliperOption = 0
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
      balanceCosts(net.i, 1)
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
    
    
    ## Major change in version 1.0.1
    my.stats1_tmp <-
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
    i = 1
    f1_sum_tmp <- c()
    f3_sum_tmp <- c()
    for (ind_tmp in names(match.list)) {
      x <- solutions[[ind_tmp]]
      f1_sum_tmp[[ind_tmp]] <- sum(paircost.v * x[1:length(paircost.v)])
      f3_sum_tmp[[ind_tmp]] <- my.stats1_tmp[1, i]
      i = i + 1
    }
    
    
    tempRhos <- as.numeric(names(match.list))
    percentageUnmatched <-
      1 - as.vector(laply(match.list, nrow) / sum(dat[, treatCol]))
    minInd <- which.min(percentageUnmatched)
    maxInd <- which.max(percentageUnmatched)
    
    
    ## Major change in version 1.0.1 
    minInd_dist <- which.min(f1_sum_tmp)
    maxInd_dist <- which.max(f1_sum_tmp)
    minInd_bal <- which.min(f3_sum_tmp)
    maxInd_bal <- which.max(f3_sum_tmp)
    
    bestRhoInd <- tempRhos[minInd]
    bestRhoInd_max <- tempRhos[maxInd]
    bestPercentageSoFar <- as.vector(percentageUnmatched)[minInd]
    bestPercentageSoFar_max <- as.vector(percentageUnmatched)[maxInd]
    bestRho1 <- rho_list[[bestRhoInd]][1]
    bestRho2 <- rho_list[[bestRhoInd]][2]
    bestRho1_max <- rho_list[[bestRhoInd_max]][1]
    bestRho2_max <- rho_list[[bestRhoInd_max]][2]
    ## Major change in version 1.0.1: do the same for other objectives as well
    
    bestRhoInd_dist <- tempRhos[minInd_dist]
    bestRhoInd_dist_max <- tempRhos[maxInd_dist]
    bestRho1_dist <- rho_list[[bestRhoInd_dist]][1]
    bestRho2_dist <- rho_list[[bestRhoInd_dist]][2]
    bestRho1_dist_max <- rho_list[[bestRhoInd_dist_max]][1]
    bestRho2_dist_max <- rho_list[[bestRhoInd_dist_max]][2]
    
    bestRhoInd_bal <- tempRhos[minInd_bal]
    bestRhoInd_bal_max <- tempRhos[maxInd_bal]
    bestRho1_bal <- rho_list[[bestRhoInd_bal]][1]
    bestRho2_bal <- rho_list[[bestRhoInd_bal]][2]
    bestRho1_bal_max <- rho_list[[bestRhoInd_bal_max]][1]
    bestRho2_bal_max <- rho_list[[bestRhoInd_bal_max]][2]
    
    ind_next = length(rho_list) + 1
    while (
      numIter <= maxIter) {
      
      
      
      rho1_search_min <- 0.1 * bestRho1 
      rho1_search_max <- 10 * bestRho1 
      rho2_search_min <- 0.1 * bestRho2 
      rho2_search_max <- 10 * bestRho2 
      r1_range <- rho1_search_max - rho1_search_min
      new_r1s <- seq(rho1_search_min,rho1_search_max, r1_range/3)
      r2_range <- rho2_search_max - rho2_search_min
      new_r2s <- seq(rho2_search_min,rho2_search_max, r2_range/3)
      
      
      rho1_search_min_max <- 0.1 * bestRho1_max 
      rho1_search_max_max <- 10 * bestRho1_max 
      rho2_search_min_max <- 0.1 * bestRho2_max 
      rho2_search_max_max <- 10 * bestRho2_max 
      r1_range_max <- rho1_search_max_max - rho1_search_min_max
      new_r1s_max <- seq(rho1_search_min_max,rho1_search_max_max, r1_range_max/3)
      r2_range_max <- rho2_search_max_max - rho2_search_min_max
      new_r2s_max <- seq(rho2_search_min_max,rho2_search_max_max, r2_range_max/3)
      
      ### Major change in version 1.0.1
      
      #### Do it for dist 
      rho1_search_min <- 0.1 * bestRho1_dist 
      rho1_search_max <- 10 * bestRho1_dist 
      rho2_search_min <- 0.1 * bestRho2_dist 
      rho2_search_max <- 10 * bestRho2_dist
      r1_range <- rho1_search_max - rho1_search_min
      new_r1s_dist <- seq(rho1_search_min,rho1_search_max, r1_range/3)
      r2_range <- rho2_search_max - rho2_search_min
      new_r2s_dist <- seq(rho2_search_min,rho2_search_max, r2_range/3)
      
      
      rho1_search_min_max <- 0.1 * bestRho1_dist_max 
      rho1_search_max_max <- 10 * bestRho1_dist_max 
      rho2_search_min_max <- 0.1 * bestRho2_dist_max 
      rho2_search_max_max <- 10 * bestRho2_dist_max 
      r1_range_max <- rho1_search_max_max - rho1_search_min_max
      new_r1s_dist_max <- seq(rho1_search_min_max,rho1_search_max_max, r1_range_max/3)
      r2_range_max <- rho2_search_max_max - rho2_search_min_max
      new_r2s_dist_max <- seq(rho2_search_min_max,rho2_search_max_max, r2_range_max/3)
      
      
      #### Do it for exclusion 
      rho1_search_min <- 0.1 * bestRho1_bal 
      rho1_search_max <- 10 * bestRho1_bal
      rho2_search_min <- 0.1 * bestRho2_bal
      rho2_search_max <- 10 * bestRho2_bal
      r1_range <- rho1_search_max - rho1_search_min
      new_r1s_bal <- seq(rho1_search_min,rho1_search_max, r1_range/3)
      r2_range <- rho2_search_max - rho2_search_min
      new_r2s_bal <- seq(rho2_search_min,rho2_search_max, r2_range/3)
      
      
      rho1_search_min_max <- 0.1 * bestRho1_bal_max 
      rho1_search_max_max <- 10 * bestRho1_bal_max 
      rho2_search_min_max <- 0.1 * bestRho2_bal_max 
      rho2_search_max_max <- 10 * bestRho2_bal_max 
      r1_range_max <- rho1_search_max_max - rho1_search_min_max
      new_r1s_bal_max <- seq(rho1_search_min_max,rho1_search_max_max, r1_range_max/3)
      r2_range_max <- rho2_search_max_max - rho2_search_min_max
      new_r2s_bal_max <- seq(rho2_search_min_max,rho2_search_max_max, r2_range_max/3)
      
      
      
      
      
      
      
      
      
      
      
      proposedNewRho1s <- unique(c(10 *max(new_r1s),10 *max(new_r1s_dist),
                                   10 *max(new_r1s_bal), 
                                   new_r1s,new_r1s_dist,new_r1s_bal,
                                   0.1 *min(new_r1s_max),new_r1s_max,
                                   0.1 *min(new_r1s_dist_max),new_r1s_dist_max,
                                   0.1 *min(new_r1s_bal_max),new_r1s_bal_max))
      proposedNewRho2s <- unique(c(10 *max(new_r2s),10 *max(new_r2s_dist),
                                   10 *max(new_r2s_bal),new_r2s,
                                   new_r2s_dist,new_r2s_bal,
                                   0.1 *min(new_r2s_max),new_r2s_max,
                                   0.1 *min(new_r2s_dist_max),new_r2s_dist_max,
                                   0.1 *min(new_r2s_bal_max),new_r2s_bal_max))
      
      for (r1 in proposedNewRho1s) {
        for(r2 in proposedNewRho2s){
          rho_list[[ind_next]] = c(r1, r2)
          temp.sol <- solveP1(
            net.i,
            f1list,
            f2list,
            f3list,
            rho1 = r1,
            rho2 = r2,
            tol = toleranceOption
          )
          
          solutions[[as.character(ind_next)]] <- temp.sol$x
          solution.nets[[as.character(ind_next)]] <- temp.sol$net
          print(paste('Matches finished for rho1 = ', 
                      r1, " and rho2 = ", r2))
          if (is.null(temp.sol$x)) {
            print(
              paste(
                'However, rho1 = ',
                r1,
                " and rho2 = ",
                r2,
                " results in empty matching."
              )
            )
          }
          ind_next = ind_next + 1}
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
      
      
      
      ## Major change in version 1.0.1
      my.stats1_tmp <-
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
      i = 1
      f1_sum_tmp <- c()
      f3_sum_tmp <- c()
      for (ind_tmp in names(match.list)) {
        x <- solutions[[ind_tmp]]
        f1_sum_tmp[[ind_tmp]] <- sum(paircost.v * x[1:length(paircost.v)])
        f3_sum_tmp[[ind_tmp]] <- my.stats1_tmp[1, i]
        i = i + 1
      }
      
      
      tempRhos <- as.numeric(names(match.list))
      percentageUnmatched <-
        1 - as.vector(laply(match.list, nrow) / sum(dat[, treatCol]))
      minInd <- which.min(percentageUnmatched)
      maxInd <- which.max(percentageUnmatched)
      
      
      ## Major change in version 1.0.1 
      minInd_dist <- which.min(f1_sum_tmp)
      maxInd_dist <- which.max(f1_sum_tmp)
      minInd_bal <- which.min(f3_sum_tmp)
      maxInd_bal <- which.max(f3_sum_tmp)
      
      bestRhoInd <- tempRhos[minInd]
      bestRhoInd_max <- tempRhos[maxInd]
      bestPercentageSoFar <- as.vector(percentageUnmatched)[minInd]
      bestPercentageSoFar_max <- as.vector(percentageUnmatched)[maxInd]
      bestRho1 <- rho_list[[bestRhoInd]][1]
      bestRho2 <- rho_list[[bestRhoInd]][2]
      bestRho1_max <- rho_list[[bestRhoInd_max]][1]
      bestRho2_max <- rho_list[[bestRhoInd_max]][2]
      ## Major change in version 1.0.1: do the same for other objectives as well
      
      bestRhoInd_dist <- tempRhos[minInd_dist]
      bestRhoInd_dist_max <- tempRhos[maxInd_dist]
      bestRho1_dist <- rho_list[[bestRhoInd_dist]][1]
      bestRho2_dist <- rho_list[[bestRhoInd_dist]][2]
      bestRho1_dist_max <- rho_list[[bestRhoInd_dist_max]][1]
      bestRho2_dist_max <- rho_list[[bestRhoInd_dist_max]][2]
      
      bestRhoInd_bal <- tempRhos[minInd_bal]
      bestRhoInd_bal_max <- tempRhos[maxInd_bal]
      bestRho1_bal <- rho_list[[bestRhoInd_bal]][1]
      bestRho2_bal <- rho_list[[bestRhoInd_bal]][2]
      bestRho1_bal_max <- rho_list[[bestRhoInd_bal_max]][1]
      bestRho2_bal_max <- rho_list[[bestRhoInd_bal_max]][2]
      
      
      
      
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
    for (ind_tmp in names(match.list)) {
      x <- solutions[[ind_tmp]]
      pair_cost_sum[[ind_tmp]] <- sum(paircost.v * x[1:length(paircost.v)])
      if (!is.null(distEuc.i)) {
        pair_cost_euc_sum[[ind_tmp]] <-
          sum(paircostEuc.v * x[1:length(paircostEuc.v)])
      }
      unmatch_ind = total_treated - length(match.list[[ind_tmp]])
      
      
      f1_sum[[ind_tmp]] <- sum(paircost.v * x[1:length(paircost.v)])
      f2_sum[[ind_tmp]] <- unmatch_ind
      f3_sum[[ind_tmp]] <- my.stats1[1, i]
      i = i + 1
      
      
    }
    
    rho_list_final = c()
    for (ind_tmp in names(match.list)) {
      rho_list_final[[ind_tmp]] <- rho_list[[as.numeric(ind_tmp)]]
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
    result$dataTable = dat
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
#' @param dist1_type One of ("euclidean", "robust_mahalanobis", "user") indicating the
#'   type of distance that are used for the first distance objective functions.
#'   NULL by default.
#' @param dist2_type One of ("euclidean", "robust_mahalanobis", "user")  charactor
#'   indicating the type of distance that are used for the second distance
#'   objective functions. NULL by default.
#' @param dist1_matrix (optional) matrix object that represents the distance
#'   matrix using the first distance measure; `dist1_type` 
#'   must be passed in as "user" if dist1_matrix is non-empty
#' @param dist2_matrix (optional) matrix object that represents the distance 
#'   matrix using the second distance measure; `dist2_type` must be passed 
#'   in as "user" if dist2_matrix is non-empty
#' @param data (optional) data frame that contain columns indicating treatment,
#'   outcome and covariates
#' @param treat_col (optional) character, name of the column indicating 
#'   treatment assignment.
#' @param dist1_col (optional) character vector names of the variables used for
#'   calculating covariate distance using first distance measure specified by
#'   dType
#' @param dist2_col (optional) character vector, names of the variables used for
#'   calculating covariate distance using second distance measure specified by
#'   dType1
#' @param exclusion_penalty (optional) numeric vector, penalty values associated
#'   with exclusion. Empty by default, where auto grid search is triggered. 
#' @param dist2_penalty (optional) numeric vector, penalty values associated 
#'   with the distance specified by `dist2_matrix` or `dist2_type`. 
#'   Empty by default, where auto grid search is tiggered.
#' @param marg_bal_col (optional) character, column name of the variable on 
#'   which to evaluate balance.
#' @param exact_col (optional) character vector, names of the variables on which
#'   to match exactly; NULL by default.
#' @param propensity_col character vector, names of columns on which to fit a
#'   propensity score model.
#' @param pscore_name (optional) character, name of the column containing fitted
#'   propensity scores; default is NULL.
#' @param ignore_col (optional) character vector of variable names that should 
#'   be ignored when constructing the internal matching. NULL by default.
#' @param max_unmatched (optional) numeric, maximum proportion of unmatched 
#'   units that can be accepted; default is 0.25.
#' @param caliper_option (optional) numeric, the propensity score caliper value
#'   in standard deviations of the estimated propensity scores; default is NULL,
#'   which is no caliper.
#' @param tol (optional) numeric, tolerance of close match distance;
#'   default is 1e-2.
#' @param max_iter (optional) integer,  maximum number of iterations to use in
#'   searching for penalty combintions that improve the matching; default is 1,
#'   where the algorithm searches for one round.
#' @param rho_max_factor (optional) numeric, the scaling factor used in proposal
#'   for penalties; default is 10.
#' @param max_pareto_search_iter (optional) numeric, the number of tries to 
#' search for the tol that yield pareto optimal solutions; default is 5.
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
#'   dist_bal_match and "Advanced" from two_dist_match. 
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
#' the `User` option in `dist1_type` and/or `dist2_type` and 
#' supplying arguments to
#' `dist1_matrix` and/or `dist2_matrix` respectively) or ask the function 
#' to create Mahalanobis or Euclidean distances on sets of covariates specified 
#' by the `dist1_col` and `dist2_col` arguments. If `dist1_type` or `dist2_type`
#' is not specified, if one of these is set to `user` and the corresponding 
#' `dist1_matrix` argument is not provided, or if one is NOT set to `User` 
#' and the corresponding `dist1_col` argument is not provided, the code would 
#' error out.
#' * User-specified distance matrices passed to `dist1_matrix` or `dist2_matrix`
#' should have row count equal to the number of treated units and column count
#' equal to the number of controls. 
#' * If the `caliper_option` argument is specified, a
#' propensity score caliper will be imposed, forbidding matches between units
#' more than a fixed distance apart on the propensity score.  The caliper will
#' be based either on a user-fit propensity score, identified in the input
#' dataframe by argument `pscore_name`, or by an internally-fit propensity score
#' based on logistic regression against the variables named in `propensity_col`.
#' If `caliper_option` is non-NULL and neither of the other arguments is 
#' specified an error will result. 
#' * `tol` controls the precision at which the
#' objective functions is evaluated. When matching problems are especially large
#' or complex it may be necessary to increase toleranceOption in order to
#' prevent integer overflows in the underlying network flow solver; generally
#' this will be suggested in appropariate warning messages. 
#' * While by default
#' tradeoffs are only assessed at penalty combinations provided by the user, the
#' user may ask for the algorithm to search over additional penalty values in
#' order to identify additional Pareto optimal solutions. `rho_max_factor` is a
#' multiplier applied to initial penalty values to discover new solutions, and
#' setting it larger leads to wider exploration; similarly, `max_iter` controls
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
#' treated_units <- idx[z==1]
#' control_units <- idx[z==0]
#' d1 <- as.matrix(d1[treated_units, control_units])
#' d2 <- as.matrix(d2[treated_units, control_units])
#' match_result_1 <- two_dist_match(data=df, treat_col="z",  dist1_matrix=d1, 
#' dist1_type= "User", dist2_matrix=d2,
#' dist2_type="User", marg_bal_col=c("x5"), exclusion_penalty=r1ss, 
#' dist2_penalty=r2ss,
#' propensity_col = c("x1"), max_iter = 0,
#' max_pareto_search_iter = 0) 
two_dist_match <-
  function(dist1_type = "user",
           dist2_type = "user",
           dist1_matrix = NULL,
           data = NULL,
           dist2_matrix = NULL,
           treat_col = NULL,
           dist1_col = NULL,
           dist2_col = NULL,
           exclusion_penalty = c(),
           dist2_penalty = c(),
           marg_bal_col = NULL,
           exact_col = NULL,
           propensity_col = NULL,
           pscore_name = NULL,
           ignore_col = NULL,
           max_unmatched = 0.25,
           caliper_option = NULL,
           tol = 1e-2,
           max_iter = 1,
           rho_max_factor = 10,
           max_pareto_search_iter=5) {
    
    
    initial_match_result <- two_dist_match_skeleton(dist1_type,
                                                    dist2_type,
                                                    dist1_matrix,
                                                    data, 
                                                    dist2_matrix,
                                                    treat_col,
                                                    dist1_col,
                                                    dist2_col,
                                                    exclusion_penalty,
                                                    dist2_penalty,
                                                    marg_bal_col,
                                                    exact_col,
                                                    propensity_col,
                                                    pscore_name,
                                                    ignore_col,
                                                    max_unmatched,
                                                    caliper_option,
                                                    tol,
                                                    max_iter,
                                                    rho_max_factor)
    ## Base case - when no iteration needed 
    tol_iter_cnt = 1
    while((!check_pareto_optimality(initial_match_result)) && 
          (max_pareto_search_iter >= 1)) {
      max_pareto_search_iter <- max_pareto_search_iter - 1
      initial_match_result <- two_dist_match_skeleton(dist1_type,
                                                      dist2_type,
                                                      dist1_matrix,
                                                      data,
                                                      dist2_matrix,
                                                      treat_col,
                                                      dist1_col,
                                                      dist2_col,
                                                      exclusion_penalty,
                                                      dist2_penalty,
                                                      marg_bal_col,
                                                      exact_col,
                                                      propensity_col,
                                                      pscore_name,
                                                      ignore_col,
                                                      max_unmatched,
                                                      caliper_option,
                                                      tol/(10^tol_iter_cnt),
                                                      max_iter,
                                                      rho_max_factor)
      tol_iter_cnt <- tol_iter_cnt + 1
    }
    if(check_pareto_optimality(initial_match_result)){
      return(initial_match_result)
    } else {
      warning("Some of the results might not be pareto optimal.
              Please decrease tol and try again.")
      return(initial_match_result)
    }
    
  }


two_dist_match_skeleton <-
  function(dist1_type = "user",
           dist2_type = "user",
           dist1_matrix = NULL,
           data = NULL,
           dist2_matrix = NULL,
           treat_col = NULL,
           dist1_col = NULL,
           dist2_col = NULL,
           exclusion_penalty = c(),
           dist2_penalty = c(),
           marg_bal_col = NULL,
           exact_col = NULL,
           propensity_col = NULL,
           pscore_name = NULL,
           ignore_col = NULL,
           max_unmatched = 0.25,
           caliper_option = NULL,
           tol = 1e-2,
           max_iter = 0,
           rho_max_factor = 10) {
    
    ## New in version 1.0.0: change naming 
    df = data
    
    rhoExclude <- exclusion_penalty
    rhoDistance <- dist2_penalty
    rho1 <- rhoExclude
    rho2 <- rhoDistance
    
    distList1 <- dist1_col
    distList2 <- dist2_col
    
    dType1 <- dist1_type
    dType2 <- dist2_type
    
    
    
    if(dType1 == 'robust_mahalanobis'){
      dType1 <- "Mahalanobis"
    } 
    
    if(dType2 == 'robust_mahalanobis'){
      dType2 <- "Mahalanobis"
    } 
    
    if(dType1 == 'euclidean'){
      dType1 <- "Euclidean"
    } 
    
    if(dType2 == 'euclidean'){
      dType2 <- "Euclidean"
    } 
    dMat1 = dist1_matrix
    dMat2 = dist2_matrix
    treatCol = treat_col

    myBalCol = marg_bal_col
    exactlist = exact_col
    propensityCols = propensity_col
    pScores = pscore_name
    ignore = ignore_col
    maxUnMatched = max_unmatched
    caliperOption = caliper_option
    toleranceOption = tol
    maxIter = max_iter
    rho.max.f = rho_max_factor
    
    
    
    
    
    
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
    #result$dataTable <- dat
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
    result$numTreat <- sum(dat$treat)
    ## 3. Define distance and cost
    
    if(is.null(distList1)&& ((dType1 == 'Euclidean') 
                             || dType1 == 'Mahalanobis')){
      distList1 <- xNames 
    }
    
    
    if(is.null(distList2)&& ((dType2 == 'Euclidean') 
                             || dType2 == 'Mahalanobis')){
      distList2 <- xNames 
    }
    
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
    
    #### Version 1.0.1 change: 
    excl.objective.edge.costs <-  excludeCosts(net.i, 1)
    bal.objective.edge.costs <-
      balanceCosts(net.i,1)
    pairEuc.objective.edge.costs <- pairCosts(distEuc.i, net.i)
    f1list = list(pair.objective.edge.costs)
    f2list = list(excl.objective.edge.costs)
    f3list = list(bal.objective.edge.costs)
    f4list = list(pairEuc.objective.edge.costs)
    
    ## Propose possible rho values for multi-objective optimization problem
    #### Major change in the auto-matic grid search and rho proposition 
    #### in version 1.0.1
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
    
    
    ## Major change in version 1.0.1
    i = 1
    f1_sum_tmp <- c()
    f3_sum_tmp <- c()
    for (ind_tmp in names(match.list)) {
      x <- solutions[[ind_tmp]]
      f1_sum_tmp[[ind_tmp]] <- sum(paircost.v * x[1:length(paircost.v)])
      f3_sum_tmp[[ind_tmp]] <- sum(paircostEuc.v * x[1:length(paircost.v)])
      i = i + 1
    }
    
    
    tempRhos <- as.numeric(names(match.list))
    percentageUnmatched <-
      1 - as.vector(laply(match.list, nrow) / sum(dat[, treatCol]))
    minInd <- which.min(percentageUnmatched)
    maxInd <- which.max(percentageUnmatched)
    
    ## Major change in version 1.0.1 
    minInd_dist <- which.min(f1_sum_tmp)
    maxInd_dist <- which.max(f1_sum_tmp)
    minInd_bal <- which.min(f3_sum_tmp)
    maxInd_bal <- which.max(f3_sum_tmp)
    
    bestRhoInd <- tempRhos[minInd]
    bestRhoInd_max <- tempRhos[maxInd]
    bestPercentageSoFar <- as.vector(percentageUnmatched)[minInd]
    bestPercentageSoFar_max <- as.vector(percentageUnmatched)[maxInd]
    bestRho1 <- rho_list[[bestRhoInd]][1]
    bestRho2 <- rho_list[[bestRhoInd]][2]
    bestRho1_max <- rho_list[[bestRhoInd_max]][1]
    bestRho2_max <- rho_list[[bestRhoInd_max]][2]
    ## Major change in version 1.0.1: do the same for other objectives as well
    
    bestRhoInd_dist <- tempRhos[minInd_dist]
    bestRhoInd_dist_max <- tempRhos[maxInd_dist]
    bestRho1_dist <- rho_list[[bestRhoInd_dist]][1]
    bestRho2_dist <- rho_list[[bestRhoInd_dist]][2]
    bestRho1_dist_max <- rho_list[[bestRhoInd_dist_max]][1]
    bestRho2_dist_max <- rho_list[[bestRhoInd_dist_max]][2]
    
    bestRhoInd_bal <- tempRhos[minInd_bal]
    bestRhoInd_bal_max <- tempRhos[maxInd_bal]
    bestRho1_bal <- rho_list[[bestRhoInd_bal]][1]
    bestRho2_bal <- rho_list[[bestRhoInd_bal]][2]
    bestRho1_bal_max <- rho_list[[bestRhoInd_bal_max]][1]
    bestRho2_bal_max <- rho_list[[bestRhoInd_bal_max]][2]
    
    ind_next = length(rho_list) + 1
    while (numIter <= maxIter) {
      
      
      rho1_search_min <- 0.1 * bestRho1 
      rho1_search_max <- 10 * bestRho1 
      rho2_search_min <- 0.1 * bestRho2 
      rho2_search_max <- 10 * bestRho2 
      r1_range <- rho1_search_max - rho1_search_min
      new_r1s <- seq(rho1_search_min,rho1_search_max, r1_range/3)
      r2_range <- rho2_search_max - rho2_search_min
      new_r2s <- seq(rho2_search_min,rho2_search_max, r2_range/3)
      
      
      rho1_search_min_max <- 0.1 * bestRho1_max 
      rho1_search_max_max <- 10 * bestRho1_max 
      rho2_search_min_max <- 0.1 * bestRho2_max 
      rho2_search_max_max <- 10 * bestRho2_max 
      r1_range_max <- rho1_search_max_max - rho1_search_min_max
      new_r1s_max <- seq(rho1_search_min_max,rho1_search_max_max, r1_range_max/3)
      r2_range_max <- rho2_search_max_max - rho2_search_min_max
      new_r2s_max <- seq(rho2_search_min_max,rho2_search_max_max, r2_range_max/3)
      
      ### Major change in version 1.0.1
      
      #### Do it for dist 
      rho1_search_min <- 0.1 * bestRho1_dist 
      rho1_search_max <- 10 * bestRho1_dist 
      rho2_search_min <- 0.1 * bestRho2_dist 
      rho2_search_max <- 10 * bestRho2_dist
      r1_range <- rho1_search_max - rho1_search_min
      new_r1s_dist <- seq(rho1_search_min,rho1_search_max, r1_range/3)
      r2_range <- rho2_search_max - rho2_search_min
      new_r2s_dist <- seq(rho2_search_min,rho2_search_max, r2_range/3)
      
      
      rho1_search_min_max <- 0.1 * bestRho1_dist_max 
      rho1_search_max_max <- 10 * bestRho1_dist_max 
      rho2_search_min_max <- 0.1 * bestRho2_dist_max 
      rho2_search_max_max <- 10 * bestRho2_dist_max 
      r1_range_max <- rho1_search_max_max - rho1_search_min_max
      new_r1s_dist_max <- seq(rho1_search_min_max,rho1_search_max_max, r1_range_max/3)
      r2_range_max <- rho2_search_max_max - rho2_search_min_max
      new_r2s_dist_max <- seq(rho2_search_min_max,rho2_search_max_max, r2_range_max/3)
      
      
      #### Do it for exclusion 
      rho1_search_min <- 0.1 * bestRho1_bal 
      rho1_search_max <- 10 * bestRho1_bal
      rho2_search_min <- 0.1 * bestRho2_bal
      rho2_search_max <- 10 * bestRho2_bal
      r1_range <- rho1_search_max - rho1_search_min
      new_r1s_bal <- seq(rho1_search_min,rho1_search_max, r1_range/3)
      r2_range <- rho2_search_max - rho2_search_min
      new_r2s_bal <- seq(rho2_search_min,rho2_search_max, r2_range/3)
      
      
      rho1_search_min_max <- 0.1 * bestRho1_bal_max 
      rho1_search_max_max <- 10 * bestRho1_bal_max 
      rho2_search_min_max <- 0.1 * bestRho2_bal_max 
      rho2_search_max_max <- 10 * bestRho2_bal_max 
      r1_range_max <- rho1_search_max_max - rho1_search_min_max
      new_r1s_bal_max <- seq(rho1_search_min_max,rho1_search_max_max, r1_range_max/3)
      r2_range_max <- rho2_search_max_max - rho2_search_min_max
      new_r2s_bal_max <- seq(rho2_search_min_max,rho2_search_max_max, r2_range_max/3)
      
      
      
      
      
      
      
      
      
      
      
      proposedNewRho1s <- unique(c(10 *max(new_r1s),10 *max(new_r1s_dist),
                                   10 *max(new_r1s_bal), 
                                   new_r1s,new_r1s_dist,new_r1s_bal,
                                   0.1 *min(new_r1s_max),new_r1s_max,
                                   0.1 *min(new_r1s_dist_max),new_r1s_dist_max,
                                   0.1 *min(new_r1s_bal_max),new_r1s_bal_max))
      proposedNewRho2s <- unique(c(10 *max(new_r2s),10 *max(new_r2s_dist),
                                   10 *max(new_r2s_bal),new_r2s,
                                   new_r2s_dist,new_r2s_bal,
                                   0.1 *min(new_r2s_max),new_r2s_max,
                                   0.1 *min(new_r2s_dist_max),new_r2s_dist_max,
                                   0.1 *min(new_r2s_bal_max),new_r2s_bal_max))
      
      for (r1 in proposedNewRho1s) {
        for(r2 in proposedNewRho2s){
          rho_list[[ind_next]] = c(r1, r2)
          temp.sol <- solveP1(
            net.i,
            f1list,
            f2list,
            f4list,
            rho1 = r1,
            rho2 = r2,
            tol = toleranceOption
          )
          
          solutions[[as.character(ind_next)]] <- temp.sol$x
          solution.nets[[as.character(ind_next)]] <- temp.sol$net
          print(paste('Matches finished for rho1 = ', 
                      r1, " and rho2 = ", r2))
          if (is.null(temp.sol$x)) {
            print(
              paste(
                'However, rho1 = ',
                r1,
                " and rho2 = ",
                r2,
                " results in empty matching."
              )
            )
          }
          ind_next = ind_next + 1}
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
      
      
      
      ## Major change in version 1.0.1
      i = 1
      f1_sum_tmp <- c()
      f3_sum_tmp <- c()
      for (ind_tmp in names(match.list)) {
        x <- solutions[[ind_tmp]]
        f1_sum_tmp[[ind_tmp]] <- sum(paircost.v * x[1:length(paircost.v)])
        f3_sum_tmp[[ind_tmp]] <- sum(paircostEuc.v * x[1:length(paircostEuc.v)])
        i = i + 1
      }
      
      
      tempRhos <- as.numeric(names(match.list))
      percentageUnmatched <-
        1 - as.vector(laply(match.list, nrow) / sum(dat[, treatCol]))
      minInd <- which.min(percentageUnmatched)
      maxInd <- which.max(percentageUnmatched)
      
      
      ## Major change in version 1.0.1 
      minInd_dist <- which.min(f1_sum_tmp)
      maxInd_dist <- which.max(f1_sum_tmp)
      minInd_bal <- which.min(f3_sum_tmp)
      maxInd_bal <- which.max(f3_sum_tmp)
      
      bestRhoInd <- tempRhos[minInd]
      bestRhoInd_max <- tempRhos[maxInd]
      bestPercentageSoFar <- as.vector(percentageUnmatched)[minInd]
      bestPercentageSoFar_max <- as.vector(percentageUnmatched)[maxInd]
      bestRho1 <- rho_list[[bestRhoInd]][1]
      bestRho2 <- rho_list[[bestRhoInd]][2]
      bestRho1_max <- rho_list[[bestRhoInd_max]][1]
      bestRho2_max <- rho_list[[bestRhoInd_max]][2]
      ## Major change in version 1.0.1: do the same for other objectives as well
      
      bestRhoInd_dist <- tempRhos[minInd_dist]
      bestRhoInd_dist_max <- tempRhos[maxInd_dist]
      bestRho1_dist <- rho_list[[bestRhoInd_dist]][1]
      bestRho2_dist <- rho_list[[bestRhoInd_dist]][2]
      bestRho1_dist_max <- rho_list[[bestRhoInd_dist_max]][1]
      bestRho2_dist_max <- rho_list[[bestRhoInd_dist_max]][2]
      
      bestRhoInd_bal <- tempRhos[minInd_bal]
      bestRhoInd_bal_max <- tempRhos[maxInd_bal]
      bestRho1_bal <- rho_list[[bestRhoInd_bal]][1]
      bestRho2_bal <- rho_list[[bestRhoInd_bal]][2]
      bestRho1_bal_max <- rho_list[[bestRhoInd_bal_max]][1]
      bestRho2_bal_max <- rho_list[[bestRhoInd_bal_max]][2]
      
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
#' matching (produced by `dist_bal_match`).
#'
#' @family numerical analysis helper functions
#' @param matching_result an object returned by the main matching function
#'   dist_bal_match
#' @param cov_list (optional) a vector of names of covariates used for evaluating
#'   covariate imbalance; NULL by default.
#' @param display_all (optional) a boolean value indicating whether or not to
#'   show the data for all possible matches; TRUE by default
#' @param stat_list (optional) a vector of statistics that are calculated for
#'   evaluating the covariate imbalance between treated and control group. The
#'   types that are supported can be found here: \link[cobalt]{bal.tab}.
#'
#' @return a named list object containing covariate balance table and statistics
#'   for numer of units being matched for each match; the names are the
#'   character of index for each match in the `matchResult`.
#'
#' @details The result can be either directly used by indexing into the list, or
#'   post-processing by the function `compare_tables` that summarizes the
#'   covariate balance information in a tidier table. Users can specify the
#'   arguments as follows: * `cov_list`: if it is set of NULL, all the covariates
#'   are included for the covariate balance table; otherwise, only the specified
#'   covariates will be included in the tabular result. * `display_all`: by
#'   default, the summary statistics for each match are included when the
#'   argument is set to TRUE. If the user only wants to see the summary
#'   statistics for matches with varying performance on three different
#'   objective values, the function would only display the matches with number
#'   of treated units being excluded at different quantiles. User can switch to
#'   the brief version by setting the parameter to FALSE. * `stat_list` is the
#'   list of statistics used for measuring balance. The argument is the same as
#'   `stats` argument in \link[cobalt]{bal.tab}, which is the function that is
#'   used for calculating the statistics. By default, only standardized
#'   difference in means is calculated.
#' @export
#'
#' @examples
#' ## Generate matches
#' data("lalonde", package="cobalt")
#' ps_cols <- c("age", "educ", "married", "nodegree", "race")
#' treat_val <- "treat"
#' response_val <- "re78"  
#' pair_dist_val <- c("age", "married", "educ", "nodegree", "race")
#' my_bal_val <- c("race")
#' r1s <- c(0.01,1,2,4,4.4,5.2,5.4,5.6,5.8,6)
#' r2s <- c(0.001)
#' match_result <- dist_bal_match(data=lalonde, treat_col= treat_val, 
#' marg_bal_col = my_bal_val, exclusion_penalty=r1s, balance_penalty=r2s, 
#' dist_col = pair_dist_val, 
#' propensity_col = ps_cols, max_iter=0)
#' 
#' ## Generate summary table for balance 
#' balance_tables <- get_balance_table(match_result)
#' balance_tables_10 <- balance_tables$'10'
get_balance_table <-
  function(matching_result,
           cov_list = NULL,
           display_all = TRUE,
           stat_list = c("mean.diffs")) {
    matchingResult = matching_result 
    covList = cov_list 
    display.all = display_all
    statList = stat_list
    
    if (matchingResult$version != "Basic") {
      stop("This graphing function only 
           works for result return by dist_bal_match")
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
#' @description This function would take the result of `get_balance_table`
#' function and combine the results in a single table. It only works for 'Basic'
#' version of the matching.
#'
#' @param balance_table a named list, which is the result from the function
#'   `get_balance_table`
#'
#' @return a dataframe with combined information
#' @export
#'
#' @examples
#' ## Generate matches
#' data("lalonde", package="cobalt")
#' ps_cols <- c("age", "educ", "married", "nodegree", "race")
#' treat_val <- "treat"
#' response_val <- "re78"  
#' pair_dist_val <- c("age", "married", "educ", "nodegree", "race")
#' my_bal_val <- c("race")
#' r1s <- c(0.01,1,2,4,4.4,5.2,5.4,5.6,5.8,6)
#' r2s <- c(0.001)
#' match_result <- dist_bal_match(data=lalonde, treat_col= treat_val, 
#' marg_bal_col = my_bal_val, exclusion_penalty=r1s, balance_penalty=r2s, 
#' dist_col = pair_dist_val, 
#' propensity_col = ps_cols, max_iter=0)
#' 
#' ## Generate summary table for comparing matches 
#' compare_tables(get_balance_table(match_result))
compare_tables <- function(balance_table) {
  balanceTable = balance_table
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
#' balance across different matches. It only works for 'Basic'
#' version of matching (using `dist_bal_match`).
#'
#'
#' @param matching_result an object returned by the main matching function
#'   dist_bal_match
#' @param cov_list (optional) factor of names of covariates that we want to
#'   evaluate covariate balance on; default is NULL. When set to NULL, the
#'   program will compare the covariates that have been used to construct a
#'   propensity model.
#' @param display_all (optional) boolean value of whether to display all the
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
#' ps_cols <- c("age", "educ", "married", "nodegree", "race")
#' treat_val <- "treat"
#' response_val <- "re78"  
#' pair_dist_val <- c("age", "married", "educ", "nodegree", "race")
#' my_bal_val <- c("race")
#' r1s <- c(0.01,1,2,4,4.4,5.2,5.4,5.6,5.8,6)
#' r2s <- c(0.001)
#' match_result <- dist_bal_match(data=lalonde, treat_col= treat_val, 
#' marg_bal_col = my_bal_val, exclusion_penalty=r1s, balance_penalty=r2s, 
#' dist_col = pair_dist_val, 
#' propensity_col = ps_cols, max_iter=0)
#' 
#' ## Generate table for comparing matches
#' compare_matching(match_result, display_all = TRUE)
compare_matching <-
  function(matching_result,
           cov_list = NULL,
           display_all = TRUE,
           stat = "mean.diff") {
    matchingResult = matching_result
    covList = cov_list
    display.all = display_all
    if (matchingResult$version != "Basic") {
      stop("This graphing function only works 
           for result return by dist_bal_match")
    }
    return(compare_tables(
      get_balance_table(
        matchingResult,
        covList,
        display_all = display.all,
        stat_list = c(stat)
      )
    ))
  }


#' An internal helper function that gives the index of matching with a wide
#' range of number of treated units left unmatched
#'
#' @param matching_result an object returned by the main matching function
#'   dist_bal_match
#'
#' @return a vector of five matching indices with the number of treated units
#'   excluded at 0th, 25th, 50th, 75th and 100th percentiles respectively.
get_five_index <- function(matching_result) {
  matchingResult = matching_result
  df <- compare_matching(matchingResult)
  matches <- colnames(df)[2:length(colnames(df))]
  return(matches)
}



#' Marginal imbalance vs. exclusion
#' @description Plotting function that visualizes the tradeoff between the total
#' variation imbalance on a specified variable and the number of unmatched
#' treated units. This function only works for the 'Basic' version of matching
#' (conducted using `dist_bal_match`).
#'
#' @family Graphical helper functions for analysis
#' @param matching_result an object returned by the main matching function
#'   dist_bal_match
#'
#' @return No return value, called for visualization of match result
#' @export
#' @examples
#' ## Generate matches 
#' data("lalonde", package="cobalt")
#' ps_cols <- c("age", "educ", "married", "nodegree", "race")
#' treat_val <- "treat"
#' response_val <- "re78"  
#' pair_dist_val <- c("age", "married", "educ", "nodegree", "race")
#' my_bal_val <- c("race")
#' r1s <- c(0.01,1,2,4,4.4,5.2,5.4,5.6,5.8,6)
#' r2s <- c(0.001)
#' match_result <- dist_bal_match(data=lalonde, treat_col= treat_val, 
#' marg_bal_col = my_bal_val, exclusion_penalty=r1s, balance_penalty=r2s, 
#' dist_col = pair_dist_val, 
#' propensity_col = ps_cols, max_iter=0)
#' ## Generate visualization of tradeoff between total variation distance and 
#' ## number of treated units left unmatched
#' get_tv_graph(match_result)
get_tv_graph <- function(matching_result) {
  matchingResult = matching_result
  if (matchingResult$version != "Basic") {
    stop("This graphing function only works for result return by dist_bal_match")
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
#' @param matching_result an object returned by the main matching function
#'   dist_bal_match
#'
#' @return No return value, called for visualization of match result
#' @export
#' @examples
#' ## Generate matches 
#' data("lalonde", package="cobalt")
#' ps_cols <- c("age", "educ", "married", "nodegree", "race")
#' treat_val <- "treat"
#' response_val <- "re78"  
#' pair_dist_val <- c("age", "married", "educ", "nodegree", "race")
#' my_bal_val <- c("race")
#' r1s <- c(0.01,1,2,4,4.4,5.2,5.4,5.6,5.8,6)
#' r2s <- c(0.001)
#' match_result <- dist_bal_match(data=lalonde, treat_col= treat_val, 
#' marg_bal_col = my_bal_val, exclusion_penalty=r1s, balance_penalty=r2s, 
#' dist_col = pair_dist_val, 
#' propensity_col = ps_cols, max_iter=0)
#' ## Generate visualization of tradeoff between pari-wise distance sum and 
#' ## number of treated units left unmatched
#' get_pairdist_graph(match_result)
get_pairdist_graph <- function(matching_result) {
  matchingResult = matching_result
  if (matchingResult$version != "Basic") {
    stop("This graphing function only works for result return by dist_bal_match")
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
#' works for 'Basic' version of matching (conducted using `dist_bal_match`).
#'
#' @family Graphical helper functions for analysis
#' @param matching_result an object returned by the main matching function
#'   dist_bal_match
#'
#' @return No return value, called for visualization of match result
#' @export
#' @examples
#' ## Generate matches 
#' data("lalonde", package="cobalt")
#' ps_cols <- c("age", "educ", "married", "nodegree", "race")
#' treat_val <- "treat"
#' response_val <- "re78"  
#' pair_dist_val <- c("age", "married", "educ", "nodegree", "race")
#' my_bal_val <- c("race")
#' r1s <- c(0.01,1,2,4,4.4,5.2,5.4,5.6,5.8,6)
#' r2s <- c(0.001)
#' match_result <- dist_bal_match(data=lalonde, treat_col= treat_val, 
#' marg_bal_col = my_bal_val, exclusion_penalty=r1s, balance_penalty=r2s, 
#' dist_col = pair_dist_val, 
#' propensity_col = ps_cols, max_iter=0)
#' ## Visualize the tradeoff between the pair-wise distance sum and 
#' ## total variation distance 
#' get_pairdist_balance_graph(match_result)
get_pairdist_balance_graph <- function(matching_result) {
  matchingResult = matching_result
  if (matchingResult$version != "Basic") {
    stop("This graphing function only works for result return by dist_bal_match")
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
#' @param matching_result an object returned by the main matching function
#'   dist_bal_match
#'
#' @return NULL
convert_index <- function(matching_result) {
  matchingResult = matching_result
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
#'   dist_bal_match
#'
#' @return NULL
matched_index <- function(matchingResult) {
  return(convert_index(matchingResult))
}


#' Get matched dataframe
#' @description A function that returns the dataframe that contains only matched
#' pairs from the original data frame with specified match index
#'
#' @param matching_result an object returned by the main matching function
#'   dist_bal_match
#' @param match_num Integer index of match that the user want to extract paired
#'   observations from
#'
#' @return dataframe that contains only matched pair data
#' @export
#' @examples
#' ## Generate Matches
#' data("lalonde", package="cobalt")
#' ps_cols <- c("age", "educ", "married", "nodegree", "race")
#' treat_val <- "treat"
#' response_val <- "re78"  
#' pair_dist_val <- c("age", "married", "educ", "nodegree", "race")
#' my_bal_val <- c("race")
#' r1s <- c(0.01,1,2,4,4.4,5.2,5.4,5.6,5.8,6)
#' r2s <- c(0.001)
#' match_result <- dist_bal_match(data=lalonde, treat_col= treat_val, 
#' marg_bal_col = my_bal_val, exclusion_penalty=r1s, balance_penalty=r2s, 
#' dist_col = pair_dist_val, 
#' propensity_col = ps_cols, max_iter=0)
#' matched_data(match_result, 1)
matched_data <- function(matching_result, match_num){
  matchingResult = matching_result
  df_tmp <- matchingResult$dataTable
  total_treated <- sum(df_tmp[matchingResult$treatmentCol])
  #df_tmp$tempSorting = df_tmp[,matchingResult$treatmentCol]
  #df_tmp = df_tmp[order(-df_tmp$tempSorting),]
  df_match_1 <- data.frame(matchingResult$matchList[[toString(match_num)]])
  df_treated <- df_tmp[as.character(rownames(df_match_1)), ] 
  df_treated$matchID <- 1:nrow(df_treated)
  df_control <- df_tmp[as.character(df_match_1$X1 + total_treated), ]
  df_control$matchID <- 1:nrow(df_control)
  res_df <- rbind(df_treated, df_control)
  if(matchingResult$version != "Basic"){
    res_df$originalID <-res_df$originalID - 
      total_treated *(1-res_df[['pseudo-treat']])
  }
  return(res_df)
}

#' Get unmatched percentage
#' @description A function that generate the percentage of unmatched units for
#' each match.
#'
#' @family numerical analysis helper functions
#' @param matching_result matchingResult object that contains information for all
#'   matches
#'
#' @return data frame with three columns, one containing the matching index, one
#'   containing the number of matched units, and one conatining the percentage
#'   of matched units (out of original treated group size).
#' @export
#' @examples
#' \dontrun{
#' get_unmatched(match_result)
#' }
get_unmatched <- function(matching_result) {
  matchingResult = matching_result
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




#' Penalty and objective values summary
#' @description Helper function to generate a dataframe with matching number,
#' penalty (rho) values, and objective function values.
#'
#' @family numerical analysis helper functions
#' @param matching_result matchingResult object that contains information for all
#'   matches.
#'
#' @return a dataframe that contains objective function values and rho values
#'   corresponding coefficients before each objective function.
#' @export
get_rho_obj <- function(matching_result){
  result <- generateRhoObj(matching_result)
  ## Convert all the columns to numeric data type 
  result[] <- lapply(result, as.numeric)
  result$match_index <- as.character(result$match_index)
  return(result)
  
  }


#' Internal helper function that converts axis name to internal variable name
#'
#' @param x the user input character for x-axis value
#' @param y the user input character for y-axis value
#' @param z the user input character for z-axis value
#'
#' @return a named list with variable names for visualization for internal use

convert_names <- function(x, y,z=NULL) {
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
  } else if(x == "distance_penalty"){
    result$x_label <- "rhoPair"
    result$x_name <- "Pair Distance Penalty"
    
  } else if(x == "exclusion_penalty"){
    result$x_label <- "rhoExclude"
    result$x_name <- "Exclusion Penalty"
  } else if(x == "dist1_penalty"){
    result$x_label <- "rhoDist1"
    result$x_name <- "Distance1 Penalty"
  } else if(x == "dist2_penalty"){
    result$x_label <- "rhoDist2"
    result$x_name <- "Distance2 Penalty"
  } else if(x == "balance_penalty"){
    result$x_label <- "rhoMarginal"
    result$x_name <- "Balance Penalty"
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
  } else if(y == "distance_penalty"){
    result$y_label <- "rhoBalance"
    result$y_name <- "Pair Distance Penalty"
    
  } else if(y == "exclusion_penalty"){
    result$y_label <- "rhoExclude"
    result$y_name <- "Exclusion Penalty"
  } else if(y == "dist1_penalty"){
    result$y_label <- "rhoDist1"
    result$y_name <- "Distance1 Penalty"
  } else if(y == "dist2_penalty"){
    result$y_label <- "rhoDist2"
    result$y_name <- "Distance2 Penalty"
  } else if(y == "balance_penalty"){
    result$y_label <- "rhoMarginal"
    result$y_name <- "Balance Penalty"
  }
  else{
    stop(
      "Wrong argument for objective function name.
      It must be one of the following:
         dist1, dist2, exclude, distPair, distMarginal."
    )
  }
  
  
  if(!is.null(z)){
    
    if (z == "dist1") {
      result$z_label <- "fDist1"
      result$z_name <- "Distance 1"
    } else if (z == "exclude") {
      result$z_label <- "fExclude"
      result$z_name <- "Number of Treated Units Unmatched"
    } else if (z == "dist2") {
      result$z_label <- "fDist2"
      result$z_name <- "Distance 2"
    } else if (z == "pair") {
      result$z_label <- "fPair"
      result$z_name <- "Pairwise Distance"
    } else if (z == "marginal") {
      result$z_label <- "fMarginal"
      result$z_name <- "Marginal Balance"
    } else {
      stop(
        "Wrong argument for objective function name.
        It must be one of the following:
           dist1, dist2, exclude, distPair, distMarginal."
      )
    }
    
  }
  
  
  
  return(result)
  
}


visualize_3d <-
  function(matching_result,
           x_axis = "dist1",
           y_axis = "dist2",
           z_axis = "exclude",
           xlab = NULL,
           ylab = NULL,
           zlab = NULL,
           main = NULL,
           display_all = FALSE,
           cond = NULL, 
           xlim = NULL,
           ylim = NULL,
           display_index = TRUE,
           average_cost = FALSE) {
    
    matchingResult = matching_result
    if (x_axis == "dist1" &&
        y_axis == "dist2" && matchingResult$version == "Basic") {
      x_axis = "pair"
      y_axis = "marginal"
    }
    if (!(x_axis %in% c("pair", "marginal", "dist1", "dist2", "exclude", 
                       "distance_penalty", "exclusion_penalty", "balance_penalty",
                       "dist1_penalty", "dist2_penalty"))) {
      stop("Wrong name for argument x_axis")
    }
    if (!(y_axis %in% c("pair", "marginal", "dist1", "dist2", "exclude",
                       "distance_penalty", "exclusion_penalty", "balance_penalty",
                       "dist1_penalty", "dist2_penalty"))) {
      stop("Wrong name for argument y_axis")
    }
    if (!(z_axis %in% c("pair", "marginal", "dist1", "dist2", "exclude"))) {
      stop("Wrong name for argument z_axis")
    }
    
    if (((!x_axis %in% c("dist1", "dist2", "exclude","distance_penalty", "exclusion_penalty", "balance_penalty",
                         "dist1_penalty", "dist2_penalty")) ||
         (!y_axis %in% c("dist1", "dist2", "exclude","distance_penalty", "exclusion_penalty", "balance_penalty",
                         "dist1_penalty", "dist2_penalty")) ||
         (!z_axis %in% c("dist1", "dist2", "exclude")))
        &&
        matchingResult$version == "Advanced") {
      stop("Please choose among 'dist1', 'dist2', 'exclude' for axis.")
    }
    
    if (((!x_axis %in% c("pair", "marginal", "exclude","distance_penalty", "exclusion_penalty", "balance_penalty",
                         "dist1_penalty", "dist2_penalty")) ||
         (!y_axis %in% c("pair", "marginal", "exclude","distance_penalty", "exclusion_penalty", "balance_penalty",
                         "dist1_penalty", "dist2_penalty")) ||
         (!z_axis %in% c("pair", "marginal", "exclude"))) &&
        matchingResult$version == "Basic") {
      stop("Please choose among 'pair', 'marginal', 'exclude'for axis.")
    }
    
    
    
    naming <- convert_names(x_axis, y_axis, z_axis)
    x_axis <- naming$x_label
    y_axis <- naming$y_label
    z_axis <- naming$z_label
    
    if (is.null(xlab)) {
      xlab <- naming$x_name
    }
    if (is.null(ylab)) {
      ylab <- naming$y_name
    }
    
    if (is.null(zlab)) {
      zlab <- naming$z_name
    }
    
    rho_obj_table <- get_rho_obj(matchingResult)
    if (!is.null(cond)) {
      rho_obj_table <- rho_obj_table[cond,]
    }
    
    
    sorted_table <- rho_obj_table[order(rho_obj_table[[naming$z_label]]), ]
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
    
    if(display_index == FALSE){
      for(i in 1:length(graph_labels)){
        graph_labels[i] = " "
      }
    }
    
    
    f_x <- as.numeric(as.vector(sorted_table[[x_axis]]))
    f_y <- as.numeric(as.vector(sorted_table[[y_axis]]))
    f_z <- as.numeric(as.vector(sorted_table[[z_axis]]))
    
    if(average_cost==TRUE){
      totalMatchedTreated = matchingResult$numTreat - 
        as.numeric(sorted_table$fExclude)
      if(x_axis!='fExclude'){
        f_x = f_x / totalMatchedTreated 
      }
      if(y_axis!='fExclude'){
        f_y = f_y / totalMatchedTreated 
      }
      if(z_axis!='fExclude'){
        f_z = f_z / totalMatchedTreated 
      }
    }
    
    if(is.null(xlim)){
      my_xlim <- c(0, max(f_x))
    } else {
      my_xlim <- xlim
    }
    
    if(is.null(ylim)){
      my_ylim <- c(0, max(f_y))
    } else {
      my_ylim <- ylim
    }

    
    my_graph_table <- sorted_table 
    my_graph_table$myx_inMOM <- f_x
    my_graph_table$myy_inMOM <- f_y
    my_graph_table$myz_inMOM <- f_z
    
    my_graph_table$mylabel_inMOM <- graph_labels
    ggplot() +
      geom_point(data = my_graph_table, aes_string(x = "myx_inMOM", y = "myy_inMOM", color = "myz_inMOM")) +
      ylim(my_ylim[1], my_ylim[2]) +
      xlim(my_xlim[1], my_xlim[2]) +
      labs(x = xlab, y = ylab) +
      labs(color = sub(" ", "\n", zlab)) +
      theme(text = element_text(size = 8)) +
      scale_color_gradient(low = "blue", high = "red") +
      geom_text(data = my_graph_table, aes_string(x = "myx_inMOM", y = "myy_inMOM", label = "mylabel_inMOM"), vjust = "inward", hjust = "inward", colour = "black")
  }


#' Visualize tradeoffs
#' @description Main visualization functions for showing the tradeoffs between
#' two of the three objective functions. A 3-d plot can be visualized 
#' where the third dimension is represented by coloring of the dots.
#'
#' @param matching_result the matching result returned by either dist_bal_match or
#'   two_dist_match.
#' @param x_axis character, naming the objective function shown on x-axis; one
#'   of ("pair", "marginal", "dist1", "dist2", "exclude", "distance_penalty", "balance_penalty",
#'   "dist1_penalty", "dist2_penalty", "exclusion_penalty"), "dist1" by default.
#' @param y_axis character, naming the objective function shown on y-axis; one
#'   of ("pair", "marginal", "dist1", "dist2", "exclude", "distance_penalty", "balance_penalty",
#'   "dist1_penalty", "dist2_penalty", "exclusion_penalty"), "dist1" by default.
#' @param z_axis character, naming the objective function for coloring; one
#'   of ("pair", "marginal", "dist1", "dist2", "exclude"), "exclude" by default.
#' @param xlab (optional) the axis label for x-axis; NULL by default.
#' @param ylab (optional) the axis label for y-axis; NULL by default.
#' @param zlab (optional) the axis label for z-axis; NULL by default.
#' @param main (optional) the title of the graph; NULL by default.
#' @param display_all (optional) whether to show all the labels for match index;
#'   FALSE by default, which indicates the visualization function only labels
#'   matches at quantiles of number of treated units being excluded.
#' @param cond (optional) NULL by default, which denotes all the matches are
#'   shown; otherwise, takes a list of boolean values indicating whether to
#'   include each match
#' @param xlim (optional) NULL by default; function automatically takes the max 
#'   of the first objective function values being plotted on x-axis; 
#'   if specified otherwise, pass in the numeric vector 
#'   c(lower_bound, upper_bound)
#' @param ylim (optional) NULL by default; function automatically takes the max 
#'   of the first objective function values being plotted on y-axis; 
#'   if specified otherwise, pass in the numeric vector 
#'   c(lower_bound, upper_bound)
#' @param display_index (optional) TRUE by default; whether to display match 
#'   index   
#' @param average_cost (optional) FALSE by default; whether to show mean cost 
#'     
#' @return No return value, called for visualization of match result
#' @details   By default, the plotting function will show the tradeoff between
#'   the first distance objective function and the marginal balance (if
#'   dist_bal_match) is used; or simply the second distance objective function, if
#'   two_dist_match is used.
#'   
#' @importFrom ggplot2 ggplot aes aes_string geom_point labs theme element_text scale_color_gradient 
#'   geom_text xlim ylim
#' @importFrom rlang .data parse_expr
#' @export
visualize <-
  function(matching_result,
           x_axis = "dist1",
           y_axis = "dist2",
           z_axis = NULL,
           xlab = NULL,
           ylab = NULL,
           zlab = NULL,
           main = NULL,
           display_all = FALSE,
           cond = NULL, 
           xlim = NULL,
           ylim = NULL,
           display_index = TRUE,
           average_cost = FALSE) {
    if(!is.null(z_axis)){
      visualize_3d(matching_result,
                   x_axis,
                   y_axis,
                   z_axis,
                   xlab,
                   ylab,
                   zlab,
                   main,
                   display_all,
                   cond, 
                   xlim,
                   ylim,
                   display_index,
                   average_cost)
    } else {
    
    
    matchingResult = matching_result
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
      rho_obj_table <- rho_obj_table[cond,]
    }
    sorted_table <- rho_obj_table %>% arrange(as.numeric(fExclude))
    
    
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
    
    if(display_index == FALSE){
      for(i in 1:length(graph_labels)){
        graph_labels[i] = " "
      }
    }
    
    
    f_x <- as.numeric(as.vector(sorted_table[[x_axis]]))
    f_y <- as.numeric(as.vector(sorted_table[[y_axis]]))
    if(average_cost==TRUE){
      totalMatchedTreated = matchingResult$numTreat - 
        as.numeric(sorted_table$fExclude)
      if(x_axis!='fExclude'){
        f_x = f_x / totalMatchedTreated 
      }
      if(y_axis!='fExlude'){
        f_y = f_y / totalMatchedTreated 
      }
    }
    
    if(is.null(xlim)){
      my_xlim <- c(0, max(f_x))
    } else {
      my_xlim <- xlim
    }
    
    if(is.null(ylim)){
      my_ylim <- c(0, max(f_y))
    } else {
      my_ylim <- ylim
    }
    
    
    plot(
      f_x,
      f_y,
      pch = 20,
      xlab = xlab,
      ylab = ylab,
      main = main,
      cex = 1.2,
      xlim = my_xlim,
      ylim = my_ylim,
      cex.lab = 1.2,
      cex.axis = 1.2,
      col = rgb(0, 0, 0, 0.3)
    )
    abline(h = 0, lty = 'dotted')
    points(f_x, f_y, col = rgb(0, 0, 1, 1), pch = 20)
    text(f_x, f_y, labels = graph_labels, pos = 2)
  }
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




#' Generate numerical summary 
#' @description Main summary functions for providing tables of numerical 
#' information in matching penalties, objective function values, and balance.
#' @param object the matching result returned by either dist_bal_match or
#'   two_dist_match.
#' @param type (optional) the type of the summary result in c("penalty", 
#'   "exclusion", "balance"). When "penalty" is passed in, the objective 
#'   function values and the penalty values are displayed for each match;
#'   when "exclusion" is passed in, the number of units being matched is  
#'   displayed for each match; when "balance" is passed in, the covariate
#'   the covariate balance table from bal.tab function in cobalt function 
#'   is displayed and user can change `covList` to specify the variables to 
#'   examine. "penalty" by default. 
#' @param cov_list (optional) factor of names of covariates that we want to
#'   evaluate covariate balance on if "balance" is passed in for `type`; 
#'   default is NULL. When set to NULL, the
#'   program will compare the covariates that have been used to construct a
#'   propensity model.
#' @param display_all (optional) boolean value of whether to display all the
#'   matches if "balance" is passed in for `type`; default is TRUE, 
#'   where all matches are displayed.
#' @param stat (optional) character of the name of the statistic used for
#'   measuring covariate balance  if "balance" is passed in for `type`; 
#'   default is "mean.diff". This argument is the
#'   same as used in "cobalt" package, see: \link[cobalt]{bal.tab}
#' @param \dots ignored.
#'
#' @return a summary dataframe of the corresponding type.
#' @exportS3Method summary multiObjMatch
summary.multiObjMatch <- function(object, type="penalty", cov_list = NULL,
                    display_all = TRUE,
                    stat = "mean.diff",...){
  matchingResult = object 
  covList = cov_list 
  display.all = display_all
  if(type == 'penalty'){
    return(get_rho_obj(matchingResult))
  }else if(type == 'exclusion'){
    return(get_unmatched(matchingResult))
  }else if(type == 'balance'){
    return(compare_matching(matchingResult, covList,
                           display.all, stat))
  }else{
    stop("Please choose the type of summary in one of the 
         following: 'penalty', exclusion', or 'balance'.")
  }

}



#' Check the representativeness of matched treated units 
#' @description Summary function to compare SMD of the key covariates in matched 
#'   and the full set of treated units.
#' @param matching_result the matching result returned by either dist_bal_match 
#'   or two_dist_match.
#' @param match_num (optional) Integer index of match that the user want to extract paired
#'   observations from. NULL by default, which will generate a table for all the matches.
#' @return a summary table of SMDs of the key covariates between the whole 
#'   treated units and the matched treated units.
#' @export
#' @importFrom dplyr bind_cols
check_representative <- function(matching_result, match_num=NULL){
  if(is.null(match_num)){
    vec_names <- as.numeric(names(matching_result$matchList))
    rnames <- rownames(check_representative(matching_result, vec_names[1]))
    func_x <- function(x){
      check_representative_helper(matching_result, x)
    }
    res <- lapply(   vec_names, func_x)
    df_final <- bind_cols(lapply(res, as.data.frame.list))
    rownames(df_final) <- rnames
    colnames(df_final) <- vec_names
    return(df_final)
    
  } else {
    return(check_representative_helper(matching_result, match_num))
    
  }
  
}




check_representative_helper <- function(matching_result, match_num) {
  match_index = match_num
  matchResult = matching_result
  original_df <- matchResult$dataTable 
  df_original_treat <- original_df[original_df[[matchResult$treatmentCol]]==1,]
  df_original_treat$is_matched_for_analysis <- 1
  res_match <- matched_data(matchResult, match_num)
  df_match_treat <- res_match[res_match[[matchResult$treatmentCol]]==1,]
  df_match_treat$is_matched_for_analysis <- 0
  df_combined <- rbind(df_original_treat[,c("is_matched_for_analysis", 
                                            matchResult$covs)],
                       df_match_treat[,c("is_matched_for_analysis",
                                         matchResult$covs)])
  
  return(bal.tab(is_matched_for_analysis~., data=df_combined)[1]$Balance[2])
  
}






#' Combine two matching result 
#'
#' @param matching_result1 the first matching result object. 
#' @param matching_result2 the second matching result object.
#'
#' @return a new matching result combining two objects. Note that 
#' the matching index for the second matching is the original name 
#' plus the maximum match index in the first matching object. 
#' @export
combine_match_result <- function(matching_result1, matching_result2){
  
  rho1_names <- names(matching_result1$rhoList) 
  rho2_names <- names(matching_result2$rhoList)
  max_rho1 <- max(as.numeric(rho1_names))
  new_rhonames_2 <- 1:length(rho2_names) + max_rho1
  new_rhonames_2 <- as.character(new_rhonames_2)
  res_new_1 <- rep(matching_result1)
  res_new_2 <- rep(matching_result2)
  
  for(n in names(res_new_2)){
    if((typeof(res_new_2[[n]]) == "list")
       && (n != 'dataTable') && (n != 'df')){
      names(res_new_2[[n]]) <- new_rhonames_2
    }
  }
  for(n in names(res_new_1)){
    if((typeof(res_new_1[[n]]) == "list") && (n != 'dataTable')){
      res_new_1[[n]] <- c(res_new_1[[n]], res_new_2[[n]])
    }
  }
  
  result <- structure(res_new_1, class = 'multiObjMatch')
  return(result)
}



#' Filter match result 
#'
#' @param matching_result the matching result object. 
#' @param filter_expr character, the filtering condition based on the
#' summary table returned by get_rho_obj. 
#'
#' @return the filtered match result
#' @export
#' @importFrom dplyr filter
#' @importFrom rlang parse_expr
filter_match_result <- function(matching_result,filter_expr){
  summ_df <- get_rho_obj(matching_result)
  summ_df_filtered <- summ_df %>% filter(eval(parse_expr(filter_expr))) 
  match_names <- summ_df_filtered$match_index
  if(length(match_names) == 0){
    warning("The filtering condition results in 0 matches. 
            The original matching object is returned.")
    return(matching_result)
  }
  filtered_result <- rep(matching_result)
  
  for(n in names(matching_result)){
    if((typeof(filtered_result[[n]]) == "list")
       && (n != 'dataTable') && (n != 'df')){
      new_list_c <- c() 
      for(m in match_names){
        new_list_c[m] <- matching_result[[n]][m]
      }
      filtered_result[[n]] <- new_list_c
      
    }
  }
  
  return(structure(filtered_result, class = 'multiObjMatch'))
}



check_pareto_optimality <- function(matching_result){
  summ_df <- get_rho_obj(matching_result)
  cnames <- colnames(summ_df)
  needed_columns <- c()
  for(c in cnames){
    if(grepl( "f",c) && c!= "fExclude"){
      needed_columns <- c(needed_columns, c)
    }
  }
  
  unique_exclude <- unique(summ_df[['fExclude']])
  optimals <- c()
  for(f_exclude in unique_exclude){
    summ_df_f <- summ_df[summ_df$fExclude == f_exclude,]
    
    
    summ_df_f <- summ_df[summ_df$fExclude == f_exclude,]
    summ_df_f1 <- summ_df_f[order(summ_df_f[[needed_columns[1]]]),]
    summ_df_f2 <- summ_df_f[order(summ_df_f[[needed_columns[2]]]),]
    pareto_optimal_f <- 
      ((!is.unsorted(summ_df_f1[[needed_columns[1]]])) && (!is.unsorted(-summ_df_f1[[needed_columns[2]]]))) || 
      ((!is.unsorted(summ_df_f2[[needed_columns[2]]])) && (!is.unsorted(-summ_df_f2[[needed_columns[1]]])))
    optimals <- c(optimals, pareto_optimal_f)
    
  }
  return(mean(optimals)==1)
}