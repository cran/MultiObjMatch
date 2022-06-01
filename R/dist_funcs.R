

#' Data precheck: Handle missing data(mean imputation) and remove redundant
#' columns; it also adds an NA column for indicating whether it's missing
#' @param X a dataframe that the user initially inputs for matching - dataframe
#'   with covariates
#'
#' @return a dataframe with modified data if necessary
data_precheck <- function(X) {
  if (is.vector(X))
    X <- matrix(X, length(X), 1)
  if (!(length(z) == (dim(X)[1]))) {
    stop("Length of z does not match row count in X")
  }
  if (!(length(exact) == length(z))) {
    stop("Length of exact does not match length of z")
  }
  if (!(all((z == 1) | (z == 0)))) {
    stop("The z argument must contain only 1s and 0s")
  }
  
  if (is.data.frame(X) || is.character(X)) {
    if (!is.data.frame(X))
      X <- as.data.frame(X)
    X.chars <-
      which(laply(X, function(y)
        'character' %in% class(y)))
    if (length(X.chars) > 0) {
      if (verbose)
        print('Character variables found in X, converting to factors.')
      for (i in X.chars) {
        X[, i] <- factor(X[, i])
        
      }
    }
    #if some variables are factors convert to dummies
    X.factors <-
      which(laply(X, function(y)
        'factor' %in% class(y)))
    
    #handle missing data
    for (i in which(laply(X, function(x)
      any(is.na(x))))) {
      if (verbose)
        print(
          paste(
            'Missing values found in column',
            i ,
            'of X; imputing and adding missingness indicators'
          )
        )
      if (i %in% X.factors) {
        #for factors, make NA a new factor level
        X[, i] <- addNA(X[, i])
      } else{
        #for numeric/logical, impute means and add a new indicator for
        #missingness
        X[[paste(colnames(X)[i], 'NA', sep = '')]] <- is.na(X[, i])
        X[which(is.na(X[, i])), i] <- mean(X[, i], na.rm = TRUE)
      }
    }
    for (i in rev(X.factors)) {
      dummyXi <- model.matrix(as.formula(paste('~', colnames(X)[i], '-1')), 
                              data = X)
      X <- cbind(X[, -i], dummyXi)
    }
    
  } else{
    #handle missing data
    for (i in c(1:ncol(X))) {
      if (any(is.na(X[, i]))) {
        X <- cbind(X, is.na(X[, i]))
        colnames(X)[ncol(X)] <- paste(colnames(X)[i], 'NA', sep = '')
        X[which(is.na(X[, i])), i] <- mean(X[, i], na.rm = TRUE)
      }
    }
    
  }
  #get rid of columns that do not vary
  varying <- apply(X, 2, function(x)
    length(unique(x)) > 1)
  if (!all(varying) &&
      verbose)
    print(
      'Constant-value columns found in X, 
      they will not be used to calculate Mahalanobis distance.'
    )
  X <- X[, which(varying), drop = FALSE]
  
  return(X)
}


#' An internal helper function that generates the data abstraction for the edge
#' weights of the main network structure.
#'
#' @param z a vector of treatment and control indicators, 1 for treatment and 0
#'   for control.
#' @param X a data frame or a numeric or logical matrix containing covariate
#'   information for treated and control units. Its row count must be equal to
#'   the length of z.
#' @param distMat a matrix of pair-wise distance specified by the user
#' @param exact an optional vector of the same length as z. If this argument is
#'   specified, treated units will only be allowed to match to control units
#'   that have equal values in the corresponding indices of the exact vector.
#'   For example, to match patients within hospitals only, one could set exact
#'   equal to a vector of hospital IDs for each patient.
#' @param dist.type one of ('propensity','user','none'). If ’propensity’ is
#'   specified (the default option), the function estimates a propensity score
#'   via logistic regression of z on X and imposes a propensity score caliper.
#'   If ’user’ is specified, the user must provide a vector of values on which a
#'   caliper will be enforced using the calip.cov argument. If ’none’ is
#'   specified no caliper is used.
#' @param calip.option a character indicating the type of caliper used
#' @param calip.cov see calip.option.
#' @param caliper a numeric value that gives the size of the caliper when the
#'   user specifies the calip.option argument as ’propensity’ or ’calip.cov’.
#' @param verbose a boolean value whether to print(cat) debug information.
#'   Default: FALSE
#'
#' @return a distance structure used for constructing the main network flow
#'   problem
build.dist.struct <-
  function(z,
           X,
           distMat,
           exact = NULL,
           dist.type = "Mahalanobis",
           calip.option = 'propensity',
           calip.cov = NULL,
           caliper = 0.2,
           verbose = FALSE) {
    cal.penalty <- 100
    
    
    if (is.null(exact))
      exact = rep(1, length(z))
    if (!(calip.option %in% c('propensity', 'user', 'none'))) {
      stop('Invalid calip.option specified.')
    }
    if (is.vector(X))
      X <- matrix(X, length(X), 1)
    if (!(length(z) == (dim(X)[1]))) {
      stop("Length of z does not match row count in X")
    }
    if (!(length(exact) == length(z))) {
      stop("Length of exact does not match length of z")
    }
    if (!(all((z == 1) | (z == 0)))) {
      stop("The z argument must contain only 1s and 0s")
    }
    
    if (is.data.frame(X) || is.character(X)) {
      if (!is.data.frame(X))
        X <- as.data.frame(X)
      X.chars <-
        which(laply(X, function(y)
          'character' %in% class(y)))
      if (length(X.chars) > 0) {
        if (verbose)
          print('Character variables found in X, converting to factors.')
        for (i in X.chars) {
          X[, i] <- factor(X[, i])
          
        }
      }
      #if some variables are factors convert to dummies
      X.factors <-
        which(laply(X, function(y)
          'factor' %in% class(y)))
      
      #handle missing data
      for (i in which(laply(X, function(x)
        any(is.na(x))))) {
        if (verbose)
          print(
            paste(
              'Missing values found in column',
              i ,
              'of X; imputing and adding missingness indicators'
            )
          )
        if (i %in% X.factors) {
          #for factors, make NA a new factor level
          X[, i] <- addNA(X[, i])
        } else{
          #for numeric/logical, impute means and add a new indicator for
          #missingness
          X[[paste(colnames(X)[i], 'NA', sep = '')]] <- is.na(X[, i])
          X[which(is.na(X[, i])), i] <- mean(X[, i], na.rm = TRUE)
        }
      }
      for (i in rev(X.factors)) {
        dummyXi <- model.matrix(as.formula(paste('~', 
                                                 colnames(X)[i], '-1')), data =
                                  X)
        X <- cbind(X[, -i], dummyXi)
      }
      
    } else{
      #handle missing data
      for (i in c(1:ncol(X))) {
        if (any(is.na(X[, i]))) {
          X <- cbind(X, is.na(X[, i]))
          colnames(X)[ncol(X)] <-
            paste(colnames(X)[i], 'NA', sep = '')
          X[which(is.na(X[, i])), i] <- mean(X[, i], na.rm = TRUE)
        }
      }
      
    }
    
    
    #get rid of columns that do not vary if the distance measure is Mahalanobis
    #distance
    if (dist.type == "Mahalanobis") {
      varying <- apply(X, 2, function(x)
        length(unique(x)) > 1)
      if (!all(varying) &&
          verbose)
        print(
          'Constant-value columns found in X, 
          they will not be used to calculate Mahalanobis distance.'
        )
      X <- X[, which(varying), drop = FALSE]
    }
    
    if (calip.option == 'propensity') {
      calip.cov <-
        glm.fit(cbind(rep(1, nrow(X)), X), z, 
                family = binomial())$linear.predictors
      cal <- sd(calip.cov) * caliper
    } else if (calip.option == 'user') {
      stopifnot(!is.null(calip.cov))
      
      cal <- sd(calip.cov) * caliper
      
    }
    
    nobs <- length(z)
    rX <- as.matrix(X)
    if (dist.type == "Mahalanobis") {
      for (j in 1:(dim(rX)[2]))
        rX[, j] <- rank(rX[, j])
      cv <- cov(rX)
      vuntied <- var(1:nobs)
      rat <- sqrt(vuntied / diag(cv))
      if (length(rat) == 1) {
        cv <- as.matrix(rat) %*% cv %*% as.matrix(rat)
      } else{
        cv <- diag(rat) %*% cv %*% diag(rat)
      }
      icov <- ginv(cv)
    }
    nums <- 1:nobs
    ctrl.nums <- 1:(sum(z == 0))
    treated <- nums[z == 1]
    
    #find distance between each treated and each control it will be connected to
    #and store in a distance structure
    dist.struct <- list()
    if (dist.type == "Euclidean") {
      distMat = as.matrix(dist(rX))
    }
    for (i in c(1:length(treated))) {
      controls <- nums[(z == 0) & (exact == exact[treated[i]])]
      control.names <- ctrl.nums[exact[z == 0] == exact[treated[i]]]
      if (dist.type == "Mahalanobis") {
        costi <-
          mahalanobis(rX[controls, , drop = FALSE], 
                      rX[treated[i],], icov, inverted = T)
      } else if (dist.type == "Euclidean") {
        costi <- as.vector(distMat[treated[i], controls])
      } else{
        costi <- as.vector(distMat[treated[i], controls])
      }
      
      if (calip.option != 'none') {
        calip.update <- rep(0, length(costi))
        calip.update[abs(calip.cov[treated[i]] - 
                           calip.cov[controls]) - cal > 0] <-
          Inf
        costi <- costi + calip.update
      }
      names(costi) <- control.names
      
      dist.struct[[i]] <- costi[is.finite(costi)]
    }
    
    if (sum(laply(dist.struct, length)) == 0)
      stop('All matches forbidden. Considering using a wider caliper?')
    return(dist.struct)
  }



#' An internal helper function that generates the data abstraction for the edge
#' weights of the main network structure using the distance matrix passed by the
#' user.
#' @param z a vector indicating whether each unit is in treatment or control
#'   group
#' @param distMat a matrix of pair-wise distance
#' @param verbose a boolean value whether to print(cat) debug information.
#'   Default: FALSE
#'
#' @return a distance structure used for constructing the main network flow
#'   problem
build.dist.struct_user <-
  function(z, distMat, verbose = FALSE) {
    if (!(all((z == 1) | (z == 0)))) {
      stop("The z argument must contain only 1s and 0s")
    }
    
    nobs <- length(z)
    
    nums <- 1:nobs
    ctrl.nums <- 1:(sum(z == 0))
    treated <- nums[z == 1]
    
    #find distance between each treated and each control it will be connected to
    #and store in a distance structure
    dist.struct <- list()
    for (i in c(1:length(treated))) {
      controls <- nums[(z == 0)]
      control.names <- ctrl.nums
      costi <- as.vector(distMat[treated[i], controls])
      names(costi) <- control.names
      dist.struct[[i]] <- costi[is.finite(costi)]
    }
    
    if (sum(laply(dist.struct, length)) == 0)
      stop('All matches forbidden. Considering using a wider caliper?')
    return(dist.struct)
  }


#' An internal helper function that combines two distance object
#'
#' @param a a distance structure object
#' @param b a distance structure object
#'
#' @return a new distance structure object whose edge weights are the sum of the
#'   corresponding edge weigths in a and b
combine_dist <- function(a, b) {
  if (length(a) != length(b))
    stop("Make sure both distance structure 
         have the same amount of treated units.")
  res = list()
  n = length(a)
  for (i in 1:n) {
    ai = a[[i]]
    bi = b[[i]]
    intersectedIndex = intersect(names(ai), names(bi))
    res[[i]] = ai[intersectedIndex] + bi[intersectedIndex]
    if ((length(intersectedIndex) != length(ai)) ||
        (length(intersectedIndex) != length(bi))) {
      stop(
        "The two distance matrices and parameters 
        do not provide the same set of edges in the network. Check your input."
      )
    }
  }
  return(res)
}
