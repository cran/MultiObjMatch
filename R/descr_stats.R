#' Generate summary statistics for matches
#'
#' @param matches One matching result from the main matching function
#' @param df the original data frame used for matching
#' @param treatCol the character of the column name for treatment vector
#' @param b.vars the vector of column names of covariates used for measuring
#'   balance
#' @param pair.vars the vector of column names used for measuring pairwise
#'   distance
#' @param extra the list of summary statistic; it must be the types that can be
#'   taken by cobalt
#'
#' @return a named vector of summary statistic
descr.stats_general <-
  function(matches,
           df,
           treatCol,
           b.vars,
           pair.vars,
           extra = FALSE) {
    treatment.status <- c(rep(1, nrow(matches)), rep(0,
                                                     nrow(matches)))
    
    #total variation imbalance
    matched.info <-
      rbind(df[df[, treatCol] == 1, ][as.numeric(rownames(matches)),
                                      , drop = FALSE], 
            df[df[, treatCol] == 0, ][as.vector(matches), , drop = FALSE])
    interact.factors.matched = as.factor(apply(matched.info[,
                          match(b.vars, colnames(matched.info)), drop = FALSE],
                                               1, function(x)
                                                 paste(x, collapse = ".")))
    fb.tab <-
      table('balance.variable' = interact.factors.matched, treatment.status)
    tv.sum <- sum(abs(fb.tab[, 1] - fb.tab[, 2]))
    tv.prop <-  sum(abs(fb.tab[, 1] - fb.tab[, 2])) / sum(fb.tab)
    names(tv.sum) <- paste(c('TV', b.vars), collapse = '.')
    names(tv.prop) <- paste(c('TV', b.vars), collapse = '.')
    if (is.null(pair.vars)) {
      return(c('TV' = tv.sum))
    }
    #chi-squared test
    
    suppressWarnings(chisq.p <-  chisq.test(fb.tab)$p.value)
    names(chisq.p) <- paste(c('Chisq', b.vars), collapse = '.')
    
    interact.pv = as.character(apply(matched.info[,
                      match(pair.vars, colnames(matched.info)), drop = FALSE],
                                     1, function(x)
                                       paste(x, collapse = ".")))
    interact.pv[is.na(interact.pv)] <- 'NA.explicit'
    exp.treated <- interact.pv[matched.info[, treatCol] == 1]
    exp.control <- interact.pv[matched.info[, treatCol] == 0]
    pair.interact.exact <-
      mean(as.character(exp.treated) == as.character(exp.control))
    names(pair.interact.exact) <-
      paste(c('Exact', pair.vars), collapse = '.')
    
    
    #do standardized diffs
    #need to generate a vector of names of non-factor variables.
    sdiff.vars <- unique(c(b.vars, pair.vars))
    factor.vars <-
      sdiff.vars[aaply(sdiff.vars, 1, function(x)
        is.factor(df[[x]]) || is.character(df[[x]]))]
    
    #augment my.slice with dummified versions of factors
    temp.df <- df
    add.vars <- c()
    for (varn in factor.vars) {
      varn.dummy <- dummy(df[[varn]], name = varn)
      for (coln in colnames(varn.dummy)) {
        temp.df[[coln]] <- varn.dummy[, which(colnames(varn.dummy) == coln)]
      }
      add.vars <- c(add.vars, colnames(varn.dummy))
    }
    sdiff.vars <- setdiff(sdiff.vars, factor.vars)
    ## SMALL Change in line 134
    ## sdiff.vars <- unique(c(sdiff.vars, add.vars, 'exp', 'comorb','age'))
    
    z <- temp.df[, treatCol]
    marg.bal <- aaply(sdiff.vars, 1, function(x) {
      sdiff(x, treatCol, orig.data = temp.df, 
            match.data = temp.df[c(which(z ==1)[as.numeric(rownames(matches))], 
                                   which(z == 0)[matches]), ])[6]
    })
    names(marg.bal) <- paste('SD', sdiff.vars, sep = '.')
    
    
    if (!extra)
      return(c(tv.prop,  chisq.p,  pair.interact.exact, marg.bal))
    
    #row-wise tests
    pval.vector <- rep(NA, nrow(fb.tab))
    for (i in 1:nrow(fb.tab)) {
      # if(nrow(fb.tab)==1){
      #   pval.vector[i] <- 1
      # }else{
      #   new.tab <- rbind(fb.tab[i, ], colSums(fb.tab[-i, ]))
      #   suppressWarnings(pval.vector[i] <- chisq.test(new.tab)$p.value)
      # }
      suppressWarnings(pval.vector[i] <-1)
    }
    names(pval.vector) <- paste('Pval', rownames(fb.tab), sep = '.')
    
    #calculate rate of exact matching on pair.vars
    exact.prop <- rep(NA, length(pair.vars))
    for (i in 1:length(exact.prop)) {
      col.idx <- match(pair.vars[i], colnames(matched.info))
      exp.treated <-
        as.character(df[, col.idx]
                     [df[, treatCol] == 1][as.numeric(rownames(matches))])
      exp.control <-
        as.character(df[, col.idx][df[, treatCol] == 0][as.vector(matches)])
      exp.treated[is.na(exp.treated)] <- 'NA.explicit'
      exp.control[is.na(exp.control)] <- 'NA.explicit'
      exact.prop[i] <-
        mean(as.character(exp.treated) == as.character(exp.control))
    }
    names(exact.prop) = paste('Exact', pair.vars, sep = '.')
    
    
    if (length(pair.vars) == 1)
      return(c(tv.prop,  chisq.p,  exact.prop, median(pval.vector)), marg.bal)
    
    return(
      c(
        'TV' = tv.sum,
        'Chi-squared pval' = chisq.p,
        pair.interact.exact,
        exact.prop,
        median(pval.vector),
        marg.bal
      )
    )
  }
