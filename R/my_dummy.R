#' This is a modified version of the function "dummy" from the R package
#' dummies. Original code Copyright (c) 2011 Decision
#' Patterns.
#'
#' Change is made to the "model.matrix" function so that the output could be
#' used for the current package.
#' 
#' @param x a data.frame, matrix or single variable or variable name
#' @param data (optional) if provided, x is the name of a column on the data
#' @param sep (optional) the separator used between variable name and the value 
#' @param drop (optional) whether to drop unused levels
#' @param fun (optional) function to coerce the value in the final matrix; 'as,integer' by default
#' @param verbose (optional) whether to print the number of variables; FALSE by default
#' @param name (optional) the column name to be selected for converting; NULL by default
dummy <-
  function (x,
            data = NULL,
            sep = "",
            drop = TRUE,
            fun = as.integer,
            verbose = FALSE,
            name = NULL)
  {
    if (is.null(data)) {
      if (is.null(name))
        name <- as.character(sys.call(1))[2]
      name <- sub("^(.*\\$)", "", name)
      name <- sub("\\[.*\\]$", "", name)
    }
    else {
      if (length(x) > 1)
        stop("More than one variable provided to produce dummy variable.")
      name <- x
      x <- data[, name]
    }
    if (drop == FALSE && inherits(x, "factor")) {
      x <- factor(x, levels = levels(x), exclude = NULL)
    }
    else {
      x <- factor(x, exclude = NULL)
    }
    if (length(levels(x)) < 2) {
      if (verbose)
        warning(name, " has only 1 level. Producing dummy variable anyway.")
      return(matrix(
        rep(1, length(x)),
        ncol = 1,
        dimnames = list(rownames(x),
                        c(paste(
                          name, sep, x[[1]], sep = ""
                        )))
      ))
    }
    mm <-
      model.matrix( ~ x - 1, model.frame( ~ x - 1))#, contrasts = FALSE)
    colnames.mm <- colnames(mm)
    if (verbose)
      cat(" ", name, ":", ncol(mm), "dummy varibles created\n")
    mm <- matrix(
      fun(mm),
      nrow = nrow(mm),
      ncol = ncol(mm),
      dimnames = list(NULL,
                      colnames.mm)
    )
    colnames(mm) <-
      sub("^x", paste(name, sep, sep = ""), colnames(mm))
    if (!is.null(row.names(data)))
      rownames(mm) <- rownames(data)
    return(mm)
  }
