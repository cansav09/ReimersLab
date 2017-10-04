function (object, test = c("Pillai", "Wilks", "Hotelling-Lawley", 
                           "Roy"), intercept = FALSE, tol = 1e-07, ...) 
{
  if (!inherits(object, "maov")) 
    stop(gettextf("object must be of class %s or %s", dQuote("manova"), 
                  dQuote("maov")), domain = NA)
  test <- match.arg(test) ##### determine which test to do
  asgn <- object$assign[object$qr$pivot[1L:object$rank]]### Rank is the number of columns? This gives the assignments of either intercepts or a variable? 
    #Pivot is the leading coefficient when a matrix is put into Row Echelon form
  uasgn <- unique(asgn) ## Narrow it down to the options for assignments
  nterms <- length(uasgn)## How many options are there? 
  effects <- object$effects ## Effects from each variable
  if (!is.null(effects)) 
    effects <- as.matrix(effects)[seq_along(asgn), , drop = FALSE]#### Just take out the effects from each group 
  rdf <- object$df.residual### Residual degrees of freedom
  nmeffect <- c("(Intercept)", attr(object$terms, "term.labels"))## labels for things
  resid <- as.matrix(object$residuals) ### matrix of residuals
  wt <- object$weights
  if (!is.null(wt)) ### if it is weighted, this next part square roots the weights
    resid <- resid * wt^0.5
  nresp <- NCOL(resid)
  if (nresp <= 1) ### If there are not more than 1 column of DV's then it will stop here. 
    stop("need multiple responses")
  if (is.null(effects)) { ### If there are no effects in the object, it sets everything to 0
    df <- nterms <- 0
    ss <- list(0)
    nmrows <- character()
  }
  else {
    df <- numeric(nterms) ### sets degrees of freedom
    ss <- list(nterms) ## creates a list with nterms
    nmrows <- character(nterms) ## creates a character object with nterms
    for (i in seq(nterms)) { ### for the number of DV's in the set, this will be repeated
      ai <- (asgn == uasgn[i]) #### Matches the assignment 0, 1 
      nmrows[i] <- nmeffect[1 + uasgn[i]]## puts the it in a respective nmrow
      df[i] <- sum(ai) ## degrees of ffreedom will be changed to the number of groups -1 
      ss[[i]] <- crossprod(effects[ai, , drop = FALSE])
    }
  }
  pm <- pmatch("(Intercept)", nmrows, 0L)#### If there's a row called Intercept which row is it
  if (!intercept && pm > 0) {
    nterms <- nterms - 1
    df <- df[-pm] ### Remove the intercepts row data
    nmrows <- nmrows[-pm] ## Remove the names for it
    ss <- ss[-pm]
  }
  names(ss) <- nmrows
  nt <- nterms
  if (rdf > 0) { #### if residual degrees of freedom is not 0
    nt <- nterms + 1
    df[nt] <- rdf
    if(object$df.residual < ncol(resid)){ ### if rank of matrix is less than the number of columns
      ss[[nt]] <- crossprod(condreg(resid,object$df.residual-1)$S)
    }else {
        ss[[nt]] <- crossprod(resid)
      } ### crossproduct of the residuals is set as an item in the ss object
    names(ss)[nt] <- nmrows[nt] <- "Residuals"
    ok <- df[-nt] > 0 ### are the degrees of freedom greater than 0? 
    eigs <- array(NA, c(nterms, nresp), dimnames = list(nmrows[-nt], ### creates a correctly sized empty array for the eigenvalues 
                                                        NULL))
    stats <- matrix(NA, nt, 5, dimnames = list(nmrows, c(test, 
                                                         "approx F", "num Df", "den Df", "Pr(>F)"))) ## makes a matrix to put the stats in 
    sc <- sqrt(sss <- diag(ss[[nt]])) ### extracts the diagonal of the residuals matrix (the variances of each variable) sets them to be the square root
    for (i in seq_len(nterms)[ok]) sss <- sss + diag(ss[[i]])
    sc[sc < sqrt(sss) * 1e-06] <- 1
    D <- diag(1/sc) ### make a matrix with the diagonals extracted before but with the inverse of the sc 
    rss.qr <- qr(D %*% ss[[nt]] %*% D, tol = tol)### multiplies by the crossproduct of the residuals by the matrix D then uses qr to decompose the matrix 
    ## tolerance here is set to 1e-07

    if (!is.null(rss.qr)) 
      for (i in seq_len(nterms)[ok]) {
        A1 <- qr.coef(rss.qr, D %*% ss[[i]] %*% D) ## takes a qr decomposed matrix and gives the coeff that result when trying to fit to the given matrix
        
        eigs[i, ] <- Re(eigen(A1, symmetric = FALSE, 
                              only.values = TRUE)$values)
        stats[i, 1L:4L] <-(test, Pillai = Pillai(eigs[i, 
                                                             ], df[i], df[nt]), Wilks = Wilks(eigs[i, ], 
                                                                                              df[i], df[nt]), `Hotelling-Lawley` = HL(eigs[i, 
                                                                                                                                           ], df[i], df[nt]), Roy = Roy(eigs[i, ], df[i], 
                                                                                                                                                                        df[nt]))
        ok <- stats[, 2L] >= 0 & stats[, 3L] > 0 & stats[, 
                                                         4L] > 0
        ok <- !is.na(ok) & ok
        stats[ok, 5L] <- pf(stats[ok, 2L], stats[ok, 
                                                 3L], stats[ok, 4L], lower.tail = FALSE)
      }
    x <- list(row.names = nmrows, SS = ss, Eigenvalues = eigs, 
              stats = cbind(Df = df, stats = stats))
  }
  else x <- list(row.names = nmrows, SS = ss, Df = df)
  class(x) <- "summary.manova"
  x
}
<bytecode: 0x123cc9aa8>
  <environment: namespace:stats>