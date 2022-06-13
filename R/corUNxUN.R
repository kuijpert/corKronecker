

corUNxUN <-
  ## Constructor for the corUNxUN class
  function(value = numeric(0), form = ~ 1, fixed = FALSE)
  {
    attr(value, "formula") <- form
    attr(value, "fixed") <- fixed
    class(value) <- c("corUNxUN", "corStruct")
    value
  }


getCovariate2 <-
  function(object, form = formula(object), data)
  {
    ## Return the primary covariate
    if (!(inherits(form, "formula"))) {
      stop("'form' must be a formula")
    }
    aux <- getCovariateFormula(form)
    if (length(all.vars(aux)) == 2) {
      out <- list(
        object[,all.vars(aux)[1]],
        object[,all.vars(aux)[2]])
      names(out) <- all.vars(aux)
      return(out)
    } else {
      ######### Should not occur: Stop #########
      stop()
      #rep(1, dim(object)[1])
    }
  }


getCovariate.corUNxUN <-
  function(object, form = formula(object), data)
  {
    if (!missing(form)) {
      form <- formula(object)
      warning("cannot change 'form'")
    }
    if (is.null(covar <- attr(object, "covariate"))) { # need to calculate it
      if (missing(data)) {
        stop("need data to calculate covariate of \"corStruct\" object")
      }
      covForm <- getCovariateFormula(form)
      grps <- if(!is.null(getGroupsFormula(form)))
        getGroups(object, data = data) ## else NULL
      if (length(all.vars(covForm)) == 2) { # covariates present
        if (is.null(grps)) {
          covar <- getCovariate2(data, covForm)
        } else {
          if (all(all.vars(covForm) == sapply(splitFormula(covForm, "+"),
                                              function(el) deparse(el[[2]])))) {
            covar1 <- split(getCovariate2(data, covForm)[[1]], grps)
            covar2 <- split(getCovariate2(data, covForm)[[2]], grps)
            names(covar1) <- levels(grps)
            names(covar2) <- levels(grps)
            covar <- list(covar1, covar2)
            names(covar) <- all.vars(covForm)

          } else {
            ######### Should not occur: stop() #########
            stop()
            #covar <- lapply(split(data, grps), getCovariate2, covForm)
          }
        }
      } else {
        ######### Should not occur: stop() #########
        stop()
        if (is.null(grps)) {
          #covar <- 1:nrow(data)
        } else {
          ######### Should not occur: stop() #########
          stop()
          #covar <- lapply(split(grps, grps), function(x) seq_along(x))
        }
      }
      if (!is.null(grps)) {
        covar <- as.list(covar)
      }
    }
    covar
  }



Initialize.corUNxUN <- function(object, data, ...)
{
  if (!is.null(attr(object, "maxCov"))) {# initialized - nothing to do
    return(object)
  }

  form <- formula(object)
  ## obtaining the groups information, if any
  if (!is.null(getGroupsFormula(form))) {
    attr(object, "groups") <- getGroups(object, form, data = data)
    attr(object, "Dim") <- Dim(object, attr(object, "groups"))
  } else {                              # no groups
    attr(object, "Dim") <- Dim(object, as.factor(rep(1, nrow(data))))
  }
  ## obtaining the covariate(s)
  attr(object, "covariate") <- getCovariate(object, data = data)

  covar <- attr(object, "covariate")
  if(!is.list(covar)) covar <- list(covar)
  if(any(c(is.numeric(unlist(covar[[1]])),is.numeric(unlist(covar[[2]])))==FALSE))
    stop("covariates must be integers for \"corUNxUN\" objects")
  if(any(unlist(mapply(function(c1,c2){
    duplicated(matrix(c(unlist(c1),unlist(c2)),ncol=2)) }
    ,covar[[1]],covar[[2]])))) {
    stop("covariates must have unique values within groups for \"corUNxUN\" objects")
  }

  covar1 <- unlist(covar[[1]]) - 1
  covar2 <- unlist(covar[[2]]) - 1
  uCov1 <- unique(covar1)
  uCov2 <- unique(covar2)
  maxCov <- c(max(uCov1) + 1, max(uCov2) + 1)

  if((length(uCov1) != maxCov[1]) | (length(uCov2) != maxCov[2])) {
    stop("unique values of the covariate  for \"corUNxUN\" objects must be a sequence of consecutive integers")
  }

  if (Dim(object)[["M"]] > 1) {
    covar <- list(
      split(covar1, getGroups(object)),
      split(covar2, getGroups(object))
    )
    names(covar) <- all.vars(getCovariateFormula(form))
    attr(object, "covariate") <- covar
  } else {
    covar <- list(covar1,covar2)
    names(covar) <- all.vars(getCovariateFormula(form))
    attr(object, "covariate") <- covar
  }
  attr(object, "maxCov") <- maxCov


  natPar <- as.vector(object)
  if (length(natPar) > 0) {
    ## parameters assumed in constrained form
    if (length(natPar) != sum(round(maxCov * (maxCov - 1) / 2))) {
      stop("initial value for \"corUNxUN\" parameters of wrong dimension")
    }
    if (max(abs(natPar)) >= 1) {
      stop("initial values for \"corUNxUN\" must be between -1 and 1")
    }

    npar <- round(maxCov * (maxCov - 1) / 2)

    natMat1 <- diag(maxCov[1])/2
    natMat1[lower.tri(natMat1)] <- natPar[1:npar[1]]
    natMat1 <- t(natMat1) + natMat1

    natMat2 <- diag(maxCov[2])/2
    natMat2[lower.tri(natMat2)] <- natPar[(npar[1]+1):(npar[1]+npar[2])]
    natMat2 <- t(natMat2) + natMat2

    natMat <- natMat2 %x% natMat1


    ## checking if positive-definite
    if (any(eigen(natMat, symmetric=TRUE, only.values=TRUE)$values <= 0)) {
      stop("initial values for \"corUNxUN\" do not define a positive-definite correlation structure")
    }
    natMat1 <- chol(natMat1)
    uncPar1 <- numeric(0)
    for(i in 2:maxCov[1]) {
      aux <- acos(natMat1[1:(i-1),i]/sqrt(cumsum(natMat1[i:1,i]^2)[i:2]))
      uncPar1 <- c(uncPar1, log(aux/(pi - aux)))
    }
    natMat2 <- chol(natMat2)
    uncPar2 <- numeric(0)
    for(i in 2:maxCov[2]) {
      aux <- acos(natMat2[1:(i-1),i]/sqrt(cumsum(natMat2[i:1,i]^2)[i:2]))
      uncPar2 <- c(uncPar2, log(aux/(pi - aux)))
    }
    uncPar <- c(uncPar1,uncPar2)
    coef(object) <- uncPar

  } else {				# initializing the parameters
    oldAttr <- attributes(object)
    object <- double(sum(round(maxCov * (maxCov - 1) / 2)))
    attributes(object) <- oldAttr
    attr(object, "factor") <- corFactor(object)
    attr(object, "logDet") <- -attr(attr(object, "factor"), "logDet")
  }
  object
}



corMatrix.corUNxUN <-
  function(object, covariate = getCovariate(object), corr = TRUE, ...)
  {
    cov1 <- covariate[[1]]
    cov2 <- covariate[[2]]

    corD <- Dim(object,
                if(is.list(cov1)) {
                  if (is.null(names(cov1)))
                    names(cov1) <- seq_along(cov1)
                  rep(names(cov1), lengths(cov1))
                } else
                  rep(1, length(cov1)))
    if (corr) {
      val <- matList.corUNxUN(
        as.vector(object),
        covariate,
        attr(object, "maxCov"),
        corD,
        mat = double(corD[["sumLenSq"]])
      )
      lD <- NULL
    } else {
      val <- corFactor(object)
      lD <- val[["logDet"]]
      val <- val[["factor"]]
    }
    if (corD[["M"]] > 1) {
      val <- split(val, rep(1:corD[["M"]], (corD[["len"]])^2))
      val <- lapply(val, function(el) {
        nel <- round(sqrt(length(el)))
        array(el, c(nel, nel))
      })
      names(val) <- names(corD[["len"]])
      val <- as.list(val)
    } else {
      val <- array(val, c(corD[["N"]], corD[["N"]]))
    }
    attr(val, "logDet") <- lD
    val
  }



getcor <- function(par,maxC) {

  n <- maxC
  work <- double(n*(n+1)/2)
  dest <- 1
  src <- 1

  par <- pi * exp(par) / (1 + exp(par))

  for(i in 1:n){
    aux <- 1
    if(i>1){
      for(j in 1:(i-1)){
        aux1 <- par[src]
        work[dest] <- aux * cos(aux1)
        aux <- aux * sin(aux1)
        dest <- dest + 1
        src <- src + 1
      }
    }
    work[dest] <- aux
    dest <- dest + 1
  }

  out <- matrix(0,nrow=n,ncol=n)
  out[upper.tri(out,diag=TRUE)] <- work

  t(out) %*% out

}



matList.corUNxUN <- function(par,covariate,maxCov,corD,mat) {

  M <- corD$M
  len <- corD$len

  stopsq <- as.integer(cumsum(len*len))
  startsq <- as.integer(stopsq - len*len +1)

  npar <- round(maxCov * (maxCov - 1) / 2)

  cormat1 <- getcor(par[1:npar[1]],maxCov[1])
  cormat2 <- getcor(par[(npar[1]+1):(npar[1]+npar[2])],maxCov[2])
  cormat <- cormat2 %x% cormat1

  d1 <- rep(0:(maxCov[1]-1),maxCov[2])
  d2 <- rep(0:(maxCov[2]-1),each=maxCov[1])

  for(i in 1:M) {

    c1 <- covariate[[1]][[i]]
    c2 <- covariate[[2]][[i]]

    index <- mapply(function(x1,x2) {  which(x1==d1 & x2==d2) },c1,c2)

    i.mat <- cormat[index,index]

    mat[startsq[i]:stopsq[i]] <- i.mat

  }

  return(mat)

}


coef.corUNxUN <- function(object, unconstrained = TRUE, ...)
{
  if (unconstrained) {
    if (attr(object, "fixed")) {
      return(numeric(0))
    } else {
      return(as.vector(object))
    }
  }

  par <- as.vector(object)
  maxCov <- attr(object,"maxCov")
  npar <- round(maxCov * (maxCov - 1) / 2)

  cormat1 <- getcor(par[1:npar[1]],maxCov[1])
  cormat2 <- getcor(par[(npar[1]+1):(npar[1]+npar[2])],maxCov[2])

  coords1 <- matrix(sapply(paste0("(",1:dim(cormat1)[2],","),
               function(el){ paste0(el,1:dim(cormat1)[2],")") }),
               nrow=dim(cormat1)[1],ncol=dim(cormat1)[2],byrow=TRUE)

  coords2 <- matrix(sapply(paste0("(",1:dim(cormat2)[2],","),
               function(el){ paste0(el,1:dim(cormat2)[2],")") }),
               nrow=dim(cormat2)[1],ncol=dim(cormat2)[2],byrow=TRUE)

  rho1 <- t(cormat1)[lower.tri(t(cormat1))]
  names(rho1) <- t(coords1)[lower.tri(t(coords1))]

  rho2 <- t(cormat2)[lower.tri(t(cormat2))]
  names(rho2) <- t(coords2)[lower.tri(t(coords2))]

  names(rho1) <- paste0("rho.1",names(rho1))
  names(rho2) <- paste0("rho.2",names(rho2))

  out <- c(rho1,rho2)
  out

}

