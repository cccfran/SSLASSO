SSLASSO <- function(X, 
                y, 
                penalty = c("separable", "adaptive"), 
                variance = c("fixed", "unknown"),
                lambda1, 
                lambda0,
                theta=0.5, 
                sigma = 1,
                a = 1, b, 
                nu = 0, xi= 0, 
                eps = 0.001, 
                max.iter=500, 
                counter=10, 
                warn = FALSE) {
  # Coersion
  penalty <- match.arg(penalty)
  variance <- match.arg(variance)

  nlambda <- length(lambda0)
  if (class(X) != "matrix") {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
  }
  if (storage.mode(X)=="integer") storage.mode(X) <- "double"
  if (class(y) != "numeric") {
    tmp <- try(y <- as.numeric(y), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("y must numeric or able to be coerced to numeric")
  }

 
  # Lambda0 should be an increasing sequence

  monotone<-sum((lambda0[-1]-lambda0[-nlambda])>0)
  if (monotone!=nlambda-1){stop("lambda0 must be a monotone increasing sequence")}

  # change lambda0 to numeric (it throws an error if integer)
  lambda0 = as.numeric(lambda0)

  # error checking
  if (lambda1 > min(lambda0) ) stop("lambda1 must be smaller than lambda0")
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to SSLASSO")

  ## Standardize
  std <- .Call("standardize", X)
  XX <- std[[1]]
  center <- std[[2]]
  scale <- std[[3]]
  nz <- which(scale > 1e-6)
  if (length(nz) != ncol(XX)) XX <- XX[ ,nz, drop=FALSE]   
  yy <- y - mean(y)
  intercept <- rep(mean(y), length(lambda0))
  
  p <- ncol(XX)
  n <- length(yy)
 
  if(missing(b)) {b = p}
  
  ## Fit
  res <- .Call("SSL_gaussian", XX, yy, penalty, variance, lambda0, eps, as.integer(max.iter), as.double(lambda1), 
                 as.double(theta), as.double(sigma),as.integer(counter), as.double(a), as.double(b), as.double(nu), as.double(xi))
  bb <- matrix(res[[1]], p, nlambda)
  iter <- res[[3]]
  thetas<-res[[4]]
  sigmas <- res[[5]]
  

  ## Warning
  if (warn & any(iter==max.iter)) {

      warning("Algorithm did not converge for the ABOVE MENTIONED values of lambda0")

      print(lambda0[iter==max.iter])
                               }
  
  if (is.nan(sigmas[length(lambda0)])) {
    warning("Algorithm did not converge: sigma is reported as NaN")
  }

  ## Unstandardize
  beta <- matrix(0, nrow=(ncol(X)), ncol=length(lambda0))
  bbb <- bb/scale[nz]
  beta[nz, ] <- bbb

  ## Names
  varnames <- if (is.null(colnames(X))) paste("V",1:ncol(X),sep="") else colnames(X)
  dimnames(beta) <- list(varnames, round(lambda0,digits=4))

  ## Select

  select<-apply(beta,2,function(x){as.numeric(x!=0)})

  ## Model

  model<-(1:p)[select[,nlambda]==1]

  ## Output
  val <- structure(list(beta = beta,
                        iter = iter,
                        lambda0 = lambda0,
                        penalty = penalty,
                        lambda1 = lambda1,
                        thetas=thetas,
                        sigmas = sigmas,
                        select=select,
                        model=model,
                        n = n),
                   class = "SSLASSO")
 
  val
}
