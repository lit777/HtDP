## check whether packages are installed
pacotes = c("invgamma", "dplyr", "MASS", "LaplacesDemon",
            "armspp", "CholWishart", "mvtnorm", "mnormt", "progress",
            "Rcpp", "mniw", "plyr")

package.check <- lapply(pacotes, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
  }
})


## set working directory
CURRENT_WORKING_DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
CURRENT_WORKING_DIR
setwd(CURRENT_WORKING_DIR)


## load CDP and HtDP functions
sourceCpp('HtDP_update.cpp')
sourceCpp('CDP_update.cpp')


## set functions
cal_by_cluster <- function(dat, indic, clus) {
  len <- length(clus)
  result <- list()
  
  # divide data set by cluster
  for (i in 1:len) {
    result[[i]] <- as.matrix(dat[which(indic==clus[i]),])
  }
  names(result) <- clus
  return(result)
  
}
dat_by_cluster <- function(dat, indic, clus) {
  
  result <- list()
  len <- length(clus)
  
  # divide data set by cluster
  for (j in 1:len) {
    result[[j]] <- as.matrix(dat[which(indic==clus[j]),])
    if (ncol(result[[j]])==1) {
      result[[j]] <- t(result[[j]])
    }
  }
  names(result) <- clus
  return(result)
}
new.dmvnorm <- function (x, mean, sigma) {
  if (is.vector(x)) 
    x <- matrix(x, ncol = length(x))
  p <- ncol(x)
  
  
  if (!missing(mean)) {
    if (!is.null(dim(mean))) 
      dim(mean) <- NULL
    if (length(mean[[1]]) != p) 
      stop("x and mean have non-conforming size")
  }
  
  dec <- Map('chol', sigma)
  centered.x <- Map('-', rep(list(t(x)), length(mean)), mean)
  tmp <- Map('backsolve', dec, centered.x, transpose=TRUE)
  tmp <- Map('*', tmp, tmp)
  rss <- Map('*', Map('colSums', tmp), 0.5)
  
  logretval <- Map('-',Map('+', Map('+', Map('sum', Map('log',Map('diag', dec))), 
                                    rep(list(0.5*p*log(2*pi)), length(mean))), 
                           rss))
  logretval <- Map('exp', logretval)
}
mat_split <- function(M, r, c){
  nr <- ceiling(nrow(M)/r)
  nc <- ceiling(ncol(M)/c)
  newM <- matrix(NA, nr*r, nc*c)
  newM[1:nrow(M), 1:ncol(M)] <- M
  
  div_k <- kronecker(matrix(seq_len(nr*nc), nr, byrow = TRUE), matrix(1, r, c))
  matlist <- split(newM, div_k)
  N <- length(matlist)
  mats <- unlist(matlist)
  dim(mats)<-c(r, c, N)
  mats <- lapply(seq(dim(mats)[3]), function(x) mats[ , , x])
  return(mats)
}
CDP <- function(dat, iter) {
  p <- ncol(dat)
  n <- nrow(dat)
  mu.y <- colMeans(dat)
  names(mu.y) <- NULL
  cov.y <- cov(dat)
  dimnames(cov.y) <- NULL
  inv.cov.y <- solve(cov.y)
  
  indicator <- list()
  xi <- list()
  rho <- NULL
  W <- list()
  beta <- NULL
  alpha <- NULL
  muj <- list()
  Sj <- list()
  cluster <- list()
  gibb.Sj <- list() ; gibb.muj <- list()
  
  # set initial hyperparameters
  W[[1]] <- mniw::rwish(1, cov.y/p, p)
  rho[1] <- rgamma(1, shape=1/4, scale = 2)
  xi[[1]] <- mvrnorm(1, mu=mu.y, Sigma=cov.y)
  alpha[1] <- LaplacesDemon::rinvgamma(1, shape = 1/2, scale = 1/2)
  beta[1] <- LaplacesDemon::rinvgamma(1, shape = 1/2, scale = p/2) + p - 1
  indicator[[1]] <- seq(1, nrow(dat))
  
  # set initial muj
  muj <- list() ; Sj <- list()
  for (i in 1:length(indicator[[1]])) {      
    S <- mniw::rwish(1, solve(beta[1]*W[[1]]), beta[1])
    muj[[i]] <- mvrnorm(1, xi[[1]], solve(rho[1]*S))
  }
  
  
  pb <- progress_bar$new(total = iter)
  for (l in 1:iter) {  
    
    
    # generate unique cluster
    cluster[[l]] <- sort(unique(indicator[[l]]))
    
    dat.by.clus <- cal_by_cluster(dat, indicator[[l]], cluster[[l]]) 
    
    nj <- unlist(Map('nrow', dat.by.clus))
    k <- length(nj)
    
    # update muj and Sj
    for (i in 1:k) {
      
      dat.quad <- dat.by.clus[[i]] - matrix(muj[[i]], byrow=TRUE, ncol=p, nrow=nj[i])
      
      Sj[[i]] <- mniw::rwish(1, 
                             solve(crossprod(dat.quad) + 
                                     rho[l]*tcrossprod(muj[[i]]-xi[[l]]) + 
                                     beta[l]*W[[l]]), 
                             nj[i]+beta[l]+1)
      
      muj[[i]] <- mvrnorm(1, 
                          (colSums(dat.by.clus[[i]]) + rho[l]*xi[[l]])/(nj[i]+rho[l]), 
                          solve(round((nj[i]+rho[l])*Sj[[i]],10)))
      
    }
    names(muj) <- names(dat.by.clus)
    names(Sj) <- names(dat.by.clus)
    
    gibb.Sj[[l]] <- Sj
    gibb.muj[[l]] <- muj
    
    # update rho
    centered_mu <- Map('-', muj, rep(list(xi[[l]]), k))
    rho_scale <- Map('t', centered_mu)
    rho_scale <- Map('%*%', rho_scale, Sj)
    rho_scale <- Map('%*%', rho_scale, centered_mu)
    rho_scale <- Reduce('+', rho_scale)
    rho[l+1] <- rgamma(1, shape=(k*p + 1/2)/2, scale=2/(rho_scale+1)) 
    
    # update xi
    sum.Sj <- Reduce('+', Sj)
    sum.Sjmuj <- Reduce('+', Map('%*%', Sj, muj))
    xi.sigma <- solve(rho[l+1]*sum.Sj + inv.cov.y)
    xi[[l+1]] <- mvrnorm(1, 
                         xi.sigma%*%(rho[l+1]*sum.Sjmuj + inv.cov.y%*%mu.y), 
                         xi.sigma)
    
    # update W
    W[[l+1]] <- mniw::rwish(1, solve(beta[l]*sum.Sj + p*inv.cov.y), beta[l]*k + p)
    
    # update beta
    sumLogDet_Sj <- Reduce('+', Map('log', Map('det', Sj)))
    tr_Sj_W <- tr(sum.Sj%*%W[[l+1]])
    
    beta[l+1] <- 1/arms(1, 
                        function(x) ((1/x+p-1)*p*k)/2*log((1/x+p-1)/2) + 
                          (1/x+p-1)*k/2*log(det(W[[l+1]])) - 
                          k*lmvgamma((1/x+p-1)/2,p) + (1/x-2)/2*sumLogDet_Sj - 
                          (1/x+p-1)/2*tr_Sj_W - 1/2*log(x) - p/2*x,
                        lower=0, upper = 1000, metropolis=TRUE, previous=1/beta[l]+p-1,
                        include_n_evaluations = FALSE) + p - 1
    
    # update indicator
    new.indicator <- indicator[[l]]
    
    part <- Map('tcrossprod',Map('-', mat_split(as.matrix(dat), 1, p), rep(list(xi[[l+1]]), n)))
    part <- Map('*', part, rho[l+1]/(rho[l+1]+1))
    part <- Map('+', part, rep(list(beta[l+1]*W[[l+1]]), n))
    det.part <- Map('det', part)
    log.det.part <- log(unlist(det.part))
    
    log_inact_prop <- c(p/2*log(rho[l+1]/(rho[l+1]+1)) - p/2*log(pi) + 
                          lgamma((beta[l+1]+1)/2) - lgamma((beta[l+1]+1-p)/2) + 
                          beta[l+1]/2*log(det(beta[l+1]*W[[l+1]]))) - (beta[l+1]+1)/2*log.det.part
    inact_prop <- exp(log_inact_prop)*alpha[l]
    
    output <- cdp_update(as.matrix(dat), new.indicator, muj, Sj, beta[[l+1]], W[[l+1]], 
                         rho[l+1], xi[[l+1]], inact_prop)
    
    muj <- output[[2]]
    Sj <- output[[3]]
    muj <- Map('as.vector', muj)
    
    indicator[[l+1]] <- as.numeric(factor(rank(output[[1]])))
    names(muj) <- sort(unique(indicator[[l+1]]))
    names(Sj) <- sort(unique(indicator[[l+1]]))
    
    
    k <- length(unique(indicator[[l+1]]))
    
    # update alpha
    alpha[l+1] <- 1/arms(1, 
                         function(x) (-k-1/2)*log(x) + lgamma(1/x) - lgamma(n+1/x) - x/2,
                         lower=0, upper = 1000, 
                         metropolis=TRUE, previous=1/alpha[l],
                         include_n_evaluations = FALSE)
    
    pb$tick()
    
  }
  
  result <- list("cluster" = cluster,
                 "indicator" = indicator,
                 "mu" = gibb.muj,
                 "S" = gibb.Sj,
                 "xi" = xi,
                 "rho" = rho,
                 "W" = W,
                 "alpha" = alpha,
                 "beta" = beta)
  
  return(result)
  
}
HtDP <- function(dat, iter) {
  p <- ncol(dat)
  n <- nrow(dat)
  mu.y <- colMeans(dat)
  names(mu.y) <- NULL
  cov.y <- cov(dat)
  dimnames(cov.y) <- NULL
  inv.cov.y <- solve(cov.y)
  
  indicator <- list()
  xi <- list()
  rho <- NULL
  alpha <- c()
  muj <- list()
  Vj <- list()
  cluster <- list()
  gibb.muj <- list()
  gibb.Vj <- list()
  
  # set initial hyperparameters
  xi[[1]] <- mvrnorm(1, mu=mu.y, Sigma=cov.y)
  rho[1] <- rgamma(1, shape=1/4, scale = 2)
  A <- 100000
  nu <- 2
  a <- matrix(rep(0, p*(iter+1)), ncol=p)
  
  a[1, ] <- LaplacesDemon::rinvgamma(p, shape = 1/2, scale = 1/(A)^2)
  diag.a <- diag(1/a[1,])
  alpha[1] <- LaplacesDemon::rinvgamma(1, shape = 1/2, scale = 1/2)
  indicator[[1]] <- seq(1, nrow(dat))
  
  # set initial Vj and muj
  Vj <- alply(mniw::riwish(length(indicator[[1]]), 
                           2*nu*diag.a, 
                           nu+p-1), 
              3)
  
  for (i in 1:length(indicator[[1]])) {      
    V <- riwish(1, 2*nu*diag.a, nu+p-1)
    muj[[i]] <- mvrnorm(1, xi[[1]], V/rho[1])
  }
  
  pb <- progress_bar$new(total = iter)
  
  for (l in 1:iter) {   
    
    # generate unique cluster
    cluster[[l]] <- sort(unique(indicator[[l]]))
    
    dat.by.clus <- cal_by_cluster(dat, indicator[[l]], cluster[[l]]) 
    
    nj <- unlist(Map('nrow', dat.by.clus))  
    k <- length(nj)  
    
    diag.a <- diag(1/a[l,])
    
    # update muj's and Vj's
    for (i in 1:k) {    
      
      dat.quad <- dat.by.clus[[i]] - matrix(muj[[i]], byrow=TRUE, ncol=p, nrow=nj[i])
      
      Vj[[i]] <- riwish(1, 
                        crossprod(dat.quad) + rho[l]*tcrossprod(muj[[i]]-xi[[l]]) + 2*nu*diag.a, 
                        nj[i]+nu+p)
      
      muj[[i]] <- mvrnorm(1, 
                          (colSums(dat.by.clus[[i]]) + xi[[l]]*rho[l]) / (nj[i]+rho[l]), 
                          Vj[[i]] / (nj[i]+rho[l]))
      
    }
    
    names(muj) <- names(dat.by.clus) 
    names(Vj) <- names(dat.by.clus)
    
    gibb.Vj[[l]] <- Vj ; gibb.muj[[l]] <- muj
    
    # update xi
    inv.Vj <- Map('solve', Vj)
    sum.inv.Vj <- Reduce('+', inv.Vj)
    inv.Vj.muj <- Map('%*%', inv.Vj, muj)
    sum.inv.Vj.muj <- Reduce('+', inv.Vj.muj)
    
    xi.sigma <- solve(sum.inv.Vj*rho[l]+inv.cov.y) 
    xi[[l+1]] <- mvrnorm(1, 
                         xi.sigma%*%(sum.inv.Vj.muj*rho[l] + inv.cov.y%*%mu.y), 
                         xi.sigma)
    
    # update rho
    centered_mu <- Map('-', muj, rep(list(xi[[l+1]]), k))
    rho_scale <- Map('t', centered_mu)
    rho_scale <- Map('%*%', rho_scale, inv.Vj)
    rho_scale <- Map('%*%', rho_scale, centered_mu)
    rho_scale <- Reduce('+', rho_scale)
    rho[l+1] <- rgamma(1, shape = p*k/2 + 1/4, scale = 1/(rho_scale/2 + 1/2))
    
    
    # update vector a
    a_scale <- diag(sum.inv.Vj)
    
    a[l+1,] <- LaplacesDemon::rinvgamma(p, 
                                        shape = (1+k*(nu+p-1))/2, 
                                        scale = nu*a_scale + 1/A^2)
    diag.a <- diag(1/a[l+1,])
    
    
    # update indicator and cluster
    new.indicator <- indicator[[l]]
    
    part <- Map('tcrossprod',Map('-', mat_split(as.matrix(dat), 1, 2), rep(list(xi[[l+1]]), n)))
    part <- Map('*', part, rho[l+1]/(rho[l+1]+1))
    part <- Map('+', part, rep(list(2*nu*diag.a), n))
    det.part <- Map('det', part)
    log.det.part <- log(unlist(det.part))
    
    inact_prop <- c(p/2*log(rho[l+1]/(rho[l+1]+1)) - p/2*log(pi) + 
                      lgamma((nu+p)/2) - lgamma((nu+p-p)/2) + 
                      (nu+p-1)/2*log(det(2*nu*diag.a))) - (nu+p)/2*log.det.part
    inact_prop <- exp(inact_prop)*alpha[l]
    
    output <- htdp_update(as.matrix(dat), new.indicator, muj, Vj, nu, 
                          diag.a, rho[l+1], xi[[l+1]], inact_prop)
    
    muj <- output[[2]]
    Vj <- output[[3]]
    muj <- Map('as.vector', muj)
    
    indicator[[l+1]] <- as.numeric(factor(rank(output[[1]])))
    names(muj) <- sort(unique(indicator[[l+1]]))
    names(Vj) <- sort(unique(indicator[[l+1]]))
    
    
    # update alpha
    k <- length(unique(indicator[[l+1]]))
    alpha[l+1] <- 1/arms(1, 
                         function(x) (-k-1/2)*log(x) + lgamma(1/x) - lgamma(n+1/x) - x/2,
                         lower=0, upper = 100000, 
                         metropolis=TRUE, previous=1/alpha[l],
                         include_n_evaluations = FALSE)
    
    pb$tick()
    
  }
  
  result <- list("cluster" = cluster,
                 "indicator" = indicator,
                 "mu" = gibb.muj,
                 "V" = gibb.Vj,
                 "xi" = xi,
                 "rho" = rho,
                 "a" = a,
                 "alpha" = alpha)
  
  return(result)
  
}


## save functions
save(cal_by_cluster,
     dat_by_cluster,
     new.dmvnorm,
     mat_split,
     CDP,
     HtDP, file="code_functions.RData")


