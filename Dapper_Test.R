library(MASS)
library(matrixsampling)
#library(MCMCpack)

source("DP_Functions.R")
source("privacy_class.R")

#' Generate MCMC samples using the dapper model
dapper_sample <- function(data_model = NULL,
                          sdp        = NULL,
                          init_par   = NULL,
                          seed       = NULL,
                          niter      = 2000,
                          warmup     = floor(niter / 2),
                          chains     = 1) {

  checkmate::assert_class(data_model, "privacy")
  checkmate::assert_count(niter)
  checkmate::assert_count(warmup)
  checkmate::assert_count(chains)
  checkmate::assert_count(niter - warmup)
  if(!is.null(seed)) checkmate::assert_count(seed)
  
  if(!is.null(seed)) set.seed(seed)
  p <- progressr::progressor(niter * chains)
  
  
  out <- dapper_chain(data_model,
                           sdp,
                           init_par,
                           niter = niter,
                           warmup = warmup,
                           prg_bar = NULL)

  f_out <- list(out)
  theta_clist <- lapply(1:chains, function(s) f_out[[s]]$sample)
  
  accept_mat <- do.call(cbind, lapply(1:chains, function(s) f_out[[s]]$accept_rate))
  
  new_dpout(theta_clist, accept_mat, data_model$varnames)

}

#' Execute a single MCMC chain
dapper_chain <- function(data_model,
                         sdp,
                         init_par,
                         niter = niter,
                         warmup = floor(niter / 2),
                         prg_bar = NULL) {
  
  checkmate::assert_class(data_model, "privacy")
  checkmate::assert_count(niter)
  checkmate::assert_count(warmup)
  checkmate::assert_count(niter - warmup)
  
  post_f <- data_model$post_f
  latent_f <- data_model$latent_f
  priv_f <- data_model$priv_f
  st_f <- data_model$st_f
  npar <- data_model$npar
  
  d <- sdp$d
  u <- sdp$u
  accept_rate <- numeric(niter)
  theta_list <- list()
  theta <- init_par
  dmat <- latent_f(init_par, sdp$n)
  st <- st_f(dmat, sdp, 1)
  nobs <- nrow(dmat)
  
  for (i in 1:niter) {
    #print(i)
    setTxtProgressBar(pb, i)
    counter <- 0
    theta_list[[i]] <- theta 
    theta <- post_f(dmat, theta, rep(0,d), 1, diag(rep(1,d)), 1, FALSE, u)
    smat <- latent_f(theta, nobs)
    
    for (j in 1:nobs) {
      xs <- smat[j, ]
      xo <- dmat[j, ]
      sn <- NULL
      
      sn <- merge_sx(st, st_f(xo, sdp, j), add = FALSE)
      sn <- merge_sx(sn, st_f(xs, sdp, j), add = TRUE)
      
      a <- priv_f(sdp, sn) - priv_f(sdp, st)
      if (log(stats::runif(1)) < a) {
        counter <- counter + 1
        dmat[j, ] <- xs
        st <- sn
      }
    }
    accept_rate[i] <- counter / nobs
    if (!is.null(prg_bar)) prg_bar(message = sprintf("Iteration %g", i))
  }
  
  theta_mat <- do.call(rbind, theta_list)  
  
  if (warmup > 0) {
    theta_mat <- theta_mat[-c(1:warmup), , drop = FALSE]
    accept_rate <- accept_rate[-c(1:warmup)]
  }

  list(sample = theta_mat, accept_rate = accept_rate)
}


check_data_model <- function(x) {
  data_model <- x$data_model
  init_par <- x$init_par
  sdp <- x$sdp
  
  if(length(init_par) != data_model$npar){
    return("Dimension of initial parameter does not match npar")
  }
  
  cmx <- data_model$latent_f(init_par, sdp$n) #latent_f changed
  if(!checkmate::test_matrix(cmx)) {
    return("latent_f() function must return a matrix")
  }
  stc <- data_model$st_f(cmx)
  t1 <- checkmate::test_numeric(sdp) & !checkmate::test_numeric(stc)
  t2 <- !checkmate::test_numeric(sdp) & checkmate::test_numeric(stc)
  t3 <- checkmate::test_matrix(sdp) & !checkmate::test_matrix(stc)
  t4 <- !checkmate::test_matrix(sdp) & checkmate::test_matrix(stc)
  if(t1 | t2 | t3 | t4) {
    return("st_f() must return the same data type as sdp")
  }
  if(!checkmate::test_number(data_model$priv_f(sdp,cmx)) | !checkmate::test_scalar(data_model$priv_f(sdp,cmx))) {
    return("priv_f() must return a scalar number")
  }
  if(!checkmate::test_class(data_model$post_f(cmx, init_par), "numeric")) {
    return("post_f() must return a numeric vector")
  }
  
  TRUE
}


post_f <- function(dmat, theta, prior_mu, prior_lambda, prior_Psi, prior_nu, 
                   rejection_sampler = FALSE, u) {
  
  n <- nrow(dmat)
  sample_mean <- colMeans(dmat)
  
  S <- cov(dmat) * (n - 1)
  
  lambda_n <- prior_lambda + n
  nu_n <- prior_nu + n
  mu_n <- (prior_lambda * prior_mu + n * sample_mean) / lambda_n
  
  mean_diff <- matrix(sample_mean - prior_mu, nrow = length(sample_mean))
  Psi_n <- prior_Psi + S + (prior_lambda * n / lambda_n) * tcrossprod(mean_diff)
  
  repeat {
    #Sigma_new <- riwish(nu_n, Psi_n)
    Sigma_new <- rinvwishart(n = 1, nu = nu_n, Omega = Psi_n)
    
    mat_Sigma_new <- matrix(Sigma_new, nrow=4, ncol=4)
    
    mu_new <- mvrnorm(1, mu = mu_n, Sigma = mat_Sigma_new / lambda_n)  
    
    if (!rejection_sampler || sqrt(sum(mu_new^2)) < u) {
      break
    }
  }
  
  return(list(mu = mu_new, Sigma = mat_Sigma_new))
}



latent_f <- function(theta, n) {

  x <- mvrnorm(n=n, mu=theta$mu, Sigma=theta$Sigma)
  
  return(x)
  
}

#' Summarise dpout object.
#'
#' @param object dp_out object
#' @param ... optional arguments to `summarise_draws()`.
#'
#' @return a summary table of MCMC statistics.
#' @export
summary.dpout <- function(object, ...) {
  posterior::summarise_draws(object$chain, ...)
}

#' Plot dpout object.
#'
#' @param x dp_out object.
#' @param ... optional arguments to `mcmc_trace()`.
#' @return trace plots.
#' @export
plot.dpout <- function(x, ...) {
  bayesplot::mcmc_trace(x$chain, ...)
}

new_dpout <- function(theta, accept_mat, varnames = NULL) {
  dp_obj <- do.call(rbind, theta)
  nr <- length(theta) * nrow(theta[[1]])
  nc <- ncol(theta[[1]])
  vn <- paste0("theta", 1:nc)
  if(!is.null(varnames)) vn <- varnames
  dl <- list(draw = as.character(1:nr),
             variable = vn)
  attr(dp_obj, 'dimnames') <- dl
  attr(dp_obj, 'nchains') <- length(theta)
  e1 <- posterior::as_draws_matrix(dp_obj)
  e2 <- accept_mat
  structure(list(chain = e1, accept_prob = accept_mat), class = c("dpout"))
}

dmod <- new_privacy(post_f = post_f,
                    latent_f = latent_f,
                    priv_f = priv_f,
                    st_f = st_f,
                    npar = 1)

niter <- 100

pb <- txtProgressBar(min = 0, max = niter, style = 3)


out <- dapper_sample(data_model = dmod,
                     sdp = sdp,
                     init_par = list(mu=rep(0,d), Sigma=diag(d)),
                     niter = niter, warmup=0)
#print('out:')
#print(out)

if (is.list(out$chain)) {
  out$chain <- do.call(rbind, out$chain[,1])
}

#summary(out)

if (round(sqrt(d), digits=0) != sqrt(d)) {
  par(mfrow = c(floor(sqrt(d)) + 1, floor(sqrt(d)) + 1))
} else {
  par(mfrow = c(sqrt(d), sqrt(d)))
}

for (i in 1:d) {
  plot(out$chain[,i], type='l')
}
