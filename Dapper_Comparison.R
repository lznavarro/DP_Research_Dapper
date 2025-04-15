args <- commandArgs(trailingOnly = TRUE)
print(args)

library(doParallel)
library(foreach)

source("DP_Functions.R")
source("privacy_class.R")
source("Dapper_Test.R")

num_cores <- 100
cl <- makeCluster(num_cores)
registerDoParallel(cl)


cat("\n", "# of cores: ", num_cores, "\n")

niter <- as.numeric(args[1])
n      <- as.numeric(args[2])
d      <- as.numeric(args[3])
p      <- as.numeric(args[4])
mu     <- as.numeric(args[5])
sa     <- as.numeric(args[6])
u <- 20
reps <- num_cores
warmup <- floor(niter/2)

start_time <- Sys.time()


# Run the simulation
results <- foreach(rep = 1:reps, .combine = rbind, .packages = c("MASS", "stats", "pracma"),
                   .export = c("dapper_sample")) %dopar% {
                     
                     set.seed(rep)
                     
                     source('DP_Functions.R')
                     source('privacy_class.R')
                     source('Dapper_Test.R')
                     
                     x <- matrix(rnorm(n * d, 
                                       mean = mu, 
                                       sd = sa), 
                                 ncol = d)
                     
                     randomVars <- rInstanceOptimal(x, n, d, u, p)
                     rInstance_result <- randomVars$noisy_mean
                     mse_rInstance <- mean((rInstance_result - mu)^2)
                     
                     dmod <- new_privacy(post_f = post_f,
                                         latent_f = latent_f,
                                         priv_f = priv_f,
                                         st_f = st_f,
                                         npar = 1)
                     
                     sdp <- list(
                       HD = randomVars$HD,
                       Cs = randomVars$Cs,
                       C_clipped = randomVars$C_clipped,
                       svect_clipped = randomVars$svect_clipped,
                       noisy_mean = randomVars$noisy_mean,
                       n = randomVars$n,
                       svects = randomVars$svects,
                       u = randomVars$u,
                       p = randomVars$p,
                       T = randomVars$T,
                       x = randomVars$x,
                       d = randomVars$d
                     )
                     
                     out <- dapper_sample(data_model = dmod, sdp = sdp,
                                          init_par = list(mu = rep(0, d), Sigma = diag(d)),
                                          niter = niter, warmup = warmup)
                     
                     chain_matrix <- do.call(rbind, out$chain[1:nrow(out$chain)])
                     chain_length <- nrow(chain_matrix)
                     
                     dapper_means <- colMeans(chain_matrix)
                     mse_Dapper <- mean((dapper_means - mu)^2)
                     
                     avg_post_means <- mean(dapper_means)
                     avg_post_vars <- mean(apply(chain_matrix, 2, var))
                     
                     data.frame(
                       n = n,
                       d = d,
                       p = p,
                       mu = mu,
                       sa = sa,
                       replicate = rep,
                       mse_rInstance = mse_rInstance,
                       mse_Dapper = mse_Dapper,
                       iterations = niter,
                       avg_post_means = avg_post_means,
                       avg_post_vars = avg_post_vars,
                       chain_length = chain_length,
                       chain = I(list(chain_matrix))
                     )
                   }

end_time <- Sys.time()
elapsed_time <- end_time - start_time

cat(sprintf("Elapsed time for parameters: n = %d, d = %d, p = %.2f, mu = %.2f, sa = %.2f: %s seconds\n",
            n, d, p, 
            mu, sa, as.numeric(elapsed_time)))

# Final results
final_results <- do.call(rbind, list(results))

outdir <- sprintf("/home/navarr72/dpwork/dapper_comparison_output_%s_%s_%s_%s_%s_%s",
                  niter, n, d, p, mu, sa)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

save(final_results, file = file.path(outdir, "sim_results.RData"))

stopCluster(cl)
