################################################################################
# fit_models.R
#
# Fits the four models (SEIR 2.1, SEIR 2.2, SEID/SER 2.1, SEID/SER 2.2) to each virus
# by Metropolis-Hastings MCMC, and computes the DIC of every fit.
#
# Workflow (one pass per virus and model)
#   cleaned/serology_model_ready.csv and cleaned/population_by_colony_year.csv  (made by prepare_data.R, read by prepare_inputs.R)
#     -> binomial likelihood of the sampled counts, given runSim() with trial parameters
#     -> fits/MCMC_<virus>_<name>_model_<m>.csv   (parameters and log-likelihood per iteration)
#     -> fits/DIC_<virus>_<name>_model_<m>.csv    (DIC of the fit)
#
# Parameters estimated (all others, e.g. the death rate, are fixed):
#   model 2.1: sigma, sigma_hib, sigma_m, beta0, beta1, beta2   (+ pi for SEID/SER)
#   model 2.2: sigma, sigma_hib, sigma_m, beta0                  (+ pi for SEID/SER)
# Priors are uniform with wide bounds. With no true values for real data, each chain
# starts at the scenario 1 values of the simulation input.
#
# DIC: deviance = -2 * log-likelihood. Dbar is its mean over the posterior draws,
# Dhat its value at the posterior mean, pD = Dbar - Dhat and DIC = Dbar + pD.
################################################################################

# setwd("C:/Users/ao571/OneDrive - University of Sussex/Desktop/Seroprevalence_Simulation_Final/RealData_v2")
dir.create("fits", showWarnings = FALSE)

library(odin)
library(readxl)

names_models <- c("SEIR", "SEID_SER")
models       <- c("2.1", "2.2")
n_iter       <- 5000
burn_in      <- 2000      # iterations discarded before summarising
set.seed(2026)

source("prepare_inputs.R")          # gives viruses, pop and obs

# To fit one virus per session, name it when running: Rscript fit_models.R WCL
# (in RStudio, or with no name, all three viruses are fitted one after the other)
if (length(commandArgs(trailingOnly = TRUE))) viruses <- viruses[commandArgs(trailingOnly = TRUE)]

# Uniform priors: the same wide bounds for every virus and model
priors <- list(sigma = c(0, 0.1), sigma_hib = c(0, 0.1), sigma_m = c(0, 0.5),
               beta0 = c(0, 0.1), beta1 = c(0, 0.1), beta2 = c(0, 1), pi = c(0, 1))


# 1. Metropolis-Hastings sampler --------------------------------------------------
# logLik: function of a named parameter vector. Proposals are log-normal random
# walks, with a step size per parameter that is adapted every 50 iterations to keep
# the acceptance rate between 20% and 40%. Every iteration prints the current
# values; every 50 iterations also prints the acceptance rates.
# Returns the chain and the log-likelihood at every iteration.

runMCMC <- function(logLik, start, lower, upper, n_iter, label = "") {
  chain <- matrix(NA, n_iter, length(start), dimnames = list(NULL, names(start)))
  logL  <- rep(NA, n_iter)
  theta <- start
  L     <- logLik(theta)
  step  <- rep(0.05, length(start))
  accepted <- rep(0, length(start))
  chain[1, ] <- theta
  logL[1]    <- L

  for (i in 2:n_iter) {
    for (j in seq_along(theta)) {
      new    <- theta
      new[j] <- theta[j] * exp(step[j] * rnorm(1))
      if (new[j] > lower[j] && new[j] < upper[j]) {
        L_new <- logLik(new)
        if (log(runif(1)) < L_new - L + log(new[j] / theta[j])) {
          theta <- new
          L     <- L_new
          accepted[j] <- accepted[j] + 1
        }
      }
    }
    cat(sprintf("%s | iteration %d of %d | %s\n", label, i, n_iter,
                paste0(names(theta), " ", signif(theta, 3), collapse = ", ")))
    if (i %% 50 == 0) {
      rate <- accepted / 50
      cat(sprintf("%s | acceptance over the last 50 iterations: %s\n", label,
                  paste0(names(start), " ", sprintf("%.2f", rate), collapse = ", ")))
      step[rate < 0.2] <- step[rate < 0.2] * 0.8
      step[rate > 0.4] <- step[rate > 0.4] * 1.2
      accepted[] <- 0
    }
    chain[i, ] <- theta
    logL[i]    <- L
  }
  list(chain = chain, logL = logL)
}


# 2. Fit every virus, compartmental model and model version -----------------------

for (name in names_models) {

  # Load this compartmental model's runSim into its own environment, because the
  # two model files use the same function names
  engine <- new.env()
  source(paste0("../Useful_Functions/model_", name, ".R"), local = engine)

  for (model in models) {
    # Fixed parameters and starting values: scenario 1 of the simulation input
    inputs <- read_xlsx(paste0("../Simulation/", name, "/input/input_for_simulation_model_", model, ".xlsx"),
                        sheet = "parameters")
    fixed  <- as.list(inputs[1, ])

    # Parameters to estimate
    estimated <- c("sigma", "sigma_hib", "sigma_m", "beta0")
    if (model == "2.1")      estimated <- c(estimated, "beta1", "beta2")
    if (name == "SEID_SER")  estimated <- c(estimated, "pi")
    start  <- unlist(fixed[estimated])
    bounds <- do.call(cbind, priors[estimated])

    for (virus in names(viruses)) {
      label <- paste(virus, name, "model", model)
      if (file.exists(paste0("fits/DIC_", virus, "_", name, "_model_", model, ".csv"))) next   # already fitted
      data_v <- obs[[virus]]

      # Only the colonies sampled for this virus, and only up to their last sampled year
      last   <- tapply(data_v$year, data_v$colony, max)
      pop_v  <- pop[pop$colony %in% names(last) & pop$year <= last[pop$colony], ]

      # Log-likelihood: binomial probability of the sampled counts, given the
      # seroprevalence that runSim gives for the trial parameters
      logLik <- function(theta) {
        par  <- fixed
        par[names(theta)] <- theta
        traj <- do.call(rbind, lapply(unique(pop_v$colony), function(colony)
          engine$runSim(par, pop_v[pop_v$colony == colony, ])))
        row  <- match(paste(data_v$colony, data_v$year, data_v$day),
                      paste(traj$colony, traj$year, traj$day))
        p    <- ifelse(data_v$age == "Adult", traj$prev_adult[row], traj$prev_juv[row])
        stopifnot(!anyNA(p))
        p    <- pmin(pmax(p, 1e-10), 1 - 1e-10)
        sum(dbinom(data_v$positive, data_v$sampling, p, log = TRUE))
      }

      fit <- runMCMC(logLik, start, bounds[1, ], bounds[2, ], n_iter, label)

      # DIC from the iterations after burn-in
      kept     <- (burn_in + 1):n_iter
      D_bar    <- mean(-2 * fit$logL[kept])
      D_hat    <- -2 * logLik(colMeans(fit$chain[kept, ]))
      dic      <- data.frame(virus = virus, model = paste(sub("_", "/", name), model),
                             Dbar = D_bar, Dhat = D_hat, pD = D_bar - D_hat, DIC = 2 * D_bar - D_hat)

      write.csv(cbind(fit$chain, logL = fit$logL),
                paste0("fits/MCMC_", virus, "_", name, "_model_", model, ".csv"), row.names = FALSE)
      write.csv(dic, paste0("fits/DIC_", virus, "_", name, "_model_", model, ".csv"), row.names = FALSE)
      cat(label, "done: DIC", round(dic$DIC, 1), "\n")
    }
  }
}
