################################################################################
# extract_trajectories.R
#
# From the MCMC chains in fits/ (no new fitting), for every virus, colony and model:
#   1. proportion of infectious bats over time, I_t / N_t
#   2. bats that died of infection over time: the number and the proportion of the
#      bats alive at that time (SEID/SER only; the SEIR model has no disease death)
#
# N_t = bats alive at day t, young and adult together, dead bats (D) excluded:
#   SEIR      S + E + I + R (+ M in the young)
#   SEID/SER  S + E1 + E2 + I + R (+ M in the young)
# Died of infection: D counts the bats that left I through death since 1 April of that
# model year (D starts at 0 each year), so
#   died_cum  = D_t, deaths since the start of the model year
#   died_new  = D_t - D_(t-1), deaths on that day
#   prop_died_cum = died_cum / N_t,   prop_died_new = died_new / N_t
#
# For each posterior draw runSim() gives the trajectories; every quantity is summarised by
# the median and the 95% band (2.5% and 97.5%) across draws. Only the colonies sampled for
# the virus are given, up to their last sampled year, as in the fits.
#
# Writes   results/trajectories/trajectories_<virus>.csv      all models, long format
#          results/trajectories/trajectories_<virus>.xlsx     one sheet per model
# day = day of the model year (1 April = 1); abs_day counts on from the first data year.
################################################################################

setwd("C:/Users/ao571/OneDrive - University of Sussex/Pierre Nouvellet's files - OneBat/Transmission_dynamics_model/Seroprevalence_Simulation_Final/RealData_v2")
dir.create("results/trajectories", recursive = TRUE, showWarnings = FALSE)

library(odin)
library(readxl)
library(writexl)

names_models <- c("SEIR", "SEID_SER")
models       <- c("2.1", "2.2")
burn_in      <- 2000
n_draws      <- 100          # posterior draws used
set.seed(21)

source("prepare_inputs.R")          # gives viruses, pop and obs
first_year <- min(data$Model_year)

engines <- lapply(setNames(names_models, names_models), function(name) {
  engine <- new.env()
  source(paste0("../Useful_Functions/model_", name, ".R"), local = engine)
  engine
})

posteriorOf <- function(virus, name, model) {
  chain <- read.csv(paste0("fits/MCMC_", virus, "_", name, "_model_", model, ".csv"))
  chain[-seq_len(burn_in), names(chain) != "logL"]
}

# Quantities of one run of runSim()
quantities <- function(r, name) {
  if (name == "SEIR") {
    alive      <- r$S_a + r$E_a + r$I_a + r$R_a + r$S_j + r$E_j + r$I_j + r$R_j + r$M_j
  } else {
    alive      <- r$S_a + r$E1_a + r$E2_a + r$I_a + r$R_a + r$S_j + r$E1_j + r$E2_j + r$I_j + r$R_j + r$M_j
  }
  out <- data.frame(alive = alive, infectious = r$I_a + r$I_j)
  out$prop_infectious <- out$infectious / out$alive
  if (name == "SEID_SER") {
    dead <- r$D_a + r$D_j
    out$died_cum      <- dead
    out$died_new      <- ave(dead, paste(r$colony, r$year), FUN = function(x) c(x[1], diff(x)))
    out$prop_died_cum <- out$died_cum / out$alive
    out$prop_died_new <- out$died_new / out$alive
  }
  out
}

trajectories <- function(posterior, name, model, pop_v) {
  inputs <- read_xlsx(paste0("../Simulation/", name, "/input/input_for_simulation_model_", model, ".xlsx"),
                      sheet = "parameters")
  fixed  <- as.list(inputs[1, ])
  draws  <- posterior[sample(nrow(posterior), n_draws), , drop = FALSE]

  grid <- NULL
  runs <- lapply(seq_len(nrow(draws)), function(i) {
    par <- fixed
    par[names(draws)] <- draws[i, ]
    r <- do.call(rbind, lapply(unique(pop_v$colony), function(colony)
      engines[[name]]$runSim(par, pop_v[pop_v$colony == colony, ])))
    if (is.null(grid)) grid <<- r[c("colony", "year", "day")]
    quantities(r, name)
  })

  out <- grid
  out$abs_day <- (out$year - first_year) * 365 + out$day
  for (metric in names(runs[[1]])) {
    values <- sapply(runs, function(r) r[[metric]])
    out[[paste0(metric, "_median")]] <- apply(values, 1, median, na.rm = TRUE)
    out[[paste0(metric, "_lower")]]  <- apply(values, 1, quantile, 0.025, na.rm = TRUE)
    out[[paste0(metric, "_upper")]]  <- apply(values, 1, quantile, 0.975, na.rm = TRUE)
  }
  out[-(1:4)] <- round(out[-(1:4)], 6)
  out
}


for (virus in names(viruses)) {
  last  <- tapply(obs[[virus]]$year, obs[[virus]]$colony, max)
  pop_v <- pop[pop$colony %in% names(last) & pop$year <= last[pop$colony], ]

  sheets <- list()
  for (name in names_models) {
    for (model in models) {
      label <- paste(name, model)
      cat("Trajectories:", virus, label, "\n")
      result <- trajectories(posteriorOf(virus, name, model), name, model, pop_v)
      sheets[[label]] <- cbind(virus = virus, model = label, result)
    }
  }
  # the SEIR sheets have no death columns; fill them with NA to stack all models in one csv
  all_columns <- unique(unlist(lapply(sheets, names)))
  long <- do.call(rbind, lapply(sheets, function(s) { s[setdiff(all_columns, names(s))] <- NA; s[all_columns] }))
  write.csv(long, paste0("results/trajectories/trajectories_", virus, ".csv"), row.names = FALSE)
  write_xlsx(sheets, paste0("results/trajectories/trajectories_", virus, ".xlsx"))
}
