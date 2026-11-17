################################################################################
# output_results.R
#
# Summarises the fits made by fit_models.R, for every virus and model.
#
#   results/DIC_table.xlsx / .csv                                  DIC of the four models for each virus
#   results/parameters/parameters_<virus>_<name>_model_<m>.png     posterior density of each parameter
#   results/seroprevalence/seroprevalence_<virus>.png              posterior seroprevalence (median and 95% band)
#                                                                  of the four models, by colony, with the observed
#                                                                  seroprevalence (95% interval) of the bats sampled
#
# The first burn_in iterations of each chain are dropped before summarising.
# Each plot is also shown on screen (in RStudio, the Plots pane) as it is saved.
################################################################################

# setwd("C:/Users/ao571/OneDrive - University of Sussex/Desktop/Seroprevalence_Simulation_Final/RealData_v2")
dir.create("results/parameters", recursive = TRUE, showWarnings = FALSE); dir.create("results/seroprevalence", showWarnings = FALSE)

library(odin)
library(readxl)
library(writexl)
library(ggplot2)

ggsave <- function(filename, plot, ...) { print(plot); ggplot2::ggsave(filename, plot, ...) }   # show each plot on screen

names_models <- c("SEIR", "SEID_SER")
models       <- c("2.1", "2.2")
burn_in      <- 2000
n_draws      <- 100          # posterior draws used for the predicted seroprevalence
set.seed(21)

source("prepare_inputs.R")          # gives viruses, pop and obs
plot_years <- min(data$Model_year):max(data$Model_year)   # the years with data

# Each compartmental model's runSim, loaded into its own environment
engines <- lapply(setNames(names_models, names_models), function(name) {
  engine <- new.env()
  source(paste0("../Useful_Functions/model_", name, ".R"), local = engine)
  engine
})

# Posterior draws (after burn-in) of one fit
posteriorOf <- function(virus, name, model) {
  chain <- read.csv(paste0("fits/MCMC_", virus, "_", name, "_model_", model, ".csv"))
  chain[-seq_len(burn_in), names(chain) != "logL"]
}


# 1. DIC table, in the layout of the workflow ---------------------------------------

dic  <- do.call(rbind, lapply(list.files("fits", pattern = "^DIC_.*_model_.*csv$", full.names = TRUE), read.csv))
rows <- c("SEIR 2.1", "SEIR 2.2", "SEID/SER 2.1", "SEID/SER 2.2")
table <- data.frame(model = rows)
for (virus in names(viruses)) {
  d <- dic[dic$virus == virus, ]
  table[[virus]] <- round(d$DIC[match(rows, d$model)], 1)
}
write.csv(table, "results/DIC_table.csv", row.names = FALSE)
write_xlsx(table, "results/DIC_table.xlsx")
print(table)


# 2. Posterior density of each parameter --------------------------------------------

plotParameters <- function(posterior, title) {
  density <- data.frame(parameter = rep(names(posterior), each = nrow(posterior)),
                        value     = unlist(posterior))
  ggplot(density, aes(value)) +
    geom_density(fill = "steelblue", colour = "steelblue", alpha = 0.4) +
    facet_wrap(~ parameter, scales = "free", ncol = 3) +
    labs(title = title, x = NULL, y = "Density") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
}


# 3. Posterior seroprevalence -------------------------------------------------------
# For each posterior draw, runSim gives the seroprevalence of every colony and day;
# the curves are the median and 95% band across draws. Only the colonies sampled
# for the virus, up to their last sampled year, are simulated.

predictSeroprevalence <- function(posterior, name, model, pop_v) {
  inputs <- read_xlsx(paste0("../Simulation/", name, "/input/input_for_simulation_model_", model, ".xlsx"),
                      sheet = "parameters")
  fixed  <- as.list(inputs[1, ])
  draws  <- posterior[sample(nrow(posterior), n_draws), , drop = FALSE]

  runs <- lapply(seq_len(nrow(draws)), function(i) {
    par <- fixed
    par[names(draws)] <- draws[i, ]
    do.call(rbind, lapply(unique(pop_v$colony), function(colony)
      engines[[name]]$runSim(par, pop_v[pop_v$colony == colony, ])))
  })
  grid <- runs[[1]][c("colony", "year", "day")]
  grid$abs_day <- (grid$year - min(plot_years)) * 365 + grid$day
  band <- function(column, group) {
    values <- sapply(runs, function(r) r[[column]])
    data.frame(grid[c("colony", "year", "abs_day")], group = group,
               median = apply(values, 1, median, na.rm = TRUE),
               lower  = apply(values, 1, quantile, 0.025, na.rm = TRUE),
               upper  = apply(values, 1, quantile, 0.975, na.rm = TRUE))
  }
  curves <- rbind(band("prev_adult", "Adult"), band("prev_juv", "Juvenile"))
  curves <- curves[curves$year %in% plot_years, ]
  curves$model <- paste(sub("_", "/", name), model)
  curves
}

plotSeroprevalence <- function(curves, obs_v, pop_v, title) {
  # Observed proportions with exact (Clopper-Pearson) 95% intervals
  obs_v$abs_day  <- (obs_v$year - min(plot_years)) * 365 + obs_v$day
  obs_v$prop     <- obs_v$positive / obs_v$sampling
  obs_v$ci_lower <- ifelse(obs_v$positive == 0, 0, qbeta(0.025, obs_v$positive, obs_v$sampling - obs_v$positive + 1))
  obs_v$ci_upper <- ifelse(obs_v$positive == obs_v$sampling, 1, qbeta(0.975, obs_v$positive + 1, obs_v$sampling - obs_v$positive))
  obs_v$group    <- obs_v$age

  p <- pop_v[pop_v$year %in% plot_years, ]
  hibernation <- data.frame(colony = p$colony,
                            xmin = (p$year - min(plot_years)) * 365 + p$hib_start,
                            xmax = (p$year - min(plot_years)) * 365 + p$hib_end)

  ggplot() +
    geom_rect(data = hibernation, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "grey85", alpha = 0.6) +
    geom_ribbon(data = curves, aes(abs_day, ymin = lower, ymax = upper, fill = group), alpha = 0.2, na.rm = TRUE) +
    geom_line(data = curves, aes(abs_day, median, colour = group), linewidth = 0.9, na.rm = TRUE) +
    geom_errorbar(data = obs_v, aes(abs_day, ymin = ci_lower, ymax = ci_upper, colour = group), width = 6) +
    geom_point(data = obs_v, aes(abs_day, prop, colour = group), size = 2) +
    facet_grid(colony ~ model) +
    scale_colour_manual(values = c(Adult = "blue", Juvenile = "red")) +
    scale_fill_manual(values = c(Adult = "blue", Juvenile = "red")) +
    scale_x_continuous(breaks = (seq_along(plot_years) - 1) * 365, labels = paste("Year", plot_years)) +
    coord_cartesian(xlim = c(0, length(plot_years) * 365), ylim = c(0, 1)) +
    labs(title = title,
         subtitle = "Lines and bands: posterior median and 95% band. Points: observed proportion with 95% interval. Grey: hibernation.",
         x = NULL, y = "Seroprevalence", colour = NULL, fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
}


# 4. Plot every virus ---------------------------------------------------------------

for (virus in names(viruses)) {
  obs_v <- obs[[virus]]
  last  <- tapply(obs_v$year, obs_v$colony, max)
  pop_v <- pop[pop$colony %in% names(last) & pop$year <= last[pop$colony], ]

  curves <- NULL
  for (name in names_models) {
    for (model in models) {
      posterior <- posteriorOf(virus, name, model)
      ggsave(paste0("results/parameters/parameters_", virus, "_", name, "_model_", model, ".png"),
             plotParameters(posterior, paste0("Posterior parameters: ", virus, ", ", sub("_", "/", name), " model ", model)),
             width = 11, height = 6, dpi = 300, bg = "white")
      curves <- rbind(curves, predictSeroprevalence(posterior, name, model, pop_v))
    }
  }
  curves$model <- factor(curves$model, levels = rows)
  ggsave(paste0("results/seroprevalence/seroprevalence_", virus, ".png"),
         plotSeroprevalence(curves, obs_v, pop_v, paste("Posterior seroprevalence:", virus)),
         width = 12, height = 2.2 * length(last) + 1.5, dpi = 200, bg = "white", limitsize = FALSE)
  cat(virus, "done\n")
}
