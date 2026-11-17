################################################################################
# prepare_inputs.R
#
# Reads the cleaned files made by prepare_data.R and gives the two inputs of the fits.
# Nothing is written. Sourced by fit_models.R and output_results.R, after their setwd().
#
#   cleaned/serology_model_ready.csv         the bats (dates, ages, serology)
#   cleaned/population_by_colony_year.csv    size, birth day, hibernation start and end
#
# Gives
#   pop   one row per colony and model year: size, birth_day, hib_start, hib_end
#   obs   a list with one table per virus (WCL, LLeidaL, Llovio):
#         colony, year, day, age, sampling (bats tested), positive (bats seropositive)
################################################################################

viruses <- c(WCL = "WCBV", LLeidaL = "LLBV", Llovio = "LLOV")   # serology column of each virus

data <- read.csv("cleaned/serology_model_ready.csv", stringsAsFactors = FALSE, fileEncoding = "UTF-8")
pop  <- read.csv("cleaned/population_by_colony_year.csv", stringsAsFactors = FALSE, fileEncoding = "UTF-8")
pop  <- pop[, c("colony", "year", "size", "birth_day", "hib_start", "hib_end")]

# One table of tested and positive bats per virus; a missing result is left out of that virus only
obs <- lapply(viruses, function(column) {
  result <- data[[column]]
  used   <- result %in% c(0, 1)
  tested <- data.frame(colony = data$Colony[used], year = data$Model_year[used], day = data$Day[used],
                       age = data$Age_group[used], sampling = 1, positive = result[used])
  out <- aggregate(cbind(sampling, positive) ~ colony + year + day + age, tested, sum)
  out[order(out$colony, out$year, out$day, out$age), ]
})

for (virus in names(obs))
  cat(sprintf("%s: %d bats tested, %d seropositive, %d colonies\n", virus,
              sum(obs[[virus]]$sampling), sum(obs[[virus]]$positive), length(unique(obs[[virus]]$colony))))
