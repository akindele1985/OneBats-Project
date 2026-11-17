################################################################################
# prepare_data.R      STEP 1 of the real-data analysis: clean and complete the data
#
# Reads    raw/Serology_complete.xlsx      the bats (never changed)
#          raw/Onebat_Data_Entry7.xlsx     sampling dates and colony sizes (Sampling Occasion),
#                                          hibernation (Colony History), country (Colony)
# Writes   cleaned/serology_model_ready.csv         one row per bat used by the models
#          cleaned/population_by_colony_year.csv    one row per colony and model year
#          cleaned/cleaned_data.xlsx                both tables as two sheets
# The cleaning follows reference/pierre_check2_original.Rmd. Run the whole script
# (Source button, or source("prepare_data.R")); it is plain R, not an Rmd.
#
# Rules (supervisor)
#   * Model year starts at the end of hibernation, 1 April (day 1). April 2024 to March 2025
#     is year 4, April 2025 to March 2026 is year 5. A March sample belongs to the year ending.
#   * Day = sampling date of the occasion (workbook), counted from 1 April.
#   * Birth day by country: France 22 June (83), Italy 1 July (92), Hungary 22 June (83),
#     Spain 22 June (83), Cyprus 23 May (53).
#   * Hibernation start and end: Colony History sheet (mid-November = day 229, early March =
#     340, mid-March = 350 ...). A colony without them gets the average over the colonies that have them.
#   * Colony size: largest size of the colony's sampling occasions in that model year (this
#     covers HU1 and HU2). Any other missing size: the colony's average, else the country's
#     average, else the average of all colonies.
#   * Age: Adult or Young; blank ages are left out. A Young bat sampled on or before the
#     colony's birth day counts as Adult. Serology NC is left out. Rows with no serology
#     result at all are left out. Recaptures stay as separate rows.
################################################################################

# setwd("C:/Users/ao571/OneDrive - University of Sussex/Desktop/Seroprevalence_Simulation_Final/RealData_v2")

library(readxl)
library(writexl)

birth_by_country <- c(FR = 83, ITA = 92, HU = 83, ES = 83, CY = 53)   # day of the model year
viruses  <- c("WCBV", "LLBV", "LLOV")
workbook <- "raw/Onebat_Data_Entry7.xlsx"

readSheet  <- function(sheet) as.data.frame(read_xlsx(workbook, sheet = sheet, col_types = "text"))[-1, ]   # row 1 is a description
colonyCode <- function(name) sub("\\s*-.*$", "", trimws(name))        # "HU1 - Mad" and "CY2a-Pegeia" -> "HU1", "CY2a"
modelYear  <- function(date) as.integer(format(date, "%Y")) - 2020 - (as.integer(format(date, "%m")) < 4)   # 2024 = year 4


# 1. OneBAT workbook: sampling occasions (date, colony size) ------------------------------
# Dates typed as text are day/month/year. Dates that Excel turned into real dates had their
# day and month swapped (when the day was 12 or less), so they are swapped back. Occasions
# named with their date (HU1-20240429) take that date.

occasions <- readSheet("Sampling Occasion")[, 1:6]
names(occasions) <- c("colony", "occasion", "date", "temperature", "humidity", "size")
occasions$occasion <- trimws(occasions$occasion)      # "HL25 7 CY1 " has a space at the end in the workbook

raw   <- occasions$date
excel <- grepl("^[0-9]+$", raw)
occasions$date        <- as.Date(raw, format = "%d/%m/%Y")
real                  <- as.Date(as.numeric(raw[excel]), origin = "1899-12-30")
swapped               <- as.Date(format(real, "%Y-%d-%m"), format = "%Y-%m-%d")   # NA when the day is above 12
occasions$date[excel] <- as.Date(ifelse(is.na(swapped), real, swapped), origin = "1970-01-01")
named <- grepl("-[0-9]{8}$", occasions$occasion)
occasions$date[named] <- as.Date(sub(".*-([0-9]{8})$", "\\1", occasions$occasion[named]), format = "%Y%m%d")
occasions$size       <- as.numeric(occasions$size)
occasions$code       <- colonyCode(occasions$colony)
occasions$model_year <- modelYear(occasions$date)
stopifnot(!anyDuplicated(occasions$occasion))


# 2. OneBAT workbook: country, birth day, hibernation start and end ------------------------

colonies <- readSheet("Colony")
country  <- setNames(trimws(colonies[, 2]), colonyCode(colonies[, 1]))

dayOfYear <- function(month, day) {                  # day of the model year (1 April = 1)
  start <- as.Date(sprintf("%d-04-01", ifelse(month >= 4, 2023, 2022)))
  as.numeric(as.Date(sprintf("2023-%02d-01", month)) + (day - 1) - start) + 1
}

# Free text such as "15-30 June", "no precise data (mid-november)", "End of November -
# early of December" -> day of the model year. "Non applicable" and empty give NA.
parseDay <- function(text) {
  text <- tolower(text)
  if (is.na(text) || grepl("non applicable|not applicable", text)) return(NA)
  parts <- strsplit(gsub("no precise data|\\(|\\)", "", text), " - ")[[1]]
  round(mean(sapply(parts, function(part) {
    months <- which(sapply(tolower(month.name), grepl, part, fixed = TRUE))
    range  <- regmatches(part, regexec("([0-9]+)-([0-9]+)", part))[[1]]
    day    <- if (length(range)) mean(as.numeric(range[2:3])) else
              if (grepl("early", part)) 5 else if (grepl("mid", part)) 15 else if (grepl("end", part)) 25 else 15
    mean(dayOfYear(months, day))
  })))
}

history <- readSheet("Colony History")[, 1:8]
names(history) <- c("year", "colony", "type", "arrival", "birth", "departure", "hib_start", "hib_end")
history$code <- colonyCode(history$colony)

firstValue <- function(column) {                      # one value per colony: the first one given
  days <- unname(sapply(history[[column]], parseDay))
  sapply(split(days, history$code), function(v) if (all(is.na(v))) NA_real_ else v[!is.na(v)][1])
}
timing <- data.frame(code = c(unique(history$code), setdiff(names(country), history$code)))
timing$country   <- unname(country[timing$code])
timing$birth_day <- unname(birth_by_country[timing$country])
for (column in c("hib_start", "hib_end")) {
  value <- unname(firstValue(column)[timing$code])
  timing[[column]] <- ifelse(is.na(value), round(mean(value, na.rm = TRUE)), value)   # no value: the average
  timing[[paste0(column, "_given")]] <- !is.na(value)
}


# 3. Serology: cleaning (as in pierre_check2) -----------------------------------------------

s <- as.data.frame(read_xlsx("raw/Serology_complete.xlsx", col_types = "text"))[, 1:15]
names(s) <- c("Id", "Colony", "Sampling", "Country", "Month", "Year", "BatTag", "Sex", "Age",
              "RStatus", "Colony_size_file", "Rphase", "WCBV", "LLBV", "LLOV")
for (i in seq_along(s)) s[[i]][s[[i]] %in% c("", "NC")] <- NA            # blank and NC are missing

s$Colony   <- sub("^([A-Z]+[0-9]+[ab]?)\\s*-\\s*", "\\1 - ", trimws(s$Colony))   # CY2a-Pegeia -> CY2a - Pegeia
s$Sampling <- trimws(s$Sampling)
s$Colony_size_file <- as.numeric(s$Colony_size_file)
for (v in viruses) s[[v]] <- as.numeric(s[[v]])                          # 0, 1 or NA

no_serology <- is.na(s$WCBV) & is.na(s$LLBV) & is.na(s$LLOV)
no_age      <- !(s$Age %in% c("Adult", "Young"))
s <- s[!no_serology & !no_age, ]
cat(sum(no_serology | no_age), "rows left out (no serology or no clear age),", nrow(s), "kept\n")


# 4. Add the OneBAT information to every bat ---------------------------------------------

s$code <- colonyCode(s$Colony)
at <- match(s$Sampling, occasions$occasion)
stopifnot(!anyNA(at))                                 # every bat belongs to a known sampling occasion

s$Date       <- occasions$date[at]
s$Model_year <- modelYear(s$Date)
s$Day        <- as.numeric(s$Date - as.Date(sprintf("%d-04-01", s$Model_year + 2020))) + 1

yearly <- aggregate(size ~ code + model_year, occasions, max)            # largest size of each colony and year
yearly$country <- country[yearly$code]
colony_mean  <- tapply(yearly$size, yearly$code, mean)
country_mean <- tapply(yearly$size, yearly$country, mean)
fillSize <- function(code) {                          # colony average, else country average, else all colonies
  sapply(code, function(k) {
    options <- c(colony_mean[k], country_mean[country[k]], mean(yearly$size))
    unname(options[!is.na(options)][1])
  })
}

size_year <- yearly$size[match(paste(s$code, s$Model_year), paste(yearly$code, yearly$model_year))]
size_occ  <- occasions$size[at]
s$Colony_size        <- ifelse(!is.na(s$Colony_size_file), s$Colony_size_file,
                        ifelse(!is.na(size_occ), size_occ,
                        ifelse(!is.na(size_year), size_year, fillSize(s$code))))
s$Colony_size_source <- ifelse(!is.na(s$Colony_size_file), "serology file",
                        ifelse(!is.na(size_occ), "OneBAT sampling occasion",
                        ifelse(!is.na(size_year), "OneBAT, largest of the year", "average (colony, country or all)")))

k <- match(s$code, timing$code)
s$Birth_day          <- timing$birth_day[k]
s$Hibernation_start  <- timing$hib_start[k]
s$Hibernation_end    <- timing$hib_end[k]
s$Hibernation_source <- ifelse(timing$hib_start_given[k] & timing$hib_end_given[k], "Colony History", "average of the colonies")

# A Young bat sampled on or before the birth day (the young are not born yet) is an Adult
s$Age_group <- ifelse(s$Age == "Young" & s$Day > s$Birth_day, "Juvenile", "Adult")

bats <- s[order(s$Colony, s$Date, s$Id), c("Id", "Colony", "Country", "Sampling", "Date", "Model_year", "Day",
          "Age", "Age_group", "WCBV", "LLBV", "LLOV", "Colony_size", "Colony_size_source",
          "Birth_day", "Hibernation_start", "Hibernation_end", "Hibernation_source")]


# 5. Population table: one row per colony and model year -----------------------------------
# Years run from 1 to the last year with data; years before the first data are the run-in.

colony_name <- tapply(s$Colony, s$code, function(x) x[1])
pop <- expand.grid(code = unique(s$code), year = 1:max(s$Model_year), stringsAsFactors = FALSE)
pop$size <- yearly$size[match(paste(pop$code, pop$year), paste(yearly$code, yearly$model_year))]
missing  <- is.na(pop$size)
pop$size[missing] <- fillSize(pop$code[missing])
pop$size_source   <- ifelse(missing, "average (colony, country or all)", "OneBAT sampling occasions")
pop <- merge(pop, timing[, c("code", "birth_day", "hib_start", "hib_end")], by = "code")
pop$colony <- colony_name[pop$code]
pop <- pop[order(pop$colony, pop$year), c("colony", "year", "size", "size_source", "birth_day", "hib_start", "hib_end")]


# 6. Write the cleaned files ----------------------------------------------------------------

dir.create("cleaned", showWarnings = FALSE)
write.csv(bats, "cleaned/serology_model_ready.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(pop,  "cleaned/population_by_colony_year.csv", row.names = FALSE, fileEncoding = "UTF-8")
write_xlsx(list(Bats = bats, Population = pop), "cleaned/cleaned_data.xlsx")

cat("\nBirth day (by country), hibernation start and end per colony; * = average used\n")
shown <- timing[timing$code %in% s$code, ]
print(data.frame(code = shown$code, birth = shown$birth_day,
                 hib_start = paste0(shown$hib_start, ifelse(shown$hib_start_given, "", "*")),
                 hib_end   = paste0(shown$hib_end,   ifelse(shown$hib_end_given,   "", "*"))), row.names = FALSE)
print(table(bats$Colony_size_source))
