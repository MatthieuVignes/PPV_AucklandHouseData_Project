#######################################################
# Joins hourly indoor data (output of aggregate_hourly()) with
# (1) the house information table (read_house_design()), and
# (2) outdoor ARC meteorology (read_arc_met()), via each house's weather station,
#     and adds the normalised-time variable trial_day.
#
# Pipeline organisation:
#   parse_house_file()  ->  aggregate_hourly()  ->  merge_outdoor()
#                                          (+ read_house_design, read_arc_met)
#
# Keys / conventions
# - Indoor `hour` (start of hour, from aggregate_hourly) is converted to 'naive' clock 
#   time with clock_time() so it lines up with ARC `hour` (hour START = ARC hour-ending 
#   stamp minus 1h).
# - house + year  -> house information (suburb, group, station, dates)
# - station + hour -> outdoor weather
#
# New columns added to the indoor hourly table
#  suburb, group, station            house information
#  state_expected                    PPV state implied by the design
#  state_mismatch                    TRUE if file's C/A state differs from design
#  week2_start                       install day (pseudo-install for controls)
#  trial_day                         whole days since week2_start
#                                    (negative = before install; 0 = install day)
#  outdoor_T, outdoor_RH, solar_rad, wind_speed, wind_dir, rain
#######################################################

#' @param hourly_long output of aggregate_hourly() (needs house, year, week,
#'   state, room, variable, hour, mean, min, max, n_obs)
#' @param design output of read_house_design()
#' @param arc output of read_arc_met() / read_arc_met_many()
merge_outdoor <- function(hourly_long, design, arc) {
  x <- hourly_long
  x$.row <- seq_len(nrow(x))
  x$hour_clock <- clock_time(x$hour)
  x$date_clock <- as.Date(format(x$hour_clock, "%Y-%m-%d"))

  # 1. house information (by house + year)
  des_cols <- c("house", "year", "suburb", "group", "station", "week2_start")
  x <- merge(x, design[, des_cols], by = c("house", "year"), all.x = TRUE, sort = FALSE)
  if (anyNA(x$group)) warning("Some rows have no matching house/year in the design table.")

  # 2. design-implied state and normalised time
  x$state_expected <- expected_state(x$group, x$year, x$week)
  x$state_mismatch <- !is.na(x$state) & x$state != x$state_expected
  x$trial_day <- as.integer(x$date_clock - x$week2_start)

  # 3. outdoor weather (by station + hour start)
  arc2 <- arc
  names(arc2)[names(arc2) == "hour"] <- "hour_clock"
  keep <- c("station", "hour_clock", "outdoor_T", "outdoor_RH", "solar_rad",
            "wind_speed", "wind_dir", "rain")
  x <- merge(x, arc2[, keep], by = c("station", "hour_clock"), all.x = TRUE, sort = FALSE)

  x <- x[order(x$.row), ]
  x$.row <- NULL
  rownames(x) <- NULL
  x
}

#' Quick coverage / consistency report for a merged table
check_merge <- function(m) {
  cat("Rows:", nrow(m), "| houses:", length(unique(m$house)),
      "| years:", paste(sort(unique(m$year)), collapse = ", "), "\n")
  cat("Rows with no outdoor temperature match:",
      sum(is.na(m$outdoor_T)), sprintf("(%.1f%%)\n", 100 * mean(is.na(m$outdoor_T))))
  if (any(is.na(m$outdoor_T))) {
    miss <- aggregate(hour_clock ~ house + year, data = m[is.na(m$outdoor_T), ],
                      FUN = function(z) sprintf("%d hours, %s to %s", length(z),
                                                format(min(z)), format(max(z))))
    cat("Unmatched hours by house-year:\n"); print(miss)
  }
  cat("Rows whose file state differs from design-expected state:",
      sum(m$state_mismatch, na.rm = TRUE), "\n")
  if (any(m$state_mismatch, na.rm = TRUE)) {
    print(unique(m[m$state_mismatch, c("house", "year", "week", "state", "state_expected")]))
  }
  invisible(m)
}

#######################################################
# Demo on one house-year file (paths from the command line):
#   Rscript R/merge_house_outdoor.R <house_file> <house_id> <design.csv> \
#           <station_lookup.csv> <ARC_2008.csv> [<ARC_2009.csv> ...]
#######################################################
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
  source(file.path(script_dir, "parse_house_file.R"))
  source(file.path(script_dir, "aggregate_hourly.R"))
  source(file.path(script_dir, "read_arc_met.R"))
  source(file.path(script_dir, "read_house_design.R"))

  if (length(args) < 5) stop("Usage: Rscript merge_house_outdoor.R <house_file> <house_id> <design.csv> <station_lookup.csv> <ARC.csv> [...]")
  house_file <- args[1]; house_id <- as.integer(args[2])

  parsed  <- parse_house_file(house_file, house_id = house_id)
  hourly  <- aggregate_hourly(parsed$measurements)
  design  <- read_house_design(args[3], args[4])
  arc     <- read_arc_met_many(args[5:length(args)])

  merged <- merge_outdoor(hourly, design, arc)
  check_merge(merged)

  cat("\nMerged table (first rows):\n")
  print(head(merged[order(merged$room, merged$variable, merged$hour_clock),
                    c("house", "year", "week", "state", "state_expected", "room", "variable",
                      "hour_clock", "mean", "suburb", "station", "trial_day",
                      "outdoor_T", "outdoor_RH", "solar_rad")], 12))
}
