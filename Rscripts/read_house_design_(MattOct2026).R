##################################################
# Reads house_dates_PPV.csv (experimental design) and the house -> weather-station 
# lookup, and returns one row per house*year ("house information" table used by the merge).
#
# house_dates_PPV.csv layout: two header rows (year labels, then column names), 
# then 30 rows: ID, Suburb, Group, then for 2008: Round.Set, Week 1 start, Week 2 start, 
# Week 2 stop, then the same four for 2009.
# Dates are written "28-Jul" (no year): the year comes from the block.
#
# Design logic encoded here:
# - Control houses: Control state in both weeks, both years.
# - Active (treatment) houses: Control in week 1 of 2008, PPV installed at the start of 
#   week 2 (= "Week 2 start"), Active in week 2 of 2008 and in both weeks of 2009.
# - Controls have no install date; their "Week 2 start" is used as a pseudo-install day 
#   so trial_day is comparable across groups.
##################################################

month_num <- c(Jan = 1, Feb = 2, Mar = 3, Apr = 4, May = 5, Jun = 6,
               Jul = 7, Aug = 8, Sep = 9, Oct = 10, Nov = 11, Dec = 12)

# "28-Jul" + 2008 -> Date  (locale-independent)
parse_dm <- function(x, year) {
  p <- strsplit(trimws(x), "-")
  d <- vapply(p, function(z) as.integer(z[1]), integer(1))
  m <- month_num[vapply(p, function(z) z[2], character(1))]
  as.Date(sprintf("%04d-%02d-%02d", year, m, d))
}

#' @param design_path path to house_dates_PPV.csv
#' @param station_path path to house_station_lookup.csv
#' @return data.frame: house, suburb, group, station, year, round_set,
#'   week1_start, week2_start, week2_stop, week1_days, week2_days
read_house_design <- function(design_path, station_path = NULL) {
  raw <- read.csv(design_path, header = FALSE, skip = 2,
                  colClasses = "character", strip.white = TRUE)
  base <- data.frame(house  = as.integer(raw[[1]]),
                     suburb = raw[[2]],
                     group  = raw[[3]],
                     stringsAsFactors = FALSE)
  blocks <- list(`2008` = 4:7, `2009` = 8:11)
  out <- lapply(names(blocks), function(yr) {
    cl <- blocks[[yr]]; y <- as.integer(yr)
    d <- cbind(base,
               year        = y,
               round_set   = raw[[cl[1]]],
               week1_start = parse_dm(raw[[cl[2]]], y),
               week2_start = parse_dm(raw[[cl[3]]], y),
               week2_stop  = parse_dm(raw[[cl[4]]], y),
               stringsAsFactors = FALSE)
    d
  })
  des <- do.call(rbind, out)
  des$week1_days <- as.integer(des$week2_start - des$week1_start)
  des$week2_days <- as.integer(des$week2_stop  - des$week2_start)

  if (!is.null(station_path)) {
    st <- read.csv(station_path, stringsAsFactors = FALSE)[, c("house", "station")]
    des <- merge(des, st, by = "house", all.x = TRUE)
    if (anyNA(des$station)) warning("Some houses have no weather station in ", station_path)
  }
  des <- des[order(des$year, des$house), ]
  rownames(des) <- NULL
  des
}

#' Expected PPV state from the design, for a given house, year and date.
#' Active houses: Control during 2008 week 1, Active afterwards.
expected_state <- function(group, year, week) {
  ifelse(group == "Control", "Control",
         ifelse(year == 2008 & week == 1, "Control", "Active"))
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  des <- read_house_design(if (length(args) >= 1) args[1] else "data/raw/house_dates_PPV.csv",
                           if (length(args) >= 2) args[2] else "data/house_station_lookup.csv")
  str(des)
  cat("\nHouses per year x group:\n"); print(table(des$year, des$group))
  cat("\nHouses per station:\n"); print(table(des[des$year == 2008, "station"]))
  cat("\nWeeks whose length is not 7 days (check against field notes):\n")
  print(des[des$week1_days != 7 | des$week2_days != 7,
            c("house", "year", "week1_start", "week2_start", "week2_stop", "week1_days", "week2_days")])
}
