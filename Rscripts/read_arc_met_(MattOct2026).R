#######################################################
# Reads one "ARC Met data <year>" file (hourly meteorology from several stations 
# laid out side by side) into ONE tidy data frame: one row per station x hour.
#
# File structure (confirmed on ARC_Met_data_2008_v1.csv):
# - a few title lines, then a row holding the STATION name above each block, then a 
#   header row ("Date","Time","Ambient Temp",...), then a units row, then data
# - each block = Date, Time, then its variables (the number of variables differs: 
#   Henderson has no Rainfall, the other two do)
# - Date is dd-mm-yy, Time is h:mm and is HOUR-ENDING: the first row is 1:00 (mean 
#   over 0:00-1:00) and the last is 24:00
# - the text "Samp<" appears in some cells (value below sampling threshold); treated 
#   as NA and flagged.
#
# Time convention:
# hour = start of the averaging hour = (hour-ending stamp) - 1h
# Times are kept as naive clock time (tz = "UTC" label, no daylight saving shifts). 
# See clock_time() below, also used by the merge step.
#######################################################

#' Convert any POSIXct to naive clock time labelled UTC
#' (keeps the printed wall-clock digits, drops the time zone, so the
#' daylight-saving change of late September cannot create NA/duplicates)
clock_time <- function(x) {
  as.POSIXct(format(x, "%Y-%m-%d %H:%M:%S"), tz = "UTC")
}

arc_var_names <- c(
  "Ambient Temp"    = "outdoor_T",  # degC
  "Wind Speed"      = "wind_speed", # m/s
  "Wind Direction"  = "wind_dir",   # deg
  "RH"              = "outdoor_RH", # %
  "Solar Radiation" = "solar_rad",  # W/m^2
  "Rainfall"        = "rain"        # mm
)

#' Read one ARC met file
#' @param path path to the CSV
#' @return data.frame(station, hour, outdoor_T, wind_speed, wind_dir,
#'   outdoor_RH, solar_rad, rain, flag_samp_below)
#'   (rain is NA for stations without a rainfall sensor)
read_arc_met <- function(path) {
  raw <- read.csv(path, header = FALSE, colClasses = "character",
                  fill = TRUE, na.strings = c("", "NA"), strip.white = TRUE)

  hdr <- which(raw[[1]] == "Date")[1]
  if (is.na(hdr)) stop("No 'Date' header row found in: ", path)
  station_row <- raw[hdr - 1, ]
  header_row  <- raw[hdr, ]
  body        <- raw[(hdr + 2):nrow(raw), ] # skip header + units rows

  date_cols <- which(header_row == "Date")
  out <- list()

  for (k in seq_along(date_cols)) {
    j0 <- date_cols[k]
    j1 <- if (k < length(date_cols)) date_cols[k + 1] - 1 else ncol(raw)
    cols <- j0:j1
    # variable columns of this block: named in the header row
    hdr_names <- unlist(header_row[cols])
    keep <- cols[!is.na(hdr_names) & hdr_names != "" & !(hdr_names %in% c("Date", "Time"))]
    st_vals <- unlist(station_row[cols]); st_vals <- st_vals[!is.na(st_vals) & st_vals != ""]
    station <- if (length(st_vals)) unname(st_vals[1]) else paste0("block", k)

    d <- body[[j0]]; tm <- body[[j0 + 1]]
    ok <- !is.na(d) & !is.na(tm) & grepl("^[0-9]{2}-[0-9]{2}-[0-9]{2}$", d)
    d <- d[ok]; tm <- tm[ok]
    end24 <- tm == "24:00"
    tm[end24] <- "0:00"
    hour_end <- as.POSIXct(paste(d, tm), format = "%d-%m-%y %H:%M", tz = "UTC")
    hour_end[end24] <- hour_end[end24] + 86400 # 24:00 -> next day 0:00

    df <- data.frame(station = station,
                     hour = hour_end - 3600, # hour START
                     stringsAsFactors = FALSE)
    samp_flag <- rep(FALSE, length(d))
    for (j in keep) {
      nm <- header_row[[j]]
      v <- body[[j]][ok]
      samp_flag <- samp_flag | (!is.na(v) & v == "Samp<")
      friendly <- arc_var_names[[nm]]
      if (is.null(friendly) || is.na(friendly)) friendly <- make.names(nm)
      df[[friendly]] <- suppressWarnings(as.numeric(v))    # "Samp<" -> NA
    }
    df$flag_samp_below <- samp_flag
    # blocks have different variable sets (e.g. no Rainfall at Henderson):
    # give every block the full set so they can be stacked
    for (v in unname(arc_var_names)) if (!v %in% names(df)) df[[v]] <- NA_real_
    out[[k]] <- df[, c("station", "hour", unname(arc_var_names), "flag_samp_below")]
  }

  res <- do.call(rbind, out)
  res <- res[!is.na(res$hour), ]
  rownames(res) <- NULL
  res
}

#' Read several ARC files (e.g. one per year) and stack them
read_arc_met_many <- function(paths) {
  do.call(rbind, lapply(paths, read_arc_met))
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  path <- if (length(args) >= 1) args[1] else "data/raw/ARC_Met_data_2008_v1.csv"
  arc <- read_arc_met(path)
  str(arc)
  cat("\nRows per station:\n"); print(table(arc$station))
  cat("\nRange of hours:\n"); print(tapply(arc$hour, arc$station, function(x) format(range(x))))
  cat("\nMissing share per variable (by station):\n")
  print(round(sapply(split(arc[, c("outdoor_T", "outdoor_RH", "solar_rad", "wind_speed")], arc$station),
                     function(d) colMeans(is.na(d))), 3))
  cat("\n'Samp<' flags per station:\n"); print(tapply(arc$flag_samp_below, arc$station, sum))
  cat("\nDuplicated station-hour keys:", sum(duplicated(arc[, c("station", "hour")])), "\n")
}
