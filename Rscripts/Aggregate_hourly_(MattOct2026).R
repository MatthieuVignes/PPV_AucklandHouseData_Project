##################################################
# Turns the tidy long output of parse_house_file() into hourly per-house/room/variable 
# summaries (mean, min, max, n_obs), ready to merge with hourly NIWA outdoor data by 
# nearest hour.
#
# Needs to be done before modelling
# 1. Raw 2-min data is heavily autocorrelated -> hourly means (consider window?) 
#    give a much better-behaved unit for the mixed model.
# 2. Collapses the "each block has its own async Time column" issue: once every block 
#    is floored to the hour, all rooms/variables share a common time axis to merge on.
##################################################

aggregate_hourly <- function(meas_df) {
  meas_df$hour <- as.POSIXct(cut(meas_df$timestamp, breaks = "hour"), tz = "Pacific/Auckland")

  agg <- aggregate(
    value ~ house + year + week + state + room + variable + hour,
    data = meas_df,
    FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                         min  = min(x, na.rm = TRUE),
                         max  = max(x, na.rm = TRUE),
                         n    = length(x))
  )
  # aggregate() packs the FUN output into a matrix column - unpack it
  stats <- as.data.frame(agg$value)
  names(stats) <- c("mean", "min", "max", "n_obs")
  out <- cbind(agg[, c("house", "year", "week", "state", "room", "variable", "hour")], stats)
  out$hour_of_day <- as.integer(format(out$hour, "%H"))
  out$date <- as.Date(out$hour)
  out
}

# Widen to one row per house/room/hour with variable columns (T, RH, CO2, CO)
widen_hourly <- function(hourly_long) {
  reshape(
    hourly_long[, c("house", "year", "week", "state", "room", "hour", "hour_of_day", "date", "variable", "mean")],
    idvar = c("house", "year", "week", "state", "room", "hour", "hour_of_day", "date"),
    timevar = "variable",
    direction = "wide"
  )
}
