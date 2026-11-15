##################################################
# Parses a raw "House_<id>_<state>YY.csv" measurement file into tidy long-format data frames.
#
# File structure note: the file is a set of column BLOCKS glued side by side. Each
# block starts with its own "TimeXXX-<code>" column followed by that block's data 
# columns, up to the next "TimeXXX" column.
# Row i in block A is NOT the same clock time as row i in block B; each block is 
# its own independently-logged time series.
# So we parse block-by-block and only merge across blocks by actual timestamp later 
# (in the aggregation/merge step), never by row position.
#
# Column suffix code, e.g. "CL18":
#   [1] state : C = Control, A = Active
#   [2] room  : L = Living room, B = Bedroom, R = Roof cavity
#   [3] week  : 1 or 2 (week of that 2-week deployment)
#   [4] year  : 8 -> 2008, 9 -> 2009
#
# Block "families" identified by the Time-column prefix:
#   TimeBW    -> indoor logger (LR or BR): CO2, CO, T, RH
#   TimeEA    -> roof cavity logger: CO2, T, RH  (no CO in RC)
#   TimePPM   -> HCHO (formaldehyde) logger, paired with a BW block
#   TimeEXT   -> external/reference channel -> IGNORED per spec
#   TimePULSE -> PPV fan pulse counter (house-level, not room-level)
##################################################

decode_code <- function(code) {
  code <- trimws(code)
  # keep only the first 4 meaningful chars (guards against stray spaces)
  code <- substr(code, 1, 4)
  list(
    state = substr(code, 1, 1), # "C" / "A"
    room  = substr(code, 2, 2), # "L" / "B" / "R"
    week  = as.integer(substr(code, 3, 3)), # 1 / 2
    year  = 2000L + as.integer(paste0("0", substr(code, 4, 4))) + 2000L - 2000L, # placeholder, fixed below
    raw_code = code
  )
}

# proper year decode: '8' -> 2008, '9' -> 2009
decode_year <- function(yr_char) {
  ifelse(yr_char == "8", 2008L,
  ifelse(yr_char == "9", 2009L, NA_integer_))
}

state_label <- function(s) ifelse(s == "C", "Control", ifelse(s == "A", "Active", NA_character_))
room_label  <- function(r) c(L = "LivingRoom", B = "Bedroom", R = "RoofCavity")[r]

#' Parse one raw house-year measurement CSV into tidy long data
#'
#' @param path path to the raw CSV
#' @param house_id integer house ID (not always recoverable from the file
#'   itself, so pass it in explicitly from the filename / design table)
#' @return list(measurements = <long df>, pulse = <long df>)
parse_house_file <- function(path, house_id) {

  raw <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                   colClasses = "character")
  col_names <- trimws(names(raw))
  names(raw) <- col_names

  is_time_col <- grepl("^Time", col_names)
  if (!is_time_col[1]) {
    stop("Expected file to start with a Time* column - check format for: ", path)
  }
  block_id <- cumsum(is_time_col)

  measurements <- list()
  pulses <- list()

  for (b in unique(block_id)) {
    cols_b <- col_names[block_id == b]
    time_col <- cols_b[1]
    data_cols <- cols_b[-1]
    if (length(data_cols) == 0) next

    prefix <- sub("^Time", "", time_col) # BW / EA / PPM / EXT / PULSE
    prefix <- sub("-.*$", "", prefix) # strip suffix if glued
    code <- sub("^[A-Za-z]+-", "", time_col) # e.g. "CL18"
    code <- trimws(code)

    if (prefix == "EXT") next # explicitly ignore EXT channels

    ts <- suppressWarnings(as.POSIXct(raw[[time_col]], format = "%Y-%m-%d %H:%M:%S", tz = "Pacific/Auckland"))

    if (prefix == "PULSE") {
      val_col <- data_cols[grepl("^PULSE", data_cols)]
      if (length(val_col) == 0) next
      # PULSE time column observed in a different format: "DD/MM/YYYY HH:MM"
      ts_pulse <- suppressWarnings(as.POSIXct(raw[[time_col]], format = "%d/%m/%Y %H:%M", tz = "Pacific/Auckland"))
      df <- data.frame(
        house    = house_id,
        code     = code,
        timestamp = ts_pulse,
        pulse_count = suppressWarnings(as.numeric(raw[[val_col[1]]])),
        stringsAsFactors = FALSE
      )
      pulses[[length(pulses) + 1]] <- df
      next
    }

    # PPM/HCHO block: usually paired alongside a BW block sharing the same
    # code but its own Time column. Treated as its own mini long block
    dcode <- decode_code(code)
    st  <- state_label(dcode$state)
    rm_ <- room_label(dcode$room)
    wk  <- dcode$week
    yr  <- decode_year(substr(code, 4, 4))

    for (v in data_cols) {
      var_name <- trimws(v)
      var_name <- sub(paste0("-", code, "\\s*$"), "", var_name)
      var_name <- trimws(var_name)
      df <- data.frame(
        house     = house_id,
        year      = yr,
        week      = wk,
        state     = st,
        room      = rm_,
        code      = code,
        instrument = prefix,
        variable  = var_name,
        timestamp = ts,
        value     = suppressWarnings(as.numeric(raw[[v]])),
        stringsAsFactors = FALSE
      )
      measurements[[length(measurements) + 1]] <- df
    }
  }

  meas_df  <- do.call(rbind, measurements)
  pulse_df <- if (length(pulses) > 0) do.call(rbind, pulses) else
    data.frame(house = integer(0), code = character(0),
               timestamp = as.POSIXct(character(0)), pulse_count = numeric(0))

  rownames(meas_df) <- NULL
  rownames(pulse_df) <- NULL

  list(measurements = meas_df, pulse = pulse_df)
}

##################################################
# Quick self-test / demo when run directly
##################################################
if (sys.nframe() == 0) {
  res <- parse_house_file("/mnt/user-data/uploads/House_1_active08.csv", house_id = 1)

  cat("=== measurements: structure ===\n")
  str(res$measurements)

  cat("\n=== unique (room, state, week, year, instrument, variable) combos ===\n")
  print(unique(res$measurements[, c("room", "state", "week", "year", "instrument", "variable")]))

  cat("\n=== sample rows ===\n")
  print(head(res$measurements[order(res$measurements$room, res$measurements$timestamp), ], 12))

  cat("\n=== pulse data ===\n")
  print(res$pulse)

  cat("\n=== timestamp ranges per block (confirms async logging) ===\n")
  agg <- aggregate(timestamp ~ room + state + week + year + instrument, data = res$measurements,
                    FUN = function(x) paste(min(x), "to", max(x)))
  print(agg)
}
