##################################################
# Paper (a) first analyses: how do outdoor conditions drive the roof cavity (RC)
# temperature, and when does the RC reach the 18 C "available energy" threshold that 
# lets a PPV system deliver useful (warm) air?
#
# Models (rstanarm / Stan; one row = one house-hour of RC temperature)
# M1  linear: RC_T ~ outdoor_T + solar + diurnal harmonics + (1 + outdoor_T | house)
# M2  logistic: 1[RC_T >= 18] ~ outdoor_T + solar + (1 | house)
# M3  sensitivity: M1 + PPV-active indicator and its interaction with outdoor_T 
# (does running PPV change the RC-outdoor relationship we are trying to characterise?)
#
# Questions answered. As posterior statements
# Q4 - how many hours a day does the RC reach 18 C (overall and by hour)
# Q5 - what outdoor temperature is needed, at a given solar level, for a 50%/80% chance 
#      that the RC reaches 18 C
#
# Usage:
# Rscript R/paper_a_model.R simulated <out_dir> (simulating weather data)
# Rscript R/paper_a_model.R real <out_dir> <raw_dir> <design.csv> \
#           <station_lookup.csv> <ARC_2008.csv> [<ARC_2009.csv>]
# Optional environment variables: PPV_ITER (default 2000), PPV_CHAINS (4)
#
# Limitations (revisit later)
# - Hourly values are autocorrelated; house random effects do not remove that. 
#   Next step: brms with an AR(1) residual and a cyclic smooth of hour of day 
#   (not available in rstanarm).
# - "Typical house" curves (Q5) set the house random effect to zero
# - PPV state is partly confounded with calendar time (2008 week 2 and 2009 are 
#   later in winter); M3 conditions on weather, and the control houses are the 
#   reference, but no explicit time term is included yet.
##################################################

library(rstanarm)

get_script_dir <- function() {
  f <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(f)) dirname(sub("^--file=", "", f[1])) else "R"
}

source_helpers <- function(dir) {
  for (f in c("parse_house_file.R", "aggregate_hourly.R", "read_arc_met.R",
              "read_house_design.R", "merge_house_outdoor.R")) {
    source(file.path(dir, f))
  }
}

##################################################
# 1. Data
##################################################

#' Build the merged hourly table from every raw house file in raw_dir.
#' File names expected like "House 12 active08.csv" or "House_12_control09.csv".
load_merged_hourly <- function(raw_dir, design_path, station_path, arc_paths) {
  files <- list.files(raw_dir, pattern = "^House[ _][0-9]+.*[0-9]{2}\\.csv$",
                      full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) stop("No House*.csv files found in ", raw_dir)
  design <- read_house_design(design_path, station_path)
  arc    <- read_arc_met_many(arc_paths)

  pieces <- list(); failed <- character(0)
  for (f in files) {
    id <- as.integer(sub("^House[ _]([0-9]+).*$", "\\1", basename(f), ignore.case = TRUE))
    res <- tryCatch({
      p <- suppressWarnings(parse_house_file(f, house_id = id))
      aggregate_hourly(p$measurements)
    }, error = function(e) { failed <<- c(failed, paste0(basename(f), ": ", conditionMessage(e))); NULL })
    if (!is.null(res)) pieces[[length(pieces) + 1]] <- res
  }
  if (length(failed)) { cat("Files that failed to parse:\n"); print(failed) }
  hourly <- do.call(rbind, pieces)
  cat("Parsed", length(pieces), "of", length(files), "files\n")
  m <- merge_outdoor(hourly, design, arc)
  check_merge(m)
  m
}

#' Reduce a merged table to RC temperature rows ready for modelling.
#' min_obs: a house-hour needs at least this many raw readings (2-min
#' logging gives ~30 per hour), so part-hours at set-up/removal are dropped.
prepare_rc_data <- function(merged, min_obs = 20) {
  d <- merged[merged$room == "RoofCavity" & merged$variable == "T", ]
  d$RC_T        <- d$mean
  d$hour_of_day <- as.integer(format(d$hour_clock, "%H"))
  d$ppv_active  <- as.integer(d$state == "Active")
  d$solar_100   <- d$solar_rad / 100                      # units of 100 W/m^2
  d$hod_sin     <- sin(2 * pi * d$hour_of_day / 24)
  d$hod_cos     <- cos(2 * pi * d$hour_of_day / 24)
  n0 <- nrow(d)
  d <- d[d$n_obs >= min_obs, ]
  d <- d[complete.cases(d[, c("RC_T", "outdoor_T", "solar_100", "house")]), ]
  d$house   <- factor(d$house)
  d$reach18 <- as.integer(d$RC_T >= 18)
  cat(sprintf("RC data: %d house-hours kept of %d (%d houses; %d with PPV active)\n",
              nrow(d), n0, nlevels(d$house), sum(d$ppv_active)))
  rownames(d) <- NULL
  d
}

#' Simulated stand-in with the SAME columns as the merged table, so the
#' whole modelling code can be tested before all real files are loaded.
#' Truth: RC_T = a_h + b_h * outdoor_T + 1.6 * solar_100 - 0.8 * ppv_active + AR(1) noise
simulate_merged <- function(n_houses = 30, days_per_house = 14, seed = 1, ppv_effect = -0.8) {
  set.seed(seed)
  controls <- c(5, 7, 12, 16, 17, 21, 23, 25, 28, 30)
  n_days <- 75
  hours <- seq(as.POSIXct("2008-07-28 00:00:00", tz = "UTC"), by = "hour", length.out = n_days * 24)
  hod <- as.integer(format(hours, "%H")); day <- rep(seq_len(n_days), each = 24)
  cloud <- runif(n_days, 0.15, 1)
  solar <- pmax(0, 450 * sin(pi * (hod - 7.5) / 10)) * (hod >= 8 & hod <= 17) * cloud[day] + rnorm(length(hod), 0, 8)
  solar <- pmax(0, solar)
  outT  <- 11 + 3.5 * sin(2 * pi * (hod - 9) / 24) + 2.5 * sin(2 * pi * day / 6.5) + rnorm(length(hod), 0, 0.5)
  outRH <- pmin(100, pmax(40, 80 - 2 * (outT - 11) + rnorm(length(hod), 0, 4)))

  a_h <- rnorm(n_houses, 1.5, 1.0); b_h <- rnorm(n_houses, 0.85, 0.07)
  out <- lapply(seq_len(n_houses), function(h) {
    start <- ((h - 1) %/% 5) * 10 * 24 + 1                 # 6 rounds of 5 houses
    idx <- start:(start + days_per_house * 24 - 1)
    week <- ifelse(seq_along(idx) <= 7 * 24, 1L, 2L)
    group <- if (h %in% controls) "Control" else "Active"
    state <- ifelse(group == "Active" & week == 2L, "Active", "Control")
    e <- as.numeric(stats::filter(rnorm(length(idx), 0, 0.6), 0.7, method = "recursive"))
    rc <- a_h[h] + b_h[h] * outT[idx] + 1.6 * solar[idx] / 100 +
      ppv_effect * (state == "Active") + e
    data.frame(house = h, year = 2008L, week = week, state = state, room = "RoofCavity",
               variable = "T", hour_clock = hours[idx], mean = rc, n_obs = 30L,
               outdoor_T = outT[idx], outdoor_RH = outRH[idx], solar_rad = solar[idx],
               group = group, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

##################################################
# 2. Models
##################################################
fit_models <- function(d, iter = 2000, chains = 4, seed = 1) {
  common <- list(chains = chains, iter = iter, cores = min(chains, 2), seed = seed, refresh = 0)

  cat("\nFitting M1 (linear, RC temperature)...\n")
  m1 <- do.call(stan_lmer, c(list(
    RC_T ~ outdoor_T + solar_100 + hod_sin + hod_cos + (1 + outdoor_T | house),
    data = d, prior = normal(0, 2.5, autoscale = TRUE), prior_intercept = normal(10, 5)), common))

  cat("Fitting M2 (logistic, RC >= 18 C)...\n")
  m2 <- do.call(stan_glmer, c(list(
    reach18 ~ outdoor_T + solar_100 + (1 | house),
    data = d, family = binomial(),
    prior = normal(0, 2.5, autoscale = TRUE), prior_intercept = normal(0, 3)), common))

  m3 <- NULL
  if (length(unique(d$ppv_active)) == 2) {
    cat("Fitting M3 (sensitivity: PPV state)...\n")
    m3 <- do.call(stan_lmer, c(list(
      RC_T ~ outdoor_T + solar_100 + hod_sin + hod_cos + ppv_active + ppv_active:outdoor_T +
        (1 + outdoor_T | house),
      data = d, prior = normal(0, 2.5, autoscale = TRUE), prior_intercept = normal(10, 5)), common))
  } else cat("Skipping M3: data contain only one PPV state.\n")

  list(m1 = m1, m2 = m2, m3 = m3)
}

diag_line <- function(fit, label) {
  s <- summary(fit)
  cat(sprintf("%s: max Rhat = %.3f, min n_eff = %d\n", label,
              max(s[, "Rhat"], na.rm = TRUE), as.integer(min(s[, "n_eff"], na.rm = TRUE))))
}

##################################################
# 3. Posterior questions
##################################################
q4_hours_above_18 <- function(m2, d) {
  pp <- posterior_epred(m2, newdata = d, re.form = NULL)        # draws x rows
  overall <- 24 * rowMeans(pp)
  cat(sprintf("\nQ4  Expected hours/day with RC >= 18 C (monitored winter days): %.1f  (95%% CrI %.1f to %.1f)\n",
              mean(overall), quantile(overall, 0.025), quantile(overall, 0.975)))
  cat(sprintf("    Observed: %.1f hours/day\n", 24 * mean(d$reach18)))
  prof <- do.call(rbind, lapply(0:23, function(h) {
    sel <- d$hour_of_day == h
    if (!any(sel)) return(NULL)
    dr <- rowMeans(pp[, sel, drop = FALSE])
    data.frame(hour = h, observed = mean(d$reach18[sel]), model_mean = mean(dr),
               lo95 = unname(quantile(dr, 0.025)), hi95 = unname(quantile(dr, 0.975)), n = sum(sel))
  }))
  cat("\n    By hour of day (observed vs model; a poor fit here flags a missing diurnal/lag term):\n")
  print(round(prof, 3), row.names = FALSE)
  list(overall = overall, profile = prof)
}

q5_thresholds <- function(m2, solar_levels = 0:5, probs = c(0.5, 0.8)) {
  dr <- as.matrix(m2, pars = c("(Intercept)", "outdoor_T", "solar_100"))
  a <- dr[, 1]; b <- dr[, 2]; cs <- dr[, 3]
  cat(sprintf("\nQ5  P(solar effect > 0) = %.3f;  P(outdoor T effect > 0) = %.3f\n", mean(cs > 0), mean(b > 0)))
  res <- do.call(rbind, lapply(solar_levels, function(S) {
    do.call(rbind, lapply(probs, function(p) {
      Tp <- (qlogis(p) - a - cs * S) / b
      data.frame(solar_W_m2 = 100 * S, prob = p, outdoor_T_mean = mean(Tp),
                 lo95 = unname(quantile(Tp, 0.025)), hi95 = unname(quantile(Tp, 0.975)))
    }))
  }))
  cat("    Outdoor temperature (C) needed for the given chance that RC >= 18 C, typical house:\n")
  print(round(res, 2), row.names = FALSE)
  res
}

q_ppv_sensitivity <- function(m3) {
  if (is.null(m3)) return(invisible(NULL))
  dr <- as.matrix(m3, pars = c("ppv_active", "outdoor_T:ppv_active"))
  cat("\nSensitivity: does PPV being active change the RC temperature at the same weather?\n")
  cat(sprintf("    shift in RC temp when PPV active (at outdoor T = 0): mean %.2f C (95%% CrI %.2f to %.2f); P(shift < 0) = %.3f\n",
              mean(dr[, 1]), quantile(dr[, 1], 0.025), quantile(dr[, 1], 0.975), mean(dr[, 1] < 0)))
  cat(sprintf("    change in the outdoor-T slope when PPV active: mean %.3f (95%% CrI %.3f to %.3f)\n",
              mean(dr[, 2]), quantile(dr[, 2], 0.025), quantile(dr[, 2], 0.975)))
  ref <- 12  # shift at a typical winter outdoor temperature
  at12 <- dr[, 1] + dr[, 2] * ref
  cat(sprintf("    net shift at outdoor T = %d C: mean %.2f C (95%% CrI %.2f to %.2f); P(shift < 0) = %.3f\n",
              ref, mean(at12), quantile(at12, 0.025), quantile(at12, 0.975), mean(at12 < 0)))
  invisible(dr)
}

##################################################
# 4. Output
##################################################
save_outputs <- function(out_dir, q4, q5, m2, d) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(q4$profile, file.path(out_dir, "q4_hourly_profile.csv"), row.names = FALSE)
  write.csv(q5, file.path(out_dir, "q5_outdoor_conditions_for_18C.csv"), row.names = FALSE)

  ok <- tryCatch({
    png(file.path(out_dir, "paper_a_q4_q5.png"), width = 1500, height = 650, res = 140, type = "cairo")
    par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
    p <- q4$profile
    plot(p$hour, p$observed, type = "n", ylim = c(0, max(p$hi95, p$observed) * 1.05),
         xlab = "Hour of day", ylab = "P(RC temperature >= 18 C)", main = "Q4: share of hours at 18 C or more")
    polygon(c(p$hour, rev(p$hour)), c(p$lo95, rev(p$hi95)), col = "grey85", border = NA)
    lines(p$hour, p$model_mean, lwd = 2); points(p$hour, p$observed, pch = 19, col = "firebrick")
    legend("topleft", c("model (95% CrI)", "observed"), lwd = c(2, NA), pch = c(NA, 19),
           col = c("black", "firebrick"), bty = "n", cex = 0.8)
    dr <- as.matrix(m2, pars = c("(Intercept)", "outdoor_T", "solar_100"))
    Tg <- seq(0, 20, 0.5); cols <- c("steelblue", "darkorange", "firebrick")
    plot(NA, xlim = c(0, 20), ylim = c(0, 1), xlab = "Outdoor temperature (C)",
         ylab = "P(RC temperature >= 18 C)", main = "Q5: needed outdoor conditions")
    for (i in seq_along(c(0, 2, 4))) {
      S <- c(0, 2, 4)[i]
      lines(Tg, sapply(Tg, function(t) mean(plogis(dr[, 1] + dr[, 2] * t + dr[, 3] * S))), col = cols[i], lwd = 2)
    }
    abline(h = c(0.5, 0.8), lty = 3, col = "grey50")
    legend("topleft", c("no sun", "200 W/m2", "400 W/m2"), col = cols, lwd = 2, bty = "n", cex = 0.8)
    dev.off(); TRUE
  }, error = function(e) { try(dev.off(), silent = TRUE); message("Figure skipped: ", conditionMessage(e)); FALSE })
  cat("\nSaved results to", out_dir, if (ok) "(including paper_a_q4_q5.png)" else "(figure skipped)", "\n")
}

run_paper_a <- function(d, out_dir, iter = 2000, chains = 4) {
  fits <- fit_models(d, iter = iter, chains = chains)

  cat("\n=== Sampler diagnostics ===\n")
  diag_line(fits$m1, "M1"); diag_line(fits$m2, "M2"); if (!is.null(fits$m3)) diag_line(fits$m3, "M3")

  cat("\n=== M1: RC temperature ~ outdoor T + solar + hour-of-day ===\n")
  print(summary(fits$m1, pars = c("(Intercept)", "outdoor_T", "solar_100", "hod_sin", "hod_cos", "sigma"),
                probs = c(0.025, 0.5, 0.975), digits = 2))
  cat("\n=== M2: P(RC >= 18 C) ~ outdoor T + solar ===\n")
  print(summary(fits$m2, pars = c("(Intercept)", "outdoor_T", "solar_100"),
                probs = c(0.025, 0.5, 0.975), digits = 2))

  q4 <- q4_hours_above_18(fits$m2, d)
  q5 <- q5_thresholds(fits$m2)
  q_ppv_sensitivity(fits$m3)

  cat("\nDescriptive check, mean outdoor conditions in hours when RC >= 18 C: ",
      sprintf("outdoor T %.1f C, solar %.0f W/m2 (n = %d hours)\n",
              mean(d$outdoor_T[d$reach18 == 1]), mean(d$solar_rad[d$reach18 == 1]), sum(d$reach18)))
  save_outputs(out_dir, q4, q5, fits$m2, d)
  invisible(list(fits = fits, q4 = q4, q5 = q5))
}

##################################################
# MAIN
##################################################
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) stop("Usage: Rscript paper_a_model.R simulated <out_dir>  |  real <out_dir> <raw_dir> <design.csv> <station.csv> <ARC.csv> [...]")
  mode <- args[1]; out_dir <- args[2]
  source_helpers(get_script_dir())
  iter   <- as.integer(Sys.getenv("PPV_ITER", "2000"))
  chains <- as.integer(Sys.getenv("PPV_CHAINS", "4"))

  if (mode == "simulated") {
    cat("*** SIMULATED DATA: results below are NOT findings, they only test the pipeline. ***\n")
    cat("    (simulated truth: slope on outdoor T ~0.85, +1.6 C per 100 W/m2 solar, PPV shift -0.8 C)\n\n")
    merged <- simulate_merged()
  } else if (mode == "real") {
    if (length(args) < 6) stop("real mode needs: <raw_dir> <design.csv> <station.csv> <ARC.csv> [...]")
    merged <- load_merged_hourly(args[3], args[4], args[5], args[6:length(args)])
  } else stop("mode must be 'simulated' or 'real'")

  d <- prepare_rc_data(merged)
  run_paper_a(d, out_dir, iter = iter, chains = chains)
}
