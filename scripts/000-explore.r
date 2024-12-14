# Notes ----------------------------------------------------------------------------------
#   Goal:   Explore the data
#   Time:   ???

# Data notes -----------------------------------------------------------------------------
#   - I drop all monitors outside of the CONUS.
#   - I drop all 1-hour samples.
#   - The "pollutant standard" is missing for 88502 monitors.

# Setup ----------------------------------------------------------------------------------
  # Load packages
  library(pacman)
  p_load(janitor, parallel, stringr, data.table, collapse, magrittr, lubridate, here)
  # nin function
  library('Hmisc', include.only = '%nin%')

# Load monitor datasets ------------------------------------------------------------------
  # Load monitor and site datasets
  site_dt = here('data', 'raw', 'aqs_sites.zip') |> fread()
  monitor_dt = here('data', 'raw', 'aqs_monitors.zip') |> fread()
  # Clean up names
  setnames(site_dt, site_dt[0, ] |> clean_names() |> names())
  setnames(monitor_dt, monitor_dt[0, ] |> clean_names() |> names())

# Load monitoring data -------------------------------------------------------------------
  # Find sets of files (daily v hourly; 88101 v 88502)
# NOTE Ignoring hourly files at the moment
  d101 = here('data', 'raw', 'daily-88101') |> dir(full.names = TRUE)
  d502 = here('data', 'raw', 'daily-88502') |> dir(full.names = TRUE)
  # Load files
  d_dt =
    mclapply(c(d101, d502), fread, showProgress = FALSE, mc.cores = detectCores()) |>
    rbindlist(use.names = TRUE, fill = TRUE)
  # Clean up names
  setnames(d_dt, d_dt[0, ] |> clean_names() |> names())
  # Drop intermediate objects
  rm(d101, d502)
  invisible(gc())

# Clean data -----------------------------------------------------------------------------
  # Focus on CONUS
  site_dt %<>% fsubset(state_code != 'CC')
  monitor_dt %<>% fsubset(state_code != 'CC')
  d_dt %<>% fsubset(state_code != 'CC')
  site_dt %<>% fsubset(as.integer(state_code) %in% funique(maps::state.fips$fips))
  monitor_dt %<>% fsubset(as.integer(state_code) %in% funique(maps::state.fips$fips))
  d_dt %<>% fsubset(as.integer(state_code) %in% funique(maps::state.fips$fips))
  # Define unique monitors
  d_dt[, `:=`(
    monitor_id = paste0(
      str_pad(state_code, 2, 'left', 0), '-',
      str_pad(county_code, 3, 'left', 0), '-',
      str_pad(site_num, 4, 'left', 0), '-',
      str_pad(parameter_code, 5, 'left', 0), '-',
      str_pad(poc, 2, 'left', 0)
    )
  )]
  # Focus on 24-hour samples
  d_dt %<>% fsubset(sample_duration != '1 HOUR')
  # Grab desired variables
  d_dt %<>%
    fselect(
      monitor_id,
      poc,
      date_local,
      event_type,
      arithmetic_mean,
      method_code
    )
  # Map monitor_id for corrected monitors back to original (POC - 20)
  d_dt[, id_raw := monitor_id]
  d_dt[method_code %in% c(736, 738), `:=`(
    monitor_id = paste0(
      str_sub(monitor_id, 1, -3),
      poc |> subtract(20) |> str_pad(2, 'left', 0)
    )
  )]
  # Drop POC
  d_dt[, poc := NULL]
  # Add indicator for 'corrected' method (i.e., 736 or 738)
  d_dt[, `:=`(
    corrected = fcase(
      method_code %in% c(736, 738), 'adj',
      default = 'unadj'
    ) |> as.factor()
  )]
  # Flag dates with 'exceptional events'
  d_dt[, has_event := 1L * (event_type != 'None')]
  d_dt[, has_event := fmax(has_event), by = .(monitor_id, date_local)]
  # Combine 'None' and 'Included' into 'raw'; 'Excluded' becomes 'excluded'
  d_dt[, `:=`(
    type = fcase(
      event_type %in% c('None', 'Included'), 'raw',
      event_type == 'Excluded', 'excluded',
      default = NA
    ) |> as.factor()
  )]
  # Pivot to wide (a row is a monitor-day)
  dw_dt =
    pivot(
      data = d_dt,
      ids = c('monitor_id', 'date_local', 'has_event'),
      values = c('arithmetic_mean'),
      names = c('type', 'corrected'),
      how = 'wider'
    )
  # Add date components
  dw_dt[, `:=`(
    y = year(date_local),
    q = quarter(date_local)
  )]
  # Start building the design values, combining relevant pieces
  # Type 1: Raw, i.e., does not exclude exceptional events or "correct" FEMs
  dw_dt[, pm_raw := raw_unadj]
  # Type 2: Exclude exceptional events; still uncorrected
  dw_dt[has_event == 0, pm_noevent := raw_unadj]
  dw_dt[has_event == 1, pm_noevent := excluded_unadj]
  # Type 3: Corrected FEMs (does not exclude exceptional events)
  dw_dt[, pm_adj := raw_unadj]
  dw_dt[!is.na(raw_adj), pm_adj := raw_adj]
  # Type 4: Exclude exceptional events and correct FEMs
  dw_dt[has_event == 0, pm_final := raw_unadj]
  dw_dt[has_event == 0 & !is.na(raw_adj), pm_final := raw_adj]
  dw_dt[has_event == 1, pm_final := excluded_unadj]
  dw_dt[has_event == 1 & !is.na(excluded_adj), pm_final := excluded_adj]
  # Drop intermediate columns
  dw_dt[, c('raw_unadj', 'excluded_unadj', 'raw_adj', 'excluded_adj') := NULL]
  # Calculate quarterly averages
  q_dt = dw_dt[, .(
    n = .N,
    n_events = fsum(has_event),
    pm_raw = fmean(pm_raw),
    pm_noevent = fmean(pm_noevent),
    pm_adj = fmean(pm_adj),
    pm_final = fmean(pm_final)
  ), by = .(monitor_id, y, q)]
  # Calculate the annual averages
  a_dt = q_dt[, .(
    n = fsum(n),
    n_events = fsum(n_events),
    pm_raw = fmean(pm_raw),
    pm_noevent = fmean(pm_noevent),
    pm_adj = fmean(pm_adj),
    pm_final = fmean(pm_final)
  ), by = .(monitor_id, y)]
  # Calculate design values (using three-year lags)
  setorder(a_dt, monitor_id, y)
  dv_dt = a_dt[, .(
    monitor_id,
    id_lag1 = shift(monitor_id, 1),
    id_lag2 = shift(monitor_id, 2),
    y,
    n = n + shift(n, 1) + shift(n, 2),
    n_events = n_events + shift(n_events, 1) + shift(n_events, 2),
    pm_raw = (pm_raw + shift(pm_raw, 1) + shift(pm_raw, 2)) / 3,
    pm_noevent = (pm_noevent + shift(pm_noevent, 1) + shift(pm_noevent, 2)) / 3,
    pm_adj = (pm_adj + shift(pm_adj, 1) + shift(pm_adj, 2)) / 3,
    pm_final = (pm_final + shift(pm_final, 1) + shift(pm_final, 2)) / 3
  )]
  # Drop years missing lags (is either NA or does not have matching lags)
  dv_dt %<>% .[!is.na(id_lag1) & !is.na(id_lag2)]
  dv_dt %<>% .[monitor_id == id_lag1 & monitor_id == id_lag2]
  # Drop intermediate columns
  dv_dt[, c('id_lag1', 'id_lag2') := NULL]
  # Drop intermediate objects
  rm(q_dt, a_dt)
  invisible(gc())

# Summarize data -------------------------------------------------------------------------
  # Violations of the new 9-unit standard by year (includes 88502)
  dv_dt[, lapply(.SD, function(x) fmean(x > 9)), by = y, .SDcols = patterns('^pm_.*')]
  # Repeat for only 88101 monitors
  dv_dt %>%
    .[str_detect(monitor_id, '88101'), ] %>%
    .[, lapply(.SD, function(x) fmean(x > 9)), by = y, .SDcols = patterns('^pm_.*')]

