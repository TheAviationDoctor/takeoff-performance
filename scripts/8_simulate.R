# ==============================================================================
#    NAME: scripts/8_simulate.R  (refactored)
#   INPUT: Climatic observations as a Parquet dataset (dir$cli_pq), one file per
#          airport, produced by 5_transform.R; calibration from cal.parquet
#          (7_calibrate.R).
# ACTIONS: For each airport, solve the operating point of each takeoff by
#          VECTORISED BISECTION instead of a one-passenger-at-a-time loop:
#            - if the aircraft fits at MTOM on TOGA  -> find the MAX thrust
#              reduction (derate) that still fits     (tom_rem = 0)
#            - else                                   -> find the MAX mass that
#              fits on TOGA, i.e. the payload to shed  (thr_red = 0)
#          Outcomes are labelled feasible / field-length-limited / thrust-
#          limited (see 6_model.R §3.3.1).
#  OUTPUT: One Parquet file of takeoff outcomes per airport in dir$tko_pq.
# RUNTIME: ~1–2 orders of magnitude faster than the former MySQL/loop version.
#  AUTHOR: Thomas D. Pellegrin <thomas@pellegr.in>
#     REV: 2026 — Parquet + arrow I/O; vectorised monotone solver.
# ==============================================================================

# ==============================================================================
# 0 Housekeeping
# ==============================================================================

rm(list = ls())

library(arrow)        # Parquet datasets (multi-threaded I/O, predicate pushdown)
library(data.table)   # In-memory vectorised compute
library(dplyr)        # open_dataset() verbs for the per-airport filter

source("scripts/0_common.R")  # adds sim$flat_*, sim$isa_lap, sim$acc_eps,
source("scripts/6_model.R")   # sim$todr_cap, dir$cli_pq, dir$tko_pq

start_time <- Sys.time()

# Cap arrow's thread pool so it composes with any outer parallelism
arrow::set_cpu_count(parallel::detectCores())

# ==============================================================================
# 1 Load the small, shared inputs once
# ==============================================================================

# 1.1 Aircraft characteristics; MTOM is the starting mass
dt_act <- fread(
  fls$act,
  colClasses = c(rep("factor", 2L), rep("integer", 5L), rep("numeric", 5L))
)[type %in% act, .(type, n, slst, bpr, s, tom_mtom = tom_max)]

# 1.2 Calibration (CLlof, CD by integer kg) from cal.parquet (7_calibrate.R).
#     Key by (type, tom) for fast joins.
dt_cal <- setDT(read_parquet(file.path(dir$cal, "cal.parquet")))
dt_cal[, type := as.factor(type)]
setkey(dt_cal, type, tom)
dt_cal[, tom_belf := min(tom), by = type]   # economic floor = break-even mass

# ==============================================================================
# 2 Helper functions
# ==============================================================================

# 2.1 Total regulatory TODR for the current state of DT (Inf if infeasible).
fn_todr <- function(DT) {
  set(DT, j = "vlof", value = fn_vlof(DT))
  (fn_dis_gnd(DT) + fn_dis_air()) * sim$tod_mul
}

# 2.2 Refresh mass-dependent aerodynamics after a trial mass changes.
fn_lookup_cal <- function(DT) {
  set(DT, j = "tom", value = as.integer(round(DT[["tom"]])))
  DT[dt_cal, `:=`(cllof = i.cllof, cd = i.cd), on = c("type", "tom")]
  invisible(DT)
}

# 2.3 Vectorised monotone solver. TODR is increasing in both `tom` and
#     `thr_red`, so for each row we seek the LARGEST value of `var` in [lo, hi]
#     for which TODR <= toda. Each row keeps its own bracket; ~n_iter passes.
fn_solve_max <- function(DT, var, lo, hi, toda, n_iter,
                         as_integer = FALSE, refresh = NULL) {
  lo <- as.numeric(rep(lo, length.out = nrow(DT)))
  hi <- as.numeric(rep(hi, length.out = nrow(DT)))
  for (k in seq_len(n_iter)) {
    mid <- (lo + hi) / 2
    if (as_integer) mid <- floor(mid)
    set(DT, j = var, value = mid)
    if (!is.null(refresh)) refresh(DT)
    feasible <- fn_todr(DT) <= toda          # Inf <= toda is FALSE, as intended
    lo <- fifelse(feasible, mid, lo)         # fits  -> answer >= mid, raise floor
    hi <- fifelse(feasible, hi,  mid)        # fails -> answer <  mid, lower ceiling
  }
  if (as_integer) floor(lo) else lo
}

# ==============================================================================
# 3 Per-airport simulation
# ==============================================================================

cli_ds <- open_dataset(dir$cli_pq)           # the whole climate dataset (lazy)
airports <- collect(distinct(cli_ds, icao))[["icao"]]

# Ensure the output directory exists before writing into it
dir.create(dir$tko_pq, showWarnings = FALSE, recursive = TRUE)

# Resume support: skip airports already written
done <- sub("\\.parquet$", "", list.files(dir$tko_pq, pattern = "\\.parquet$"))
airports <- setdiff(airports, done)

fn_simulate <- function(this_icao) {

  # 3.1 Read this airport's observations (corrected ps & rho from 5_transform)
  dt_cli <- cli_ds |>
    filter(icao == this_icao) |>
    select(year, obs, icao, lat, lon, zone, ssp, huss, ps, tas,
           rho, hdw, rwy, toda) |>
    collect() |>
    setDT()
  if (!nrow(dt_cli)) return(invisible(NULL))

  # 3.2 Cross with aircraft types and attach MTOM calibration
  dt <- dt_cli[, as.list(dt_act), by = dt_cli]
  set(dt, j = "type", value = as.factor(dt[["type"]]))
  set(dt, j = "tom",  value = dt[["tom_mtom"]])
  fn_lookup_cal(dt)

  # 3.3 Regime split on the MTOM / TOGA takeoff
  set(dt, j = "thr_red", value = 0L)
  todr0 <- fn_todr(dt)
  fits  <- todr0 <= dt[["toda"]]

  dtA <- dt[fits]                             # derate regime  (tom_rem = 0)
  dtB <- dt[!fits]                            # offload regime (thr_red = 0)

  # 3.4 Derate regime: max thrust reduction that still fits, at MTOM
  if (nrow(dtA)) {
    set(dtA, j = "tom", value = dtA[["tom_mtom"]]); fn_lookup_cal(dtA)
    thr_star <- fn_solve_max(dtA, "thr_red", 0L, sim$thr_ini,
                             toda = dtA[["toda"]], n_iter = 6L, as_integer = TRUE)
    set(dtA, j = "thr_red", value = pmin(thr_star, sim$thr_ini))
    set(dtA, j = "tom_rem", value = 0L)
    set(dtA, j = "todr",    value = ceiling(fn_todr(dtA)))
    set(dtA, j = "outcome", value = "feasible")
  }

  # 3.5 Offload regime: max mass that fits on TOGA; payload shed = MTOM - mass
  if (nrow(dtB)) {
    set(dtB, j = "thr_red", value = 0L)
    mass_star <- fn_solve_max(dtB, "tom", dtB[["tom_belf"]], dtB[["tom_mtom"]],
                              toda = dtB[["toda"]], n_iter = 12L,
                              refresh = fn_lookup_cal)
    set(dtB, j = "tom", value = as.integer(floor(mass_star))); fn_lookup_cal(dtB)
    todrB <- fn_todr(dtB)
    set(dtB, j = "tom_rem", value = dtB[["tom_mtom"]] - dtB[["tom"]])
    set(dtB, j = "todr",    value = ifelse(is.finite(todrB), ceiling(todrB), NA_integer_))
    set(dtB, j = "outcome", value = fifelse(
      todrB <= dtB[["toda"]],          "feasible",
      fifelse(is.finite(todrB),        "infeasible_field",
                                       "infeasible_thrust")))
    set(dtB, i = dtB[, .I[todr > sim$todr_cap & outcome == "feasible"]],
        j = "outcome", value = "infeasible_field")
  }

  out <- rbind(dtA, dtB, fill = TRUE)

  # 3.6 Write one Parquet file for this airport (partition-friendly, resumable)
  cols <- c("year","obs","icao","lat","lon","zone","ssp","type","huss","ps",
            "tas","rho","hdw","rwy","toda","todr","vlof","thr_red","tom_rem",
            "outcome")
  write_parquet(out[, ..cols],
                file.path(dir$tko_pq, paste0(this_icao, ".parquet")),
                compression = "zstd")
  invisible(NULL)
}

# ==============================================================================
# 4 Run. arrow is internally threaded, so a sequential airport loop already uses
#    all cores for I/O + compute. To add coarse parallelism, wrap in
#    parallel::mclapply (fork) — but give each worker its own arrow thread budget.
# ==============================================================================

invisible(lapply(airports, fn_simulate))

# ==============================================================================
# 5 Housekeeping
# ==============================================================================

Sys.time() - start_time

# EOF
