# ==============================================================================
#    NAME: scripts/6_model.R
#   INPUT: Aircraft characteristics & takeoff conditions from the calling script
# ACTIONS: Simulate takeoffs and return the ground distance required in m
#  OUTPUT: A vector of ground distances in m (Inf where the takeoff is infeasible)
# RUNTIME: Variable based on input data table size
#  AUTHOR: Thomas D. Pellegrin <thomas@pellegr.in>
#    YEAR: 2023
#     REV: 2026 — (a) engine flat-rating thrust-temperature lapse (§3.2.1)
#                 (b) explicit infeasibility handling, replaces the acc clamp
#                 (c) runway slope converted to radians
# ==============================================================================

# ==============================================================================
# 1 Define a function to calculate the liftoff speed based on:
# cllof = dimensionless lift coefficient at liftoff
# g     = gravitational acceleration constant in m/s²
# tom   = takeoff mass in kg
# rho   = air density in kg/m³
# s     = wing surface area in m²
# Adapted from Blake (2009).
# ==============================================================================

fn_vlof <- function(DT) {

  vlof <- sqrt(
    (2L * DT[, tom] * sim$g) /
    (DT[, rho] * DT[, s] * DT[, cllof])
  )

} # End of fn_vlof function

# ==============================================================================
# 2 Define a function to calculate the horizontal dist. covered by the aircraft
# during the first-segment climb, in m. Adapted from Filippone (2012, p. 258).
# ==============================================================================

fn_dis_air <- function() {

  dis_air <- sim$scr_hgt / tan(sim$clb_ang * pi / 180L) * sim$ft_to_m

} # End of fn_dis_air function

# ==============================================================================
# 3 Define a function to calculate the ground distance required.
# Returns Inf for any takeoff that cannot accelerate to liftoff (thrust-limited).
# ==============================================================================

fn_dis_gnd <- function(DT) {

  # ============================================================================
  # 3.1 Calculate the airspeed and groundspeed intervals in m/s
  # Groundspeed is airspeed plus headwind.
  # ============================================================================

  vtas <- as.vector(
    mapply(FUN = seq, from = DT[, hdw], to = DT[, vlof],
      length.out = sim$int_stp)
  )

  vgnd <- as.vector(
    mapply(FUN = seq, from = 0L, to = DT[, vlof] - DT[, hdw],
      length.out = sim$int_stp)
  )

  # ============================================================================
  # 3.2 Calculate the propulsive force in N. Adapted from Sun et al. (2020).
  # ============================================================================

  vsnd  <- rep(sqrt(sim$adb_idx * sim$rsp_air * DT[, tas]), each = sim$int_stp)
  vmach <- vtas / vsnd
  dp    <- rep(DT[, ps] / sim$isa_ps, each = sim$int_stp)
  bpr   <- rep(DT[, bpr], each = sim$int_stp)

  g0 <-  .0606 * bpr  + .6337
  y  <- -.4327 * dp^2L + 1.3855 * dp    + .0472
  x  <-  .1377 * dp^3L - .4374  * dp^2L + 1.3003 * dp
  z  <-  .9106 * dp^3L - 1.7736 * dp^2L + 1.8697 * dp

  tr <- y - .377 * (1L + bpr) / sqrt((1L + .82 * bpr) * g0) * z * vmach +
    (.23 + .19 * sqrt(bpr)) * x * vmach^2L

  slst <- rep(DT[, slst], each = sim$int_stp)
  n    <- rep(DT[, n], each = sim$int_stp)
  fmax <- tr * slst * n

  # ============================================================================
  # 3.2.1 Apply the engine flat-rating (thrust-temperature lapse) [REV. 2026]
  # Max takeoff thrust is held to a corner temperature (local ISA + a per-engine
  # deviation, sim$flat_dev) and decays above it. The multiplier is exactly 1 at
  # or below the corner, so the ISA calibration in 7_calibrate.R is undisturbed.
  # For per-engine values, add flat_dev / flat_exp columns to aircraft.csv and
  # replace the sim$ scalars with rep(DT[, flat_dev], each = sim$int_stp) etc.
  # ============================================================================

  tas         <- rep(DT[, tas], each = sim$int_stp)
  t_isa_local <- sim$isa_tas * dp^(sim$rsp_air * sim$isa_lap / sim$g)
  t_corner    <- t_isa_local + sim$flat_dev
  flat_mult   <- data.table::fifelse(
    tas <= t_corner, 1, (t_corner / tas)^sim$flat_exp
  )
  fmax <- fmax * flat_mult

  thr_red <- rep(DT[, thr_red], each = sim$int_stp)
  frto    <- fmax * (100L - thr_red) / 100L

  # ============================================================================
  # 3.3 Calculate the acceleration in m/s² up to liftoff. Adapted from Blake
  # (2009). Runway slope is converted from degrees to radians [REV. 2026].
  # ============================================================================

  tom   <- rep(DT[, tom], each = sim$int_stp)
  w     <- tom * sim$g
  s     <- rep(DT[, s], each = sim$int_stp)
  cllof <- rep(DT[, cllof], each = sim$int_stp)
  cd    <- rep(DT[, cd], each = sim$int_stp)
  rho   <- rep(DT[, rho], each = sim$int_stp)
  q     <- .5 * rho * vtas^2L

  acc <- (sim$g / w) * (frto - (sim$rwy_frc * w) -
    (cd - sim$rwy_frc * cllof) * q * s - (w * sin(sim$rwy_slp * pi / 180L)))

  # ============================================================================
  # 3.3.1 Flag thrust-/energy-limited takeoffs [REV. 2026]
  #
  # If the net accelerating force falls to (or below) zero anywhere before
  # liftoff, the aircraft cannot reach Vlof: the takeoff is thrust-limited and
  # physically infeasible at this mass / derate / condition, irrespective of
  # runway length. These are real outcomes (and become more common once the
  # flat-rating in §3.2.1 cuts hot-day thrust), so we flag them and return Inf
  # for the affected takeoffs.
  #
  # This replaces the previous `acc[acc < 0] <- 1e-3` clamp, which fabricated a
  # finite-but-meaningless distance, conflated thrust-limited with field-length-
  # limited cases, and could write a garbage TODR to the results for the most
  # extreme (and most policy-relevant) hot-and-high observations.
  #
  # The integrand v/a legitimately diverges as a -> 0, so near-equilibrium
  # *feasible* cases yield very large finite TODR; the calling script caps these
  # at sim$todr_cap and labels them field-length-infeasible (a distinct mode).
  # ============================================================================

  # Per-takeoff minimum acceleration (column-wise min over each int_stp block)
  acc_min  <- do.call(pmin, asplit(matrix(acc, nrow = sim$int_stp), 1L))
  feasible <- acc_min > sim$acc_eps

  # ============================================================================
  # 3.4 Increment the horizontal takeoff distances in m. Adapted from Blake
  # (2009). Adaptive window widths reset per takeoff, so an Inf/NaN arising in an
  # infeasible block cannot leak into a neighbouring (feasible) takeoff.
  # ============================================================================

  bar_width <- rep(
    x = c(seq.int(2L), rep(x = 2L, each = sim$int_stp - 2L)),
    times = nrow(DT)
  )

  acc_bar  <- frollmean(x = acc,  n = bar_width, adaptive = TRUE)
  vgnd_bar <- frollmean(x = vgnd, n = bar_width, adaptive = TRUE)

  vlof     <- rep(DT[, vlof], each = sim$int_stp)
  hdw      <- rep(DT[, hdw],  each = sim$int_stp)
  vgnd_int <- (vlof - hdw) / (sim$int_stp - 1L)

  inc <- vgnd_bar * vgnd_int / acc_bar
  cum <- frollsum(
    x = inc, n = rep(x = seq(1L:sim$int_stp), times = nrow(DT)), adaptive = TRUE
  )

  # ============================================================================
  # 3.5 Assemble the ground distance in m; set infeasible takeoffs to Inf.
  # Adapted from Blake (2009) and Gratton et al (2020).
  # ============================================================================

  dis_gnd <- cum[seq(sim$int_stp, length(cum), sim$int_stp)]
  dis_gnd[!feasible] <- Inf
  return(dis_gnd)

} # End of fn_dis_gnd function

# EOF
