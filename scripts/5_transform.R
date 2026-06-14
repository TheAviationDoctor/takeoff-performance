# ==============================================================================
#    NAME: scripts/5_transform.R
#   INPUT: ~2.2e9 long climate observations read from the dir$imp_pq Parquet set
# ACTIONS: Pivot the data
#          Calculate the air density, wind vector, and active runway
#  OUTPUT: ~4.4e8 wide climate observations as one Parquet file per airport in
#          dir$cli_pq (consumed by 8_simulate.R and 9_analyze.R)
# RUNTIME: ~4.25 hours (3.8 GHz CPU / 128 GB DDR4 RAM / SSD)
#  AUTHOR: Thomas D. Pellegrin <thomas@pellegr.in>
#    YEAR: 2023
# ==============================================================================

# ==============================================================================
# 0 Housekeeping
# ==============================================================================

# Clear the environment
rm(list = ls())

# Load the required libraries
library(arrow)
library(data.table)
library(parallel)
library(stringr)

# Import the common settings
source("scripts/0_common.R")

# Start a script timer
start_time <- Sys.time()

# Clear the console
cat("\014")

# Set the number of CPU cores for parallel processing (leave one core free)
crs <- parallel::detectCores() - 1L

# ==============================================================================
# 1 Fetch the data that we need
# ==============================================================================

# Fetch the list of airports and runways in the sample from the population
# Parquet (orog was filled in by 4_import.R)
dt_smp <- setDT(arrow::read_parquet(fls$pop))[
  traffic > sim$pop_thr, .(icao, lat, lon, zone, elev, orog, rwy, toda)
]

# Recast column types
set(x = dt_smp, j = "icao", value = as.factor(dt_smp[, icao]))
set(x = dt_smp, j = "zone", value = as.factor(dt_smp[, zone]))

# Convert the runway's name (e.g., RW26R) to its magnetic heading in degrees
# (e.g., 260) for later headwind calculation
dt_smp[, hdg := as.numeric(substr(rwy, 3L, 4L)) * 10L]

# For two runways with the same magnetic heading at a given airport (e.g. RWY26R
# and RWY26L), keep the one with the longest TODA (i.e. the most favorable case)
dt_smp <- dt_smp[dt_smp[, .I[which.max(toda)], by = .(icao, hdg)]$V1]

# Return the resulting count of runways to the console
nrow(dt_smp)

# ==============================================================================
# 2 Open the imported climate dataset (lazy; per-airport filter pushed down)
# ==============================================================================

# Ensure the output directory exists before the workers write into it
dir.create(dir$cli_pq, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 3 Transform the climate data for each airport
# ==============================================================================

fn_transform <- function(apt) {

  # The airport key arrives as a factor element; arrow's predicate and the file
  # name both need a plain character scalar
  apt <- as.character(apt)

  # Offset the start of each worker by a random duration to spread disk I/O load
  Sys.sleep(time = sample(x = 1L:(crs * 10L), size = 1L))

  # Inform the log file
  fn_log(apt, "(1/6) Fetching climate observations...")

  # ============================================================================
  # 3.1 Fetch and prepare the climate data for the current airport
  # ============================================================================

  # Fetch the climate data for the current airport from the Parquet dataset
  # (arrow pushes the icao filter down and prunes row groups by icao)
  dt_nc <- arrow::open_dataset(dir$imp_pq) |>
    dplyr::filter(icao == apt) |>
    dplyr::select(obs, icao, ssp, var, val) |>
    dplyr::collect() |>
    setDT()

  # Recast column types
  set(x = dt_nc, j = "icao", value = as.factor(dt_nc[, icao]))
  set(x = dt_nc, j = "ssp",  value = as.factor(dt_nc[, ssp]))
  set(x = dt_nc, j = "var",  value = as.factor(dt_nc[, var]))

  # Inform the log file
  fn_log(apt, "(2/6) Pivoting", format(x = nrow(dt_nc), big.mark = ","),
    "observations...")

  # Pivot the dataset from long to wide format
  dt_nc <- dcast.data.table(
    data      = dt_nc,
    formula   = obs + icao + ssp ~ var,
    value.var = "val"
  )

  # ============================================================================
  # 3.2 Correct surface pressure to the airport's field elevation [REV. 2026]
  # The model reports ps at its grid-cell mean elevation (z_model = orog), which
  # for coarse ~100 km cells can differ substantially from an airport's true
  # field elevation (z_field = elev). Reduce ps hypsometrically from z_model to
  # z_field so density (and the thrust model's pressure ratio) reflect the
  # airport, not the grid cell. This corrects a one-directional bias that is
  # largest at elevated airports (JNB, MEX, DEN, BOG, ADD, ...). Both elevations
  # are in metres; elev is set in 1_population.R, orog in 4_import.R.
  # ============================================================================

  # Field and model elevations for this airport in m (constant across its rows)
  z_field <- dt_smp[icao == apt, elev][1L]
  z_model <- dt_smp[icao == apt, orog][1L]

  # Hypsometric reduction of ps from the model elevation to the field elevation
  set(
    x     = dt_nc,
    j     = "ps",
    value = dt_nc[, ps] *
      ((dt_nc[, tas] - sim$isa_lap * (z_field - z_model)) / dt_nc[, tas])^
        (sim$g / (sim$rsp_air * sim$isa_lap))
  )

  # ============================================================================
  # 3.3 Calculate the air density of moist air at the current airport
  # ============================================================================

  # Inform the log file
  fn_log(apt, "(3/6) Calculating air density for",
    format(x = nrow(dt_nc), big.mark = ","), "observations...")

  # ============================================================================
  # 3.3.1 Air density from specific humidity (huss) [REV. 2026]
  # huss is co-sampled with ps/tas/uas/vas in the 6hrPt table, so no temporal
  # realignment is needed and the water-vapour partial pressure follows directly
  # from huss and ps — no relative humidity, no saturation-pressure polynomial,
  # and no supersaturation to cap. Adapted from the partial-pressure form of
  # Dalton's law with eps = Mw / Md = sim$mwr.
  # ============================================================================

  # Water-vapour partial pressure in Pa:  e = p * q / (eps + (1 - eps) * q)
  pv <- dt_nc[, ps] * dt_nc[, huss] /
    (sim$mwr + (1 - sim$mwr) * dt_nc[, huss])

  # Dry-air partial pressure in Pa
  pd <- dt_nc[, ps] - pv

  # Moist-air density in kg/m3 (ideal gas, partial pressures)
  set(
    x     = dt_nc,
    j     = "rho",
    value = pd / (sim$rsp_air * dt_nc[, tas]) +
            pv / (sim$rsp_h2o * dt_nc[, tas])
  )

  # ============================================================================
  # 3.4 Merge with the list of runways
  # ============================================================================

  # Extract the list of runways for the current airport
  dt_rwys <- dt_smp[icao == apt, ]

  # Return the Cartesian product of observations times runway headings
  dt_nc <- merge(x = dt_nc, y = dt_rwys, by = "icao", allow.cartesian = TRUE)

  # Inform the log file
  fn_log(apt, "(4/6) Calculating headwinds for",
    format(x = nrow(dt_nc), big.mark = ","), "observations...")

  # ============================================================================
  # 3.5 Calculate the wind vector for each runway
  # ============================================================================

  # Calculate the airport's wind speed in m/s
  set(x = dt_nc, j = "wnd_spd", value = sqrt(dt_nc[, uas]^2L + dt_nc[, vas]^2L))

  # Calculate the airport's wind direction in °
  set(
    x = dt_nc,
    j = "wnd_dir",
    value = (180L + (180L / pi) * atan2(dt_nc[, uas], dt_nc[, vas])) %% 360L
  )
  
  # Calculate each runway's headwind speed in m/s
  set(
    x = dt_nc,
    j = "hdw",
    value = dt_nc[, wnd_spd] * cos(abs(dt_nc[, hdg] - dt_nc[, wnd_dir])
      * pi / 180L)
  )

  # Keep only the runway with the maximum headwind speed (presumed to be
  #  the active runway) for each observation and experiment (SSP)
  dt_nc <- dt_nc[dt_nc[, .I[which.max(hdw)], by = .(obs, ssp)]$V1]

  # ============================================================================
  # 3.6 Write the data in wide format to a per-airport Parquet file
  # ============================================================================

  # Inform the log file
  fn_log(apt, "(5/6) Writing", format(x = nrow(dt_nc), big.mark = ","),
    "observations to Parquet...")

  # Create the year column
  set(
    x = dt_nc,
    j = "year",
    value = format.Date(x = dt_nc[, obs], format = "%Y")
  )

  # Select which columns to write and in which order
  cols <- c(
    "year", "obs", "icao", "lat", "lon", "zone", "ssp",
    "huss", "ps", "tas", "rho",
    "hdw", "rwy", "toda"
  )

  # Write one Parquet file for this airport (one airport per worker => no
  # contention); this is the dir$cli_pq dataset that 8_simulate.R reads
  arrow::write_parquet(
    x    = dt_nc[, ..cols],
    sink = file.path(dir$cli_pq, paste0(apt, ".parquet")),
    compression = "zstd"
  )

  # Inform the log file
  fn_log(apt, "(6/6) Written", format(x = nrow(dt_nc), big.mark = ","),
    "observations to Parquet.")

} # End of the fn_transform function

# ==============================================================================
# 4 Handle the parallel computation
# ==============================================================================

# Distribute the sample airports across the CPU cores
fn_par_lapply(
  crs = crs,
  pkg = c("arrow", "data.table", "dplyr", "stringr"),
  lst = unique(dt_smp[, icao], by = "icao"),
  fun = fn_transform
)

# ==============================================================================
# 5 Housekeeping
# ==============================================================================

# Stop the script timer
Sys.time() - start_time

# EOF
