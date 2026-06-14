# ==============================================================================
#    NAME: scripts/3_download.R
#   INPUT: Search criteria for the climate models, defined in section 1 below
# ACTIONS: Query the Earth System Grid Federation (ESGF) for fitting NetCDF data
#  OUTPUT: CSV file listing the matching NetCDF results
# RUNTIME: ~4 seconds (3.8 GHz CPU / 128 GB DDR4 RAM / SSD)
#  AUTHOR: Thomas D. Pellegrin <thomas@pellegr.in>
#    YEAR: 2023
# ==============================================================================

# ==============================================================================
# 0 Housekeeping
# ==============================================================================

# Clear the environment
rm(list = ls())

# Load the required libraries
library(dplyr)
library(epwshiftr)

# Import the common settings
source("scripts/0_common.R")

# Start a script timer
start_time <- Sys.time()

# Clear the console
cat("\014")

# ==============================================================================
# 1 Perform the query
# ==============================================================================

# Query the ESGF server
nc_files <- rbind(
  # Single query: all five variables share the 6hrPt (instantaneous) table.
  # huss (specific humidity) replaces hurs (relative humidity): huss is
  # co-published with ps/tas/uas/vas at 6hrPt so no temporal realignment needed.
  esgf_query(
    activity   = "ScenarioMIP",
    variable   = c("ps", "tas", "uas", "vas", "huss"),
    frequency  = c("6hrPt"),
    experiment = c("ssp126", "ssp245", "ssp370", "ssp585"),
    source     = "MPI-ESM1-2-HR",
    variant    = "r1i1p1f1",
    replica    = FALSE,
    latest     = TRUE,
    resolution = "100 km",
    type       = "File",
    limit      = 10000L,
    data_node  = NULL
  ),
  # Fixed orography (orog, fx table) for the field-elevation surface-pressure
  # correction in 5_transform.R. orography is time-invariant, so the field is
  # not re-published per scenario: take the single CMIP/historical fx field.
  esgf_query(
    activity   = "CMIP",
    variable   = c("orog"),
    frequency  = c("fx"),
    experiment = c("historical"),
    source     = "MPI-ESM1-2-HR",
    variant    = "r1i1p1f1",
    replica    = FALSE,
    latest     = TRUE,
    resolution = "100 km",
    type       = "File",
    limit      = 10000L,
    data_node  = NULL
  )
)

# ==============================================================================
# 2 Parse the results
# ==============================================================================

# Remove duplicate tracking ids (there is a weird ESGF server-side issue with a
# few identical hurs files appearing twice with the same tracking_id).
nc_files <- nc_files[!rev(duplicated(rev(nc_files$tracking_id))), ]

# Count number of unique datasets
length(unique(nc_files$dataset_id))

# Sum size of combined dataset in GB
sum(nc_files$file_size) / 10^9

# Average size per file
mean(nc_files$file_size)

# Count number of files
nrow(nc_files)

# Display the number of files by experiment (SSP) and variable
nc_files %>%
  group_by(experiment_id) %>%
  summarize(
    huss = sum(variable_id == "huss"),
    orog = sum(variable_id == "orog"),
    ps   = sum(variable_id == "ps"),
    tas  = sum(variable_id == "tas"),
    uas  = sum(variable_id == "uas"),
    vas  = sum(variable_id == "vas")
  )

# Save the query results to a file for later reference
write.csv(x = nc_files, file = fls$net)

# ==============================================================================
# 3 Housekeeping
# ==============================================================================

# Stop the script timer
Sys.time() - start_time

# EOF
