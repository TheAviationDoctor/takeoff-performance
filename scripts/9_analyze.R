# ==============================================================================
#    NAME: scripts/9_analyze.R
#   INPUT: Climate and takeoff Parquet datasets output by earlier scripts
# ACTIONS: Summarise the climate and takeoff data with arrow and plot the results
#  OUTPUT: Plot and summary data files saved to disk
# RUNTIME: Dominated by the arrow aggregations over dir$cli_pq and dir$tko_pq.
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
library(ggplot2)
library(rnaturalearth)
library(scales)
library(viridis)

# Import the common settings
source("scripts/0_common.R")

# Start a script timer
start_time <- Sys.time()

# Clear the console
cat("\014")

# ==============================================================================
# 1. Climate change summary
# ==============================================================================

# ==============================================================================
# 1.1 Read and summarise the climate model data from dir$cli_pq. Variables:
# tas  = Near-surface air temperature in °C
# huss = Near-surface specific humidity in g/kg [REV. 2026, was hurs]
# ps   = Near-surface air pressure in Pa
# rho  = Near-surface air density in kg/m³ [REV. 2026, single huss-based rho]
# hdw  = Near-surface headwind in m/s
# ==============================================================================

# Create the summary table (runtime: ~15 minutes)
# Aggregate the per-observation climate to max/avg/min per airport-year-SSP
# (arrow pushes the group-by down to its C++ engine; only the summary collects)
dt_cli <- arrow::open_dataset(dir$cli_pq) |>
  dplyr::group_by(year, ssp, zone, icao, lat, lon) |>
  dplyr::summarise(
    max_tas  = max(tas),  avg_tas  = mean(tas),  min_tas  = min(tas),
    max_huss = max(huss), avg_huss = mean(huss), min_huss = min(huss),
    max_ps   = max(ps),   avg_ps   = mean(ps),   min_ps   = min(ps),
    max_rho  = max(rho),  avg_rho  = mean(rho),  min_rho  = min(rho),
    max_hdw  = max(hdw),  avg_hdw  = mean(hdw),  min_hdw  = min(hdw),
    .groups  = "drop"
  ) |>
  dplyr::collect() |>
  setDT()

# Recast column types
set(x = dt_cli, j = "year", value = as.integer(dt_cli[, year]))
set(x = dt_cli, j = "zone", value = as.factor(dt_cli[, zone]))
set(x = dt_cli, j = "ssp",  value = as.factor(dt_cli[, ssp]))
set(x = dt_cli, j = "icao", value = as.factor(dt_cli[, icao]))

# Convert temperatures from °K to °C
cols <- c("max_tas", "avg_tas", "min_tas")
dt_cli[, (cols) := lapply(X = .SD, FUN = "-", sim$k_to_c), .SDcols = cols]

# Convert near-surface specific humidity from kg/kg to g/kg
cols <- c("max_huss", "avg_huss", "min_huss")
dt_cli[, (cols) := lapply(X = .SD, FUN = "*", 10^3), .SDcols = cols]

# Convert near-surface air pressure from Pa to hPa
cols <- c("max_ps", "avg_ps", "min_ps")
dt_cli[, (cols) := lapply(X = .SD, FUN = "/", 10^2), .SDcols = cols]

# Convert near-surface air density from kg/m³ to g/m³
cols <- c("max_rho", "avg_rho", "min_rho")
dt_cli[, (cols) := lapply(X = .SD, FUN = "*", 10^3), .SDcols = cols]

# Recode frigid airports to temperate
dt_cli[zone == "Frigid", zone := "Temperate"]

# Save the base values to disk
fwrite(
  x    = dt_cli,
  file = paste(dir$res, "dt_cli_base_values_by_airport.csv", sep = "/")
)

# Declare climatic variables and their statistics
cols <- list()
# Base values
cols$max <- paste("max", names(cli), sep = "_")
cols$avg <- paste("avg", names(cli), sep = "_")
cols$min <- paste("min", names(cli), sep = "_")
cols$all <- c(cols$max, cols$avg, cols$min)
# Locally-estimated scatterplot smoothing (LOESS) of base values
cols$max_loe <- paste(cols$max, "loe", sep = "_")
cols$avg_loe <- paste(cols$avg, "loe", sep = "_")
cols$min_loe <- paste(cols$min, "loe", sep = "_")
cols$all_loe <- paste(cols$all, "loe", sep = "_")
# Absolute changes to base values
cols$max_abs <- paste(cols$max, "abs", sep = "_")
cols$avg_abs <- paste(cols$avg, "abs", sep = "_")
cols$min_abs <- paste(cols$min, "abs", sep = "_")
cols$all_abs <- paste(cols$all, "abs", sep = "_")
# Relative changes to base values
cols$max_rel <- paste(cols$max, "rel", sep = "_")
cols$avg_rel <- paste(cols$avg, "rel", sep = "_")
cols$min_rel <- paste(cols$min, "rel", sep = "_")
cols$all_rel <- paste(cols$all, "rel", sep = "_")
# Absolute changes to LOESS values
cols$max_loe_abs <- paste(cols$max_loe, "abs", sep = "_")
cols$avg_loe_abs <- paste(cols$avg_loe, "abs", sep = "_")
cols$min_loe_abs <- paste(cols$min_loe, "abs", sep = "_")
cols$all_loe_abs <- paste(cols$all_loe, "abs", sep = "_")
# Relative changes to LOESS values
cols$max_loe_rel <- paste(cols$max_loe, "rel", sep = "_")
cols$avg_loe_rel <- paste(cols$avg_loe, "rel", sep = "_")
cols$min_loe_rel <- paste(cols$min_loe, "rel", sep = "_")
cols$all_loe_rel <- paste(cols$all_loe, "rel", sep = "_")

# ==============================================================================
# 1.2 Summarize climate change by airport
# ==============================================================================

# ==============================================================================
# 1.2.1 Calculate LOESS values by airport
# ==============================================================================

# Create a new data table for summarizing by airport
dt_cli_apt <- copy(dt_cli)

# Add locally-estimated scatterplot smoothing (LOESS) to base values
dt_cli_apt[,
  (cols$all_loe) := lapply(
    X   = .SD,
    FUN = function(x) {
      predict(loess(formula = x ~ year, span = .75, model = TRUE))
    }
  ),
  by      = c("ssp", "icao"),
  .SDcols = cols$all
]

# Save the LOESS values to disk
fwrite(
  x    = dt_cli_apt[, !cols$all, with = FALSE],
  file = paste(dir$res, "dt_cli_loess_values_by_airport.csv", sep = "/")
)

# ==============================================================================
# 1.2.2 Summarize the changes in LOESS values by airport
# ==============================================================================

# Calculate the absolute difference between each year and the first, by airport
dt_cli_apt[,
  (cols$all_loe_abs) := lapply(
    X   = .SD,
    FUN = function(x) {
      (x - x[1:1])
    }
  ),
  by      = c("ssp", "zone", "icao"),
  .SDcols = cols$all_loe
]

# Save the final-year changes in LOESS values by airport to disk
fwrite(
  x = dt_cli_apt[year == dt_cli_apt[which.max(year), year]
    ][, !c(cols$all, cols$all_loe), with = FALSE
    ][, (cols$all_loe_abs) := round(.SD, 1L), .SDcols = cols$all_loe_abs
    ][, melt(.SD, id.vars = c("year", "zone", "ssp", "icao", "lat", "lon"))
    ][, dcast(.SD, formula = year + icao + zone + lat + lon + variable ~ ssp)
    ][, variable := gsub("_abs", "", variable)
    ],
  file = paste(dir$res, "dt_cli_loess_changes_by_airport.csv", sep = "/")
)

# ==============================================================================
# 1.2.3 Plot the changes in LOESS values by airport onto a choropleth map
# ==============================================================================

# Transform the data
dt_plt <- dt_cli_apt[year == dt_cli_apt[which.max(year), year]
  ][, !c(cols$all, cols$all_loe), with = FALSE]

# Define the world object from the Natural Earth package
world <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")

# Create a function to plot results
fn_plot <- function(col) {

  ggplot() +
  geom_sf(data = world, fill = "gray") +
  coord_sf(expand = FALSE) +
  # Define the scales
  scale_x_continuous(breaks = c(-180L, 180L)) +
  scale_y_continuous(
    breaks = unique(unlist(x = geo, use.names = FALSE)),
    limits = c(-90L, 90L)
  ) +
  scale_color_viridis(
    direction = -1L,
    name      = cli[[gsub("max|avg|min|_|loe|abs", "", col)]],
    option    = "magma"
  ) +
  facet_wrap(facets = vars(toupper(ssp))) +
  geom_point(
    data     = dt_plt,
    mapping  = aes(x = lon, y = lat, color = get(col)),
    shape    = 20L,
    size     = 1L
  ) +
  # Add parallel labels
  geom_text(
    data  = world,
    color = "gray",
    hjust = 0L,
    label = "Arctic circle",
    size  = 1.5,
    x     = -179L,
    y     = unique(unlist(x = geo, use.names = FALSE))[5] - 2L
  ) +
  geom_text(
    data  = world,
    color = "gray",
    hjust = 0L,
    label = "Tropic of Cancer",
    size  = 1.5,
    x     = -179L,
    y     = unique(unlist(x = geo, use.names = FALSE))[4] - 2L
  ) +
  geom_text(
    data  = world,
    color = "gray",
    hjust = 0L,
    label = "Tropic of Capricorn",
    size  = 1.5,
    x     = -179L,
    y     = unique(unlist(x = geo, use.names = FALSE))[3] - 2L
  ) +
  geom_text(
    data  = world,
    color = "gray",
    hjust = 0L,
    label = "Antarctic circle",
    size  = 1.5,
    x     = -179L,
    y     = unique(unlist(x = geo, use.names = FALSE))[2] - 2L
  ) +
  # Add zonal labels
  geom_text(
    angle = 90L,
    data  = world,
    color = "gray",
    hjust = .5,
    label = unique(names(geo))[1],
    size  = 2L,
    x     = -175L,
    y     = mean(
      c(
        unique(unlist(x = geo, use.names = FALSE))[1],
        unique(unlist(x = geo, use.names = FALSE))[2]
      )
    )
  ) +
  geom_text(
    angle = 90L,
    data  = world,
    color = "gray",
    hjust = .5,
    label = unique(names(geo))[2],
    size  = 2L,
    x     = -175L,
    y     = mean(
      c(
        unique(unlist(x = geo, use.names = FALSE))[2],
        unique(unlist(x = geo, use.names = FALSE))[3]
      )
    )
  ) +
  geom_text(
    angle = 90L,
    data  = world,
    color = "gray",
    hjust = .5,
    label = unique(names(geo))[3],
    size  = 2L,
    x     = -175L,
    y     = mean(
      c(
        unique(unlist(x = geo, use.names = FALSE))[3],
        unique(unlist(x = geo, use.names = FALSE))[4]
      )
    )
  ) +
  geom_text(
    angle = 90L,
    data  = world,
    color = "gray",
    hjust = .5,
    label = unique(names(geo))[2],
    size  = 2L,
    x     = -175L,
    y     = mean(
      c(
        unique(unlist(x = geo, use.names = FALSE))[4],
        unique(unlist(x = geo, use.names = FALSE))[5]
      )
    )
  ) +
  geom_text(
    angle = 90L,
    data  = world,
    color = "gray",
    hjust = .5,
    label = unique(names(geo))[1],
    size  = 2L,
    x     = -175L,
    y     = mean(
      c(
        unique(unlist(x = geo, use.names = FALSE))[5],
        unique(unlist(x = geo, use.names = FALSE))[6]
      )
    )
  ) +
  theme_light() +
  theme(
    axis.title      = element_blank(),
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    legend.key.size = unit(.2, plt$units),
    text            = element_text(size = plt$text)
  )

  # Save the plot
  ggsave(
    filename = paste(
      "9_cli_map_of_",
      gsub("_abs", "", col),
      ".png",
      sep = ""
    ),
    plot   = last_plot(),
    device = plt$device,
    path   = plt$path,
    scale  = plt$scale,
    height = plt$height,
    width  = plt$width,
    units  = plt$units,
    dpi    = plt$dpi
  )

} # End of the fn_plot function

# Generate the plots
mapply(
  FUN = fn_plot,
  col = cols$all_loe_abs
)

# ==============================================================================
# 1.2.4 Plot the correlation between absolute latitude & changes in LOESS values
# ==============================================================================

# Create a function to plot results
fn_plot <- function(col) {

  ggplot(
    data    = dt_plt,
    mapping = aes(x = abs(lat), y = get(col))
  ) +
  # Define the scales
  scale_x_continuous(
    breaks = seq(from = 0L, to = 70L, length.out = 8L),
    limits = c(0L, 70L),
    name   = "Latitude (absolute)"
  ) +
  scale_y_continuous(
    name = paste(
      "Changes in",
      gsub("_", ". ", gsub("_loe_abs", "", col)),
      cli[[gsub("max|avg|min|_|loe|abs", "", col)]],
      sep = " "
    )
  ) +
  facet_wrap(facets = vars(toupper(ssp))) +
  geom_point(
    shape    = 20L,
    size     = 1L
  ) +
  geom_smooth(formula = y ~ x, method = "loess", linewidth = 1) +
  theme_light() +
  theme(
    axis.text        = element_text(size = plt$axis.text),
    axis.title       = element_text(size = plt$axis.title),
    strip.text       = element_text(size = plt$strip.text),
    panel.grid.minor = element_blank()
  )

  # Save the plot
  ggsave(
    filename = paste(
      "9_cli_correlation_plot_of_",
      gsub("_abs", "", col),
      "_with_latitude.png",
      sep = ""
    ),
    plot   = last_plot(),
    device = plt$device,
    path   = plt$path,
    scale  = plt$scale,
    height = plt$height,
    width  = plt$width,
    units  = plt$units,
    dpi    = plt$dpi
  )

} # End of the fn_plot function

# Generate the plots
mapply(
  FUN = fn_plot,
  col = cols$all_loe_abs
)

# ==============================================================================
# 1.3 Summarize climate change by zone
# ==============================================================================

# Create a new data table for summarizing by zone
dt_cli_zon <- copy(dt_cli)

# ==============================================================================
# 1.3.1 Summarize base and LOESS values by zone
# ==============================================================================

# Declare independent variables for grouping
grp <- c("year", "ssp", "zone")

# Average the max, mean, and min annual airport values by zone and globally
dt_cli_zon <- rbind(
  dt_cli_zon[, lapply(X = .SD, FUN = mean), by = grp, .SDcols = cols$all],
  dt_cli_zon[, zone := "Global"][, lapply(X = .SD, FUN = mean),
    by      = grp,
    .SDcols = cols$all
  ]
)

# Calculate LOESS values from the base values summarized by zone
dt_cli_zon[,
  (cols$all_loe) := lapply(
    X   = .SD,
    FUN = function(x) {
      predict(loess(formula = x ~ year, span = .75, model = TRUE))
    }
  ),
  by      = c("ssp", "zone"),
  .SDcols = cols$all
]

# ==============================================================================
# 1.3.2 Plot a global overview of the mean of each climate variable by SSP
# ==============================================================================

# Select base values for the global zone
dt_plt <- dt_cli_zon[zone == "Global"
][, c("year", "ssp", cols$avg), with = FALSE
][, melt(.SD, id.vars = c("year", "ssp"))
][, variable := gsub("avg_", "", variable)
][, variable := factor(variable, gsub("avg_", "", cols$avg))
]

# Select the starting label data (LOESS instead of base values)
labs_start <- dt_cli_zon[
  zone == "Global", .SD[which.min(year)],
  by = "ssp", .SDcols = c("year", cols$avg_loe)
][, melt(.SD, id.vars = c("year", "ssp"))
][, variable := gsub("avg_", "", variable)
][, variable := gsub("_loe", "", variable)
][, variable := factor(
  variable,
  gsub("avg_", "", gsub("_loe", "", cols$avg_loe))
)
]

# Select the ending label data (LOESS instead of base values)
labs_end <- dt_cli_zon[
  zone == "Global", .SD[which.max(year)],
  by = "ssp", .SDcols = c("year", cols$avg_loe)
][, melt(.SD, id.vars = c("year", "ssp"))
][, variable := gsub("avg_", "", variable)
][, variable := gsub("_loe", "", variable)
][, variable := factor(
  variable,
  gsub("avg_", "", gsub("_loe", "", cols$avg_loe))
)
]

# Rename facet labels to include units
labs <- paste(names(cli), unlist(unname(cli)), sep = " ")
names(labs) <- levels(dt_plt$variable)

# Build the plot
ggplot(data    = dt_plt) +
  geom_line(
    linewidth    = .2,
    mapping      = aes(
      x          = year,
      y          = value
    )
  ) +
  geom_smooth(
    formula      = y ~ x,
    method       = "loess",
    linewidth    = .5,
    mapping      = aes(
      x          = year,
      y          = value
    )
  ) +
  # Starting value labels
  geom_label(
    data = labs_start,
    aes(
      x          = year,
      y          = value,
      label      = sprintf(fmt = "%.1f", value)
    ),
    alpha        = .5,
    fill         = "white",
    label.r      = unit(0L, "lines"),
    label.size   = 0L,
    nudge_x      = 1.5,
    size         = plt$label.text
  ) +
  # Ending value labels
  geom_label(
    data = labs_end,
    aes(
      x          = year,
      y          = value,
      label      = sprintf(fmt = "%.1f", value)
    ),
    alpha        = .5,
    fill         = "white",
    label.r      = unit(0L, "lines"),
    label.size   = 0L,
    nudge_x      = -1.5,
    size         = plt$label.text
  ) +
  scale_x_continuous(name = "Year",  n.breaks = 5L) +
  scale_y_continuous(name = "Value") +
  facet_grid(
    rows = vars(variable),
    cols = vars(toupper(ssp)),
    scales       = "free_y",
    labeller     = labeller(variable = labs)
  ) +
  theme_light() +
  theme(
    axis.title.y = element_blank(),
    text         = element_text(size = plt$text)
  )

# Save the plot
ggsave(
  filename = "9_cli_global_overview.png",
  plot     = last_plot(),
  device   = plt$device,
  path     = plt$path,
  scale    = plt$scale,
  width    = plt$width,
  height   = plt$height,
  units    = plt$units,
  dpi      = plt$dpi
)

# ==============================================================================
# 1.3.3 Plot the annual base and LOESS values as line plots for all zones
# ==============================================================================

# Order the zones so the facets display in alphabetical order
dt_cli_zon[, zone := factor(
    zone,
    levels = sort(unique(levels(dt_cli_zon[, zone])))
  )
]

# Define independent variables for grouping
grp <- c("ssp", "zone")

# Create a function to plot results
fn_plot <- function(col) {

  # Build the plot
  ggplot(
    data = dt_cli_zon,
    mapping = aes(
      x     = year,
      y     = dt_cli_zon[[as.character(col)]]
    )
  ) +
    # Plot the base values
    geom_line(linewidth = .2) +
    # Plot the LOESS values
    geom_smooth(formula = y ~ x, method = "loess", linewidth = .5) +
    # Starting value labels
    geom_label(
      data = dt_cli_zon[, .SD[which.min(year)], by  = grp],
      aes(
        x = year,
        y = dt_cli_zon[, .SD[which.min(year)], by   = grp]
        [[paste(as.character(col), "loe", sep   = "_")]],
        label = sprintf(fmt = "%.1f",
          x = dt_cli_zon[, .SD[which.min(year)], by = grp]
          [[paste(as.character(col), "loe", sep = "_")]]
        )
      ),
      alpha      = .5,
      fill       = "white",
      label.r    = unit(0L, "lines"),
      label.size = 0L,
      nudge_x    = 5L,
      size       = plt$label.text
    ) +
    # Ending value labels
    geom_label(
      data = dt_cli_zon[, .SD[which.max(year)], by  = grp],
      aes(
        x = year,
        y = dt_cli_zon[, .SD[which.max(year)], by   = grp]
        [[paste(as.character(col), "loe", sep   = "_")]],
        label = sprintf(fmt = "%.1f",
          x = dt_cli_zon[, .SD[which.max(year)], by = grp]
          [[paste(as.character(col), "loe", sep = "_")]]
        )
      ),
      alpha      = .5,
      fill       = "white",
      label.r    = unit(0L, "lines"),
      label.size = 0L,
      nudge_x    = -5L,
      size       = plt$label.text
    ) +
    # Define the scales
    scale_x_continuous(name = "Year",  n.breaks = 5L) +
    scale_y_continuous(name = "Value", labels   = label_comma(accuracy = .1)) +
    facet_grid(
      rows = vars(zone),
      cols = vars(toupper(ssp)),
      scales = "free_y"
    ) +
    theme_light() +
    theme(
      axis.title.y = element_blank(),
      text         = element_text(size = plt$text)
    )

  # Save the plot
  ggsave(
    filename = tolower(paste("9_cli_lineplot_of_", col, ".png", sep = "")),
    plot     = last_plot(),
    device   = plt$device,
    path     = plt$path,
    scale    = plt$scale,
    width    = plt$width,
    height   = plt$height,
    units    = plt$units,
    dpi      = plt$dpi
  )

} # End of the fn_plot function

# Generate the plots
mapply(
  FUN = fn_plot,
  col = cols$all
)

# ==============================================================================
# 1.3.4 Summarize the changes in LOESS values by zone
# ==============================================================================

# Calculate the absolute difference between each year and the first, by zone
dt_cli_zon[,
  (cols$all_loe_abs) := lapply(
    X   = .SD,
    FUN = function(x) {
      (x - x[1:1])
    }
  ),
  by      = c("ssp", "zone"),
  .SDcols = cols$all_loe
]

# Save the final-year changes in LOESS values by zone to disk
fwrite(
  x = dt_cli_zon[year == dt_cli_zon[which.max(year), year]
  ][, !c(cols$all, cols$all_loe), with = FALSE
  ][, (cols$all_loe_abs) := round(.SD, 1L), .SDcols = cols$all_loe_abs
  ][, melt(.SD, id.vars = c("year", "zone", "ssp"))
  ][, dcast(.SD, formula = year + zone + variable ~ ssp)
  ][, zone := factor(zone, levels = sort(unique(levels(dt_cli_zon[, zone]))))
  ][, variable := gsub("_abs", "", variable)
  ][, variable := factor(
        variable,
        c(rbind(cols$max_loe, cols$avg_loe, cols$min_loe))
      )
  ][order(variable, zone)
  ],
  file = paste(dir$res, "dt_cli_loess_changes_by_zone.csv", sep = "/")
)

# ==============================================================================
# 1.4 Sensitivity analysis of rho to tas, ps, and huss [REV. 2026]
# Density is now computed with the pipeline's own moist-air formulation from
# specific humidity (huss), identical to 5_transform.R §3.3, instead of the
# former masscor/relative-humidity call.
# ==============================================================================

# Set the number of data points to plot
res <- 11L

# Moist-air density (kg/m3) from huss, ps (Pa) and tas (K). Matches 5_transform.
fn_rho <- function(tas, ps, huss) {
  pv <- ps * huss / (sim$mwr + (1 - sim$mwr) * huss)
  pd <- ps - pv
  pd / (sim$rsp_air * tas) + pv / (sim$rsp_h2o * tas)
}

# Representative sea-level specific humidity span for the humidity sensitivity:
# ISA air is dry, so vary huss from 0 up to ~saturation at ISA sea level (kg/kg)
huss_max <- 0.01

# Build a data table for the sensitivity analysis (SI units: K, Pa, kg/kg)
dt_cli_sa <- data.table(
  isa_tas  = rep(x = sim$isa_tas, times = res),
  isa_ps   = rep(x = sim$isa_ps,  times = res),
  isa_huss = rep(x = 0,           times = res),
  # Fractional change applied to the independent variables
  scale    = seq(from = 0L, to = .1, length.out = res)
)

# Flex one independent variable at a time
dt_cli_sa[,
          var_tas  := isa_tas  * (1L + scale)            # Increase tas by 10%
][,
  var_ps   := isa_ps   * (1L - scale)            # Decrease ps by 10%
][,
  var_huss := isa_huss + huss_max / max(scale) * scale   # 0 -> huss_max
]

# Calculate the sensitivity of air density to each independent variable
dt_cli_sa[, rho_tas  := fn_rho(var_tas, isa_ps,  isa_huss)]
dt_cli_sa[, rho_ps   := fn_rho(isa_tas, var_ps,  isa_huss)]
dt_cli_sa[, rho_huss := fn_rho(isa_tas, isa_ps,  var_huss)]

# Calculate relative changes in the air density
dt_cli_sa[,
          rho_tas_rel  := abs(rho_tas  / rho_tas[1:1] - 1L)
][,
  rho_ps_rel   := abs(rho_ps   / rho_ps[1:1] - 1L)
][,
  rho_huss_rel := abs(rho_huss / rho_huss[1:1] - 1L)
]

# Save the data to disk
fwrite(
  x    = dt_cli_sa,
  file = paste(dir$res, "dt_cli_rho_sensitivity_analysis.csv", sep = "/")
)

# Plot the relative changes in air density (DV) based on changes to the IVs
ggplot(
  data = dt_cli_sa,
  mapping = aes(x = scale)
) +
  # Add lines
  geom_line(mapping = aes(y = rho_tas_rel),  linewidth = 1L) +
  geom_line(mapping = aes(y = rho_ps_rel),   linewidth = 1L) +
  geom_line(mapping = aes(y = rho_huss_rel), linewidth = 1L) +
  # Add labels
  geom_label(mapping = aes(x = .1, y = max(rho_tas_rel),  label = "tas")) +
  geom_label(mapping = aes(x = .1, y = max(rho_ps_rel),   label = "ps")) +
  geom_label(mapping = aes(x = .1, y = max(rho_huss_rel), label = "huss")) +
  # Define scales
  scale_x_continuous(
    name   = "Absolute percentage of change in tas, ps, or huss",
    labels = scales::label_percent()
  ) +
  scale_y_continuous(
    name   = "Absolute percentage of change in rho",
    labels = scales::label_percent()
  ) +
  theme_light() +
  theme(
    text = element_text(size = plt$text)
  )

# Save the plot
ggsave(
  filename = paste("9_cli_rho_sensitivity_analysis.png", sep = ""),
  plot     = last_plot(),
  device   = plt$device,
  path     = plt$path,
  scale    = plt$scale,
  height   = plt$height,
  width    = plt$width,
  units    = plt$units,
  dpi      = plt$dpi
)

# ==============================================================================
# 2. Takeoff outcomes summary
# ==============================================================================

# ==============================================================================
# 2.1 Read and summarise the takeoff simulation data [REV. 2026]
# Source: the per-airport Parquet dataset written by 8_simulate.R (dir$tko_pq),
# read and aggregated with arrow instead of the former MySQL `tko` round-trip.
# Each row carries an `outcome` of feasible / infeasible_field / infeasible_thrust
# (see 6_model.R §3.3.1 and 8_simulate.R), so feasibility is read straight from
# `outcome` rather than re-derived from todr <= toda. The iteration count (itr)
# no longer exists: the vectorised bisection solver replaced the per-passenger
# loop. Counts per (year, ssp, zone, icao, lat, lon, type):
# tko_ok_thr_min       = feasible takeoffs at maximum derate (thr_red = thr_ini)
# tko_ok_thr_mid       = feasible takeoffs at partial derate (0 < thr_red < ini)
# tko_ok_thr_max_no_rm = feasible at full thrust (thr_red = 0), no payload removal
# tko_ok_thr_max_rm    = feasible at full thrust (thr_red = 0), with payload removal
# tko_ok_thr_max       = all feasible takeoffs at full thrust (thr_red = 0)
# tko_ok               = all feasible takeoffs (any derate or payload removal)
# tko_ko               = all infeasible takeoffs (field-length- or thrust-limited)
# tko                  = count of all takeoffs (feasible or not)
# ==============================================================================

# Aggregate the outcomes over each airport's observations (arrow pushes the
# group-by down to its C++ engine, so only the small summary is collected)
dt_tko <- arrow::open_dataset(dir$tko_pq) |>
  dplyr::group_by(year, ssp, zone, icao, lat, lon, type) |>
  dplyr::summarise(
    tko_ok_thr_min       = sum(outcome == "feasible" & thr_red == sim$thr_ini),
    tko_ok_thr_mid       = sum(outcome == "feasible" &
                                 thr_red >= 1L & thr_red <= sim$thr_ini - 1L),
    tko_ok_thr_max_no_rm = sum(outcome == "feasible" &
                                 thr_red == 0L & tom_rem == 0L),
    tko_ok_thr_max_rm    = sum(outcome == "feasible" &
                                 thr_red == 0L & tom_rem > 0L),
    tko_ok_thr_max       = sum(outcome == "feasible" & thr_red == 0L),
    tko_ok               = sum(outcome == "feasible"),
    tko_ko               = sum(outcome != "feasible"),
    tko                  = dplyr::n(),
    .groups              = "drop"
  ) |>
  dplyr::collect() |>
  setDT()

# Coerce arrow's 64-bit integer counts to base numeric for downstream maths
cnt <- grep("^tko", names(dt_tko), value = TRUE)
dt_tko[, (cnt) := lapply(.SD, as.numeric), .SDcols = cnt]

# Recast column types
set(x = dt_tko, j = "year", value = as.integer(dt_tko[, year]))
set(x = dt_tko, j = "zone", value = as.factor(dt_tko[, zone]))
set(x = dt_tko, j = "ssp",  value = as.factor(dt_tko[, ssp]))
set(x = dt_tko, j = "icao", value = as.factor(dt_tko[, icao]))
set(x = dt_tko, j = "type", value = as.factor(dt_tko[, type]))

# Recode frigid airports to temperate
dt_tko[zone == "Frigid", zone := "Temperate"]

# Replace the aircraft types by their body types
levels(dt_tko$type) <- bod

# Summarize the data by body type
dt_tko <- dt_tko[,
  lapply(X = .SD, FUN = sum),
  by = c("year", "ssp", "zone", "icao", "lat", "lon", "type")
]

# Save the base values to disk
fwrite(
  x    = dt_tko,
  file = paste(dir$res, "dt_tko_base_values_by_airport.csv", sep = "/")
)

# Declare output variables
cols     <- list()
cols$bas <- grep("tko", names(dt_tko), value = TRUE) # Absolute base values
cols$rel <- paste(cols$bas, "rel", sep = "_")        # Relative base values
cols$loe <- paste(cols$rel, "loe", sep = "_")        # LOESS values
cols$dif <- paste(cols$loe, "dif", sep = "_")        # Changes in LOESS values

# ==============================================================================
# 2.1.1 Infeasibility analysis by zone, SSP, and year [REV. 2026]
# The share of takeoffs that cannot be flown — split into field-length-limited
# (finite TODR exceeding TODA or todr_cap) and thrust-limited (cannot reach
# Vlof) — is itself a publishable result, so summarise and plot it on its own.
# ==============================================================================

# Aggregate the three outcome categories by zone, SSP, and year
dt_inf <- arrow::open_dataset(dir$tko_pq) |>
  dplyr::group_by(year, ssp, zone) |>
  dplyr::summarise(
    feasible          = sum(outcome == "feasible"),
    infeasible_field  = sum(outcome == "infeasible_field"),
    infeasible_thrust = sum(outcome == "infeasible_thrust"),
    tko               = dplyr::n(),
    .groups           = "drop"
  ) |>
  dplyr::collect() |>
  setDT()

# Coerce arrow's 64-bit integer counts to base numeric
cnt <- c("feasible", "infeasible_field", "infeasible_thrust", "tko")
dt_inf[, (cnt) := lapply(.SD, as.numeric), .SDcols = cnt]

# Recast column types
set(x = dt_inf, j = "year", value = as.integer(dt_inf[, year]))
set(x = dt_inf, j = "zone", value = as.factor(dt_inf[, zone]))
set(x = dt_inf, j = "ssp",  value = as.factor(dt_inf[, ssp]))

# Recode frigid airports to temperate, then re-aggregate the merged zone
dt_inf[zone == "Frigid", zone := "Temperate"]
dt_inf <- dt_inf[,
  lapply(X = .SD, FUN = sum),
  by      = c("year", "ssp", "zone"),
  .SDcols = cnt
]

# Compute the infeasible shares as fractions of all takeoffs
dt_inf[, `:=`(
  share_field  = infeasible_field  / tko,
  share_thrust = infeasible_thrust / tko,
  share_infeas = (infeasible_field + infeasible_thrust) / tko
)]

# Save the infeasibility summary to disk
fwrite(
  x    = dt_inf,
  file = paste(dir$res, "dt_tko_infeasible_share_by_zone.csv", sep = "/")
)

# Plot the infeasible share over time, by climate zone and SSP
ggplot(
  data    = dt_inf,
  mapping = aes(x = year, y = share_infeas, color = ssp)
) +
  geom_line(linewidth = .5) +
  scale_x_continuous("Year") +
  scale_y_continuous("Share of infeasible takeoffs", labels = scales::percent) +
  scale_color_viridis(discrete = TRUE, name = "SSP") +
  facet_wrap(~zone, ncol = 1L, scales = "free_y") +
  theme_light()

# Save the plot
ggsave(
  filename = "9_infeasible_share_by_zone.png",
  plot     = last_plot(),
  device   = "png",
  path     = "plots",
  scale    = 1L,
  width    = 6L,
  height   = NA,
  units    = "in",
  dpi      = "retina"
)

# ==============================================================================
# 2.2 Summarize takeoff outcomes by airport
# ==============================================================================

# ==============================================================================
# 2.2.1 Calculate LOESS values by airport
# ==============================================================================

# Create a new data table for summarizing by airport
dt_tko_apt <- copy(dt_tko)

# Add locally-estimated scatterplot smoothing (LOESS) to base values
dt_tko_apt[,
  # Convert absolute values to relative (percentage)
  (cols$rel) := lapply(
    X   = .SD,
    FUN = function(x) {
      x / tko * 100L
    }
  ),
  .SDcols = cols$bas
  # Add LOESS values
][, (cols$loe) := lapply(
    X   = .SD,
    FUN = function(x) {
      predict(loess(formula = x ~ year, span = .75, model = TRUE))
    }
  ),
  by      = c("ssp", "zone", "icao", "type"),
  .SDcols = cols$rel
][,
  # Remove unneeded columns
  c(cols$bas, cols$rel) := NULL
]

# Save the base values to disk
fwrite(
  x    = dt_tko_apt,
  file = paste(dir$res, "dt_tko_loess_values_by_airport.csv", sep = "/")
)

# ==============================================================================
# 2.2.2 Summarize the changes in LOESS values by airport
# ==============================================================================

# Calculate the absolute difference between each year and the first, by airport
dt_tko_apt[,
  (cols$dif) := lapply(
    X   = .SD,
    FUN = function(x) {
      (x - x[1:1])
    }
  ),
  by      = c("ssp", "zone", "icao", "type"),
  .SDcols = cols$loe
][, (cols$loe) := NULL]

# Save the final-year changes in LOESS values by airport to disk
fwrite(
  x = dt_tko_apt[year == dt_tko_apt[which.max(year), year]
    ][, (cols$dif) := round(.SD, 1L), .SDcols = cols$dif
    ][, melt(
          data    = .SD,
          id.vars = c("year", "ssp", "zone", "icao", "lat", "lon", "type")
        )
    ][, dcast(
          data    = .SD,
          formula = year + zone + icao + lat + lon + variable ~ type + ssp)
    ],
  file = paste(dir$res, "dt_tko_loess_changes_by_airport.csv", sep = "/")
)

# ==============================================================================
# 2.2.3 Plot the changes in LOESS values by airport onto a choropleth map
# ==============================================================================

# Define the world object from the Natural Earth package
world <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")

# Create a function to plot results
fn_plot <- function(body, cols) {

  ggplot() +
  geom_sf(data = world, fill = "gray") +
  coord_sf(expand = FALSE) +
  # Define the scales
  scale_x_continuous(breaks = c(-180L, 180L)) +
  scale_y_continuous(
    breaks = unique(unlist(x = geo, use.names = FALSE)),
    limits = c(-90L, 90L)
  ) +
  scale_color_viridis(
    direction = -1L,
    name      = "in p. p.",
    option    = "magma"
  ) +
  facet_wrap(facets = vars(toupper(ssp))) +
  geom_point(
    data = dt_tko_apt[year == dt_tko_apt[which.max(year), year] & type == body],
    mapping  = aes(
      x      = lon,
      y      = lat,
      color  = .data[[as.character(cols)]]
    ),
    shape    = 20L,
    size     = 1L
  ) +
  # Add parallel labels
  geom_text(
    data  = world,
    color = "gray",
    hjust = 0L,
    label = "Arctic circle",
    size  = 1.5,
    x     = -179L,
    y     = unique(unlist(x = geo, use.names = FALSE))[5] - 2L
  ) +
  geom_text(
    data  = world,
    color = "gray",
    hjust = 0L,
    label = "Tropic of Cancer",
    size  = 1.5,
    x     = -179L,
    y     = unique(unlist(x = geo, use.names = FALSE))[4] - 2L
  ) +
  geom_text(
    data  = world,
    color = "gray",
    hjust = 0L,
    label = "Tropic of Capricorn",
    size  = 1.5,
    x     = -179L,
    y     = unique(unlist(x = geo, use.names = FALSE))[3] - 2L
  ) +
  geom_text(
    data  = world,
    color = "gray",
    hjust = 0L,
    label = "Antarctic circle",
    size  = 1.5,
    x     = -179L,
    y     = unique(unlist(x = geo, use.names = FALSE))[2] - 2L
  ) +
  # Add zonal labels
  geom_text(
    angle = 90L,
    data  = world,
    color = "gray",
    hjust = .5,
    label = unique(names(geo))[1],
    size  = 2L,
    x     = -175L,
    y     = mean(
      c(
        unique(unlist(x = geo, use.names = FALSE))[1],
        unique(unlist(x = geo, use.names = FALSE))[2]
      )
    )
  ) +
  geom_text(
    angle = 90L,
    data  = world,
    color = "gray",
    hjust = .5,
    label = unique(names(geo))[2],
    size  = 2L,
    x     = -175L,
    y     = mean(
      c(
        unique(unlist(x = geo, use.names = FALSE))[2],
        unique(unlist(x = geo, use.names = FALSE))[3]
      )
    )
  ) +
  geom_text(
    angle = 90L,
    data  = world,
    color = "gray",
    hjust = .5,
    label = unique(names(geo))[3],
    size  = 2L,
    x     = -175L,
    y     = mean(
      c(
        unique(unlist(x = geo, use.names = FALSE))[3],
        unique(unlist(x = geo, use.names = FALSE))[4]
      )
    )
  ) +
  geom_text(
    angle = 90L,
    data  = world,
    color = "gray",
    hjust = .5,
    label = unique(names(geo))[2],
    size  = 2L,
    x     = -175L,
    y     = mean(
      c(
        unique(unlist(x = geo, use.names = FALSE))[4],
        unique(unlist(x = geo, use.names = FALSE))[5]
      )
    )
  ) +
  geom_text(
    angle = 90L,
    data  = world,
    color = "gray",
    hjust = .5,
    label = unique(names(geo))[1],
    size  = 2L,
    x     = -175L,
    y     = mean(
      c(
        unique(unlist(x = geo, use.names = FALSE))[5],
        unique(unlist(x = geo, use.names = FALSE))[6]
      )
    )
  ) +
  theme_light() +
  theme(
    axis.title      = element_blank(),
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    legend.key.size = unit(.2, plt$units),
    text            = element_text(size = plt$text)
  )

  # Save the plot
  ggsave(
    filename = tolower(
      paste(
        "9_tko_map_of_",
        body,
        "_",
        cols,
        ".png",
        sep = ""
      )
    ),
    plot   = last_plot(),
    device = plt$device,
    path   = plt$path,
    scale  = plt$scale,
    height = plt$height,
    width  = plt$width,
    units  = plt$units,
    dpi    = plt$dpi
  )

} # End of the fn_plot function

# Combine the aircraft bodies and output variables to be plotted
mix <- expand.grid(
  body = names(bod),
  cols = cols$dif
)

mapply(
  FUN  = fn_plot,
  body = mix$body,
  col  = mix$cols
)

# ==============================================================================
# 2.3 Summarize takeoff outcomes by climate zone
# ==============================================================================

# ==============================================================================
# 2.3.1 Calculate LOESS values by climate zone
# ==============================================================================

# Create a new data table for summarizing by airport
dt_tko_zon <- copy(dt_tko)

# Summarize the data by group, then combine them
dt_tko_zon <- rbind(
  # Zonal summary by group
  dt_tko_zon[, lapply(X = .SD, FUN = sum),
    by      = c("year", "ssp", "zone", "type"),
    .SDcols = cols$bas
  ],
  # Global summary by group
  dt_tko_zon[, zone := "Global"][, lapply(X = .SD, FUN = sum),
    by      = c("year", "ssp", "zone", "type"),
    .SDcols = cols$bas
  ]
)

# Create relative and LOESS output variables as a percentage of all takeoffs
dt_tko_zon[,
  (cols$rel) := lapply(
    X   = .SD,
    FUN = function(x) {
      x / tko
    }
  ),
  .SDcols = cols$bas
][, (cols$loe) := lapply(
    X   = .SD,
    FUN = function(x) {
      predict(loess(formula = x ~ year, span = .75, model = TRUE))
    }
  ),
  by      = c("ssp", "zone", "type"),
  .SDcols = cols$rel
]

# Save the data to disk
fwrite(
  x    = dt_tko_zon,
  file = paste(dir$res, "dt_tko_loess_values_by_zone.csv", sep = "/")
)

# ==============================================================================
# 2.3.2 Plot the results onto chronological lineplots
# ==============================================================================

# Order the zones so the facets display in alphabetical order
dt_tko_zon[, zone := factor(
  zone,
  levels = sort(unique(levels(dt_tko_zon[, zone])))
)
]

# Declare independent variables for grouping
grp <- c("ssp", "zone")

# Create a function to plot results
fn_plot <- function(body, cols) {

  # Build the plot
  ggplot(
    data    = dt_tko_zon[type == body],
    mapping = aes(
      x     = year,
      y     = dt_tko_zon[type == body][[as.character(cols)]]
    )
  ) +
    geom_line(linewidth = .2) +
    geom_smooth(formula = y ~ x, method = "loess", linewidth = .5) +
    # Starting value labels
    geom_label(
      data = dt_tko_zon[type == body][, .SD[which.min(year)], by = grp],
      aes(
        x = year,
        y = dt_tko_zon[type == body][, .SD[which.min(year)], by = grp]
          [[paste(as.character(cols), "loe", sep = "_")]],
          label = sprintf(fmt = "%1.1f%%", dt_tko_zon[type == body]
          [, .SD[which.min(year)], by = grp]
          [[paste(as.character(cols), "loe", sep = "_")]] * 100L)
      ),
      alpha      = .5,
      fill       = "white",
      label.r    = unit(0L, "lines"),
      label.size = 0L,
      nudge_x    = 4L,
      size       = 2L
    ) +
    # Ending value labels
    geom_label(
      data = dt_tko_zon[type == body][, .SD[which.max(year)], by = grp],
      aes(
        x = year,
        y = dt_tko_zon[type == body][, .SD[which.max(year)], by = grp]
          [[paste(as.character(cols), "loe", sep = "_")]],
          label = sprintf(fmt = "%1.1f%%", dt_tko_zon[type == body]
          [, .SD[which.max(year)], by = grp]
          [[paste(as.character(cols), "loe", sep = "_")]] * 100L)
      ),
      alpha      = .5,
      fill       = "white",
      label.r    = unit(0L, "lines"),
      label.size = 0L,
      nudge_x    = -4L,
      size       = 2L
    ) +
    scale_x_continuous(name = "Year", n.breaks = 3L) +
    scale_y_continuous(
      name   = "Value",
      labels = scales::label_percent(accuracy = .1)
    ) +
    facet_grid(
      rows = vars(zone),
      cols = vars(toupper(ssp)),
      scales = "free_y"
    ) +
    theme_light() +
    theme(
      axis.title.y = element_blank(),
      text         = element_text(size = plt$text)
    )

# Save the plot
ggsave(
  filename = tolower(
    paste(
      "9_tko_lineplot_of_",
      body,
      "_",
      cols,
      ".png",
      sep = ""
    )
  ),
  plot     = last_plot(),
  device   = plt$device,
  path     = plt$path,
  scale    = plt$scale,
  height   = plt$height,
  width    = plt$width,
  units    = plt$units,
  dpi      = plt$dpi
)

} # End of the fn_plot function

# Combine the aircraft bodies and output variables to be plotted
mix <- expand.grid(
  body = names(bod),
  cols = cols$rel
)

# Generate the plots
mapply(
  FUN  = fn_plot,
  body = mix$body,
  cols = mix$cols
)

# ==============================================================================
# 2.3.3 Summarize the changes in LOESS values by climate zone
# ==============================================================================

# Calculate the absolute difference between each year and the first, by zone
dt_tko_zon[,
  (cols$dif) := lapply(
    X     = .SD,
    FUN   = function(x) {
      (x - x[1:1])
    }
  ),
  by      = c("ssp", "zone", "type"),
  .SDcols = cols$loe
]

# Save the final-year changes in LOESS values by zone to disk
fwrite(
  x = dt_tko_zon[year == dt_tko_zon[which.max(year), year]
  ][, c("year", "ssp", "zone", "type", cols$dif), with = FALSE
  ][, (cols$dif) := round(.SD * 100L, 1L), .SDcols = cols$dif
  ][, melt(data = .SD, id.vars = c("year", "zone", "ssp", "type"))
  ][, dcast(data = .SD, formula = year + zone + variable ~ type + ssp)
  ][, zone := factor(zone, levels = sort(unique(levels(dt_tko_zon[, zone]))))
  ][, variable := factor(variable, cols$dif)
  ][order(variable, zone)
  ],
  file = paste(dir$res, "dt_tko_loess_changes_by_zone.csv", sep = "/")
)

# ==============================================================================
# 3. Research questions summary
# ==============================================================================

# ==============================================================================
# 3.1 Create, fetch, and cleanse the data. Variables:
# avg_todr    = Mean takeoff distance required in m
# avg_thr_red = Mean thrust reduction in percentage points of TOGA
# avg_tom_rem = Mean takeoff mass reduction in kg
# ==============================================================================

# Read and average the feasible takeoffs per airport [REV. 2026]
# Source: dir$tko_pq via arrow (was the MySQL `res` summary). Only feasible
# takeoffs contribute, matching the former WHERE todr <= toda; their todr is
# finite, while infeasible rows (excluded here) carry todr = NA.
dt_res <- arrow::open_dataset(dir$tko_pq) |>
  dplyr::filter(outcome == "feasible") |>
  dplyr::group_by(year, ssp, zone, icao, lat, lon, type) |>
  dplyr::summarise(
    avg_todr    = mean(todr),
    avg_thr_red = mean(thr_red),
    avg_tom_rem = mean(tom_rem),
    .groups     = "drop"
  ) |>
  dplyr::collect() |>
  setDT()

# Recast column types
set(x = dt_res, j = "year", value = as.integer(dt_res[, year]))
set(x = dt_res, j = "zone", value = as.factor(dt_res[, zone]))
set(x = dt_res, j = "ssp",  value = as.factor(dt_res[, ssp]))
set(x = dt_res, j = "icao", value = as.factor(dt_res[, icao]))
set(x = dt_res, j = "type", value = as.factor(dt_res[, type]))

# Recode frigid airports to temperate
dt_res[zone == "Frigid", zone := "Temperate"]

# Combine the aircraft types to narrow/widebody
levels(dt_res$type) <- bod

# Convert thrust reduction below TOGA to thrust as a percentage of TOGA
dt_res[, avg_thr := (100L - avg_thr_red) / 100L][, avg_thr_red := NULL]

# Convert payload removal in kg to passengers based on pax mass assumptions
dt_res[, avg_pax_rem := avg_tom_rem / sim$pax_avg][, avg_tom_rem := NULL]


# Summarize the data by body type
dt_res <- dt_res[,
  lapply(X = .SD, FUN = mean),
  by = c("year", "ssp", "zone", "icao", "lat", "lon", "type")
]

# Save the base values to disk
fwrite(
  x    = dt_res,
  file = paste(dir$res, "dt_res_base_values_by_airport.csv", sep = "/")
)

# Declare output variables
cols     <- list()
cols$bas <- grep("todr|thr|pax", names(dt_res), value = TRUE) # Base values
cols$loe <- paste(cols$bas, "loe", sep = "_")                 # LOESS values

# ==============================================================================
# 3.2 Summarize results by airport
# Some airports do not have enough successful takeoffs to use the LOESS method
# ==============================================================================

# ==============================================================================
# 3.2.1 Calculate changes in base values by airport
# ==============================================================================

# Create a new data table for summarizing by airport
dt_res_apt <- copy(dt_res)

# Define change variable
cols$dif <- paste(cols$bas, "dif", sep = "_")

# Calculate the absolute difference between each year and the first, by airport
dt_res_apt[,
  (cols$dif) := lapply(
    X   = .SD,
    FUN = function(x) {
      (x - x[1:1])
    }
  ),
  by      = c("ssp", "icao", "type"),
  .SDcols = cols$bas
][, (cols$bas) := NULL]

# Save the final-year changes in base values by airport to disk
fwrite(
  x = dt_res_apt[year == dt_res_apt[which.max(year), year]
  ][, avg_thr_dif := avg_thr_dif * 100L # Change to percentage
  ][, (cols$dif) := round(.SD, 1L), .SDcols = cols$dif
  ][, melt(
    data    = .SD,
    id.vars = c("year", "ssp", "zone", "icao", "lat", "lon", "type")
  )
  ][, dcast(
    data    = .SD,
    formula = year + zone + icao + lat + lon + variable ~ type + ssp)
  ],
  file = paste(dir$res, "dt_res_base_changes_by_airport.csv", sep = "/")
)

# ==============================================================================
# 3.2.2 Boxplot the changes in base values by airport
# ==============================================================================

# Add a global zone to the plot data
dt_plt <- rbind(
  dt_res_apt[year == dt_res_apt[which.max(year), year]],
  dt_res_apt[year == dt_res_apt[which.max(year), year]][, zone := "Global"]
)

# Order the climate zones so they display in alphabetical order
dt_plt[, zone := factor(zone, levels = sort(unique(levels(dt_plt[, zone]))))]

# Save the quantile values by SSP and aircraft type
fwrite(
  x = dt_plt[
    zone == "Global",
    lapply(X = .SD, FUN = quantile, probs = sim$quantiles),
    .SDcols = cols$dif, by = c("ssp", "type")
  ][, avg_thr_dif := avg_thr_dif * 100L
  ][, (cols$dif) := round(.SD, 1L), .SDcols = cols$dif
  ][, quantile := rep(x = sim$quantiles, times = 8L)
  ],
  file = paste(dir$res, "dt_res_base_changes_quantiles.csv", sep = "/")
)

# Create a function to plot results
fn_plot <- function(body, cols) {

  # Build the plot
  ggplot(
    data = dt_plt,
    mapping = aes(x = zone, y = .data[[as.character(cols)]])
  ) +
    stat_boxplot(geom = "errorbar", linewidth = .3) + # Add whisker ends
    geom_boxplot(outlier.size = .2, linewidth = .3) +
    scale_x_discrete(name = NULL) +
    scale_y_continuous(
      name   = NULL,
      labels = ifelse(
        cols == "avg_thr_dif",
        scales::label_percent(accuracy = 1L),
        scales::label_comma(accuracy = 1L)
      )
    ) +
    # Zoom into the canvas with different limits for each variable
    coord_cartesian(
      ylim = if (cols == "avg_pax_rem_dif") { c(-10L, 10L) }
        else if (cols == "avg_thr_dif")     { c(-.02, .02) }
        else if (cols == "avg_todr_dif")    { c(-100L, 100L) }
    ) +
    stat_summary(fun = mean, size = .01) + # Display the mean onto the boxplot
    facet_grid(
      rows = vars(type),
      cols = vars(toupper(ssp))
    ) +
    theme_light() +
    theme(
      axis.text  = element_text(size = plt$axis.text),
      axis.title = element_text(size = plt$axis.title),
      strip.text = element_text(size = plt$strip.text)
    )
  
  # Save the plot
  ggsave(
    filename = paste(
      "9_res_boxplot_of_",
      cols,
      ".png",
      sep = ""
    ),
    plot     = last_plot(),
    device   = plt$device,
    path     = plt$path,
    scale    = plt$scale,
    height   = plt$height,
    width    = plt$width,
    units    = plt$units,
    dpi      = plt$dpi
  )
  
} # End of the fn_plot function

# Combine the aircraft bodies and output variables to be plotted
mix <- expand.grid(
  body = names(bod),
  cols = cols$dif
)

# Generate the plots
mapply(
  FUN  = fn_plot,
  body = mix$body,
  cols = mix$cols
)

# ==============================================================================
# 3.2.3 Plot the changes in base values by airport onto a choropleth map
# ==============================================================================

# Define the world object from the Natural Earth package
world <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")

# Create a function to plot results
fn_plot <- function(body, cols, labs) {
print(cols)
  ggplot() +
  geom_sf(data = world, fill = "gray") +
  coord_sf(expand = FALSE) +
  # Define the scales
  scale_x_continuous(breaks = c(-180L, 180L)) +
  scale_y_continuous(
    breaks = unique(unlist(x = geo, use.names = FALSE)),
    limits = c(-90L, 90L)
  ) +
  scale_color_viridis(
    direction = -1L,
    name      = labs,
    option    = "magma"
  ) +
  # facet_wrap(~toupper(ssp)) +
  facet_wrap(facets = vars(toupper(ssp))) +
  geom_point(
    data = dt_res_apt[year == dt_res_apt[which.max(year), year] & type == body],
    mapping  = aes(
      x      = lon,
      y      = lat,
      color  = .data[[as.character(cols)]]
    ),
    shape    = 20L,
    size     = 1L
  ) +
  # Add parallel labels
  geom_text(
    data  = world,
    color = "gray",
    hjust = 0L,
    label = "Arctic circle",
    size  = 1.5,
    x     = -179L,
    y     = unique(unlist(x = geo, use.names = FALSE))[5] - 2L
  ) +
  geom_text(
    data  = world,
    color = "gray",
    hjust = 0L,
    label = "Tropic of Cancer",
    size  = 1.5,
    x     = -179L,
    y     = unique(unlist(x = geo, use.names = FALSE))[4] - 2L
  ) +
  geom_text(
    data  = world,
    color = "gray",
    hjust = 0L,
    label = "Tropic of Capricorn",
    size  = 1.5,
    x     = -179L,
    y     = unique(unlist(x = geo, use.names = FALSE))[3] - 2L
  ) +
  geom_text(
    data  = world,
    color = "gray",
    hjust = 0L,
    label = "Antarctic circle",
    size  = 1.5,
    x     = -179L,
    y     = unique(unlist(x = geo, use.names = FALSE))[2] - 2L
  ) +
  # Add zonal labels
  geom_text(
    angle = 90L,
    data  = world,
    color = "gray",
    hjust = .5,
    label = unique(names(geo))[1],
    size  = 2L,
    x     = -175L,
    y     = mean(
      c(
        unique(unlist(x = geo, use.names = FALSE))[1],
        unique(unlist(x = geo, use.names = FALSE))[2]
      )
    )
  ) +
  geom_text(
    angle = 90L,
    data  = world,
    color = "gray",
    hjust = .5,
    label = unique(names(geo))[2],
    size  = 2L,
    x     = -175L,
    y     = mean(
      c(
        unique(unlist(x = geo, use.names = FALSE))[2],
        unique(unlist(x = geo, use.names = FALSE))[3]
      )
    )
  ) +
  geom_text(
    angle = 90L,
    data  = world,
    color = "gray",
    hjust = .5,
    label = unique(names(geo))[3],
    size  = 2L,
    x     = -175L,
    y     = mean(
      c(
        unique(unlist(x = geo, use.names = FALSE))[3],
        unique(unlist(x = geo, use.names = FALSE))[4]
      )
    )
  ) +
  geom_text(
    angle = 90L,
    data  = world,
    color = "gray",
    hjust = .5,
    label = unique(names(geo))[2],
    size  = 2L,
    x     = -175L,
    y     = mean(
      c(
        unique(unlist(x = geo, use.names = FALSE))[4],
        unique(unlist(x = geo, use.names = FALSE))[5]
      )
    )
  ) +
  geom_text(
    angle = 90L,
    data  = world,
    color = "gray",
    hjust = .5,
    label = unique(names(geo))[1],
    size  = 2L,
    x     = -175L,
    y     = mean(
      c(
        unique(unlist(x = geo, use.names = FALSE))[5],
        unique(unlist(x = geo, use.names = FALSE))[6]
      )
    )
  ) +
  theme_light() +
  theme(
    axis.title      = element_blank(),
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    legend.key.size = unit(.2, plt$units),
    text            = element_text(size = plt$text)
  )

  # Save the plot
  ggsave(
    filename = tolower(
      paste(
        "9_res_map_of_",
        body,
        "_",
        cols,
        ".png",
        sep = ""
      )
    ),
    plot   = last_plot(),
    device = plt$device,
    path   = plt$path,
    scale  = plt$scale,
    height = plt$height,
    width  = plt$width,
    units  = plt$units,
    dpi    = plt$dpi
  )

} # End of the fn_plot function

# Combine the aircraft bodies and output variables to be plotted
mix <- expand.grid(
  body = names(bod),
  cols = cols$dif
)

mapply(
  FUN  = fn_plot,
  body = mix$body,
  cols = mix$cols,
  labs = rep(x = c("in m", "in p. p.", "in pax"), each = 2L)
)

# ==============================================================================
# 3.3 Summarize results by climate zone
# ==============================================================================

# ==============================================================================
# 3.3.1 Calculate LOESS values by climate zone
# ==============================================================================

# Create a new data table for summarizing by airport
dt_res_zon <- copy(dt_res)

# Summarize the data by group, then combine them
dt_res_zon <- rbind(
  # Zonal summary by group
  dt_res_zon[, lapply(X = .SD, FUN = mean),
    by      = c("year", "ssp", "zone", "type"),
    .SDcols = cols$bas
  ],
  # Global summary by group
  dt_res_zon[, zone := "Global"][, lapply(X = .SD, FUN = mean),
    by      = c("year", "ssp", "zone", "type"),
    .SDcols = cols$bas
  ]
)

# Create LOESS output variables from base values
dt_res_zon[, (cols$loe) := lapply(
    X   = .SD,
    FUN = function(x) {
      predict(loess(formula = x ~ year, span = .75, model = TRUE))
    }
  ),
  by      = c("ssp", "zone", "type"),
  .SDcols = cols$bas
]

# Save the data to disk
fwrite(
  x    = dt_res_zon,
  file = paste(dir$res, "dt_res_loess_values_by_zone.csv", sep = "/")
)

# ==============================================================================
# 3.3.2 Plot the results onto chronological lineplots
# ==============================================================================

# Order the zones so they display in alphabetical order
dt_res_zon[, zone := factor(zone,
    levels = sort(unique(levels(dt_res_zon[, zone])))
  )
]

# Declare independent variables for grouping
grp <- c("ssp", "zone")

# Create a function to plot results
fn_plot <- function(body, cols) {

  # Build the plot
  ggplot(
    data    = dt_res_zon[type == body],
    mapping = aes(
      x     = year,
      y     = dt_res_zon[type == body][[as.character(cols)]]
    )
  ) +
    geom_line(linewidth = .2) +
    geom_smooth(formula = y ~ x, method = "loess", linewidth = .5) +
    # Starting value labels
    geom_label(
      data    = dt_res_zon[type == body][, .SD[which.min(year)], by = grp],
      aes(
        x     = year,
        y     = dt_res_zon[type == body][, .SD[which.min(year)], by = grp]
                  [[paste(as.character(cols), "loe", sep = "_")]],
        label = sprintf(
                  fmt = ifelse(
                    cols == "avg_thr",
                    "%1.1f%%",
                    ifelse(cols == "avg_todr", "%1.0f", "%.1f")
                  ),
                  dt_res_zon[type == body][, .SD[which.min(year)], by = grp]
                    [[paste(as.character(cols), "loe", sep = "_")]] *
                    ifelse(cols == "avg_thr", 100L, 1L)
                )
      ),
      alpha      = .5,
      fill       = "white",
      label.r    = unit(0L, "lines"),
      label.size = 0L,
      nudge_x    = 4L,
      size       = 2L
    ) +
    # Ending value labels
    geom_label(
      data    = dt_res_zon[type == body][, .SD[which.max(year)], by = grp],
      aes(
        x     = year,
        y     = dt_res_zon[type == body][, .SD[which.max(year)], by = grp]
                  [[paste(as.character(cols), "loe", sep = "_")]],
        label = sprintf(
                  fmt = ifelse(
                    cols == "avg_thr",
                    "%1.1f%%",
                    ifelse(cols == "avg_todr", "%1.0f", "%.1f")
                  ),
                  dt_res_zon[type == body][, .SD[which.max(year)], by = grp]
                    [[paste(as.character(cols), "loe", sep = "_")]] *
                    ifelse(cols == "avg_thr", 100L, 1L)
                )
      ),
      alpha      = .5,
      fill       = "white",
      label.r    = unit(0L, "lines"),
      label.size = 0L,
      nudge_x    = -4L,
      size       = 2L
    ) +
    scale_x_continuous(name = "Year", n.breaks = 3L) +
    scale_y_continuous(
      name   = "Value",
      labels = ifelse(
        cols == "avg_thr",
        scales::label_percent(accuracy = .1),
        ifelse(
          cols == "avg_todr",
          scales::label_comma(accuracy = 1L),
          scales::label_comma(accuracy = .1)
        )
      )
    ) +
    facet_grid(
      rows   = vars(zone),
      cols   = vars(toupper(ssp)),
      scales = "free_y"
    ) +
    theme_light() +
    theme(
      axis.title.y = element_blank(),
      text         = element_text(size = plt$text)
    )

  # Save the plot
  ggsave(
    filename = tolower(
      paste(
        "9_res_lineplot_of_",
        body,
        "_",
        cols,
        ".png",
        sep = ""
      )
    ),
    plot     = last_plot(),
    device   = plt$device,
    path     = plt$path,
    scale    = plt$scale,
    height   = plt$height,
    width    = plt$width,
    units    = plt$units,
    dpi      = plt$dpi
  )

} # End of the fn_plot function

# Combine the aircraft bodies and output variables to be plotted
mix <- expand.grid(
  body = names(bod),
  cols = cols$bas
)

# Generate the plots
mapply(
  FUN  = fn_plot,
  body = mix$body,
  cols = mix$cols
)

# ==============================================================================
# 3.3.3 Summarize the changes in LOESS values by climate zone
# ==============================================================================

# Define change variable
cols$dif <- paste(cols$loe, "dif", sep = "_")

# Calculate the absolute difference between each year and the first, by zone
dt_res_zon[,
  (cols$dif) := lapply(
    X     = .SD,
    FUN   = function(x) {
      (x - x[1:1])
    }
  ),
  by      = c("ssp", "zone", "type"),
  .SDcols = cols$loe
]

# Save the final-year changes in LOESS values by zone to disk
fwrite(
  x = dt_res_zon[year == dt_res_zon[which.max(year), year]
  ][, c("year", "ssp", "zone", "type", cols$dif), with = FALSE
  ][, avg_thr_loe_dif := avg_thr_loe_dif * 100L # Change to percentage
  ][, (cols$dif) := round(.SD, 1L), .SDcols = cols$dif
  ][, melt(data = .SD, id.vars = c("year", "zone", "ssp", "type"))
  ][, dcast(data = .SD, formula = year + zone + variable ~ type + ssp)
  ][, zone := factor(zone, levels = sort(unique(levels(dt_res_zon[, zone]))))
  ][, variable := factor(variable, cols$dif)
  ][order(variable, zone)
  ],
  file = paste(dir$res, "dt_res_loess_changes_by_zone.csv", sep = "/")
)

# ==============================================================================
# 6 Housekeeping
# ==============================================================================

# Stop the script timer
Sys.time() - start_time

# EOF