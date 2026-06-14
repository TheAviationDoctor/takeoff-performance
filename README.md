# Climatic simulation of takeoff performance
A fork of my previously-developed code for my doctoral research on simulating the impacts of climate change on civil aircraft takeoff performance. Now with enhanced features, including CMIP7 readiness, engine flat-rating thrust–temperature lapse, proper
infeasibility handling replacing the old `acc[acc < 0] <- 1e-3` clamp, a runway slope radians fix, replacing `hurs` with `huss` for the humidity variable, and using a parquet file for the climate data.

## Reproducibility (renv)

R package versions are pinned with [`renv`](https://rstudio.github.io/renv/). After
cloning, restore the exact library recorded in `renv.lock`:

```r
install.packages("renv")   # if renv is not already available
renv::restore()            # installs all pinned packages into renv/library/
```

`renv` activates automatically for new R sessions via `.Rprofile`. The project
library (`renv/library/`) is not committed; it is rebuilt from `renv.lock`.

Four packages are no longer on CRAN for current R versions and are restored from
the [CRAN archive](https://cran.r-project.org/src/contrib/Archive/): `PCICt`,
`ncdf4.helpers`, `epwshiftr`, and `rgeos`. `renv::restore()` fetches them from
the archive automatically; if that fails on your platform, install them manually,
e.g. `renv::install("rgeos@0.6-4")`. Note that `ncdf4` requires the system NetCDF
library and `RMySQL` requires the MySQL/MariaDB client library to be present.
