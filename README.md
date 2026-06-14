# PhD
A fork of my previously-developed code for my doctoral research on simulating the impacts of climate change on civil aircraft takeoff performance. Now with enhanced features, including CMIP7 readiness, engine flat-rating thrust–temperature lapse, proper
infeasibility handling replacing the old `acc[acc < 0] <- 1e-3` clamp, a runway slope radians fix, replacing `hurs` with `huss` for the humidity variable, and using a parquet file for the climate data.
