### First Passage Time ====
# Authors: Quinn Webber, Michel Laforge, Maegwin Bonar, Alec Robitaille,
#           Chris Hart, Sana Zabihi-Seissan, Eric Vander Wal
# Inputs: Collar data
# Outputs: FPT

### Packages ---
libs <- c('data.table', 'adehabitatLT','segmented', 'ggplot2', 'gridExtra')
lapply(libs, require, character.only = TRUE)

### Set variables ----
source('R/variables.R')

### Input data ----
locs <- readRDS(paste0(derived, 'cleaned-locs.Rds'))

locs[, hour := hour(itime)]

### Sensitivity Analysis ----
# Select point where standard deviation overlaps with max point.
locs.traj <- as.ltraj(xy = data.frame(locs[, .(EASTING, NORTHING)]),
                      date = locs[, datetime],
                      id = locs[, ANIMAL_ID])

sensitivity.fpt <- fpt(locs.traj, seq(0, 10000, length = 50))

sensitivity.var <- varlogfpt(sensitivity.fpt, graph = FALSE)

sensitivity <- data.table(
  r = colnames(sensitivity.var),
  var = colMeans(sensitivity.var, na.rm = TRUE),
  seq = attributes(sensitivity.var)$radii,
  sd = apply(sensitivity.var, 2, sd, na.rm = T)
)

### Plotting ----
## Segmented Linear model
# run model on all points < 4500m
sensLT5400 <- sensitivity[seq <= 4500]

m1 <- lm(var ~ seq, data = sensLT5400[2:23])
mod1 <- segmented(m1, seg.Z = ~seq, data = sensLT5400)
summary(mod1)

mod1$psi.history[[5]][[1]]

saveRDS(mod1, "data/derived-data/FPT sensitivity/FPTbrokenStickMod1.RDS")
saveRDS(sensitivity, "data/derived-data/FPT sensitivity/FPTsensitivity.RDS")


### ARS ----
ARS <- function(in.dt, radii) {
  # inputs: in.dt - dataframe to work with
  #         radii - radius of ARS zone
  locs.traj <-
    as.ltraj(
      xy = data.frame(in.dt[, .(EASTING, NORTHING)]),
      date = in.dt[, datetime],
      id = in.dt[, ANIMAL_ID],
      typeII = TRUE
    )

  fpt <- fpt(locs.traj, radii, unit = 'hours')
  mn.fpt <- meanfpt(fpt, graph = FALSE)

  list(meanFPT = (mn.fpt$r1),
       ANIMAL_ID = (rownames(mn.fpt)))
}

radii <- mod1$psi.history[[5]][[1]]
mn.fpt.by.yearblock <- locs[, ARS(.SD, radii), by = .(year, block),
                            .SDcols = c('EASTING', 'NORTHING',
                                        'datetime', 'ANIMAL_ID')]

# Merge mean FPT back on
setkey(locs, 'ANIMAL_ID')
setkey(mn.fpt.by.yearblock, 'ANIMAL_ID')

locs <- merge(locs, mn.fpt.by.yearblock,
              by = c('ANIMAL_ID', 'block', 'year'))

# Need to convert NA values to 240 hours b/c these individuals never left radii.
locs[is.na(meanFPT), meanFPT := 240]
locs[meanFPT > 240, meanFPT := 240]


### Output ----
# Sub only unique rows of meanFPT and blockidyr
out <- unique(locs[, .(meanFPT, blockidyr)])

out[, fptScaled := scale(meanFPT)]

saveRDS(out, paste0(derived, 'first-passage-time.Rds'))

####
message('=== FPT COMPLETE ===')
