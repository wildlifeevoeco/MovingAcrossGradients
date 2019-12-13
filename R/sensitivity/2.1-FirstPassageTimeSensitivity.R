
### First Passage Time Sensitivity ====
# Authors: Quinn Webber, Michel Laforge, Maegwin Bonar, Chris Hart,
#           Alec Robitaille, Sana Zabihi-Seissan, Eric Vander Wal
# Inputs: Collar data
# Outputs: FPT

libs <- c('data.table', 'adehabitatLT','segmented', 'ggplot2', 'gridExtra')
lapply(libs, require, character.only = TRUE)

### Input raw data ----
locs <- readRDS('data/derived-data/cleaned-locs.Rds')
locs[,hour := hour(itime)]


### load first passage time 2 hour intervals
fpt.dt2hr <- readRDS('data/derived-data/first-passage-time.Rds')


colnames(fpt.dt2hr)[1] = "meanFPT2hr"

########### 4 hour intervals #############

locs4hr <- locs[hour == 0 | hour == 4 | hour == 8 | hour == 12 | hour == 16 | hour == 20]

### Sensitivity Analysis ----
# Supposed to select the 'peak' of the curve.
# We went with point where standard deviation overlaps with max point.
locs.traj4hr <- as.ltraj(xy = data.frame(locs4hr[, .(EASTING, NORTHING)]),
                      date = locs4hr[, datetime],
                      id = locs4hr[, ANIMAL_ID])

sensitivity.fpt4hr <- fpt(locs.traj4hr, seq(0,10000, length = 50))

sensitivity.var4hr <- varlogfpt(sensitivity.fpt4hr, graph = FALSE)

sensitivity4hr <- data.table(r = colnames(sensitivity.var4hr),
                          var = colMeans(sensitivity.var4hr, na.rm=TRUE),
                          seq = attributes(sensitivity.var4hr)$radii,
                          sd = apply(sensitivity.var4hr, 2, sd, na.rm=T))

saveRDS(sensitivity4hr, "data/derived-data/sensitivity/FPTsensitivity4hr.RDS")

### Plotting ----
## Segmented Linear model
sensLT5400_4hr <- sensitivity4hr[seq <= 4700] ## run model on all points < 4500m

m1_4hr <- lm(var ~ seq, data = sensLT5400_4hr)
mod1_4hr <- segmented.lm(m1_4hr, seg.Z = ~seq, data = sensLT5400_4hr)
summary(mod1_4hr)

saveRDS(mod1_4hr, "data/derived-data/sensitivity/FPTbrokenStickMod4hr.RDS")


### ARS ----
ARS <- function(in.dt,radii){
  # inputs: in.dt - dataframe to work with
  #         radii - radius of ARS zone
  locs.traj <- as.ltraj(xy = data.frame(in.dt[, .(EASTING, NORTHING)]),
                        date = in.dt[, datetime],
                        id = in.dt[, ANIMAL_ID],
                        typeII = TRUE)

  fpt <- fpt(locs.traj, radii, unit = 'hours')
  mn.fpt <- meanfpt(fpt, graph = FALSE)

  list(meanFPT = (mn.fpt$r1),
       ANIMAL_ID = (rownames(mn.fpt)))
}

radii_4hr <- mod1_4hr$psi.history[[5]][[1]]
mn.fpt.by.yearblock4hr <- locs4hr[, ARS(.SD, radii_4hr), by = .(year, block),
                            .SDcols = c('EASTING', 'NORTHING',
                                        'datetime', 'ANIMAL_ID')]

# Merge mean FPT back on
setkey(locs4hr, 'ANIMAL_ID')
setkey(mn.fpt.by.yearblock4hr, 'ANIMAL_ID')

locs4hr <- merge(locs4hr, mn.fpt.by.yearblock4hr,
                 by = c('ANIMAL_ID', 'block', 'year'))

# Need to convert NA values to 360 hours b/c these individuals never left radii.
locs4hr[is.na(meanFPT), meanFPT := 240]
locs4hr[meanFPT > 240, meanFPT := 240]

### Export ----
# Sub only unique rows of meanFPT and blockByIDYear
fpt.dt4hr <- unique(locs4hr[, .(meanFPT, blockidyr)])

colnames(fpt.dt4hr)[1] = "meanFPT4hr"

saveRDS(fpt.dt4hr, 'data/derived-data/sensitivity/first-passage-time4hr.Rds')


####### 8 hour intervals ##########

locs8hr <- locs[hour == 0 | hour == 8 | hour == 16]

### Sensitivity Analysis ----
# Supposed to select the 'peak' of the curve.
# We went with point where standard deviation overlaps with max point.
locs.traj8hr <- as.ltraj(xy = data.frame(locs8hr[, .(EASTING, NORTHING)]),
                      date = locs8hr[, datetime],
                      id = locs8hr[, ANIMAL_ID])

sensitivity.fpt8hr <- fpt(locs.traj8hr, seq(0,10000, length = 50))

sensitivity.var8hr <- varlogfpt(sensitivity.fpt8hr, graph = FALSE)

sensitivity8hr <- data.table(r = colnames(sensitivity.var8hr),
                          var = colMeans(sensitivity.var8hr, na.rm=TRUE),
                          seq = attributes(sensitivity.var8hr)$radii,
                          sd = apply(sensitivity.var8hr, 2, sd, na.rm=T))

saveRDS(sensitivity8hr, "data/derived-data/sensitivity/FPTsensitivity8hr.RDS")
### Plotting ----
## Segmented Linear model
sensLT5400_8hr <- sensitivity8hr[seq <= 4900] ## run model on all points < 4500m

m1_8hr <- lm(var ~ seq, data = sensLT5400_8hr)
mod1_8hr <- segmented.lm(m1_8hr, seg.Z = ~seq, data = sensLT5400_8hr)
summary(mod1_8hr)
saveRDS(mod1_8hr, "data/derived-data/sensitivity/FPTbrokenStickMod8hr.RDS")

radii_8hr <- mod1_8hr$psi.history[[5]][[1]]
mn.fpt.by.yearblock8hr <- locs8hr[, ARS(.SD, radii_8hr), by = .(year, block),
                               .SDcols = c('EASTING', 'NORTHING',
                                           'datetime', 'ANIMAL_ID')]

# Merge mean FPT back on
setkey(locs8hr, 'ANIMAL_ID')
setkey(mn.fpt.by.yearblock8hr, 'ANIMAL_ID')

locs8hr <- merge(locs8hr, mn.fpt.by.yearblock8hr,
                  by = c('ANIMAL_ID', 'block', 'year'))

# Need to convert NA values to 360 hours b/c these individuals never left radii.
locs8hr[is.na(meanFPT), meanFPT := 240]
locs8hr[meanFPT > 240, meanFPT := 240]

### Export ----
# Sub only unique rows of meanFPT and blockByIDYear
fpt.dt8hr <- unique(locs8hr[, .(meanFPT, blockidyr)])

colnames(fpt.dt8hr)[1] = "meanFPT8hr"

saveRDS(fpt.dt8hr, 'data/derived-data/sensitivity/first-passage-time8hr.Rds')


####### 12 hour intervals ##########

locs12hr <- locs[hour == 0 | hour == 12]

### Sensitivity Analysis ----
# Supposed to select the 'peak' of the curve.
# We went with point where standard deviation overlaps with max point.
locs.traj12hr <- as.ltraj(xy = data.frame(locs12hr[, .(EASTING, NORTHING)]),
                      date = locs12hr[, datetime],
                      id = locs12hr[, ANIMAL_ID])

sensitivity.fpt12hr <- fpt(locs.traj12hr, seq(0,10000, length = 50))

sensitivity.var12hr <- varlogfpt(sensitivity.fpt12hr, graph = FALSE)

sensitivity12hr <- data.table(r = colnames(sensitivity.var12hr),
                          var = colMeans(sensitivity.var12hr, na.rm=TRUE),
                          seq = attributes(sensitivity.var12hr)$radii,
                          sd = apply(sensitivity.var12hr, 2, sd, na.rm=T))
saveRDS(sensitivity12hr, "data/derived-data/sensitivity/FPTsensitivity12hr.RDS")

### Plotting ----
## Segmented Linear model
sensLT5400_12hr <- sensitivity12hr[seq <= 4900] ## run model on all points < 4500m

m1_12hr <- lm(var ~ seq, data = sensLT5400_12hr[2:23])
mod1_12hr <- segmented.lm(m1_12hr, seg.Z = ~seq, data = sensLT5400_12hr)
summary(mod1_12hr)
saveRDS(mod1_12hr, "data/derived-data/sensitivity/FPTbrokenStickMod12hr.RDS")


radii_12hr <- mod1_12hr$psi.history[[5]][[1]]
mn.fpt.by.yearblock12hr <- locs12hr[, ARS(.SD, radii_12hr), by = .(year, block),
                               .SDcols = c('EASTING', 'NORTHING',
                                           'datetime', 'ANIMAL_ID')]

# Merge mean FPT back on
setkey(locs12hr, 'ANIMAL_ID')
setkey(mn.fpt.by.yearblock12hr, 'ANIMAL_ID')

locs12hr <- merge(locs12hr, mn.fpt.by.yearblock12hr,
              by = c('ANIMAL_ID', 'block', 'year'))

# Need to convert NA values to 360 hours b/c these individuals never left radii.
locs12hr[is.na(meanFPT), meanFPT := 240]
locs12hr[meanFPT > 240, meanFPT := 240]

### Export ----
# Sub only unique rows of meanFPT and blockByIDYear
fpt.dt12hr <- unique(locs12hr[, .(meanFPT, blockidyr)])

colnames(fpt.dt12hr)[1] = "meanFPT12hr"

saveRDS(fpt.dt12hr, 'data/derived-data/sensitivity/first-passage-time12hr.Rds')


