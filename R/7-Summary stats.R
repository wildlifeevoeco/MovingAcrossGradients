## Summary Stats ====
# Authors: Quinn Webber, Michel Laforge, Maegwin Bonar, Alec Robitaille,
#           Chris Hart, Sana Zabihi-Seissan, Eric Vander Wal

### Packages ----
libs <- c('data.table')
lapply(libs, require, character.only = TRUE)

# Final DT
DT <- readRDS("data/derived-data/DT_final_brns.RDS")

# Cleaned locs
clean = readRDS("data/derived-data/cleaned-locs.Rds")


### Summary numbers for supplementary material ----
# number of unique individuals
length(unique(DT$ANIMAL_ID))

# unique individuals for each year by herd combination
DT[, uniqueN(ANIMAL_ID), by = .(HERD, year)]

# number of unique individuals in each herd
sum_herd <- DT[, uniqueN(ANIMAL_ID), by = .(HERD)]

# number of unique individuals in each year
sum_yr <- DT[, uniqueN(ANIMAL_ID), by = .(year)]


# mean proportion of 'good' pixels in HR
mean(DT$propGood)
sd(DT$propGood)
range(DT$propGood)


# count number of years each individual was tracked
sum_id_year  <- DT[, uniqueN(year), by = .(ANIMAL_ID)]

# mean +/- sd, min and max number of years tracked
mean(sum_id_year$V1)
sd(sum_id_year$V1)
min(sum_id_year$V1)
max(sum_id_year$V1)


# number of observations (in each block*ID*year) for each individual
count_id <- DT[, .N, by = .(ANIMAL_ID)]

# average +/-, min and max for
# number of observations per individual
mean(count_id$N)
sd(count_id$N)
min(count_id$N)
max(count_id$N)

# Herd/year specific Moran's I values
DT[, mean(moran), by = .(HERD)]
DT[, sd(moran), by = .(HERD)]


# Mean +/-, min and max FPT values
mean(DT$meanFPT)
sd(DT$meanFPT)
min(DT$meanFPT)
max(DT$meanFPT)

# Mean +/-, min and max KDE values
mean(DT$areaRatio)
sd(DT$areaRatio)
min(DT$areaRatio)
max(DT$areaRatio)


# herd/year specific
# mean +/- FPT & KDE values
DT[, mean(meanFPT), by = .(HERD)]
DT[, sd(meanFPT), by = .(HERD)]
DT[, mean(areaRatio), by = .(HERD)]
DT[, sd(areaRatio), by = .(HERD)]


# number of fixes for each herd*year combination
clean[, .N, by = .(HERD, year)]
clean[, uniqueN(ANIMAL_ID), by = .(HERD, year)]
clean[, uniqueN(ANIMAL_ID), by = .(HERD)]
clean[, uniqueN(ANIMAL_ID), by = .(year)]

# number of fixes for each unique herd
clean[, .N, by = .(HERD)]

# number of fixes for each unique year
clean[, .N, by = .(year)]

id_block <- clean[, .N, by = .(block, year, ANIMAL_ID)]
median(id_block$N)
range(id_block$N)


# fixes per day
id_day <- clean[, .N, by = .(JDate, year, ANIMAL_ID)]
median(id_day$N)
range(id_day$N)
