### All Locs - Cleaning ====
# Authors: Quinn Webber, Michel Laforge, Maegwin Bonar, Alec Robitaille,
#           Chris Hart, Sana Zabihi-Seissan, Eric Vander Wal
# Inputs: Collar data
# Outputs: Prepped data, EE asset

### Packages ----
libs <- c('data.table', 'ggplot2', 'rgdal', 'spatsoc', 'toast')
lapply(libs, require, character.only = TRUE)

### Set variables ----
source('R/variables.R')


### Input raw data ----
dropCols <- c('SPECIES', 'EPSG_CODE', 'Map_Quality',
							'COLLAR_FILE_ID', 'EXCLUDE', 'VENDOR_CL', 'DOP',
							'NAV', 'VALIDATED', 'LOCQUAL', 'COLLAR_ID')

locs <- fread(paste0(raw, 'AllCaribouDataRaw.csv'), drop = dropCols)


### Date columns ----
locs[, idate := as.IDate(FIX_DATE)]
locs[, itime := as.ITime(FIX_TIME)]

locs[, year := year(idate)]
locs[, month := month(idate)]

locs[, JDate := yday(idate)]

### Subset ----
message(paste('=== SUBSETTING === \n', 'sex:', sex, '\n',
							'collar type:', collarType, '\n', 'year >=', lowestYear, '\n',
							'HERD in:', paste(herdList, collapse = ', '), '\n',
							'DOY between:', paste(doy, collapse = ', '),
							'\n ================='))

locs <- locs[SEX == sex & COLLAR_TYPE_CL == collarType &
						 	year >= lowestYear & HERD %in% herdList &
						 	JDate %between% doy]

### More date ----
# Datetime (after filter)
locs[, datetime := as.POSIXct(paste(idate, itime, sep = ' '))]

# UTM zone 21N
locs[, (projCols) := as.data.table(project(cbind(X_COORD, Y_COORD), utm21N))]

message(paste('=== BOUNDING BOX === \n',
							'min EASTING:', lowEast, '\n',
							'max EASTING:', highEast, '\n',
							'min NORTHING:', lowNorth, '\n',
							'max NORTHING:', highNorth, '\n ================='))

locs <- locs[(lowEast < EASTING & EASTING < highEast) &
						 	(lowNorth < NORTHING & EASTING < highNorth)]


### toast::step_length ----
step_length(
	locs,
	coords = projCols,
	time = 'datetime',
	splitBy = c('ANIMAL_ID', 'year'),
	moverate = TRUE
)

# Drop more than 30km/hr movements
message(paste('=== MOVEMENT RATE === \n',
							'max movement rate:', maxMoveRate, '\n ================='))
locs <- locs[moveRate < maxMoveRate]
# locs[moveRate > 15000]

### Assign Block IDs ----
# How many unique individual IDs by year and herd?
# locs[, .('numbIndividuals' = uniqueN(ANIMAL_ID)),
# 		 by = .(HERD, year)]

# Split up the season into blocks
message(paste('=== ASSIGN BLOCK === \n',
							'length of study period:', doy[2] - doy[1], '\n',
							'block length:', blockLength, '\n ================='))

group_times(locs, 'datetime', paste0(blockLength, ' days'))

# Assign a non-numerical descriptor for each block * ID * Year
locs[, blockidyr :=
		 	paste(.BY[[1]], .BY[[2]], .BY[[3]], sep = '_'),
		 by = .(block, ANIMAL_ID, year)]

# Determine the number of locs by animal*block
locs[, NumbByblockidyr := .N, by = .(blockidyr)]

# Determine the number of locs for each day
locs[, NumbByDayYear := .N, by = .(blockidyr, JDate)]

# Drop insufficient locs in a block and block days
message(paste('=== REQUIRE N LOCS === \n',
							'cutoff n by blockidyear:', nBlock, '\n',
							'median n blockidyr:', median(locs$NumbByblockidyr), '\n',
							'range n blockidyr:',
							paste(range(locs$NumbByblockidyr), collapse = '-'), '\n',
							'cutoff n by day for each blockidyr:', nDayBlock, '\n',
							'median n by day block:', median(locs$NumbByDayYear), '\n',
							'range n by day block:',
							paste(range(locs$NumbByDayYear), collapse = '-'), '\n ================='
	)
)



locs <- locs[NumbByblockidyr > nBlock]
locs <- locs[!(blockidyr %in%
							 	locs[NumbByDayYear < nDayBlock,
							 			 unique(blockidyr)])]

### Export it ----
locs[, ptId := .I]

info <- unique(locs[, .(blockidyr, block, ANIMAL_ID, year, HERD)])
blockLimits <- locs[, .(minDOY = min(JDate), maxDOY = max(JDate)),
										by = block]

saveRDS(locs, paste0(derived, 'cleaned-locs.Rds'))
saveRDS(info, paste0(derived, 'info-blockidyr.Rds'))
saveRDS(blockLimits, paste0(derived, 'info-block-limits.Rds'))

message('=== PREP COMPLETE ===')
