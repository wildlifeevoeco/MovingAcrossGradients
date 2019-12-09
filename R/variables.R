### Set options ----
options(scipen = 999)


### Set paths ----
raw <- 'data/raw-data/'
derived <- 'data/derived-data/'

### Set variables ----
# Variables for subsetting
sex <- 'F'
collarType <- 'GPS'
lowestYear <- 2002
herdList <- c('BUCHANS', 'GREY', 'LAPOILE',
							'MIDRIDGE', 'POTHILL', 'TOPSAILS')
doy <- c(161, 211)
tLims <- irg:::julseq[jul == doy[1] | jul == doy[2]]$t

# Filter by bounds
lowEast <- 0; highEast <- 800000
lowNorth <- 5200000; highNorth <- 6000000

# Max moverate
maxMoveRate <- 30000

# Block length
blockLength <- 10

# Filter by N locs by blockidyear and by day
nBlock <- 90
nDayBlock <- 6

### Set variables ----
# Kernel areas to return
kernelAreas <- c(50, 95)

# Percent for vertices
percent <- 95

# Regular grid cellsize
regularGridSize <- c(250, 250)

### Projection ----
projCols <- c('EASTING', 'NORTHING')

utm21N <- '+proj=utm +zone=21 ellps=WGS84'
