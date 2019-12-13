### Range use ====
# Authors: Quinn Webber, Michel Laforge, Maegwin Bonar, Alec Robitaille,
#           Chris Hart, Sana Zabihi-Seissan, Eric Vander Wal
# Inputs: Collar data
# Outputs: Vertices, range use

### Packages ----
libs <- c('data.table', 'sp', 'rgeos', 'rgdal')
lapply(libs, require, character.only = TRUE)

### Set variables ----
source('R/variables.R')
message(paste('=== RANGE USE PARAMETERS === \n',
							'kernel areas compared:',
							paste(kernelAreas, collapse = ' - '), '\n',
							'percent used for vertices:', percent, '\n',
							'grid size to determine patchiness of hr:',
							paste(regularGridSize, collapse = ' by '),
							'\n ================='))


### Input data ----
locs <- readRDS(paste0(derived, 'cleaned-locs.Rds'))


### HR by block ----
# For each herd, make each blockidyr's kernel
kerns <- lapply(locs[, unique(HERD)], function(h) {
	adehabitatHR::kernelUD(
		SpatialPointsDataFrame(
			coords = locs[HERD == h, .(EASTING, NORTHING)],
			data = locs[HERD == h, .(blockidyr)],
			proj4string = CRS(utm21N)
		)
	)
})

# For each herd's kernels, get the kernel area for specified areas
lsareas <- lapply(seq_along(kerns), function(k) {
	lapply(seq_along(kerns[[k]]), function(x) {
		tryCatch(
			data.table(t(
				adehabitatHR::kernel.area(kerns[[k]][[x]], kernelAreas)
			),
			blockidyr = names(kerns[[k]][x])),
			warning = function(w)
				data.table(warn = w[['message']],
									 blockidyr = names(kerns[[k]][x]))
		)
	})
})
areas <- rbindlist(unlist(lsareas, recursive = FALSE), fill = TRUE)
setnames(areas, as.character(kernelAreas), c('area50', 'area95'))

### Calculate Range Use ----
areas[, areaRatio := area50 / area95]
areas[, areaRatioScaled := scale(areaRatio)]

### Vertices ----
lsverts <- lapply(seq_along(kerns), function(k) {
	lapply(seq_along(kerns[[k]]), function(x) {
		tryCatch(
			adehabitatHR::getverticeshr(kerns[[k]][[x]], percent,
																	ida = names(kerns[[k]][x])),
			error = function(e) {
				list(err = e[['message']], nm = names(kerns[[k]][x]))
			}
		)
	})
})


## Regularly sample polygons
# Union vertices
verts <- lapply(lsverts, function(v) {
	noerrs <- sapply(v, function(i) is.null(i[['err']]))
	do.call(rbind, v[noerrs])
})

together <- gUnaryUnion(do.call(rbind, verts))


# Make a grid
grd <- SpatialPoints(makegrid(together, cellsize = regularGridSize),
										 proj4string = CRS(utm21N))

# Intersect the grid with the unioned verts
o <- over(grd, together)
reg <- grd[!is.na(o)]
reg <- SpatialPointsDataFrame(reg, data = data.frame(id = seq_along(reg)))
proj4string(reg) <- utm21N

# Which points match which verts
whichverts <- lapply(seq_along(verts), function(i) {
	gContains(verts[[i]], reg, byid = TRUE,
						returnDense = FALSE)
})

### Output ====
saveRDS(reg, paste0(derived, 'regular-pts.Rds'))
saveRDS(whichverts, paste0(derived, 'match-grid-verts.Rds'))
saveRDS(verts, paste0(derived, 'blockidyr-verts.Rds'))
saveRDS(areas, paste0(derived, 'areas.Rds'))

# 200K chunks
chunks <- data.table(from = c(1, 2e5, 4e5),
										 to = c(2e5, 4e5, nrow(reg)),
										 chunk = c(1, 2, 3))
chunks[, {
	writeOGR(reg[from:to,],
					 paste0('data/derived-data/reg-pts-', .BY[[1]]),
					 paste0('reg-pts-', .BY[[1]]),
					 driver = 'ESRI Shapefile', overwrite_layer = TRUE)
	1
}, by = chunk]

###
message('=== RANGE USE COMPLETE ===')
