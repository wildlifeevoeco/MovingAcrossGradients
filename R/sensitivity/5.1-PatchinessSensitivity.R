### Patchiness Sensitivity ====
# Alec Robitaille
# December 2018

### Packages ----
libs <- c('data.table', 'sp', 'raster', 'irg')
lapply(libs, require, character.only = TRUE)

### Set variables ----
source('R/variables.R')

### Input ----
locs <- readRDS(paste0(derived, 'cleaned-locs.Rds'))

fits <- readRDS(paste0(derived, 'patch-irg.Rds'))

matchverts <- readRDS(paste0(derived, 'match-grid-verts.Rds'))
blocklimits <- readRDS(paste0(derived, 'info-block-limits.Rds'))
blockids <- readRDS(paste0(derived, 'info-blockidyr.Rds'))
reg <- readRDS(paste0(derived, 'regular-pts.Rds'))

### Merge nearest points ----
# Add x, y back on
regDT <- as.data.table(reg)
setnames(regDT, c('x1', 'x2'), c('x', 'y'))
fits[regDT, c('x', 'y') := .(x, y), on = 'id']


### Moran's ----
# Reduce full regpts IRG data down to blockidyr values
trans <- Reduce(c, matchverts)

# Rescale doy in blocklimits
alloc.col(blocklimits)
scale_doy(blocklimits, doy = 'minDOY')[, tMin := t]
scale_doy(blocklimits, doy = 'maxDOY')[, tMax := t][, t := NULL]

# List points within each blockherdyear
uniquepts <- blockids[, .(pts = {
	nms <- names(trans)[names(trans) %in% blockidyr]
	pts <- trans[names(trans) %in% nms]
	list(unique(unlist(pts)))
}), by = .(block, HERD, year)]

uniquepts[, blockherdyear := paste(block, HERD, year, sep = '-')]

blockids[, blockherdyear := paste(block, HERD, year, sep = '-')]

allQuantiles <- lapply(seq(0.50, 0.95, by = 0.05), function(probs) {
	# Generate "good/bad" pixels for each blockherdyear
	pixels <- uniquepts[, {
		bhy <- .BY[[1]]
		blk <- block
		blocklims <- blocklimits[block == blk]
		selpts <- pts[[1]]
		DT <- fits[id %in% selpts][yr == year &
															 	between(t, blocklims$tMin, blocklims$tMax)]
		if (nrow(DT) > 0) {
			# Calculate mean IRG for each pixel
			DT[, meanIRG := mean(irg), id]

			# Calculate 0.8 quantile
			DT[, q := quantile(meanIRG, probs = probs, na.rm = TRUE)]

			# Binarize IRG compared to quantile
			DT[, good := meanIRG > q]

			unique(DT[, .(id, x, y, meanIRG, good, q)])
		}

	}, by = blockherdyear]

	# Calculate Morans
	blockids[, blockherdyear := paste(block, HERD, year, sep = '-')]

	morans <- blockids[, {
		biy <- .BY[[1]]
		bhy <- blockherdyear

		p <- pixels[blockherdyear == bhy &
									id %in% unlist(trans[names(trans) == biy])][!is.na(good)]

		out <- list(npixels = p[, .N],
								blockherdyear = blockherdyear,
								propGood = p[, mean(good, na.rm = TRUE)],
								probs = probs)


		if (out[['npixels']] == 0) {
			append(out,
						 list(
						 	numbNA = NA_integer_,
						 	numbNonNA = NA_integer_,
						 	tInPatch = NA_real_,
						 	moran = NA_real_,
						 	quantile = NA_real_,
						 	error = 'no pixels'
						 ))
		} else {
			r <- raster(
				SpatialPixelsDataFrame(p[, .(x, y)],
															 p[, .(good)],
															 proj4string = CRS(utm21N))
			)

			# Observed
			obs <- locs[blockidyr == biy, .(blockidyr, good = extract(r, matrix(c(EASTING, NORTHING), ncol = 2)))]
			numbNA <- as.integer(obs[, sum(is.na(good), na.rm = TRUE)])
			numbNonNA <- as.integer(obs[, .N - numbNA])
			tInPatch <- as.double(obs[, sum(good, na.rm = TRUE) / numbNonNA])


			out <- append(out,
										list(
											numbNA = numbNA,
											numbNonNA = numbNonNA,
											tInPatch = tInPatch
										))

			if (nrow(p) < 9) {
				append(out,
							 list(
							 	moran = NA_real_,
							 	quantile = NA_real_,
							 	error = 'less than 9'
							 ))
			} else if (nrow(p) >= 9  &
								 p[meanIRG < q, .N] == out[['npixels']]) {
				append(out,
							 list(
							 	moran = NA_real_,
							 	quantile = NA_real_,
							 	error = 'nothing in good patch'
							 ))
			} else if (nrow(p) >= 9 & p[meanIRG < q, .N] < out[['npixels']]) {
				append(out, list(
					moran = as.double(Moran(r)),
					quantile = as.double(p[, unique(q)]),
					error = NA_character_
				))
			}
		}
	}, by = blockidyr]
})

out <- rbindlist(allQuantiles)


### Output ----
saveRDS(out, paste0(derived, 'quantile-sensitivity.Rds'))


####
message('=== PATCHINESS COMPLETE ===')
