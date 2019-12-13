### Heatmap Sensitivity ====
# Authors: Quinn Webber, Michel Laforge, Maegwin Bonar, Chris Hart,
#           Alec Robitaille, Sana Zabihi-Seissan, Eric Vander Wal


### Packages ----
libs <- c('raster', 'lme4', 'piecewiseSEM',
          'data.table', 'ggplot2')
lapply(libs, require, character.only = TRUE)


### Set variables ----
source('R/variables.R')


### Input ----
Sens <- readRDS(paste0(derived, 'quantile-sensitivity.Rds'))

fpt <- readRDS(paste0(derived, 'first-passage-time.Rds'))
range <- readRDS(paste0(derived, 'areas.Rds'))
info <- readRDS(paste0(derived, 'info-blockidyr.Rds'))
patchiness <- readRDS(paste0(derived, 'patchiness.Rds'))

# Merge the fpt, moran and RSF scores
DT <- Reduce(function(x, y) x[y, on = "blockidyr"],
						 list(fpt, patchiness, info, range))


### Prep ----
# Drop where moran is NA
DT <- DT[!is.na(moran)]

# Set year as factor
DT[, year := factor(year)]

# Cast block as factor
DT[, block := factor(block)]

# Scale moran
DT[, moranScaled := scale(moran)]

### By percent ----
byPercent <- lapply(seq(0.50, 0.95, by = 0.05), function(prb) {
	sub <- na.omit(Sens[probs == prb],
											col = 'moran')

	probDT <- Reduce(function(x, y) x[y, on = 'blockidyr'],
									 list(fpt, range, info, sub))

	probDT[, moranScaled := scale(moran)]

	heat <-
		lmer(tInPatch ~ fptScaled * moranScaled + areaRatioScaled * moranScaled +
				 	(1 | ANIMAL_ID),
				 data = probDT)

	fpt <- data.frame(
		x = rep(seq(min(probDT$moranScaled, na.rm = TRUE),
								max(probDT$moranScaled, na.rm = TRUE),
								by = ((max(probDT$moranScaled, na.rm = TRUE) - min(probDT$moranScaled, na.rm = TRUE))/200)), 201),
		y = rep(seq(min(probDT$fptScaled, na.rm = TRUE),
								max(probDT$fptScaled, na.rm = TRUE),
								by = ((max(probDT$fptScaled, na.rm = TRUE) - min(probDT$fptScaled, na.rm = TRUE))/200)), each = 201))

	kde <- data.frame(
		x = rep(seq(min(probDT$moranScaled, na.rm = TRUE),
								max(probDT$moranScaled, na.rm = TRUE),
								by = ((max(probDT$moranScaled, na.rm = TRUE) - min(probDT$moranScaled, na.rm = TRUE))/200)), 201),
		y = rep(seq(min(probDT$areaRatioScaled, na.rm = TRUE),
								max(probDT$areaRatioScaled, na.rm = TRUE),
								by = ((max(probDT$areaRatioScaled, na.rm = TRUE) - min(probDT$areaRatioScaled, na.rm = TRUE))/200)), each = 201))

	fpt$z <- ((fixef(heat)[1]) +
							(fixef(heat)[3] * fpt$x) +
							(fixef(heat)[4] * mean(probDT$areaRatioScaled)) +
							(fixef(heat)[2] * fpt$y) +
							(fixef(heat)[5] * (fpt$x) * (fpt$y)) +
							(fixef(heat)[6] * (fpt$x) * mean(probDT$areaRatioScaled)))

	kde$z <- ((fixef(heat)[1]) +
							(fixef(heat)[3] * kde$x) +
							(fixef(heat)[2] * mean(probDT$fptScaled)) +
							(fixef(heat)[4] * kde$y) +
							(fixef(heat)[6] * (kde$x) * (kde$y)) +
							(fixef(heat)[5] * (kde$x) * mean(probDT$fptScaled)))

	fpt$xnew <-
		(fpt$x - (min(fpt$x))) / (max(fpt$x) - min(fpt$x))
	fpt$ynew <-
		(fpt$y - (min(fpt$y))) / (max(fpt$y) - min(fpt$y))

	kde$xnew <-
		(kde$x - (min(kde$x))) / (max(kde$x) - min(kde$x))
	kde$ynew <-
		(kde$y - (min(kde$y))) / (max(kde$y) - min(kde$y))

	rastFPT <- cbind(fpt$xnew, fpt$ynew, fpt$z)
	rastKDE <- cbind(kde$xnew, kde$ynew, kde$z)

	heatFPT <- rasterFromXYZ(rastFPT)
	heatKDE <- rasterFromXYZ(rastKDE)

	contKDE <- rasterToContour(heatKDE, nlevels = 5)
	contFPT <- rasterToContour(heatFPT, nlevels = 5)

	rbPal <- colorRampPalette(c('#ffffe5', '#f7fcb9', '#d9f0a3',
															'#addd8e', '#78c679', '#41ab5d',
															'#238443', '#006837', '#004529'))(75)

	png(paste0('graphics/Supplement3/HeatSensPieces/FPT', prb * 100, '.png'),
		width = 3.15,
		height = 3.15,
		units = 'cm',
		res = 600
	)
	par(mar = c(0, 0, 0, 0))
	image(heatFPT, col = rbPal,
				zlim = c(min(c(minValue(heatKDE), minValue(heatFPT))),
								 max(c(maxValue(heatKDE), maxValue(heatFPT)))))
	lines(contFPT)
	dev.off()

	png(paste0('graphics/Supplement3/HeatSensPieces/KDE', prb * 100, '.png'),
		width = 3.15,
		height = 3.15,
		units = 'cm',
		res = 600
	)
	par(mar = c(0, 0, 0, 0))
	image(heatKDE, col = rbPal,
				zlim = c(min(c(minValue(heatKDE), minValue(heatFPT))),
								 max(c(maxValue(heatKDE), maxValue(heatFPT)))))
	lines(contKDE)
	dev.off()
})

### Legend ----
rbPal <- colorRampPalette(c('#ffffe5', '#f7fcb9', '#d9f0a3',
														'#addd8e', '#78c679', '#41ab5d',
														'#238443', '#006837', '#004529'))

legend_im <- as.raster(matrix(rev(rbPal(20)), ncol = 1))

png(
	'graphics/Supplement3/HeatSensPieces/Legend.png',
	height = 7.3,
	width = 1,
	units = 'cm',
	res = 600
)

par(mar = c(0, 0, 0, 0))
plot(legend_im)
dev.off()
