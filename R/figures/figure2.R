### Figure 2 - Heatmap ====
# Authors: Quinn Webber, Michel Laforge, Maegwin Bonar, Chris Hart,
#           Alec Robitaille, Sana Zabihi-Seissan, Eric Vander Wal


### Packages ----
libs <- c('raster', 'lme4', 'piecewiseSEM',
					'data.table', 'ggplot2')
lapply(libs, require, character.only = TRUE)


### Set variables ----
source('R/variables.R')


### Load data ----
fpt <- readRDS(paste0(derived, 'first-passage-time.Rds'))
range <- readRDS(paste0(derived, 'areas.Rds'))
info <- readRDS(paste0(derived, 'info-blockidyr.Rds'))
patchiness <- readRDS(paste0(derived, 'patchiness.Rds'))

# Merge the fpt, moran and RSF scores
DT <- Reduce(function(x, y) x[y, on = "blockidyr"],
						 list(fpt, patchiness, info, range))

### Prep ----
# Drop NAs
DT <- DT[!is.na(moran)]

# Scale morans
DT[, moranScaled := scale(moran)]

# Cast year as factor
DT[, year := factor(year)]

# Cast block as factor
DT[, block := factor(block)]


### Models ----
m1 <- lmer(tInPatch ~ fptScaled * moranScaled + areaRatioScaled * moranScaled +
			 	(1 | ANIMAL_ID),
			 data = DT)

rsquared(m1)


# Make heatmap dataframes
HeatmapdataFPT <- data.frame(
	x = rep(seq(min(DT$moranScaled, na.rm = TRUE),
							max(DT$moranScaled, na.rm = TRUE),
							by = ((max(DT$moranScaled, na.rm = TRUE) - min(DT$moranScaled, na.rm = TRUE))/200)), 201),
	y = rep(seq(min(DT$fptScaled, na.rm = TRUE),
							max(DT$fptScaled, na.rm = TRUE),
							by = ((max(DT$fptScaled, na.rm = TRUE) - min(DT$fptScaled, na.rm = TRUE))/200)), each = 201))

HeatmapdataKDE <- data.frame(
	x = rep(seq(min(DT$moranScaled, na.rm = TRUE),
							max(DT$moranScaled, na.rm = TRUE),
							by = ((max(DT$moranScaled, na.rm = TRUE) - min(DT$moranScaled, na.rm = TRUE))/200)), 201),
	y = rep(seq(min(DT$areaRatioScaled, na.rm = TRUE),
							max(DT$areaRatioScaled, na.rm = TRUE),
							by = ((max(DT$areaRatioScaled, na.rm = TRUE) - min(DT$areaRatioScaled, na.rm = TRUE))/200)), each = 201))


# Model
HeatmapdataFPT$z <- ((fixef(m1)[1]) + (fixef(m1)[3] * HeatmapdataFPT$x) +
										 	(fixef(m1)[4] * mean(DT$areaRatioScaled)) +
										 	(fixef(m1)[2] * HeatmapdataFPT$y) +
										 	(fixef(m1)[5] * (HeatmapdataFPT$x) * (HeatmapdataFPT$y)) +
										 	(fixef(m1)[6] * (HeatmapdataFPT$x) * mean(DT$areaRatioScaled)))

HeatmapdataKDE$z <- ((fixef(m1)[1]) + (fixef(m1)[3] * HeatmapdataKDE$x) +
										 	(fixef(m1)[2] * mean(DT$fptScaled)) +
										 	(fixef(m1)[4] * HeatmapdataKDE$y) +
										 	(fixef(m1)[6] * (HeatmapdataKDE$x) * (HeatmapdataKDE$y)) +
										 	(fixef(m1)[5] * (HeatmapdataKDE$x) * mean(DT$fptScaled)))


HeatmapdataFPT$xnew <- (HeatmapdataFPT$x - (min(HeatmapdataFPT$x)))/(max(HeatmapdataFPT$x) - min(HeatmapdataFPT$x))
HeatmapdataFPT$ynew <- (HeatmapdataFPT$y - (min(HeatmapdataFPT$y)))/(max(HeatmapdataFPT$y) - min(HeatmapdataFPT$y))

HeatmapdataKDE$xnew <- (HeatmapdataKDE$x - (min(HeatmapdataKDE$x)))/(max(HeatmapdataKDE$x) - min(HeatmapdataKDE$x))
HeatmapdataKDE$ynew <- (HeatmapdataKDE$y - (min(HeatmapdataKDE$y)))/(max(HeatmapdataKDE$y) - min(HeatmapdataKDE$y))

RastdataFPT <- cbind(HeatmapdataFPT$xnew, HeatmapdataFPT$ynew, HeatmapdataFPT$z)
RastdataKDE <- cbind(HeatmapdataKDE$xnew, HeatmapdataKDE$ynew, HeatmapdataKDE$z)

HeatmapFPT <- rasterFromXYZ(RastdataFPT)
HeatmapKDE <- rasterFromXYZ(RastdataKDE)


# Labels
labels<-seq(-3,3,1)
moranLoc<-(labels-(min(DT$moranScaled)))/(max(DT$moranScaled)-min(DT$moranScaled))

FPTLoc<-(labels-(min(DT$fptScaled)))/(max(DT$fptScaled)-min(DT$fptScaled))

KDELoc<-(labels-(min(DT$areaRatioScaled)))/(max(DT$areaRatioScaled)-min(DT$areaRatioScaled))



# Contour lines
ContKDE <- rasterToContour(HeatmapKDE, nlevels = 5)
ContFPT <- rasterToContour(HeatmapFPT, nlevels = 5)

# Palette
rbPal <- colorRampPalette(c('#ffffe5','#f7fcb9','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#006837','#004529'))(75)

### Output ----
pdf("graphics/Figures/Fig2_Heatmap.pdf",
		useDingbats = FALSE, height = 3.5, width = 7)
par(mar = c(3.5, 3, 1, 0), mfrow = c(1, 2))
legend_im <- as.raster(matrix(rev(rbPal), ncol = 1))
image(HeatmapFPT, col = rbPal, axes = FALSE, xlim = c(0, 1.15),xlab=NA,
			zlim = c(min(c(minValue(HeatmapKDE), minValue(HeatmapFPT))),
							 max(c(maxValue(HeatmapKDE), maxValue(HeatmapFPT)))))
lines(ContFPT)
axis(1,at=moranLoc,labels=labels,cex.axis=0.8)
segments(0,0,1,0)
axis(2,at=FPTLoc,labels=labels,cex.axis=0.8, las = 1)
axis(2,at=c(-100,100),labels=NA,cex.axis=0.8)
mtext("Forage patch aggregation", side = 1, padj = 3, adj = 0.4)
mtext("First-passage time", side = 2, padj = 2, outer = T)
mtext("A", side = 3,  at = 0.05)
par(mar = c(3.5, 1.5, 1, 0))
image(HeatmapKDE, col = rbPal, axes = FALSE, xlim = c(0, 1.25), xlab=NA,
			zlim = c(min(c(minValue(HeatmapKDE), minValue(HeatmapFPT))),
							 max(c(maxValue(HeatmapKDE), maxValue(HeatmapFPT)))))
lines(ContKDE)
axis(1,at=moranLoc,labels=labels,cex.axis=0.8)
axis(2,at=KDELoc,labels=labels,cex.axis=0.8, las = 1)
axis(2,at=c(-100,100),labels=NA,cex.axis=0.8)
segments(0,0,1,0)
mtext("Forage patch aggregation", side = 1, padj = 3, adj = 0.4)
mtext("Range-use ratio", side = 2, padj = 28.4, outer = T)
rasterImage(legend_im, 1.02, 0.2, 1.09, 0.8)
mtext("B", side = 3,  at = 0.05)
text(x = 1.18, y = 0.8, round(max(c(maxValue(HeatmapKDE), maxValue(HeatmapFPT))), digits = 2))
text(x = 1.18, y = 0.2, round(min(c(minValue(HeatmapKDE), minValue(HeatmapFPT))), digits = 2))
dev.off()


saveRDS(m1, 'data/derived-data/figure2-m1.Rds')
