### Calculate IRG ====
# Alec Robitaille
# December 2018

### Packages ----
libs <- c('data.table', 'irg')
lapply(libs, require, character.only = TRUE)

### Set variables ----
source('R/variables.R')

### Input ----
DT <- rbindlist(lapply(1:3, function(i){
	fread(paste0(derived, 'reg-pts-NDVI-', i, '.csv'))
}))


### Calculate IRG for each point within each HR ----
# Filter NDVI time series
filter_qa(DT, qa = 'SummaryQA', good = c(0, 1))

filter_winter(DT, probs = 0.025, limits = c(60L, 300L),
							doy = 'DayOfYear', id = 'id')

filter_roll(DT, window = 3L, id = 'id', method = 'median')

filter_top(DT, probs = 0.925, id = 'id')


# Scale NDVI time series
scale_doy(DT)
scale_ndvi(DT)

## Drop NAs
DT <- DT[!is.na(scaled)]

# Guess starting parameters
model_start(DT)

# Double logistic model parameters given starting parameters for nls
mods <- model_params(
	DT,
	returns = 'models',
	xmidS = 'xmidS_start',
	xmidA = 'xmidA_start',
	scalS = 0.05,
	scalA = 0.01
)

## Fit double logistic curve to NDVI time series
modsNA <- na.omit(mods, invert = TRUE)

# Drop NAs and rep each row
modsDOY <- na.omit(mods)[rep(1:.N, each = diff(doy))]
#modsDOY[, jul := seq(doy[1], doy[2]-1, by = 1), by = yr]
modsDOY[, t := rep(irg:::julseq[doy[1]:(doy[2])]$t, length.out = .N)]

# Rescale the julian days
#scale_doy(modsDOY, 'jul')

model_ndvi(modsDOY, observed = TRUE)

# Calculate IRG for each day of the year
calc_irg(modsDOY)

### Output ----
saveRDS(fit, paste0(derived, 'patch-irg.Rds'))

###
message('=== CALC IRG COMPLETE ===')
