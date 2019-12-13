### MCMC Reaction Norms ====
# Authors: Quinn Webber, Michel Laforge, Maegwin Bonar, Alec Robitaille,
#           Chris Hart, Sana Zabihi-Seissan, Eric Vander Wal

### Packages ----
libs <- c('data.table',
          'lme4','MCMCglmm',
          'tidyr', 'sqldf', 'gridExtra',
          'ggplot2', 'lmtest','dplyr')
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

# Set year as factor
DT[, year := factor(year)]


DT <- DT[!is.na(moran)]
DT[, moranScaled := scale(moran)]

saveRDS(DT, "data/derived-data/DT_final_brns.Rds")


################### MODEL SELECTION BIVARIATE MCMCglmm ###################

### Model 1: Herd
prior0 <- list(R=list(V=diag(12), nu=3),
               G=list(G1=list(V=diag(2), nu=3,
                              alpha.V=diag(var(DT$meanFPT),2,2))))

mcmc1 <- MCMCglmm(cbind(fptScaled, areaRatioScaled) ~ trait-1,
                         rcov =~ idh(trait:HERD):units,
                         family = c("gaussian","gaussian"),
                         #prior = prior0,
                         nitt=420000,
                         burnin=20000,
                         thin=100,
                         verbose = TRUE,
                         data = DT,
                         pr=T,saveX = TRUE,saveZ = TRUE)

saveRDS(mcmc1, "data/derived-data/MCMC-models/mod1.RDS")

### Model 2: Herd + Ind
prior1 <- list(R=list(V=diag(12), nu=3),
               G=list(G1=list(V=diag(2), nu=3,
                              alpha.V=diag(var(DT$meanFPT),2,2))))

mcmc2 <- MCMCglmm(cbind(fptScaled, areaRatioScaled) ~ trait-1,
                         random =~ us(trait):ANIMAL_ID,
                         rcov =~ idh(trait:HERD):units,
                         family = c("gaussian","gaussian"),
                         prior = prior1,
                         nitt=420000,
                         burnin=20000,
                         thin=100,
                         verbose = TRUE,
                         data = DT,
                         pr=T,saveX = TRUE,saveZ = TRUE)

saveRDS(mcmc2, "data/derived-data/MCMC-models/mod2.RDS")


### Model 3: Herd + Ind (with fixed effects)
mcmc3 <- MCMCglmm(cbind(fptScaled, areaRatioScaled) ~ trait-1 +
                           trait:moranScaled +
                           trait:year +
                           trait:HERD +
                           trait:block,
                         random =~ us(trait):ANIMAL_ID,
                         rcov =~ idh(trait:HERD):units,
                         family = c("gaussian","gaussian"),
                         prior = prior1,
                         nitt=420000,
                         burnin=20000,
                         thin=100,
                         verbose = TRUE,
                         data = DT,
                         pr=T,saveX = TRUE,saveZ = TRUE)

saveRDS(mcmc3, "data/derived-data/MCMC-models/mod3.RDS")

### Model 4: Herd + Ind (with fixed effects interaction)
mcmc4 <- MCMCglmm(cbind(fptScaled, areaRatioScaled) ~ trait-1 +
                           trait:moranScaled +
                           trait:year +
                           trait:HERD*trait:block,
                         random =~ us(trait):ANIMAL_ID,
                         rcov =~ idh(trait:HERD):units,
                         family = c("gaussian","gaussian"),
                         prior = prior1,
                         nitt=420000,
                         burnin=20000,
                         thin=100,
                         verbose = TRUE,
                         data = DT,
                         pr=T,saveX = TRUE,saveZ = TRUE)

saveRDS(mcmc4, "data/derived-data/MCMC-models/mod4.RDS")

### Model 5: Herd + Ind x Env (with fixed effects)
prior2 <- list(R=list(V=diag(12), nu=4),
               G=list(G1=list(V=diag(4), nu=4,
                              alpha.V=diag(var(DT$meanFPT),4,4))))

mcmc5 <- MCMCglmm(cbind(fptScaled, areaRatioScaled) ~ trait-1 +
                           trait:moranScaled +
                           trait:year +
                           trait:HERD +
                           trait:block,
                         random =~ us(trait + moranScaled:trait):ANIMAL_ID,
                         rcov =~ idh(trait:HERD):units,
                         family = c("gaussian","gaussian"),
                         prior = prior2,
												 nitt=420000,
												 burnin=20000,
												 thin=100,
                         verbose = TRUE,
                         data = DT,
                         pr=T,saveX = TRUE,saveZ = TRUE)

saveRDS(mcmc5, "data/derived-data/MCMC-models/mod5.RDS")

mcmc1 <- readRDS("data/derived-data/MCMC-models/mod1.RDS")
mcmc2 <- readRDS("data/derived-data/MCMC-models/mod2.RDS")
mcmc3 <- readRDS("data/derived-data/MCMC-models/mod3.RDS")
mcmc4 <- readRDS("data/derived-data/MCMC-models/mod4.RDS")
mcmc5 <- readRDS("data/derived-data/MCMC-models/mod5.RDS")



### Calculate Deviance Information Critera (DIC) ####
DIC <- data.table(DIC = c(mcmc1$DIC,
                       mcmc2$DIC, mcmc3$DIC,
                       mcmc4$DIC, mcmc5$DIC))

DIC$deltaDIC <- DIC$DIC - min(DIC$DIC)
DIC


### TABLE 2 ###
param_table <- data.table(model = c(1, 2, 3, 4, 5),
													deltaDIC = DIC$deltaDIC,
	      FPT = c(0, median(mcmc2$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"]),
	      				median(mcmc3$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"]),
	      				median(mcmc4$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"]),
	      				median(mcmc5$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"])),
				lowerFPT = c(0, HPDinterval(mcmc2$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"])[1],
										 HPDinterval(mcmc3$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"])[1],
										 HPDinterval(mcmc4$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"])[1],
										 HPDinterval(mcmc5$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"])[1]),
				upperFPT = c(0, HPDinterval(mcmc2$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"])[2],
										 HPDinterval(mcmc3$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"])[2],
										 HPDinterval(mcmc4$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"])[2],
										 HPDinterval(mcmc5$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"])[2]),
				KDE = c(0, median(mcmc2$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"]),
								median(mcmc3$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"]),
								median(mcmc4$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"]),
								median(mcmc5$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"])),
				lowerKDE = c(0, HPDinterval(mcmc2$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"])[1],
										 HPDinterval(mcmc3$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"])[1],
										 HPDinterval(mcmc4$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"])[1],
										 HPDinterval(mcmc5$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"])[1]),
				upperKDE = c(0, HPDinterval(mcmc2$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"])[2],
										 HPDinterval(mcmc3$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"])[2],
										 HPDinterval(mcmc4$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"])[2],
										 HPDinterval(mcmc5$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"])[2]),
				FPT_IxE = c(0, 0, 0, 0,
										median(mcmc5$VCV[,"traitfptScaled:moranScaled:traitfptScaled:moranScaled.ANIMAL_ID"])),
				FPT_IxE_lower = c(0, 0, 0, 0,
										HPDinterval(mcmc5$VCV[,"traitfptScaled:moranScaled:traitfptScaled:moranScaled.ANIMAL_ID"])[1]),
				FPT_IxE_upper = c(0, 0,0,0,
										HPDinterval(mcmc5$VCV[,"traitfptScaled:moranScaled:traitfptScaled:moranScaled.ANIMAL_ID"])[2]),
        KDE_IxE = c(0, 0, 0, 0,
						median(mcmc5$VCV[,"traitareaRatioScaled:moranScaled:traitareaRatioScaled:moranScaled.ANIMAL_ID"])),
				KDE_IxE_lower = c(0, 0, 0, 0,
									HPDinterval(mcmc5$VCV[,"traitareaRatioScaled:moranScaled:traitareaRatioScaled:moranScaled.ANIMAL_ID"])[1]),
				KDE_IxE_upper = c(0, 0,0,0,
									HPDinterval(mcmc5$VCV[,"traitareaRatioScaled:moranScaled:traitareaRatioScaled:moranScaled.ANIMAL_ID"])[2]))

param_table

## Intercept Repeatability estimates for each herd (FPT):
rep_FPT_buchans <- mcmc5$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"]/(
  mcmc5$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"] +
    mcmc5$VCV[,"traitfptScaled:HERDBUCHANS.units"])

rep_FPT_grey <- mcmc5$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"]/(
  mcmc5$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"] +
    mcmc5$VCV[,"traitfptScaled:HERDGREY.units"])

rep_FPT_lapoile <- mcmc5$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"]/(
  mcmc5$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"] +
    mcmc5$VCV[,"traitfptScaled:HERDLAPOILE.units"])

rep_FPT_middleridge <- mcmc5$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"]/(
  mcmc5$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"] +
    mcmc5$VCV[,"traitfptScaled:HERDMIDRIDGE.units"])

rep_FPT_pothill <- mcmc5$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"]/(
  mcmc5$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"] +
    mcmc5$VCV[,"traitfptScaled:HERDPOTHILL.units"])

rep_FPT_topsails <- mcmc5$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"]/(
  mcmc5$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"] +
    mcmc5$VCV[,"traitfptScaled:HERDTOPSAILS.units"])

rep_KDE_buchans <- mcmc5$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"]/(
  mcmc5$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"] +
    mcmc5$VCV[,"traitareaRatioScaled:HERDBUCHANS.units"])

rep_KDE_grey <- mcmc5$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"]/(
  mcmc5$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"] +
    mcmc5$VCV[,"traitareaRatioScaled:HERDGREY.units"])

rep_KDE_lapoile <- mcmc5$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"]/(
  mcmc5$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"] +
    mcmc5$VCV[,"traitareaRatioScaled:HERDLAPOILE.units"])

rep_KDE_middleridge <- mcmc5$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"]/(
  mcmc5$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"] +
    mcmc5$VCV[,"traitareaRatioScaled:HERDMIDRIDGE.units"])

rep_KDE_pothill <- mcmc5$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"]/(
  mcmc5$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"] +
    mcmc5$VCV[,"traitareaRatioScaled:HERDPOTHILL.units"])

rep_KDE_topsails <- mcmc5$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"]/(
  mcmc5$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"] +
    mcmc5$VCV[,"traitareaRatioScaled:HERDTOPSAILS.units"])

repeats <- data_frame(Traits = c("Pothill",
                                 "Topsails",
                                 "Lapoile",
                                 "Grey River",
                                 "Buchans",
                                 "Middleridge"),
                      Estimate_FPT = c(median(rep_FPT_pothill),median(rep_FPT_topsails),
                                       median(rep_FPT_lapoile),median(rep_FPT_grey),
                                       median(rep_FPT_buchans),median(rep_FPT_middleridge)),
                      Lower_FPT = c(HPDinterval(rep_FPT_pothill)[,"lower"],HPDinterval(rep_FPT_topsails)[,"lower"],
                                    HPDinterval(rep_FPT_lapoile)[,"lower"],HPDinterval(rep_FPT_grey)[,"lower"],
                                    HPDinterval(rep_FPT_buchans)[,"lower"],HPDinterval(rep_FPT_middleridge)[,"lower"]),
                      Upper_FPT = c(HPDinterval(rep_FPT_pothill)[,"upper"],HPDinterval(rep_FPT_topsails)[,"upper"],
                                    HPDinterval(rep_FPT_lapoile)[,"upper"],HPDinterval(rep_FPT_grey)[,"upper"],
                                    HPDinterval(rep_FPT_buchans)[,"upper"],HPDinterval(rep_FPT_middleridge)[,"upper"]),
                      Estimate_KDE = c(median(rep_KDE_pothill),median(rep_KDE_topsails),
                                       median(rep_KDE_lapoile),median(rep_KDE_grey),
                                       median(rep_KDE_buchans),median(rep_KDE_middleridge)),
                      Lower_KDE = c(HPDinterval(rep_KDE_pothill)[,"lower"],HPDinterval(rep_KDE_topsails)[,"lower"],
                                    HPDinterval(rep_KDE_lapoile)[,"lower"],HPDinterval(rep_KDE_grey)[,"lower"],
                                    HPDinterval(rep_KDE_buchans)[,"lower"],HPDinterval(rep_KDE_middleridge)[,"lower"]),
                      Upper_KDE = c(HPDinterval(rep_KDE_pothill)[,"upper"],HPDinterval(rep_KDE_topsails)[,"upper"],
                                    HPDinterval(rep_KDE_lapoile)[,"upper"],HPDinterval(rep_KDE_grey)[,"upper"],
                                    HPDinterval(rep_KDE_buchans)[,"upper"],HPDinterval(rep_KDE_middleridge)[,"upper"]),
											Vres_FPT = c(median(mcmc5$VCV[,"traitfptScaled:HERDBUCHANS.units"]),
																	 median(mcmc5$VCV[,"traitfptScaled:HERDGREY.units"]),
																	 median(mcmc5$VCV[,"traitfptScaled:HERDLAPOILE.units"]),
																	 median(mcmc5$VCV[,"traitfptScaled:HERDMIDRIDGE.units"]),
																	 median(mcmc5$VCV[,"traitfptScaled:HERDPOTHILL.units"]),
																	 median(mcmc5$VCV[,"traitfptScaled:HERDTOPSAILS.units"])),
											Vres_KDE = c(median(mcmc5$VCV[,"traitareaRatioScaled:HERDBUCHANS.units"]),
																   median(mcmc5$VCV[,"traitareaRatioScaled:HERDGREY.units"]),
																   median(mcmc5$VCV[,"traitareaRatioScaled:HERDLAPOILE.units"]),
																   median(mcmc5$VCV[,"traitareaRatioScaled:HERDMIDRIDGE.units"]),
																   median(mcmc5$VCV[,"traitareaRatioScaled:HERDPOTHILL.units"]),
																   median(mcmc5$VCV[,"traitareaRatioScaled:HERDTOPSAILS.units"])))

repeats

### average repeats across herds
mean(repeats$Estimate_FPT)
sd(repeats$Estimate_FPT)

mean(repeats$Estimate_KDE)
sd(repeats$Estimate_KDE)

##################  CORRELATIONS ##################
## CORRELATIONS between Intercept FPT and Intercept KDE:
mcmc_FPT_KDE_ints <- mcmc5$VCV[,"traitfptScaled:traitareaRatioScaled.ANIMAL_ID"]/
  (sqrt(mcmc5$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"])*
     sqrt(mcmc5$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"]))

## CORRELATIONS between Slope FPT and Intercept FPT:
mcmc_FPT_int_slope <- mcmc5$VCV[,"traitfptScaled:moranScaled:traitfptScaled.ANIMAL_ID"]/
  (sqrt(mcmc5$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"])*
     sqrt(mcmc5$VCV[,"traitfptScaled:moranScaled:traitfptScaled:moranScaled.ANIMAL_ID"]))

## CORRELATIONS between Slope KDE and Intercept KDE:
mcmc_KDE_int_slope <- mcmc5$VCV[,"traitareaRatioScaled:moranScaled:traitareaRatioScaled.ANIMAL_ID"]/
  (sqrt(mcmc5$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"])*
     sqrt(mcmc5$VCV[,"traitareaRatioScaled:moranScaled:traitareaRatioScaled:moranScaled.ANIMAL_ID"]))

## CORRELATIONS between Slope FPT and Slope KDE:
mcmc_KDE_FPT_slope <- mcmc5$VCV[,"traitareaRatioScaled:moranScaled:traitfptScaled:moranScaled.ANIMAL_ID"]/
  (sqrt(mcmc5$VCV[,"traitareaRatioScaled:moranScaled:traitareaRatioScaled:moranScaled.ANIMAL_ID"])*
     sqrt(mcmc5$VCV[,"traitfptScaled:moranScaled:traitfptScaled:moranScaled.ANIMAL_ID"]))

## CORRELATIONS between Slope FPT and Int KDE:
mcmc_KDEint_FPTslope <- mcmc5$VCV[,"traitareaRatioScaled:traitfptScaled:moranScaled.ANIMAL_ID"]/
  (sqrt(mcmc5$VCV[,"traitareaRatioScaled:traitareaRatioScaled.ANIMAL_ID"])*
     sqrt(mcmc5$VCV[,"traitfptScaled:moranScaled:traitfptScaled:moranScaled.ANIMAL_ID"]))

## CORRELATIONS between Slope KDE and Int FPT:
mcmc_KDEslope_FPTint <- mcmc5$VCV[,"traitfptScaled:traitareaRatioScaled:moranScaled.ANIMAL_ID"]/
  (sqrt(mcmc5$VCV[,"traitfptScaled:traitfptScaled.ANIMAL_ID"])*
     sqrt(mcmc5$VCV[,"traitareaRatioScaled:moranScaled:traitareaRatioScaled:moranScaled.ANIMAL_ID"]))

mcmc_cor <- data_frame(Traits = c("FPT Intercept, HR Intercept",
                                  "FPT Intercept, FPT Slope",
                                  "HR Intercept, HR Slope",
                                  "FPT Slope, HR Slope"),
                       Estimate = c(median(mcmc_FPT_KDE_ints),median(mcmc_FPT_int_slope),
                                    median(mcmc_KDE_int_slope),median(mcmc_KDE_FPT_slope)),
                       Lower = c(HPDinterval(mcmc_FPT_KDE_ints)[,"lower"],HPDinterval(mcmc_FPT_int_slope)[,"lower"],
                                 HPDinterval(mcmc_KDE_int_slope)[,"lower"],HPDinterval(mcmc_KDE_FPT_slope)[,"lower"]),
                       Upper = c(HPDinterval(mcmc_FPT_KDE_ints)[,"upper"],HPDinterval(mcmc_FPT_int_slope)[,"upper"],
                                 HPDinterval(mcmc_KDE_int_slope)[,"upper"],HPDinterval(mcmc_KDE_FPT_slope)[,"upper"]))

mcmc_cor
