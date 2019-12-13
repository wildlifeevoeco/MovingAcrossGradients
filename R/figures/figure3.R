### Figure 3 - Reaction Norm Figures ====
# Authors: Quinn Webber, Michel Laforge, Maegwin Bonar, Chris Hart,
#           Alec Robitaille, Sana Zabihi-Seissan, Eric Vander Wal

### Packages ----
libs <- c('data.table',
					'MCMCglmm',
					'tidyr', 'gridExtra',
					'ggplot2', 'dplyr')
lapply(libs, require, character.only = TRUE)


### Data ----
DT <- readRDS("data/derived-data/DT_final_brns.Rds")


### Reaction Norm Figures ----
df_brns <- tibble(ID = DT$ANIMAL_ID,
											KDE = DT$areaRatio,
											FPT = DT$meanFPT,
											moranScaled = DT$moranScaled,
											Year = DT$year,
											block = DT$block,
											herd = DT$HERD)


## run models for Figure 3 and export for figure script
p.var_FPT <- var(df_brns$FPT,na.rm=TRUE)

prior_FPT <- list(G=list(G1=list(V=diag(2)*(p.var_FPT/2), nu=1,
															 alpha.V=diag(2)*p.var_FPT/2)),
							  	R=list(V=diag(1)*(p.var_FPT/2), nu=1))

mcmcFPT <- MCMCglmm(scale(FPT) ~ scale(KDE) + moranScaled +
											Year + herd + block,
										random =~ us(1 + moranScaled):ID,
										rcov = ~units,
										family = "gaussian",
										prior = prior_FPT,
										nitt=420000,
										burnin=20000,
										thin=100,
										verbose = TRUE,
										data = df_brns,
										pr=TRUE,
										saveX = TRUE,
										saveZ = TRUE)

saveRDS(mcmcFPT, "data/derived-data/mcmc/mod_FPT.RDS")

### Reaction Norm for KDE ----
p.var_KDE <- var(df_brns$KDE, na.rm = TRUE)

prior_KDE <- list(G=list(G1=list(V=diag(2)*(p.var_KDE/2), nu=1,
															 alpha.V=diag(2)*p.var_KDE/2)),
							  	R=list(V=diag(1)*(p.var_KDE/2), nu=1))

mcmcKDE <- MCMCglmm(scale(KDE) ~ scale(FPT) + moranScaled +
											Year + herd + block,
										random =~ us(1 + moranScaled):ID,
										rcov = ~units,
										family = "gaussian",
										prior = prior_KDE,
										nitt=420000,
										burnin=20000,
										thin=100,
										verbose = TRUE,
										data = df_brns,
										pr=TRUE,
										saveX = TRUE,
										saveZ = TRUE)

saveRDS(mcmcKDE, "data/derived-data/mcmc/mod_KDE.RDS")

### load MCMC files to generate figure
mcmcFPT <- readRDS("data/derived-data/mcmc/mod_FPT.RDS")
mcmcKDE <- readRDS("data/derived-data/mcmc/mod_KDE.RDS")

# Bind predictions of the random intercepts model
# and random regression model to the original data;
# Select variables of interest;
# Get mean values at each opp_size for each individual
# (within each 'type' of data), so that we are
# essentially averaging across block effects;
# Convert to a 'long' format for plotting.
## FPT
df_fpt <- cbind(df_brns,
								fit = predict(mcmcFPT, marginal = NULL)) %>%
	group_by(ID, herd, moranScaled) %>%
	summarise(fit = mean(fit),
						FPT = mean(FPT)) %>%
	gather(Type, Value,
				 fit:FPT)

df_fit_fpt = setDT(df_fpt)[Type == "fit"]

saveRDS(df_fit_fpt, "data/derived-data/mcmc/df_model_fit_fpt.RDS")

## KDE
df_kde <- cbind(df_brns,
								fit = predict(mcmcKDE, marginal = NULL)) %>%
	group_by(ID, herd, moranScaled) %>%
	summarise(fit = mean(fit),
						KDE = mean(KDE)) %>%
	gather(Type, Value,
				 fit:KDE)

df_fit_kde = setDT(df_kde)[Type == "fit"]

saveRDS(df_fit_kde, "data/derived-data/mcmc/df_model_fit_kde.RDS")

# extract slopes and intsercepts for population level
intFPT = lm(Value ~ moranScaled, data = df_fit_fpt)$coefficients[1]
slopeFPT = lm(Value ~ moranScaled, data = df_fit_fpt)$coefficients[2]
intKDE = lm(Value ~ moranScaled, data = df_fit_kde)$coefficients[1]
slopeKDE = lm(Value ~ moranScaled, data = df_fit_kde)$coefficients[2]

### Figure 1 ----
pdf("graphics/Figures/Fig3_BRNs.pdf", height = 4, useDingbats = FALSE)
a = ggplot(df_fit_fpt, aes(x = moranScaled, y = Value, group = factor(ID))) +
	geom_smooth(
		aes(moranScaled, Value, group = ID),
		color = "darkgrey",
		size = 0.25,
		method = lm,
		se = FALSE
	) +
	#geom_abline(intercept = intFPT, slope = slopeFPT, size = 1) +
	xlab("Forage patch aggregation") +
	ylab("First-passage time") +
	ggtitle('A') +
	theme(
		legend.position = 'none',
		plot.title = element_text(size = 9),
		axis.text = element_text(size = 8, color = "black"),
		axis.title = element_text(size = 9),
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),
		panel.border = element_rect(
			colour = "black",
			fill = NA,
			size = 0.5
		)
	)

b = ggplot(df_fit_kde, aes(x = moranScaled, y = Value, group = factor(ID))) +
	geom_smooth(
		aes(moranScaled, Value, group = ID),
		colour = "darkgrey",
		size = 0.25,
		method = lm,
		se = FALSE
	) +
	#geom_abline(intercept = intKDE, slope = slopeKDE, size = 1) +
	xlab("Forage patch aggregation") +
	ylab("Range-use ratio") +
	ggtitle('B') +
	theme(
		legend.position = 'none',
		plot.title = element_text(size = 9),
		axis.text = element_text(size = 8, color = "black"),
		axis.title = element_text(size = 9),
		panel.grid.minor = element_blank(),
		panel.background = element_blank(),
		panel.border = element_rect(
			colour = "black",
			fill = NA,
			size = 0.5
		)
	)
grid.arrange(a, b, ncol = 2)
dev.off()
