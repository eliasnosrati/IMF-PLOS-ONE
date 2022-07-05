
###################################################################################################
###################################################################################################

# Structural adjustment and infectious disease
# Elias Nosrati

###################################################################################################
###################################################################################################

# Clear environment and set working directory
ls()
rm(list = ls())
setwd("")

# Packages
library(haven)
library(stargazer)
library(plm)
library(lmtest)
library(viridis)
library(tidyverse)

###################################################################################################
###################################################################################################

# Data
data <- read_csv("data.csv")
names(data)

###################################################################################################
###################################################################################################

# Baseline mortality model: lag 1 year
FE <- plm(resp_deaths ~ lagBA2PRI +
			factor(year) | . - lagBA2PRI + lagZPRI,
		data = data, 
		model = "within",
		index = c("country", "year"))
summary(FE, vcov = vcovSCC(FE, type = "HC3", cluster = "group"))

# Baseline mortality model: lag 5 years
FE <- plm(resp_deaths ~ lag5BA2PRI +
			factor(year) | . - lag5BA2PRI + lag5ZPRI,
		data = data, 
		model = "within",
		index = c("country", "year"))
summary(FE, vcov = vcovSCC(FE, type = "HC3", cluster = "group"))

# Baseline mortality model: lag 10 years
FE <- plm(resp_deaths ~ lag10BA2PRI +
			factor(year) | . - lag10BA2PRI + lag10ZPRI,
		data = data, 
		model = "within",
		index = c("country", "year"))
summary(FE, vcov = vcovSCC(FE, type = "HC3", cluster = "group"))

###################################################################################################

# Control: GDP
FE <- plm(resp_deaths ~ lagBA2PRI +
			log(wdi_gdpcapcon2010) + factor(year) | . - lagBA2PRI + lagZPRI,
		data = data, 
		model = "within",
		index = c("country", "year"))
summary(FE, vcov = vcovSCC(FE, type = "HC3", cluster = "group"))

# Control: financial crisis
FE <- plm(resp_deaths ~ lagBA2PRI +
			crisis_LV + factor(year) | . - lagBA2PRI + lagZPRI,
		data = data, 
		model = "within",
		index = c("country", "year"))
summary(FE, vcov = vcovSCC(FE, type = "HC3", cluster = "group"))

# Control: reserves
FE <- plm(resp_deaths ~ lagBA2PRI +
			reserves_WDI + factor(year) | . - lagBA2PRI + lagZPRI,
		data = data, 
		model = "within",
		index = c("country", "year"))
summary(FE, vcov = vcovSCC(FE, type = "HC3", cluster = "group"))

# Control: democracy
FE <- plm(resp_deaths ~ lagBA2PRI +
			p_ipolity2 + factor(year) | . - lagBA2PRI + lagZPRI,
		data = data, 
		model = "within",
		index = c("country", "year"))
summary(FE, vcov = vcovSCC(FE, type = "HC3", cluster = "group"))

# Control: egalitarian democracy
FE <- plm(resp_deaths ~ lagBA2PRI +
			egaldem + factor(year) | . - lagBA2PRI + lagZPRI,
		data = data, 
		model = "within",
		index = c("country", "year"))
summary(FE, vcov = vcovSCC(FE, type = "HC3", cluster = "group"))

# Control: coup
FE <- plm(resp_deaths ~ lagBA2PRI +
			coup_any + factor(year) | . - lagBA2PRI + lagZPRI,
		data = data, 
		model = "within",
		index = c("country", "year"))
summary(FE, vcov = vcovSCC(FE, type = "HC3", cluster = "group"))

# Control: UNGA voting alignment
FE <- plm(resp_deaths ~ lagBA2PRI +
			s_unga3g7 + factor(year) | . - lagBA2PRI + lagZPRI,
		data = data, 
		model = "within",
		index = c("country", "year"))
summary(FE, vcov = vcovSCC(FE, type = "HC3", cluster = "group"))

# Control: average female education
FE <- plm(resp_deaths ~ lagBA2PRI +
			meduf + factor(year) | . - lagBA2PRI + lagZPRI,
		data = data, 
		model = "within",
		index = c("country", "year"))
summary(FE, vcov = vcovSCC(FE, type = "HC3", cluster = "group"))

# Control: hospital beds
FE <- plm(resp_deaths ~ lagBA2PRI +
			hosp_beds + factor(year) | . - lagBA2PRI + lagZPRI,
		data = data, 
		model = "within",
		index = c("country", "year"))
summary(FE, vcov = vcovSCC(FE, type = "HC3", cluster = "group"))

###################################################################################################
###################################################################################################

# Visualise
# Lag 1
FE1 <- plm(resp_deaths ~ lagBA2PRI +
			factor(year) | . - lagBA2PRI + lagZPRI,
		data = data, 
		model = "within",
		index = c("country", "year"))

# Extract parameter estimates and covariance matrix
P1 <- matrix(summary(FE1)$coefficients[, 1])
V1 <- matrix(vcovSCC(FE1, type = "HC3", cluster = "group"), 27, 27)

# To account for estimation uncertainty, draw parameter estimates from multivariate Normal distribution
set.seed(959455)
sim.par <- mvtnorm::rmvnorm(n = 100000, mean = P1, sigma = V1)

# Set covariate values
IMF0 <- c(0, rep(0, 25), 1)
IMF1 <- c(1, rep(0, 25), 1)

# Simulate expected values for IMF == 0
E0_1 <- NULL
for (i in 1:100000) {
	E0_1[i] <- IMF0 %*% sim.par[i, ]
}

# Simulate expected values for IMF == 1
E1_1 <- NULL
for (i in 1:100000) {
	E1_1[i] <- IMF1 %*% sim.par[i, ]
}

# Calculate first difference
FD1 <- (E1_1 - E0_1)
FD1 <- tibble(FD = FD1, lag = "1 year")

###################################################################################################

# Lag 5
FE2 <- plm(resp_deaths ~ lag5BA2PRI +
			factor(year) | . - lag5BA2PRI + lag5ZPRI,
		data = data, 
		model = "within",
		index = c("country", "year"))

# Extract parameter estimates and covariance matrix
P2 <- matrix(summary(FE2)$coefficients[, 1])
V2 <- matrix(vcovSCC(FE2, type = "HC3", cluster = "group"), 23, 23)

# Estimation uncertainty
set.seed(959455)
sim.par <- mvtnorm::rmvnorm(n = 100000, mean = P2, sigma = V2)

# Set covariate values
IMF0 <- c(0, rep(0, 21), 1)
IMF1 <- c(1, rep(0, 21), 1)

# Simulate expected values for IMF == 0
E0_2 <- NULL
for (i in 1:100000) {
	E0_2[i] <- IMF0 %*% sim.par[i, ]
}

# Simulate expected values for IMF == 1
E1_2 <- NULL
for (i in 1:100000) {
	E1_2[i] <- IMF1 %*% sim.par[i, ]
}

# Calculate first difference
FD2 <- (E1_2 - E0_2)
FD2 <- tibble(FD = FD2, lag = "5 years")

###################################################################################################

# Lag 10
FE3 <- plm(resp_deaths ~ lag10BA2PRI +
			factor(year) | . - lag10BA2PRI + lag10ZPRI,
		data = data, 
		model = "within",
		index = c("country", "year"))

# Extract parameter estimates and covariance matrix
P3 <- matrix(summary(FE3)$coefficients[, 1])
V3 <- matrix(vcovSCC(FE3, type = "HC3", cluster = "group"), 18, 18)

# Estimation uncertainty
set.seed(959455)
sim.par <- mvtnorm::rmvnorm(n = 100000, mean = P3, sigma = V3)

# Set covariate values
IMF0 <- c(0, rep(0, 16), 1)
IMF1 <- c(1, rep(0, 16), 1)

# Simulate expected values for IMF == 0
E0_3 <- NULL
for (i in 1:100000) {
	E0_3[i] <- IMF0 %*% sim.par[i, ]
}

# Simulate expected values for IMF == 1
E1_3 <- NULL
for (i in 1:100000) {
	E1_3[i] <- IMF1 %*% sim.par[i, ]
}

# Calculate first difference
FD3 <- (E1_3 - E0_3)
FD3 <- tibble(FD = FD3, lag = "10 years")

###################################################################################################

# Gather
FD <- bind_rows(FD1, FD2, FD3)
names(FD) <- c("fd", "lag")
FD$lag <- factor(FD$lag, 
	levels = c("1 year", "5 years", "10 years"))

write_csv(FD, "FD_pri.csv")

# Plot
g <- ggplot(FD, aes(x = fd, fill = lag)) +
		geom_density(adjust = 3, alpha = 0.6) +
		geom_vline(aes(xintercept = 94.6438), 
			colour = "red", linetype = "dashed", size = 1) +
		geom_vline(aes(xintercept = 72.3388), 
			colour = "darkgreen", linetype = "dashed", size = 1) +
		geom_vline(aes(xintercept = 46.4370), 
			colour = "blue", linetype = "dashed", size = 1) +
		labs(
			fill = "Time lag", 
			x = "Excess deaths per 100,000 population caused by IMF-mandated privatisation reforms", 
			y = "Density") +
		scale_x_continuous(breaks = seq(-40, 220, 20)) +
		coord_fixed(4000) +
		theme_bw() +
		theme(legend.position = "bottom")

###################################################################################################
###################################################################################################

# Instrument validity
FE <- plm(BA2PRI ~ ZPRI + log(wdi_gdpcapcon2010) + reserves_WDI + 
			crisis_LV + coup_any + p_ipolity2 + s_unga3g7 + factor(year),
		data = data, 
		model = "within",
		index = c("country", "year"))
summary(FE, vcov = vcovHC(FE, method = "arellano", cluster = "group"))

# F-test
pwaldtest(FE, test = "F", vcov = vcovHC(FE, method = "arellano", cluster = "group"))
car::linearHypothesis(FE, c("ZPRI = 0")) # Chi-square = 40.464 (1 df), p < 0.001

###################################################################################################
###################################################################################################
