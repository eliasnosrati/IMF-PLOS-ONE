
###################################################################################################
###################################################################################################

# Structural adjustment and infectious disease burden
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

# Sensitivity plot for all-cause mortality
b <- 76	
fun1 <- function(x) b / x

f <- ggplot(data.frame(x = seq(0.1, 1, 0.01)), aes(x = x)) +
	stat_function(fun = fun1)
p1 <- layer_data(f)
g1 <- ggplot(p1, aes(x, y)) + 
	geom_line(colour = "violet", size = 1) +
	labs(
		x = "Difference in prevalence of confounder", 
		y = "Excess mortality rate caused by confounder") +
	scale_x_continuous(breaks = seq(0, 1, 0.1)) +
	scale_y_continuous(breaks = seq(0, 1000, 50)) +
	theme_minimal()

ggsave("sensit1.jpg", g1)

###################################################################################################
###################################################################################################

# Sensitivity plot for all-cause DALYs
b <- 95
fun2 <- function(x) b / x

f <- ggplot(data.frame(x = seq(0.1, 1, 0.01)), aes(x = x)) +
	stat_function(fun = fun2)
p2 <- layer_data(f)
g2 <- ggplot(p2, aes(x, y)) + 
	geom_line(colour = "darkblue", size = 1) +
	labs(
		x = "Difference in prevalence of confounder", 
		y = "Excess mortality rate caused by confounder") +
	scale_x_continuous(breaks = seq(0, 1, 0.1)) +
	scale_y_continuous(breaks = seq(0, 1000, 50)) +
	theme_minimal()

ggsave("sensit2.jpg", g2)

###################################################################################################
###################################################################################################
