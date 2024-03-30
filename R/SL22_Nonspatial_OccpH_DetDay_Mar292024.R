
#first, loading the packages, note that SpOccupancy and others must be installed
library(spOccupancy)
library(coda)
library(stars)
library(ggplot2)
library(MCMCvis)
library(cowplot)
library(dbplyr)
library(dplyr)
library(lubridate)
library(sf)
library(tibble)
library(tidyr)

#now loading in the SL22 data
setwd("C:/Users/natha/OneDrive/Documents/IBIO6070")
SL22 <- read.csv("C:/Users/natha/OneDrive/Documents/IBIO6070/SL22_Nate2.csv")
#Note that replicate 1 = summer, replicate 2 = fall, and replicate 3 = winter sampling campaigns

#We should first clean the data to remove outliers and NA's 

SL22 <- SL22 %>%
  mutate(PH = ifelse(PH < 4 | PH > 10, NA, PH))


#Before we fit a model, we must make an appropriate data frame with our master sheet


#Now organising the detection data into an array

detection_matrix <- xtabs(Detect ~ Station_ID + Replicate, data = SL22)
y <- as.data.frame.matrix(detection_matrix)

#organising detection covariets into an array (only date)

SL22$Day_of_Year <- yday(ymd(SL22$Date))
Date_data <- xtabs(Day_of_Year ~ Station_ID + Replicate, data = SL22)
Date <- as.data.frame.matrix(Date_data)


det.covs <- list(day = Date)
str(det.covs)


#now making occupancy formula

occ.covs <- aggregate(PH ~ Station_ID, data = SL22, FUN = function(x) mean(x, na.rm = TRUE))
str(occ.covs)

#structuring coords

sf_data <- st_as_sf(SL22, coords = c("Long", "Lat"), crs = 4326)
sf_data_transformed <- st_transform(sf_data, crs = "+proj=utm +zone=17 +datum=WGS84")
SL22$Easting <- st_coordinates(sf_data_transformed)[, 1]
SL22$Northing <- st_coordinates(sf_data_transformed)[, 2]

#Now I have Easting and Northing column :)
#making a matrix of coords

coords <- SL22 %>%
  select(Station_ID, Northing, Easting) %>%
  group_by(Station_ID) %>%
  summarise(Northing = first(Northing), Easting = first(Easting)) %>%
  as.data.frame() %>%
  column_to_rownames(var = "Station_ID") %>%
  as.matrix()
str(coords)

#Now packing everyhing into one final list

data.Pack <- list(y = y, 
                  occ.covs = occ.covs, 
                  det.covs = det.covs, 
                  coords = coords)
str(data.Pack)

#This completes the data prep portion of the model

# now we will try fitting the data to a simple single species occupancy model

#first, we make formulas for occ and det

occ.formula <- ~ scale(pH)
det.formula <- ~ scale(day)

#specifying inits

# Format with explicit specification of inits for alpha and beta
# with four detection parameters and three occurrence parameters 
# (including the intercept).
NZ.inits <- list(alpha = c(0, 0, 0, 0), 
                   beta = c(0, 0, 0), 
                   z = apply(data.Pack$y, 1, max, na.rm = TRUE))
# Format with abbreviated specification of inits for alpha and beta.
NZ.inits <- list(alpha = 0, 
                   beta = 0, 
                   z = apply(data.Pack$y, 1, max, na.rm = TRUE))

#Specifying Priors

NZ.priors <- list(alpha.normal = list(mean = 0, var = 2.72), 
                    beta.normal = list(mean = 0, var = 2.72))

# specificing MCMC 

n.samples <- 5000
n.burn <- 3000
n.thin <- 2
n.chains <- 3

# Run it

out <- PGOcc(occ.formula = occ.formula, 
             det.formula = det.formula, 
             data = data.Pack, 
             inits = NZ.inits, 
             n.samples = n.samples, 
             priors = NZ.priors, 
             n.omp.threads = 1, 
             verbose = TRUE, 
             n.report = 1000, 
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains)

#Now checking results

names(out)
summary(out)
plot(out$beta.samples, density = FALSE) # Occupancy parameters.
plot(out$alpha.samples, density = FALSE) # Detection parameters.

#Finding Baysian P value - this needs fix

ppc.out <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out)

ppc.df <- data.frame(fit = ppc.out$fit.y, 
                     fit.rep = ppc.out$fit.y.rep, 
                     color = 'lightskyblue1')
ppc.df$color[ppc.df$fit.rep > ppc.df$fit] <- 'lightsalmon'
plot(ppc.df$fit, ppc.df$fit.rep, bg = ppc.df$color, pch = 21, 
     ylab = 'Fit', xlab = 'True')
lines(ppc.df$fit, ppc.df$fit, col = 'black')

#WAIC - smaller is better

waicOcc(out)


