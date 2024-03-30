#Spatial Occupancy Model

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

#structure coords

#structuring coords

sf_data <- st_as_sf(SL22, coords = c("Long", "Lat"), crs = 4326)
sf_data_transformed <- st_transform(sf_data, crs = "+proj=utm +zone=17 +datum=WGS84")
SL22$Easting <- st_coordinates(sf_data_transformed)[, 1]
SL22$Northing <- st_coordinates(sf_data_transformed)[, 2]

#Now I have Easting and Northing column :)
#making the matrix

coords <- SL22 %>%
  select(Station_ID, Northing, Easting) %>%
  group_by(Station_ID) %>%
  summarise(Northing = first(Northing), Easting = first(Easting)) %>%
  as.data.frame() %>%
  column_to_rownames(var = "Station_ID") %>%
  as.matrix()
str(coords)


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

occ.covs <- aggregate(cbind(Conductivity, TDS, PH, Easting, Northing) ~ Station_ID, data = SL22, FUN = function(x) mean(x, na.rm = TRUE))
str(occ.covs)
pH <- occ.covs$PH
TDS <- occ.covs$TDS
Conductivity <- occ.covs$Conductivity


#Now packing everyhing into one final list

data.Pack <- list(y = y, 
                  occ.covs = occ.covs, 
                  det.covs = det.covs, 
                  coords = coords)
str(data.Pack)

#Now Lets Set the specs

occ.formula <- ~ scale(pH) + scale(Conductivity) + scale(TDS)
det.formula <- ~ scale(day)

# Pair-wise distances between all sites
dist.SL <- dist(data.Pack$coords)
# Exponential covariance model
cov.model <- "exponential"
# Specify list of inits
NZ.inits <- list(alpha = 0, 
                 beta = 0, 
                 z = apply(data.Pack$y, 1, max, na.rm = TRUE), 
                 sigma.sq = 2, 
                 phi = 3 / mean(dist.SL), 
                 w = rep(0, nrow(data.Pack$y)))


# Batch inits

batch.length <- 25
n.batch <- 400
n.burn <- 2000
n.thin <- 20
n.chains <- 3

SL.tuning <- list(phi = 1)

min.dist <- min(dist.SL)
max.dist <- max(dist.SL)
NZ.priors <- list(beta.normal = list(mean = 0, var = 2.72), 
                  alpha.normal = list(mean = 0, var = 2.72), 
                  sigma.sq.ig = c(2, 1), 
                  phi.unif = c(3/max.dist, 3/min.dist))

n.omp.threads <- 1
verbose <- TRUE
n.report <- 100 # Report progress at every 100th batch.

# RUN IT

out.sp <- spPGOcc(occ.formula = occ.formula, 
                  det.formula = det.formula, 
                  data = data.Pack, 
                  inits = NZ.inits, 
                  n.batch = n.batch, 
                  batch.length = batch.length, 
                  priors = NZ.priors, 
                  cov.model = cov.model, 
                  NNGP = TRUE, 
                  n.neighbors = 5,
                  tuning = SL.tuning, 
                  n.report = n.report, 
                  n.burn = n.burn, 
                  n.thin = n.thin, 
                  n.chains = n.chains)

summary(out.sp)

ppc.sp.out <- ppcOcc(out.sp, fit.stat = 'freeman-tukey', group = 2)
summary(ppc.sp.out)

waicOcc(out.sp)

