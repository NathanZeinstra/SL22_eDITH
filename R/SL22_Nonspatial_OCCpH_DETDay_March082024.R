
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
  filter(Temp >= -40 & Temp <= 40, PH >= 4 & PH <= 10)

SL22 <- SL22 %>%
  filter(!is.na(Date) & !is.na(pH))

#Before we fit a model, we must make an appropriate data frame with our master sheet


#Now organising the detection data into an array

y <- matrix(SL22$Detect, nrow = 141, ncol = 3, byrow = TRUE)
y <- apply(y, c(1, 2), as.numeric)
str(y)

#organising detection covariets into an array (only date)

#converting date coloum to numeric

SL22$Date <- dmy(SL22$Date)
SL22$Date <- ymd(SL22$Date)
SL22$Date <- as.numeric(yday(SL22$Date))
class(SL22$Date)

Date <- SL22 %>%
  select(Station_ID, Replicate, Date) %>%
  rename(Day = Date)

Date <- Date %>%
  pivot_wider(names_from = Replicate, values_from = Day) %>%
  data.matrix()

Date <- Date[complete.cases(Date), ]

str(Date)

#Now I've made a new variable, Date, which can be used as a detection covariate.

#Now making a formula, starting simple with one det cov

det.covs <- list(day = Date)
str(det.covs)

#organising occupancy variables into array (only pH for now)

pH <- SL22$PH
str(pH)

#now making occupancy formula

occ.covs <- data.frame(pH)
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

#Attempting Prediction

pH.pred <- (SL22$PH - mean(data.Pack$occ.covs[, 1])) / sd(data.Pack$occ.covs[, 1])
X.0 <- cbind(1, pH.pred, pH.pred^2) 
out.pred <- predict(out, X.0)

#This needs work too

