#!/usr/bin/env Rscript
# Author: David Scott
# Contact: 
# Date:  November 2018
# Description: Wrangles biotraits data and estimate starting parameetr values

rm(list=ls()) # clears workspace
graphics.off() # closes all open graphics

#### Packages ####
require(dplyr)

#### The Data ####

cat("\nLoading the Data. \n")

# load dataset
BioData <- read.csv("../data/BioTraits.csv", header = TRUE)
 
# Inspect the dataset 
#dplyr::tbl_df(BioData) # displays top rows 
#dplyr::dim_desc(BioData) # prints dimensions
#dplyr::glimpse(BioD) # prints structure
#utils::View(BioData) #displays in spreadsheet like form

# Create a DataFrame (tbl_df)
BioDF <- tbl_df(BioData)

# keep only desired columns
BioDF <- select(BioDF, FinalID, OriginalTraitName, OriginalTraitDef,
                OriginalTraitValue, OriginalTraitUnit, ConTemp, ConTempUnit)

############################################################################################
#..........................................................................................#
#                                 Wrangling Raw Data                                       #
#..........................................................................................#
############################################################################################

cat("\nWrangling the Data. \n")

# replace any absent data with NA
BioDF[BioDF == ""] = NA

# remove rows with NA from originial trait value and ConTemp
BioDF <- filter(BioDF, !is.na(OriginalTraitValue))
BioDF <- filter(BioDF, !is.na(ConTemp))
#BioDF <- BioDF[!is.na(BioDF$OriginalTraitValue), ]

# remove zero or negative OriginalTraitValues 
BioDF <- filter(BioDF, OriginalTraitValue > 0)
#BioDF <- BioDF[as.numeric(as.character(BioDF$OriginalTraitValue)) > 0,]

# number of points to be filteres
number_ofpoints = 8
# filter response curve with < 8 points (schoolfield has 6 params)
BioDF <- BioDF %>% group_by(FinalID) %>% filter(n() >= number_ofpoints) #%>% droplevels() %>% ungroup()
# save number of points filtered so it can be read in LaTeX
cat(number_ofpoints, file = "../output/Points.tex")

# calculate log of OriginalTraitValue and add as new column
BioDF <- mutate(BioDF, OriginalTraitValue_log = log(OriginalTraitValue))

# convert temperature from celsius to kelvin. add as new column
BioDF <- mutate(BioDF, ConTemp_K = ConTemp + 273.15)

#names(BioDF)

# species = length(unique(BioDF$Consumer))
# species
# 
# species = length(unique(BioDF$ConSpecies))
# species
# 
# habitat = length(unique(BioDF$Habitat))
# 
# trait = length(unique(BioDF$OriginalTraitName))
# 
# 
# climate = length(unique(BioDF$Climate))

# summarise the data
#BioDF <- group_by(BioDF, FinalID)

# BioDF_Summary <- summarize(BioDF, Count = n(), 
#                            av_trait_log = mean(OriginalTraitValue_log, na.rm = T),
#                            av_temp = mean(ConTemp_K, na.rm = T))
# BioDF_Summary

# Create unique ID for each response curve 
UniqueID = unique(BioDF$FinalID)
# get number of curves to be plotted 
number_ofcurves <- length(UniqueID)
# save them to be read into LaTeX document
cat(number_ofcurves, file = "../output/Curves.tex")
############################################################################################
#..........................................................................................#
#                              Parameter Starting Values                                   #
#..........................................................................................#
############################################################################################

## constants 
# set reference temperature 
ref <- 283.19 #(10 C)

K <- 8.617e-5 # boltzman constant

#............................................................................#
## Functions

#ResponseSub <- BioDF[BioDF$FinalID == "MTD2079", ]
# 
# plot(ResponseSub$OriginalTraitValue ~ ResponseSub$ConTemp_K)
# 
# plot(ResponseSub$OriginalTraitValue_log ~ X_K)

# fucntion to estimates schoolfield startign params
schoolfield <- function(ResponseSub, FinalID){
  # set ConTemp_K to X and OriginalTraitValue_log to Y
  X <- ResponseSub$ConTemp_K
  Y <- ResponseSub$OriginalTraitValue

  # maximum trait value (Y) observed
  Ypeak<-max(Y) 
  # index of peak 
  peak_index <- which.max(Y)
  # temp (in K) with highest trait value (the peak)
  Xpeak <- X[peak_index] 
  
  # index each point before and after peak, to get associated temperature values
  pre_index <- which(X <= Xpeak)
  post_index <- which(X >= Xpeak)
  pre <- X[pre_index]
  post <- X[post_index]
  
  # index of min and max y values
  Y_min <- which.min(Y)
  #Ymax <- which.max(Y)
  #............................................................................#
  # calculate B0 as the temo closest to ref temperature
  if (min(X) > ref){
    B0 <- min(Y) 
    #B0 <- min(Y[which(X)], na.rm = TRUE)
  } else {
    B0 <- max(Y[which(X <= ref)])
  }
  #............................................................................#
  ## transform the curve (this flips it horizontally)
  # convert temp values (all temps and pre and post peak temps)
  X_K = 1/(K * X)
  pre_K = 1/(K * pre)
  post_K = 1/(K * post)
  # convert Y to log scale
  Y_log <- ResponseSub$OriginalTraitValue_log
  # get the mean of Y values (logged)
  Y_mean <- mean(Y_log)
  #............................................................................#
  if (length(pre_K) < 3){
    # get linear model of post peak
    m <- lm(Y_log[post_index] ~ post_K)
    # high temp deactivation energy (post peak slope) - this is transformed 
    Eh <- summary(m)$coefficients[2]
    E = Eh * 0.5
    El = E * 0.5
    Tl = Xpeak - 10
    K_Th <- (Y_mean - (summary(m)$coefficients[1]))/Eh
    Th <- 1/(K * K_Th)
    curve = "post"
  } 
  else if (length(post_K) < 3){
    # slope of rise of the curve (now on the right)
    m <- lm(Y_log[pre_index] ~ pre_K)
    # high temp deactivation energy (post peak slope) - this is transformed 
    E <- summary(m)$coefficients[2]
    Eh = E * 2
    El = E * 0.5
    Th = Xpeak + 10
    K_Tl <- (Y_mean - (summary(m)$coefficients[1]))/E
    Tl <- 1/(K * K_Tl)
    curve = "pre"
  } 
  else {
    # slope of rise of the curve (now on the right)
    m1 <- lm(Y_log[pre_index] ~ pre_K)
    # slope of fall of the curve (now on the left)
    m2 <- lm(Y_log[post_index] ~ post_K)
    # activation energy (pre peak slope) - this is transformed 
    E <- summary(m1)$coefficients[2]
    # high temp deactivation energy (post peak slope) - this is transformed 
    Eh <- summary(m2)$coefficients[2]
    # low temp deactivation energy - this is transformed 
    El <- E * 0.5
    #................
    # mean of log trait values minus y intercept of E, divided by E (E is still transformed)
    K_Tl <- (Y_mean - (summary(m1)$coefficients[1]))/E
    K_Th <- (Y_mean - (summary(m2)$coefficients[1]))/Eh
    # convert back to standard temp (Kelvin) 
    Tl <- 1/(K * K_Tl)
    Th <- 1/(K * K_Th)
    curve = "uni"
  }
  # get slope of incline and decline (of now transformed curve)
  # if peak is in the first two points
  # if (peak_index == Y_min || peak_index == Y_min + 1){
  #   # get linear model of whole curve
  #   m1 <- lm(Y_log ~ X_K)
  #   m2 <- lm(Y_log ~ X_K)
  #   curve = "post"
  # } # or last two points or
  # else if (peak_index == length(Y) || peak_index + 1 == length(Y)){
  #   # get linear model of whole curve
  #   m1 <- lm(Y_log ~ X_K)
  #   m2 <- lm(Y_log ~ X_K)
  #   curve = "pre"
  # }
  # else {
  #   # slope of rise of the curve (now on the right)
  #   m1 <- lm(Y_log[pre_index] ~ pre_K)
  #   # slope of fall of the curve (now on the left)
  #   m2 <- lm(Y_log[post_index] ~ post_K)
  #   curve = "uni"
  # }
  # #............................................................................#
  # ### assign parameter values
  # # activation energy (pre peak slope) - this is transformed 
  # E <- summary(m1)$coefficients[2]
  # # high temp deactivation energy (post peak slope) - this is transformed 
  # Eh <- summary(m2)$coefficients[2]
  # # low temp deactivation energy - this is transformed 
  # El <- E * 0.5
  # #............................................................................#
  # # get the mean of Y values (logged)
  # Y_mean <- mean(Y_log)
  # # get midpoint 
  # 
  # # mean of log trait values minus y intercept of E, divided by E (E is still transformed)
  # K_Tl <- (Y_mean - (summary(m1)$coefficients[1]))/E
  # K_Th <- (Y_mean - (summary(m2)$coefficients[1]))/Eh
  # # convert back to standard temp (Kelvin) 
  # Tl <- 1/(K * K_Tl)
  # Th <- 1/(K * K_Th)
  # # if tl or th = infinity 
  # if (Tl == Inf || is.na(Tl) || Tl > Xpeak*3){
  #   Tl = Xpeak - 10 
  # }
  # if (Th == Inf || is.na(Th) || Th > Xpeak*3){
  #   Th = Xpeak + 10
  # }
  # 
  #............................................................................#
  # return absolute values of E, Eh, and El
  params <- data.frame("FinalID" = FinalID, "B0" = B0,
                       "E" = abs(E), "El" = abs(El),
                       "Eh" = Eh, "Tl" =Tl, "Th" = Th, "Tpeak" = Xpeak, "Curve" = curve)
  return(params)
}
 

#............................................................................#
## function for estimatign briere starting params

briere <- function(ResponseSub, FinalID){
  # assign variables
  X <- ResponseSub$ConTemp
  Y <- ResponseSub$OriginalTraitValue
  # maximum feasibale temperature for trait (Tm)
  Tm <- X[which.max(X)]
  # minimum feasibale temperature for trait (T0)
  T0 <- X[which.min(X)]
  #
  #b_B0 <- 
  params <- data.frame("FinalID" = FinalID, 
                       "b_Tm" = Tm, 
                       "b_T0" = T0)
  return(params)
}

#............................................................................#
## create epty data frame to store params
school <- data.frame("FinalID" = character(0), "B0"= numeric(0), "E"= numeric(0), 
                     "El"= numeric(0), "Eh"= numeric(0), "Tl"= numeric(0), 
                     "Th"= numeric(0), "Tpeak"= numeric(0),  "Curve" = character(0))

brie <- data.frame("FinalID" = character(0), 
                   "b_Tm"= numeric(0), 
                   "b_T0"= numeric(0))

#............................................................................#
## create new parameters and add to dataframes

cat("\nCalculating Parameter Starting Values Schoolfield and Briere models\n")

# initial parameter values for cubic polynomial 
BioDF <- mutate(BioDF, c_B0 = 1)
BioDF <- mutate(BioDF, c_B1 = 1)
BioDF <- mutate(BioDF, c_B2 = 1)
BioDF <- mutate(BioDF, c_B3 = 1)

# briere normalisation constant
BioDF <- mutate(BioDF, b_B0 = 0.1)

# loop through each unique response curves and call functions
for (i in UniqueID){
  # subset each unique response curve
  ResponseSub <- BioDF[BioDF$FinalID == i, ]
  s <- schoolfield(ResponseSub, i)
  b <- briere(ResponseSub, i)
  # append params to data frame
  school <- rbind(school, s)
  brie <- rbind(brie, b)
}

#............................................................................#
## Merge new dataframes 

cat("\nMerging Dataframes and Saving to Data Directory\n")

BioDF_params <- merge(BioDF, brie, by = "FinalID")
BioDF_params <- merge(BioDF_params, school, by = "FinalID")

BioDF_params$Curve
final <- BioDF_params[!duplicated(BioDF_params$FinalID),]
table(final$Curve)

#............................................................................#
## Save modfied data to csv file 
write.csv(BioDF_params, "../data/BioTraits_Params.csv", row.names = FALSE)

cat("\nScript Complete. \n")


# to get E 
# this si the slope of linear regression to the right
# peak of log trait values
# we use temp in K, trait value, botzman, and peak)
# get temp that less than or equal to peak temp
# get lm of it with coefficients

# Eh is the slope to the elft hand side

# if peak is on far left and right can set E and Eh to the be the same slope
# i.e the one that goes through all the points 

# El is 0.5 x E

# to get B0
# this is the trait value at the ref temp
# we use our refernce temp, temp in K and trait values 
# if min tem greater than ref temp 
#     return log of the minimum temp 
# else return log of max temp below or equal to ref temp

# Th is temp with highest trait value 

# Tl was lowest temp with a trait value