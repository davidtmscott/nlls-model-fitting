#!/usr/bin/env Rscript
# Author: David Scott
# Contact: 
# Date:  November 2018
# Description: Plotting BioDF data (wrangled)

rm(list=ls()) # clears workspace
graphics.off() # closes all open graphics

#### Packages ####
require(ggplot2)
require(dplyr)

library(reshape2)
library(xtable)

############# Load the dataset ###############
# raw data and estimated parameter values. Calculated in Wrangle.R script
Bio_DF <- read.csv("../data/BioTraits_Params.csv", header = TRUE, stringsAsFactors = FALSE)
# optimised parameters from pythons lmfit minimize function
BioParams <- read.csv("../data/BioTraits_FinalParams.csv", header = TRUE)

############## constants ###################
K <- 8.617e-5 # boltzman constant

############### functions #####################

briere <- function(NewTemp, ParamSub){
  # Runs briere model in Celsius
  # B0 is a normalisation consatnt
  # Tm is the maximum temp recorded 
  # T0 is the minimum tem recorded
  # Returns predicted values 
  
  # variable 
  Temp = NewTemp
  # params
  B0 = ParamSub$b_B0
  Tm = ParamSub$b_Tm
  T0 = ParamSub$b_T0
  # the model
  B = B0 * Temp * (Temp - T0) * (Tm - Temp)**1/2
  return(B)
}


cubic <- function(NewTemp, ParamSub){
  # variable 
  Temp = NewTemp
  # params
  B0 = ParamSub$c_B0
  B1 = ParamSub$c_B1
  B2 = ParamSub$c_B2
  B3 = ParamSub$c_B3
  # th model
  B = B0 + (B1 * Temp) + (B2 * ((Temp)**2)) + (B3 * ((Temp)**3))
  return(B)
}


schoolfield <- function(NewTemp_K, ParamSub){
  # variable 
  Temp = NewTemp_K
  # params
  B0 = ParamSub$s_B0
  E = ParamSub$s_E
  Eh = ParamSub$s_Eh
  El = ParamSub$s_El
  Tl = ParamSub$s_Tl
  Th = ParamSub$s_Th
  # the model
  B = (B0 * exp(((0-E)/K)*((1/Temp)-(1/283.15)))) / (1 + exp((El/K)*((1/Tl) - (1/Temp))) + exp((Eh/K)*((1/Th)-(1/Temp))))
  return(B)
}


schoolfield_high <- function(NewTemp_K, ParamSub){
  # variable 
  Temp = NewTemp_K
  # params
  B0 = ParamSub$sh_B0
  E = ParamSub$sh_E
  Eh = ParamSub$sh_Eh
  Th = ParamSub$sh_Th
  # the model
  B = (B0 * exp(((0-E)/K)*((1/Temp)-(1/283.15)))) / (1 + exp((Eh/K)*((1/Th)-(1/Temp))))
  return(B)
}


schoolfield_low <- function(NewTemp_K, ParamSub){
  # variable 
  Temp = NewTemp_K
  # params
  B0 = ParamSub$sl_B0
  E = ParamSub$sl_E
  El = ParamSub$sl_El
  Tl = ParamSub$sl_Tl
  # the model
  B = (B0 * exp(((0-E)/K)*((1/Temp)-(1/283.15)))) / (1 + exp((El/K)*((1/Tl)-(1/Temp))))
  return(B)
}

#..............................................................................
# plot all five models with optimised parameters
#..............................................................................

plot_it <- function(i){
  # create new df
  New_DF <- as.data.frame(matrix(nrow=50, ncol=7))
  
  colnames(New_DF) <- c("NewTemp_C", "NewTemp_K", "Briere", "Cubic", 
                        "Schoolfield", "Schoolfield_High", "Schoolfield_Low") 
  # subset data
  PlotSub <- Bio_DF[Bio_DF$FinalID == i, ] # observed data
  ParamsSub <- BioParams[BioParams$FinalID == i, ] # optimised parameters

  # create new data for smooth lines
  New_DF$NewTemp_C <- seq(min(PlotSub$ConTemp), max(PlotSub$ConTemp), length=50)
  New_DF$NewTemp_K <- seq(min(PlotSub$ConTemp_K), max(PlotSub$ConTemp_K), length=50)
  
  # predicted trait values
  New_DF$Briere <- briere(New_DF$NewTemp_C, ParamsSub)
  New_DF$Cubic <- cubic(New_DF$NewTemp_C, ParamsSub)
  New_DF$Schoolfield <- schoolfield(New_DF$NewTemp_K, ParamsSub)
  New_DF$Schoolfield_Low <- schoolfield_low(New_DF$NewTemp_K, ParamsSub)
  New_DF$Schoolfield_High <- schoolfield_high(New_DF$NewTemp_K, ParamsSub)
  
  # put data into long format for ggplot
  New_DF = melt(New_DF, id.vars = c("NewTemp_C", "NewTemp_K"), 
                variable.name = "Model", 
                value.name = "PredictedTrait")
  
  # make the plot. Use observed temperature (celsisu) versus observed trait values for points
  a = ggplot(data = PlotSub, aes(x = ConTemp, y = OriginalTraitValue)) + 
    geom_point(shape = 21, colour = "black", fill = "white", size = 1.3, stroke = 1.2) + 
    # plot new sequenced temperature (in celsius) versus predicted trait values from each model for lines. 
    geom_line(data = New_DF, aes(x = NewTemp_C, y = PredictedTrait, colour = Model, linetype=Model)) +
    # dash the phenomological models witha  dashed line
    scale_linetype_manual(values=c("longdash", "longdash", "solid", "solid", "solid")) +
    scale_color_manual(values=c("rosybrown4", "lightsteelblue3", "seagreen3", "salmon3", "steelblue3")) +
    # labels, theme and legend
    xlab("Temperature - Celsius") + ylab("Trait Value") +
    theme_classic() +
    theme(legend.position="bottom")
        
  # save plot 
  suppressMessages(ggsave(paste0("../output/Model_Plots/", i, ".pdf"), a))
}


###########################################
cat("\nWrangling Data, Generating new TempData and Estimating Trait Values with Optimised Parameters.\n")

# Create unique ID for each response curve  
UniqueID = unique(BioParams$FinalID)

# for (i in UniqueID){
#   plot_it(i)
# }

#..............................................................................
# estimated parameters and optimised parameters as a table
#..............................................................................
### estimated startign parameters for briere and schoolfield models
param_est <- aggregate(Bio_DF[, 10:22], list(Bio_DF$FinalID), mean) %>% ungroup()
# get the mean across all thermal peformance curves
param_est_mean <- sapply(param_est[, 2:14], median, na.rm = T)
# extract those used in schoolfield common to schoolfield low and schoolfield high (reduced versions of model)
Sh_est_mean <- param_est_mean[c("B0", "E", "Eh", "Th")]
Sl_est_mean <- param_est_mean[c("B0", "E", "El", "Tl")]
# join them all together. 
param_est_mean <- c(param_est_mean, Sh_est_mean, Sl_est_mean)   

# create a vetor of variable names to retriee (params values)
p_names <- c("c_B0", "c_B1", "c_B2", "c_B3", 
            "b_B0", "b_Tm", "b_T0", 
            "s_B0", "s_E", "s_El", "s_Eh", "s_Tl", "s_Th", 
            "sh_B0", "sh_E", "sh_Eh", "sh_Th", 
            "sl_B0", "sl_E", "sl_El", "sl_Tl")

### extract the optimised parameters for all models
param_opt <- select(BioParams, p_names)

# get the mean across all TPC's
param_opt_mean <- sapply(param_opt, median, na.rm = T)

# create vector of parameter names
param_names <- c("$B_0$", "$B_1$", "$B_2$", "$B_3$", 
                 "$B_0$", "$T_m$", "$T_0$", 
                 "$B_0$", "$E$", "$E_l$", "$E_h$", "$T_l$", "$T_h$", 
                 "$B_0$", "$E$", "$E_h$", "$T_h$", 
                 "$B_0$", "$E$", "$E_l$", "$T_l$")
# create a vector of model names
model_names <- c(rep("Cubic", 4), rep("Briere", 3), rep("Schoolfield", 6), 
                 rep("Schoolfield High", 4), rep("Schoolfield Low", 4))

# make dataframe with model names, parameter names, mean estimated parameter values and mean optimised parameter values
p <- data.frame(model_names, param_names, param_est_mean, param_opt_mean)

# update the column names
colnames(p) <- c("Model", "Parameter", "Estimated Value", "Optimised Value")

p$Model <- as.factor(p$Model)

# make dataframe into a LaTeX table 
ptab <- xtable(p, type = "latex",
               align = c("c", "c", "c", "c", "c"))
# save the table to LaTeX document
print(ptab, file = "../output/Params.tex", append = FALSE,
      table.placement = "h", caption.placement="top",
      sanitize.text.function = function(x){x},
      floating=FALSE, latex.environments=NULL,booktabs=TRUE,
      include.rownames = FALSE, hline.after = c(-1,0, 4, 7, 13, 17, nrow(p)))

#..............................................................................
## Model Comparison - table and plots
#..............................................................................

# create seperate dataframes for AIC, BIC, adjusted R^2 and R^2 values of each response curve
AIC_df <- select(BioParams, FinalID, c_aic, b_aic, s_aic, sl_aic, sh_aic)
AICc_df <- select(BioParams, FinalID, c_aicc, b_aicc, s_aicc, sl_aicc, sh_aicc)
BIC_df <- select(BioParams, FinalID, c_bic, b_bic, s_bic, sl_bic, sh_bic)
RSq_df <- select(BioParams, FinalID, c_RSq, b_RSq, s_RSq, sl_RSq, sh_RSq)
RSq_adj_df <- select(BioParams, FinalID, c_RSq_adj, b_RSq_adj, s_RSq_adj, sl_RSq_adj, sh_RSq_adj)
colnames(AIC_df) <- c("FinalID", "Cubic", "Briere", "Schoolfield", "Schoolfield Low", "Schoolfield High")
colnames(AICc_df) <- c("FinalID", "Cubic", "Briere", "Schoolfield", "Schoolfield Low", "Schoolfield High")
colnames(BIC_df) <- c("FinalID", "Cubic", "Briere", "Schoolfield", "Schoolfield Low", "Schoolfield High")
colnames(RSq_df) <- c("FinalID", "Cubic", "Briere", "Schoolfield", "Schoolfield Low", "Schoolfield High")
colnames(RSq_adj_df) <- c("FinalID", "Cubic", "Briere", "Schoolfield", "Schoolfield Low", "Schoolfield High")

# AICc_boxplot <- melt(BIC_df, id.vars = c("FinalID"), 
#                      variable.name = "Model", 
#                      value.name = "Score")
# 
# ggplot(AICc_boxplot, aes(Model, Score)) +
# geom_boxplot(outlier.colour="black", outlier.shape=16,
#              outlier.size=2, notch=FALSE)

# get number of converged, not converged (DNC), median AIC, BIC and R^2
# subset all columns except FinalID with [, 2:6]
Converged <- sapply(AIC_df[,2:6], function(x) sum(!is.na(x)))
DNC <- sapply(AIC_df[,2:6], function(x) sum(is.na(x)))
Median_RSq <-sapply(RSq_df[,2:6], median, na.rm = T)
Median_RSq_adj <-sapply(RSq_adj_df[,2:6], median, na.rm = T)
Median_AIC <-sapply(AIC_df[,2:6], median, na.rm = T)
Median_AICc <-sapply(AICc_df[,2:6], median, na.rm = T)
Median_BIC <-sapply(BIC_df[,2:6], median, na.rm = T)

# store them together in new dataframe
s <- data.frame(Converged, DNC, Median_RSq, Median_RSq_adj, Median_AIC, Median_AIC, Median_BIC)
colnames(s) <- c("Converged", "No Converge", "$R^2$", "$adjusted R^2$", "AIC", "AICc", "BIC")

# create LaTeX table from dataframe
stab <- xtable(s, type = "latex", 
               align = c("c", "c", "c", "c", "c", "c", "c", "c"))

# save table to LaTeX document
print(stab, file = "../output/Scores.tex", append = FALSE,
      table.placement = "h", caption.placement="top", 
      sanitize.text.function = function(x){x},
      floating=FALSE, latex.environments=NULL,booktabs=TRUE,
      include.rownames = TRUE, hline.after = c(-1,0, nrow(s)))

#........................................................................................
## Criteria scores

compare <- function(x_df, i){
  # subset dataframe
  a <- x_df[x_df$FinalID == i, ] 
  # find column with the minimum value
  Min <- a[which(colnames(a) == a$Min)]
  # its value
  Minvalue <- Min[1,]
  # upper bounds of acceptable difference
  upper <- Minvalue + 2
  # loop through and replace value with na if noth within distance of 2
  for (i in a[,2:6]){
    if ((!(is.na(i))) && (!(i <= upper))){
        index = which(a == i)
        a[index] = NA
    }
  }
  return(a)
}

#...........................................................................................
# create new dataframe to store results
min_AIC_df <- data.frame("FinalID" = character(0), "Cubic" = numeric(0), "Briere" = numeric(0), "Schoolfield" = numeric(0), 
                    "Schoolfield Low" = numeric(0), "Schoolfield High" = numeric(0), "Min" = character(0))

min_AICc_df <- data.frame("FinalID" = character(0), "Cubic" = numeric(0), "Briere" = numeric(0), "Schoolfield" = numeric(0), 
                         "Schoolfield Low" = numeric(0), "Schoolfield High" = numeric(0), "Min" = character(0))

min_BIC_df <- data.frame("FinalID" = character(0), "Cubic" = numeric(0), "Briere" = numeric(0), "Schoolfield" = numeric(0), 
                    "Schoolfield Low" = numeric(0), "Schoolfield High" = numeric(0), "Min" = character(0))

#..............................................................................................
# convert FinalID to character class
AIC_df$FinalID <- as.character(AIC_df$FinalID)
AICc_df$FinalID <- as.character(AICc_df$FinalID)
BIC_df$FinalID <- as.character(BIC_df$FinalID)

# find minimum in each row and assign name to new Min column
AIC_df$Min <- colnames(AIC_df)[apply(AIC_df,1,which.min)]
AICc_df$Min <- colnames(AICc_df)[apply(AICc_df,1,which.min)]
BIC_df$Min <- colnames(BIC_df)[apply(BIC_df,1,which.min)]

# find maximum
RSq_adj_df$Max <- colnames(RSq_adj_df)[apply(RSq_adj_df,1,which.max)]
RSq_df$Max <- colnames(RSq_df)[apply(RSq_df,1,which.max)]

# get table of frequency of each was max Rsquared or adjusted R squared
max_adj_RSq <- table(RSq_adj_df$Max)
max_RSq <- table(RSq_df$Max)

## these are for plots below
# melt them
max_adj_RSq <- melt(max_adj_RSq, variable.name = "Model", 
                value.name = "adj_Rsquared")
# make dataframe
max_adj_RSq <- data.frame(max_adj_RSq)
# dubplicate for each criteria and add column of names
#max_adj_RSq <- rbind(max_adj_RSq, max_adj_RSq, max_adj_RSq)
# add criteria names
#max_adj_RSq$Criteria <- c(rep("AIC", 5), rep("AICc", 5), rep("BIC", 5))
#colnames(max_adj_RSq) <- c("Model", "Freq", "Criteria") 
colnames(max_adj_RSq) <- c("Model", "adj_Rsquared") 

rownames(max_adj_RSq) <- max_adj_RSq[,1]
max_adj_RSq <- select(max_adj_RSq, adj_Rsquared)

# apply function 
for (i in UniqueID){
  a <- compare(AIC_df, i)
  min_AIC_df <- rbind(min_AIC_df, a)
}

for (i in UniqueID){
  ac <- compare(AICc_df, i)
  min_AICc_df <- rbind(min_AICc_df, ac)
}

for (i in UniqueID){
  b <- compare(BIC_df, i)
  min_BIC_df <- rbind(min_BIC_df, b)
}

# get number fo times each model was best according to each criteria
p_aic <- sapply(min_AIC_df[,2:6], function(x) sum(!is.na(x)))
p_aicc <- sapply(min_AICc_df[,2:6], function(x) sum(!is.na(x)))
p_bic <- sapply(min_BIC_df[,2:6], function(x) sum(!is.na(x)))

# melt them for ggplot (long format)
p_aic = melt(p_aic, variable.name = "Model", 
                 value.name = "AIC")

p_aicc = melt(p_aicc, variable.name = "Model", 
             value.name = "AICc")

p_bic = melt(p_bic, variable.name = "Model", 
             value.name = "BIC")

# combine them and make index names a column
mm <-cbind(p_aic, p_aicc, p_bic)
mm <- tibble::rownames_to_column(mm, "Model")

# melt once more
mm<-melt(mm, id.vars = c("Model"), 
     variable.name = "Criteria", 
     value.name = "Freq")


# plot frequency of lowest AIC scores per model with model type
b <- ggplot(data=mm, aes(x=factor(Criteria), y=Freq, fill = Model)) +
  geom_bar(stat="identity", position = "dodge", width = 0.8, colour = "black") +
  scale_fill_manual(values=c("steelblue4", "slategray3", "black",  "steelblue3",  "slateblue4")) +
  geom_text(aes(label=Freq), position = position_dodge(width= 0.8),  vjust=-0.3, size=3.5) +
  xlab("Criteria") + ylab("Count") +
  theme_classic() +
  theme(legend.position="bottom") 
  
#b <- b + geom_point(data = max_adj_RSq, position = position_dodge(0.8), size = 2, color = "red")

# save plot
ggsave("../output/AIC_Plots/Lowest_AICscores.pdf", b) 


c <- ggplot(max_adj_RSq, aes(Model, Freq)) +
  geom_linerange(aes(x = Model, ymin = 0, ymax = Freq, group = Model),
                 color = "lightgray", size = 1.5, position = position_dodge(0.8)) +
  #geom_text(aes(label=Freq), position = position_dodge(width= 0.8), vjust=-0.9,   size=3.5) +
  geom_point(aes(color = Model), position = position_dodge(0.8), size = 5) +
  scale_color_manual(values = c("red", "purple3", "black",  "steelblue3",  "seagreen4"))+
  xlab("Model") + ylab("Count") +
  theme_classic() +
  theme(legend.position="none")

print(c)

#...........................................................................................
#...............................................................................
# barplot


# # create a dataframe of how many times each model had the lowest AIC score
# min_freq <- as.data.frame(table(AIC_df$Min))
# # add a column with the type fo model
# min_freq$Type <- c(rep("Phenomenological", 2), rep("Mechanistic", 3))
# 
# # plot frequency of lowest AIC scores per model with model type
# b <- ggplot(data=min_freq, aes(x=Var1, y=Freq, fill = Type)) +
#   geom_bar(stat="identity", width = 0.8, colour = "black") +
#   #scale_fill_brewer(palette="Blues") +
#   scale_fill_manual(values=c("steelblue4", "slategray3")) +
#   geom_text(aes(label=Freq), vjust=-0.3, size=3.5) +
#   xlab("Model") + ylab("Frequency") +
#   theme_classic() +
#   theme(legend.position="bottom")
# print(b)

# save plot
#ggsave("../output/AIC_Plots/Lowest_AICscores.pdf", b)  









