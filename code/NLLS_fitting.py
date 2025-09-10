#!/usr/bin/env python3
# Date: January 2019

""" Uses NLLS method to fit models to BioTraits data in Python3""" 

__appname__ = '[NLLS_fitting.py]'
__author__ = 'David Scott'
__version__ = '0.0.1'
__license__ = "License for this code/program"

## imports 
import pandas as pd
import numpy as np
import time
#import random

from lmfit import minimize, Parameters, report_fit
from model_functions import *

import matplotlib.pyplot as plt

## constants 
#from scipy.constants import Boltzmann as K
from math import e
K = 8.617e-5

## data 
# import data as a dataframe using pandas
BioDF = pd.read_csv("../data/BioTraits_Params.csv", sep = ",")

## explore dataframe
#BioDF.head
#BioDF.tail 
#BioDF.describe 
#BioDF.values
#BioDF.columns
#BioDF.shape
#BioDF.dtype

# set seed
#random.seed(5)

# start timer
start = time.time()

# create a test data subset of one TPC 
#testdf = BioDF[BioDF["FinalID"] == "MTD2071"]

#.........................................................
# get list of unique response curves by FinalID
uniqueID = BioDF.FinalID.unique()

# get the length to determine size of dataframes for results
frame_size = len(uniqueID)

print("There is a total of {} individual response curves to be fitted!\
    \nHold on tight!".format(frame_size))

#..............................................................................
# The main loop
#..............................................................................

# counters for failed and succcesful converges. 
cubic_fail = 0
cubic_converge = 0

briere_fail = 0
briere_converge = 0

school_fail = 0
school_converge = 0

school_high_fail = 0
school_high_converge = 0

school_low_fail = 0
school_low_converge = 0

# index for inserting data to predefined dataframes. 
indexx = 0
# list to stores lists
all_res = []

for i in uniqueID:
    # subset dataframe 
    Subset = BioDF.loc[BioDF["FinalID"] == i]
    # get temp and Trait values
    X = np.array(Subset.ConTemp)
    Y = np.array(Subset.OriginalTraitValue)
    # kelvin
    X_K = np.array(Subset.ConTemp_K)
    #..............................................
    # cubic 
    # call function 
    c_res = cubic(Subset, i, Temp = X, Trait = Y, n = 3)
    # if returns None (did not converge)
    if c_res == None:
        c_res = [i, np.nan, np.nan, np.nan, np.nan, 
                np.nan, np.nan, np.nan, np.nan, np.nan]
        cubic_fail += 1
    # else add results (final params and stats) to dataframe
    else:
        cubic_converge +=1
    #...........................................
    # briere
    b_res = briere(Subset, Temp = X, Trait = Y, n = 3)
    if b_res == None:
        b_res = [np.nan, np.nan, np.nan, np.nan, np.nan, 
                np.nan, np.nan, np.nan]
        briere_fail += 1
    else:
        briere_converge +=1
    #..............................................
    #schoolfield
    s_res = schoolfield(Subset, Temp = X_K, Trait = Y, n = 5)
    if s_res == None:
        s_res = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 
                np.nan, np.nan, np.nan, np.nan, np.nan]
        school_fail += 1
    else:
        school_converge +=1
    # #..............................................
    # #schoolfield_high (reduced model)
    sh_res = schoolfield_high(Subset, Temp = X_K, Trait = Y, n = 5)
    if sh_res == None:
        sh_res = [np.nan, np.nan, np.nan, np.nan, np.nan, 
                np.nan, np.nan, np.nan, np.nan]
        school_high_fail += 1
    else:
        school_high_converge +=1
    # #..............................................
    # # schoolfield low (reduced model)
    sl_res = schoolfield_low(Subset, Temp = X_K, Trait = Y, n = 5)
    if sl_res == None:
        sl_res = [np.nan, np.nan, np.nan, np.nan, np.nan, 
                np.nan, np.nan, np.nan, np.nan]
        school_low_fail += 1
    else:
        school_low_converge +=1
    
    # concat results into one big list
    a = c_res + b_res + s_res + sl_res + sh_res
    # append to list of lists
    all_res.append(a)

    # print progress
    if indexx == 300:
        print("Gone through 300 curves! \nA third of the way there!!")
    if indexx == 600:
        print("Gone through 600 curves! \nTwo thirds of the way there!!")

    indexx += 1


#.........................................................................................................

# success and failure of the cubic polynomial 
print("\nThe Cubic polynomial           :       {} SUCCESFULLY converged and {} FAILED\
    ".format(cubic_converge, cubic_fail))

print("\nThe Briere model               :       {} SUCCESFULLY converged and {} FAILED\
    ".format(briere_converge, briere_fail))

print("\nThe Schoolfield model          :       {} SUCCESFULLY converged and {} FAILED\
    ".format(school_converge, school_fail))

print("\nThe Schoolfield (high) model   :       {} SUCCESFULLY converged and {} FAILED\
    ".format(school_high_converge, school_high_fail))

print("\nThe Schoolfield (low) model    :       {} SUCCESFULLY converged and {} FAILED\
    ".format(school_low_converge, school_low_fail))


#.........................................................................................................

# list of column names in correct order c for cubic, s for schoolfiel, sl for schoolfield low etc etc
all_names = ['FinalID', 
            'c_B0', 'c_B1', 'c_B2', 'c_B3', 'c_bic', 'c_aic', 'c_aicc', 'c_RSq', 'c_RSq_adj', 
            'b_B0', 'b_T0', 'b_Tm', 'b_bic', 'b_aic', 'b_aicc', 'b_RSq', 'b_RSq_adj', 
            's_B0', 's_E', 's_El', 's_Eh', 's_Tl', 's_Th', 's_bic', 's_aic', 's_aicc', 's_RSq', 's_RSq_adj', 
            'sh_B0', 'sh_E', 'sh_Eh', 'sh_Th', 'sh_bic', 'sh_aic', 'sh_aicc', 'sh_RSq', 'sh_RSq_adj', 
            'sl_B0', 'sl_E', 'sl_El', 'sl_Tl', 'sl_bic', 'sl_aic', 'sl_aicc', 'sl_RSq', 'sl_RSq_adj']

# create data frame with column names
final_params_df = pd.DataFrame(all_res, columns=all_names)

#...........................................................
# # save the dataframes to a csv 
#...........................................................
final_params_df.to_csv("../data/BioTraits_FinalParams.csv", sep=',')    

print("\nMerged dataframes and saved as single csv (BioTraits_FinalParams.csv) to Data directory.")

# final run time
end = time.time()
run_time = (end - start) / 60
print("The NLLS model fitting took %.2f minutes to complete." %run_time)
