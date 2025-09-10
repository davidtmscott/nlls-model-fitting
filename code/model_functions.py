#!/usr/bin/env python3
# Date: January 2019

""" Defines non-linear models as functions to be 
fitted to TPC using Python3
""" 

__appname__ = '[model_functions.py]'
__author__ = 'David Scott'
__version__ = '0.0.1'
__license__ = "License for this code/program"

###############################################################################

########## SET A SEED

###############################################################################

# imports 
import pandas as pd
import numpy as np
from lmfit import minimize, Parameters, report_fit

#import matplotlib.pyplot as plt

# constants 
#from scipy.constants import Boltzmann as K
from math import e
K = 8.617e-5

# data 
#BioDF = pd.read_csv("../data/BioTraits_Params.csv", sep = ",")

# create a test data subset of one TPC 
#testdf = BioDF[BioDF["FinalID"] == "MTD3767"]

#.............................................................................
# goodness of fit measures
#..............................................................................

def fit_measure(resid_func, out, Temp, Trait):
     """
     Calculate goodness of fit measures 
     Rsquared
     adjusted R squared
     """
     # sample size
     n = len(Temp)
     # number of parameters
     p = len(out.params)

     # residual sum of squares
     RSS = sum(resid_func(out.params, Temp, Trait)**2)
     # total sum of squares
     TSS = sum((Trait - np.mean(Trait))**2)
     # R squared
     RSq = 1 - (RSS/TSS)
     # adjusted R squared 
     # nn = samle size, p = number of parameters
     #RSq_adj = 1 - ((RSS/(n - p - 1))/(TSS/(n - 1)))
     RSq_adj = 1 - (((1 - RSq) * (n - 1))/(n - p - 1))
     output = [RSq, RSq_adj]
     return output

def calc_AICc(out, Temp):
     """ 
     Calculates AICc score 
     AIC with an extra term to correct fo rsmall sample sizes
     """
     n = len(Temp)
     k = len(out.params)

     AICc = out.aic + ((2*(k**2) + 2*k)/(n-k-1))
     return AICc

#............................................................#
#                                                            #
### Phenomenological Models                                  #
#                                                            #
#............................................................#

def cubic_resids(params, Temp, Trait):
     """
     Calculate residuals for cubic polynomial model.

     Variable
     T    -    Temperature (Celcius)

     Parameter
     B0   -    No Mechanistic interpretation
     B1   - 
     B2   -
     B3   -
 
     """
     # variable
     T = Temp

     # param values
     B0 = params["B0"].value
     B1 = params["B1"].value
     B2 = params["B2"].value
     B3 = params["B3"].value

     # the model 
     B = B0 + (B1 * T) + (B2 * ((T)**2)) + (B3 * ((T)**3))

     # the residual differenc between observed data and estimated data (B) 
     return B - Trait

def cubic(Subset, FinalID, Temp, Trait, n):
     """
     General Cubic Polynomial model (Phenomenological).
     
     Optimises paramter values using minimize from lmfit package

     Calls Functions:
     fit_measure(resid_func, out, Temp, Trait)
     calc_AICc(out, Temp)

     Returns a list of:
     Optimised Parameters: B0, B1, B2, B3
     BIC
     AIC
     AICc
     Rsquared
     adjusted Rsquared
     """
     
     # variable values
     # Temp = np.asarray(Subset.ConTemp)
     # Trait = np.asarray(Subset.OriginalTraitValue)

     # estimated parameter values - can change
     B0 = np.array(Subset.c_B0)[0] 
     B1 = np.array(Subset.c_B1)[0]
     B2 = np.array(Subset.c_B2)[0]
     B3 = np.array(Subset.c_B3)[0]

     # estimated parameter values - cannot change
     B0_orig = B0
     B1_orig = B1
     B2_orig = B2
     B3_orig = B3

     # an initial bestfit list with an arbitarily large AIC
     #         [FinalID, B0, B1, B2, B3, Chisqr, BIC, AIC]  
     bestfit = [FinalID, 0, 0, 0, 0, 0, 100000, 0]

     # DNC - Did Not Converge flag
     # this ensures the above "best" does not get returned if none converge
     DNC = True

     for i in range(n):
          try:
               if i != 0:
                    # resample param values
                    B0 = np.random.normal(B0_orig) 
                    B1 = np.random.normal(B1_orig)
                    B2 = np.random.normal(B2_orig)
                    B3 = np.random.normal(B3_orig)
               # create dictinary of params
               params = Parameters()
               # add with tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)
               params.add_many(("B0",  B0, True, None, None, None, None),
                              ("B1", B1, True, None, None, None, None),
                              ("B2", B2, True, None, None, None, None),
                              ("B3", B3, True, None, None, None, None))
               # minimize residuals
               out = minimize(cubic_resids, params, args = (Temp, Trait))
               #...............................................................
               # write error report
               #A = report_fit(out.params)
               #...............................................................
               ## store results of best fit (based on aic score)
               if out.aic < bestfit[6]:
                    # if try gets to this point, it has converged at least once
                    # so set DNC to False
                    DNC = False  
                    # calculate AICc
                    AICc = calc_AICc(out, Temp)
                    # calculate goodness of fit measures 
                    goodness_of_fit = fit_measure(cubic_resids, out, Temp, Trait)
                    # bestfit takes final params and measures of fit
                    bestfit = [FinalID, 
                    out.params["B0"].value, out.params["B1"].value, 
                    out.params["B2"].value, out.params["B3"].value, 
                    out.bic, out.aic, AICc]
                    # merge best fit and goodness fo fit 
                    bestfit = bestfit + goodness_of_fit            
               # calculate final result
               #final = Trait + out.residual  
          except:
               print("Error")
     # print(final)
     # plt.plot(Temp, Trait, 'o')
     # plt.plot(Temp, final, 'r')
     # plt.show()
     if not DNC:
          return bestfit
     else:
          return None  

# cubic(testdf, "MTD2071", n = 5)
# # create a test data subset of one TPC 
# testdf = BioDF[BioDF["FinalID"] == "MTD2071"]


#...........................................................................................
#  Briere Model
#...........................................................................................

def briere_resids(params, Temp, Trait):
     """
     Briere model (Phenomenological)

     Variables: 
     T    -    Temperature (Celcius) 
     
     Parameters:
     B0   -    Normalisation Constant
     T0   -    Minimum feasible temperature for trait
     Tm   -    Maximum feasible temperature for trait
     """
     
     # variables 
     T = Temp

     # param values
     B0 = params["B0"].value
     T0 = params["T0"].value
     Tm = params["Tm"].value

     # the model 
     B = B0 * T * (T - T0) * (Tm - T)**1/2

     # return the residuals
     return B - Trait

def briere(Subset, Temp, Trait, n):
     """
     Briere model (Phenomenological)

     Optimises paramter values using minimize from lmfit package

     Calls Functions:
     fit_measure(resid_func, out, Temp, Trait)
     calc_AICc(out, Temp)

     Returns a list of:
     Optimised Parameters: B0, Tm, T0
     BIC
     AIC
     AICc
     Rsquared
     adjusted Rsquared

     """
     #** look at np.array vs np.asarray

     # variable values
     # Temp = np.array(Subset.ConTemp)
     # Trait = np.array(Subset.OriginalTraitValue)

     # estimated parameter values 
     B0 = np.array(Subset.b_B0)[0]
     #B0 = Subset.loc[Subset.index[0], "b_B0"]
     T0 = np.array(Subset.b_T0)[0]
     Tm = np.array(Subset.b_Tm)[0]

     # save orginal estimated param values
     B0_orig = B0
     T0_orig = T0
     Tm_orig = Tm 

     # temp peak - using as a bound
     Tpeak = (np.array(Subset.Tpeak)[0] - 273.15)

     # an initial bestfit list with an arbitarily large AIC
     #         [B0, T0, Tm, Chisqr, BIC, AIC] 
     bestfit = [0, 0, 0, 0, 0, 100000, 0]

     # DNC - Did Not Converge flag
     # this ensures the above "best" does not get returned if none converge
     DNC = True

     for i in range(n):
          try:
               # for every try after the first one
               if i != 0:
                    # resample param values
                    B0 = np.random.normal(B0_orig) 
                    T0 = np.random.normal(T0_orig)
                    Tm = np.random.normal(Tm_orig)

               # create dictinary of params
               params = Parameters()
               # add with tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)
               params.add_many(("B0",  B0, True, 0, 1,  None, None),
                              ("T0", T0, True, -10, Tpeak, None, None),
                              ("Tm", Tm, True, Tpeak, 100, None, None))

               # minimize residuals
               out = minimize(briere_resids, params, args = (Temp, Trait))
               #...............................................................
               # write error report
               #A = report_fit(out.params)
               #...............................................................
               ## store results of best fit (based on aic score)
               if out.aic < bestfit[5]:
                    # if try gets to this point, it has converged at least once
                    # so set DNC to False
                    DNC = False  
                    # calculate AICc
                    AICc = calc_AICc(out, Temp)
                    # calculate goodness of fit measures 
                    goodness_of_fit = fit_measure(briere_resids, out, Temp, Trait)
                    # bestfit takes final params and measures of fit
                    bestfit = [ 
                    out.params["B0"].value, out.params["T0"].value, 
                    out.params["Tm"].value, out.bic, out.aic, AICc]
                    # merge best fit and goodness fo fit 
                    bestfit = bestfit + goodness_of_fit            
               # calculate final result to test plot
               #final = Trait + out.residual
          # except Exception as e:
          #      print(e)
          except IOError:
               pass  
          #except:
               #print("Error")
          
     # print(final)
     # plt.plot(Temp, Trait, 'o')
     # plt.plot(Temp, final, 'r')
     # plt.show()

     if not DNC:
          return bestfit
     else:
          return None  




#............................................................#
#                                                            #
#                                                            #
### Mechanistic Models                                       #
#                                                            #
#............................................................#

#...........................................................................................
#   Full schoolfield model
#...........................................................................................

def schoolfield_resids(params, Temp, Trait):
     """ 
     Calculates Residuals for full Schoolfield model 
     Mechanistic - Uses thermodynamics and Enzyme kinetics.

     Output:
     B   -     Trait value at given temperature (in kelvin).

     Constant:
     K   -     Boltzmann constant.

     Parameters:
     B0  -     Trait value at ref temp of 283.15 K (10 C)
     E   -     Activation energy controlling rise of curve up to the peak
     El  -     Low temp de-activation energy
     Tl  -     Enzyme 50% low-temp deactivated 
     Eh  -     Enzye high temp deactivation energy
     Th  -     Enzyme 50% high-temp deactivated
     """

     #variables 
     T = Temp

     # params
     B0 = params["B0"].value
     E = params["E"].value
     El = params["El"].value
     Tl = params["Tl"].value
     Eh = params["Eh"].value
     Th = params["Th"].value
   
     # the model 
     B = (B0 * np.exp(((0-E)/K)*((1/T)-(1/283.15)))) / (1 + np.exp((El/K)*((1/Tl) - (1/T))) + np.exp((Eh/K)*((1/Th)-(1/T))))
     # this should be wrapped with np.log to get the log of the function

     # return residuals
     return B - Trait


def schoolfield(Subset, Temp, Trait, n):
     """
     Runs FUll Schoolfield model (Mechanistic)

     Optimises paramter values using minimize from lmfit package

     Calls Functions:
     fit_measure(resid_func, out, Temp, Trait)
     calc_AICc(out, Temp)
     
     Returns a list of:
     Optimised Parameters: B0, E, El, Tl, Eh, Th
     BIC
     AIC
     AICc
     Rsquared
     adjusted Rsquared
     """

     # variable values
     # Temp = np.array(Subset.ConTemp_K)
     # Trait = np.array(Subset.OriginalTraitValue)

     # estimated parameters - can change
     B0 = np.array(Subset.B0)[0]
     E = np.array(Subset.E)[0]
     El = np.array(Subset.El)[0]
     Eh = np.array(Subset.Eh)[0]
     Tl = np.array(Subset.Tl)[0]
     Th = np.array(Subset.Th)[0]
     
     # estimated params - cannot change
     B0_orig = B0
     E_orig = E
     El_orig = El
     Eh_orig = Eh
     Th_orig = Th
     Tl_orig = Tl

     # temp peak - using as a bound
     Tpeak = np.array(Subset.Tpeak)[0]

     # an initial bestfit list with an arbitarily large AIC
     #         [B0, E, El, Eh, Th, Tl, BIC, AIC ] 
     bestfit = [0, 0,  0,  0,  0,  0, 0, 100000, 0]

     # DNC - Did Not Converge flag
     # this ensures the above "best" does not get returned if none converge
     DNC = True
     #.............................................................................
     # repeat multiple times to get the best converge 
     for i in range(n):
          # this try and except block handles error (being our estimated params dont converge)
          # this ensures the code runs for n times without stoppign even if its hits an error
          try:
               if i != 0:
                    B0 = np.random.normal(B0_orig) 
                    E =  abs(np.random.normal(E_orig)) 
                    El = abs(np.random.normal(El_orig))
                    Eh = abs(np.random.normal(Eh_orig))
                    Th = np.random.normal(Th_orig) 
                    Tl = np.random.normal(Tl_orig)

               # create dictinary of parameters. Can modify attributes of each.
               params = Parameters()
               # add with tuples:(NAME, VALUE, VARY, MIN, MAX, EXPR, BRUTE_STEP)
               params.add_many(("B0", B0, True, 0, 10, None, None),
                              ("E",  E,  True, 0, 3, None, None), 
                              ("El", El, True, 0, 3, None, None),
                              ("Eh", Eh, True, 0, 6, None, None),
                              ("Th", Th, True, Tpeak, 400, None, None),
                              ("Tl", Tl, True, 270, Tpeak, None, None))

               # e and el should be between zero and minus infinity 
               # minimize residuals
               out = minimize(schoolfield_resids, params, args = (Temp, Trait))
               #...............................................................
               # write error report
               #A = report_fit(out.params)
               #...............................................................
               ## store results of best fit (based on aic score)
               if out.aic < bestfit[7]:
                    # if try gets to this point, it has converged at least once
                    DNC = False
                    # calculate AICc
                    AICc = calc_AICc(out, Temp)
                    # calculate goodness of fit measures 
                    goodness_of_fit = fit_measure(schoolfield_resids, out, Temp, Trait)
                    # bestfit takes final params and measures of fit
                    bestfit = [
                    out.params["B0"].value, out.params["E"].value, 
                    out.params["El"].value, out.params["Eh"].value, 
                    out.params["Tl"].value, out.params["Th"].value,
                    out.bic, out.aic, AICc]  
                    # merge best fit and goodness fo fit 
                    bestfit = bestfit + goodness_of_fit       
               # calculate final result to test plot
               #final = Trait + out.residual  
          except Exception as e:
               pass
               #print(e)
          #except IOError:
               #pass
          
     # print(final)
     # plt.plot(Temp, Trait, 'o')
     # plt.plot(Temp, final, 'r')
     # plt.show()
     # print(out.params)

     if not DNC:
          return bestfit
     else:
          return None 

#...........................................................................................
#   High temperature inactivation - reduced schoolfield model
#...........................................................................................

def school_high_resids(params, Temp, Trait):
     """ 
     Calculates Residuals for
     High temperature inactivation reduced schoolfield model 
     Mechanistic - Uses thermodynamics and Enzyme kinetics.

     Output:
     B   -     Trait value at given temperature (in kelvin).

     Constant:
     K   -     Boltzmann constant.

     Parameters:
     B0  -     Trait value at ref temp of 283.15 K (10 C)
     E   -     Activation energy controlling rise of curve up to the peak
     Eh  -     Enzye high temp deactivation energy
     Th  -     Enzyme 50% high-temp deactivated
     """

     # variable 
     T = Temp

     # params
     B0 = params["B0"].value
     E = params["E"].value
     Eh = params["Eh"].value
     Th = params["Th"].value

     # the model 
     B = (B0 * np.exp(((0-E)/K)*((1/T)-(1/283.15)))) / (1 + np.exp((Eh/K)*((1/Th)-(1/T))))
     # this should be wrapped with np.log to get the log of the function

     return B - Trait
 

def schoolfield_high(Subset, Temp, Trait, n):
     """
     High temperature inactivation reduced schoolfield model
     (Mechanistic)

     Optimises paramter values using minimize from lmfit package

     Calls Functions:
     fit_measure(resid_func, out, Temp, Trait)
     calc_AICc(out, Temp)
     
     Returnsa list of:
     Optimised Parameters: B0, E, Eh, Th
     BIC
     AIC
     AICc
     Rsquared
     adjusted Rsquared
     """

     # variable values
     # Temp = np.array(Subset.ConTemp_K)
     # Trait = np.array(Subset.OriginalTraitValue)
     Temp = Temp
     Trait = Trait

     # estimated parameters - can change
     B0 = np.array(Subset.B0)[0]
     E = np.array(Subset.E)[0]
     Eh = np.array(Subset.Eh)[0]
     Th = np.array(Subset.Th)[0]

     # estimated params - cannot change
     B0_orig = B0
     E_orig = E
     Eh_orig = Eh
     Th_orig = Th

     # temp peak - using as a bound
     Tpeak = np.array(Subset.Tpeak)[0]

     # an initial bestfit list with an arbitarily large AIC 
     #          B0, E, Eh, Th, BIC, AIC
     bestfit = [0, 0, 0, 0, 0, 100000, 0]

     # DNC - Did Not Converge flag
     # this ensures the above "best" does not get returned if none converge
     DNC = True
     #.............................................................................
     # repeat multiple times to get the best converge 
     for i in range(n):
          # this try and except block handles error (being our estimated params dont converge)
          # this ensures the code runs for n times without stoppign even if its hits an error
          try:
               if i != 0:
                    B0 = np.random.normal(B0_orig) 
                    E =  abs(np.random.normal(E_orig)) 
                    Eh = abs(np.random.normal(Eh_orig)) 
                    Th = np.random.normal(Th_orig) 

               # create dictinary of parameters. Can modify attributes of each.
               params = Parameters()
               # add with tuples:(NAME, VALUE, VARY, MIN, MAX, EXPR, BRUTE_STEP)
               params.add_many(("B0", B0, True, 0, 10, None, None),
                              ("E",  E,  True, 0, 3, None, None),
                              ("Eh", Eh, True, 0, 6, None, None),
                              ("Th", Th, True, Tpeak, 400, None, None))

               # minimize residuals
               out = minimize(school_high_resids, params, args = (Temp, Trait))
               #...............................................................
               # write error report
               #A = report_fit(out.params)
               #...............................................................
               ## store results of best fit (based on aic score)
               if out.aic < bestfit[5]:
                    # if try gets to this point, it has converged at least once
                    DNC = False  
                    # calculate goodness of fit measures 
                    goodness_of_fit = fit_measure(school_high_resids, out, Temp, Trait)
                    # calculate AICc
                    AICc = calc_AICc(out, Temp)
                    # bestfit takes final params and measures of fit
                    bestfit = [ 
                    out.params["B0"].value, out.params["E"].value, 
                    out.params["Eh"].value, out.params["Th"].value, 
                    out.bic, out.aic, AICc] 
                    # merge best fit and goodness fo fit 
                    bestfit = bestfit + goodness_of_fit                 
               # calculate final result to test plot
               #final = Trait + out.residual  
          except Exception as e:
               pass
               #print(e)
          #except IOError:
          #     pass
          
     # print(final)
     # plt.plot(Temp, Trait, 'o')
     # plt.plot(Temp, final, 'r')
     # plt.show()

     if not DNC:
          return bestfit
     else:
          return None 

#...........................................................................................
#     Low temperature inactivation - reduced schoolfield model
#...........................................................................................


def school_low_resids(params, Temp, Trait):
     """ 
     Calculates Residuals for reduced Schoolfield model 
     Low temperature inactivation reduced schoolfield model
     Mechanistic - Uses thermodynamics and Enzyme kinetics.

     Output:
     B   -     Trait value at given temperature (in kelvin).

     Constant:
     K   -     Boltzmann constant.

     Parameters:
     B0  -     Trait value at ref temp of 283.15 K (10 C)
     E   -     Activation energy controlling rise of curve up to the peak
     El  -     Low temp de-activation energy
     Tl  -     Enzyme 50% low-temp deactivated      
     """

     # values
     T = Temp

     # params
     B0 = params["B0"].value
     E = params["E"].value
     El = params["El"].value
     Tl = params["Tl"].value

     # the model 
     B1 = (B0 * np.exp(((0-E)/K)*((1/T)-(1/283.15)))) 
     B2 = (1 + np.exp((El/K)*((1/Tl)-(1/T))))
     B = B1 / B2
     # this should be wrapped with np.log to get the log of the function
     return B - Trait

# while rsquared is less than 0.5

def schoolfield_low(Subset, Temp, Trait, n):
     """
     Low temperature inactivation reduced schoolfield model 
     (Mechanistic)

     Optimises paramter values using minimize from lmfit package

     Calls Functions:
     fit_measure(resid_func, out, Temp, Trait)
     calc_AICc(out, Temp)
     
     Returns a list of:
     Optimised Parameters: B0, E, El, Tl
     BIC
     AIC
     AICc
     Rsquared
     adjusted Rsquared
     """

     # variable values
     # Temp = np.array(Subset.ConTemp_K)
     # Trait = np.array(Subset.OriginalTraitValue)

     # estimated parameters - can change
     B0 = np.array(Subset.B0)[0]
     E = np.array(Subset.E)[0]
     El = np.array(Subset.El)[0]
     Tl = np.array(Subset.Tl)[0]

     # estimated params - cannot change
     B0_orig = B0
     E_orig = E
     El_orig = El
     Tl_orig = Tl

     # temp peak - using as a bound
     Tpeak = np.array(Subset.Tpeak)[0]

     # an initial bestfit list with an arbitarily large AIC 
     #         [B0, E, El, Tl, BIC, AIC]
     bestfit = [0, 0, 0, 0, 0, 100000, 0]

     # DNC - Did Not Converge flag
     # this ensures the above "best" does not get returned if none converge
     DNC = True
     #.............................................................................
     # repeat multiple times to get the best converge 
     for i in range(n):
          # this try and except block handles error (being our estimated params dont converge)
          # this ensures the code runs for n times without stoppign even if its hits an error
          try:
               if i != 0:
                    B0 = np.random.normal(B0_orig) 
                    E =  abs(np.random.normal(E_orig))
                    El = abs(np.random.normal(El_orig))
                    Tl = np.random.normal(Tl_orig)
                    
               # create dictinary of parameters. Can modify attributes of each.
               params = Parameters()
               # add with tuples:(NAME, VALUE, VARY, MIN, MAX, EXPR, BRUTE_STEP)
               params.add_many(("B0", B0, True, 0, 10, None, None),
                              ("E",  E,  True, 0, 3, None, None), 
                              ("El", El, True, 0, 3, None, None),
                              ("Tl", Tl, True, 270, Tpeak, None, None))
               # minimize residuals
               out = minimize(school_low_resids, params, args = (Temp, Trait))
               #...............................................................
               # write error report
               #A = report_fit(out.params)
               #..............................................................
               #...............................................................
               ## store results of best fit (based on aic score)
               if out.aic < bestfit[5]:
                    # if try gets to this point, it has converged at least once
                    DNC = False
                    # calculate goodness of fit measures 
                    goodness_of_fit = fit_measure(school_low_resids, out, Temp, Trait)
                    # calculate AICc
                    AICc = calc_AICc(out, Temp)
                    # bestfit takes final params and measures of fit
                    bestfit = [
                    out.params["B0"].value, out.params["E"].value, 
                    out.params["El"].value, out.params["Tl"].value, 
                    out.bic, out.aic, AICc]
                    # merge best fit and goodness fo fit 
                    bestfit = bestfit + goodness_of_fit                
               # calculate final result to test plot
               #final = Trait + out.residual  
          except Exception as e:
               pass
               #print(e)
          #except IOError:
               #pass
          
     # print(final)
     # print(bestfit)
     # plt.plot(Temp, Trait, 'o')
     # plt.plot(Temp, final, 'r')
     # plt.show()

     if not DNC:
          return bestfit
     else:
          return None 

