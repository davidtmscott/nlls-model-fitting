# Non-linear least squares model fitting
## Author: David Scott 
## Date: _Mar - 2019_

### Computing mini project completed as part of MSc 

#### Description: 
All scripts are stored in the code directory, all results produced are stored in the output directory and all data files used initially or created were stored in the data directory. LaTeX document for compiling the report (pdf format) and bibtex for references are stored in the Writeup directory. To run the project, run the bash script found in the code direcoty. This runs the entire project workflow from the intial data to the final LaTeX document, providing meaningfull output throughout. 

__Languages__ : Python, R, Bash, LaTeX

__Models__ : Cubic Polynomial, Briere, Full Schoolfield, Schoolfield Low, Schoolfield High

A total of five models were fitted to the data, of which two were phenomenological and three were mechanistic.  

#### Phenomenological Models
The general **Cubic Polynomial model** is phenomenological and thus it's parameters $(B_{0}, B_{1}, B_{2}, B_{3})$, do not describe any underlying biological processes. Here, $B$ is the predicted trait value at a given temperature $T$ in Celsius. 

$$
B = B_{0} + B_{1}T + B_{2}T + B_{3}T
$$

The **Briere model** is also phenominological. Temperature $T$  was taken in Celsius. $B$ is the predicted trait value. The parameter $B_{0}$ is a normalisation constant. The minimum and maximum feasible temperature beyond which the trait value goes to zero (is no longer active) are represented by parameters $T_{0}$ and $T_{m}$ respectively. 

$$
B = B_{0}T(T - T_0)\sqrt{T_{m} - T}
$$

#### Mechanistic Models
The **Schoolfield model**, models the effect of thermodynamics and enzyme kinetics on metabolic rates as a result of temerature ([Schoolfield et al., 1981](https://doi.org/10.1016/0022-5193(81)90246-0)). Temperature was taken in Kelvin. 

$$
B = \frac{B_{0}e^{\frac{-E}{k}(\frac{1}{T} - \frac{1}{283.15})}} {1 + e^{\frac{E_{l}}{k}(\frac{1}{T_{l}} - \frac{1}{T}} + e^{\frac{E_{h}}{k}(\frac{1}{T_{h}} - \frac{1}{T})}}
$$


The model returns a value of metabolism ($B$), at a given temperature $T$ based on the parameter values: Activation energy ($E$), Low temperature deactivation energy ($E_l$), Temperature for which the enzyme is half active and half low-temperature suppressed ($T_l$), High temperature deactivation energy ($E_h$), Temperature at which the enzyme is half active and half high-temperature suppressed ($T_h$). $B_0$ is the trait value at a reference temperature of 283.15 Kelvin ($10^{\circ}C$). This represents the growth rate at low temperature, controlling the offset of the curve. $K$ is the Boltzmann constant ($8.617 \times 10^{-5}$ eV $\cdot K^{-1}$).

The Schoolfield high temperature inactivation energy model (**Schoolfield High**) is a simplified version of the Schoolfield model, used when low temperature deactivation energy is weak. 

$$
B = \frac{B_{0}e^{\frac{-E}{k}(\frac{1}{T} - \frac{1}{283.15})}} {1 + e^{\frac{E_{h}}{k}(\frac{1}{T_{h}} - \frac{1}{T})}}
$$

Likewise, a further simplified version of the Schoolfield model can be used in cases where high temperature inactivation of the trait is weak or not recorded. This is the Schoolfield low temperature inactivation energy model (**Schoolfield Low**). 

$$
B = \frac{B_{0}e^{\frac{-E}{k}(\frac{1}{T} - \frac{1}{283.15})}} {1 + e^{\frac{E_{l}}{k}(\frac{1}{T_{l}} - \frac{1}{T})}}
$$
	
__Parameters__ : 

* Estimated starting parameters are stored in ../data/BioTraits_Params.csv
* Optimised parameters (post minimisation) are stored in ../data/BioTraits_Params.csv

__Approximate Total Run Time:__

#### Workflow Overview:
The project runs in the following sequence:

##### 1. wrangle.R
R script that imports raw data (BioTraits.csv) from the data directory and filters data. It then estimates starting parameters for the Briere model and all Schoolfield models. Saves data as BioTraits_Params.csv to data directory. 

* Packages used: dplyr
* Imported Data: ../data/BioTraits.csv
* Outputted Data: ../data/BioTraits_Params.csv

##### 2. model_functions.py
Python script that defines functions to fit models to TPC's. Defines functions to return residual between observed and predicted values for each model. Also defines a function that uses lmfit minimise function to call the residual function and adjusts the parameter values to minimise the residuals. 

##### 3. NLLS_fitting.py 
Python script, loops through each unique TPC and uses NLLS to fit each of the five models.
Saves data as BioTraits_FinalParams.csv to data directory. 

imports functions from model_functions.py

* Packages used: lmfit, pandas, numpy, time
* Scripts used: model_functions.py (also in code Directory, see #2)
* Imported Data: ../data/BioTraits_Params.csv
* Outputted Data: ../data/BioTraits_FinalParams.csv 

##### 4. plotting.R
R script that imports BioTraits_FinalParams.csv from data directory and Plots all models that converged, create tables for LaTeX (saved to output directory) and produces some overall statistics. 

* Packages used: dplyr, reshape2, ggplot2, xtable
* Imported Data: ../data/BioTraits_FinalParams.csv 
                 ../data/BioTraits.csv
* Outputted Data: ../output/Model_Plots/MTD####.pdf 

##### 5. Writeup.tex
LaTeX document

#### Tree map
```
.
├── code
│   ├── model_functions.py :              Python        Defines models as functions to be  fitted to TPC using Python3
│   ├── NLLS_fitting.py :                 Python        Uses NLLS method to fit models to BioTraits data in Python3
│   ├── plotting.R :                      R             Plotting, Tables and Results Script
│   ├── run_project.sh :                  Bash          Bash to run each component of project 
│   └── wrangle.R :                       R             Wrangles biotraits data and estimate starting parameter values 
├── data
│   ├── BioTraits.csv                     -     Original data, used wrangle.R
│   ├── BioTraits_FinalParams.csv         -     Produced by NLLS_fitting.py, contains optimised parameters post minimisation. Is used in Plotting.R
│   └── BioTraits_Params.csv              -     Produced by wrangle.R, contains estimated parameters. Is used in NLLS_fitting.py
├── README.tmp
├── Report
│   ├── Writeup.bib                       -     Bibtex Refences
│   └── Writeup.tex                       -     LaTeX script
└── output
    ├── AIC_Plots
    │   └── Lowest_AICscores.pdf
    ├── Curves.tex
    ├── Model_Plots
    │   ├── MTD2079.pdf
    │   ├── ....
    │   ├── Many many more plots
    │   ├── ....
    │   └── MTD608.pdf
    ├── Params.tex
    ├── Points.tex
    ├── Scores.tex
    └── Writeup.pdf

8 directories, 1006 files

```
