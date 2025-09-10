# Non-linear least squares model fitting
## Author: David Scott 
## Date: _Mar - 2019_

### Computing mini project completed as part of MSc 

#### Description: 
All scripts are stored in the Code directory, all results produced are stored in the Results directory and all data files used initially or created were stored in the Data directory. LaTeX document for compiling the report (pdf format) and bibtex for references are stored in the Writeup directory. To run the project, run the bash script found in the Code direcoty. This runs the entire project workflow from the intial data to the final LaTeX document, providing meaningfull output throughout. 

__Languages__ : Python, R, Bash, LaTeX

__Models__ : Cubic Polynomial, Briere, Full Schoolfield, Schoolfield Low, Schoolfield High

__Parameters__ : 

* Estimated starting parameters are stored in ../data/BioTraits_Params.csv
* Optimised parameters (post minimisation) are stored in ../data/BioTraits_Params.csv

__Approximate Total Run Time: __

#### Workflow Overview:
The project runs in the following sequence:

##### 1. Wrangle.R
R script that imports raw data (BioTraits.csv) from the Data directoryand filters data. It then estimates starting parameters for the Briere model and all Schoolfield models. Saves data as BioTraits_Params.csv to Data directory. 

* Packages used: dplyr
* Imported Data: ../data/BioTraits.csv
* Outputted Data: ../data/BioTraits_Params.csv

##### 2. model_functions.py
Python script that defines functions to fit models to TPC's. Defines functions to return residual between observed and predicted values for each model. Also defines a function that uses lmfit minimise function to call the residual function and adjusts the parameter values to minimise the residuals. 

##### 3. NLLS_fitting.py 
Python script, loops through each unique TPC and uses NLLS to fit each of the five models.
Saves data as BioTraits_FinalParams.csv to Data directory. 

imports functions from model_functions.py

* Packages used: lmfit, pandas, numpy, time
* Scripts used: model_functions.py (also in Code Directory, see #2)
* Imported Data: ../data/BioTraits_Params.csv
* Outputted Data: ../data/BioTraits_FinalParams.csv 

##### 4. Plotting.R
R script that imports BioTraits_FinalParams.csv from Data directory and Plots all models that converged, create tables for LaTeX (saved to Results directory) and produces some overall statistics. 

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
│   ├── run_project.sh :              Bash          Bash to run each component of MiniProject 
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
