# #!/bin/bash
# # Author: David Scott 
# # Script: run_project.sh
# # Desc: script to run each component of MiniProject
# # Arguments: none
# # Date: February 2019

echo -e "\n$(tput setaf 3)Hello World$(tput sgr 0)"

STARTTIME=$(date +%s)

echo -e "\n$(tput setaf 4)Running wrangle.R\n" 
Rscript wrangle.R
echo -e "\n$(tput setaf 4)Finished Running wrangle.R$(tput sgr 0)\n" 

# might have to make the python scripts executable
echo -e "\n$(tput setaf 2)Running model_functions.py$(tput sgr 0)\n" 
python3 model_functions.py 
echo -e "\n$(tput setaf 2)Finished Running model_functions.py$(tput sgr 0)\n" 

echo -e "\n$(tput setaf 2)Running NLLS_fitting.py \n" 
python3 NLLS_fitting.py 
echo -e "\n$(tput setaf 2)Finished Running NLLS_fitting.py$(tput sgr 0)\n" 

echo -e "\n$(tput setaf 4)Running plotting.R \n" 
Rscript plotting.R
echo -e "\n$(tput setaf 4)Finished Running plotting.R$(tput sgr 0)\n"

echo -e "\n$(tput setaf 5)Compiling LaTeX$(tput sgr 0)\n" 
pdflatex ../Report/Writeup.tex   
pdflatex ../Report/Writeup.tex 
bibtex Writeup                  # bibtex cannot have file with .tex
pdflatex ../Report/Writeup.tex 
pdflatex ../Report/Writeup.tex 
mv Writeup.pdf ../Results       # moves output (pdf version of texfile) to Results
evince ../Results/Writeup.pdf & # displays pdf output on screen 

## Cleanup - Removing all other files created other than .bib .pdf and .tex
#remove as these files are not required
rm *~
rm *.aux
rm *.bbl   
rm *.blg    
rm *.dvi
rm *.log
rm *.nav 
rm *.out
rm *.snm
rm *.toc

echo -e "\n$(tput setaf 5)LaTeX Compiled and unwanted files removed$(tput sgr 0).\n" 


ENDTIME=$(date +%s)

echo "It takes $(($ENDTIME - $STARTTIME)) seconds to run the miniproject"

echo -e "$(tput setaf 3)Goodbye World$(tput sgr 0)\n"
