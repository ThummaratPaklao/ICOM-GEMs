# ICOM-GEMs
Integration of gene co-expression network and metabolic network based on flux balance anaysis

This document covers instruction on how to run the integration of gene co-expression network and metabolic model method in MATLAB.

(or you can access in the cobra toolbox : https://opencobra.github.io/cobratoolbox/latest/modules/analysis/ICONGEMs/index.html)

## REQUIREMENTS ##
1. Matlab (version 2018a or better)
2. Cobra Toolbox
2. Gurobi solver (version 9.0.1 or better, free academic)
3. Gene expression profile 
  Note that the first column of gene expression data should have gene symbols/names used in the GPR association of the genome scale metabolic model. First row of gene expression data should have condition names.

## Usage ##
Step 1. Open Matlab program and enter:

    >> initCobraToolbox;

at Matlab command line.

Step 2. Read SBML model

    >> model = readCbModel('model_name');

Step 3. Read gene expression data

    >> [exp txt] = xlsread('gene_expression_profile_file.csv');

  Where exp is the numeric value of gene expression data in gene_expression_profile_file.csv file and txt is text data in gene_expression_profile_file.csv. Numeric values in inner spreadsheet rows and columns appear as empty character vectors in txt.

Step 4. Perform analysis

    >> [solICONGEMs,boundEf]  = ICONGEMs(model, exp, txt, condition, threshold,alpha, numericalFlag);

  The algorithm of integration of co-expression network and metabolic model is completed by using function ICONGEMs where the input is model file in step2, exp and txt in step 3 and row vector of condition that are wanted to calculate flux distribution (default is all conditions) and threshold for constructing co-expression network (default value is 0.9). The alpha value is the proportion of biomass (value in range (0,1]).Parameter numericFlag is 1 if using Human Recon  (Default = 0).

  After the algorithm is finished, solICONGEMs for the predicted metabolic fluxes will be added to the Workspace. Numerical flux values can be examined in more detail by double-clicking solICONGEMs. The boundEf is upperbound of E-flux method. Moreover, the output of this algorithm is reported in result.csv file.

## Citation ##
Paklao, T., Suratanee, A. & Plaimas, K. ICON-GEMs: integration of co-expression network in genome-scale metabolic models, shedding light through systems biology. BMC Bioinformatics 24, 492 (2023). https://doi.org/10.1186/s12859-023-05599-0


