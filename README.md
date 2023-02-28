# MD_glasso
For the subject "Modelling under dependence" WS22/23 at LMU. 
Contains code and visualisations used for the Graphical Lasso section of the project "Graphical lasso and Mixed Graphical Model".

Authors: <br />
Patricio Massaro <br />
Dusan Urosevic

# Contents
- glassoProtein.R - main file that contains all visualisations from the presentation and final report.
- gaussian_experiment.R - syntehtic data experiments where we try to reconstruct a gaussian model using Friedmans implementation of the Graphical Lasso (also called Covariance Lasso).
- test.R and sample_mixed_model.R - not used for the final report. The first is just a test of visualisation packages, the second was an attempt at recovering the mixed model using the MGM package

# Installation
The R files can be run as is in R studio. Each file starts with installations of the necessary packages. **Note** - when constructing the dataframes we assume that the folder structure is unchanged (all data files are in the dataGLasso or data folder).
