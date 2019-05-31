# Partitioned-Survival-Analysis
This repository contains functions useful for conducting partitioned survival analysis in R and provides a few examples.

The file 'SurvFunctions_Final' contains 2 functions:

1. fit.fun is the function that fits a number of parametric survival models to the data, estimates the Akaike and bayesian information criterion and stores all the output on a list of survival models and their goodness of fit. 

2. partsurv is the function that partitions the survival probability estimates into progression-free, progressed, and overall survival.

The file 'Partitioned_Report.html' illustrates the use of the functions through an artificial example.

