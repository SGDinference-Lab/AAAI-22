# SGD Inference Simulation Experiments
This repository contains R codes to replicate the results in the following paper:

**Fast and Robust Online Inference with Stochastic Gradient Descent via Random Scaling**

Last update: 2021-10-02

## Notes
1. The **linear** directory contains all source codes for the linear regression model.
2. The **logistic** directory contains all source codes for the logistic regression model.
3. The **logistic_single_parameter** directory contains all source codes for the logistic regression model, where we update only a single element of the variance-covariance matrix.
4. In each simulation design, 1000 replications are divided into 100 codes with 10 replications. In each code, the Seed number is set to be from 15673 to 15772. 


## How to replicate the results of linear regression
1. Run **sgd_linear.R** with appropriate input arguments. 
	1. If you have access a cluster system, use appropriate scripts to run all 100 codes of each design.
	2. Otherwise, type in the terminal: **Rscript sgd_linear.R 15673 05 0.505 0.5 01 TRUE 1**
	3. Please see the comments section in the code **sgd_linear.R** for the input argument dictionary.
	4. In each run, increase the seed number(15673) by one up to 15772.
2. The results of simulations will be collected as .RData files in the subfolder. For example, ../d-05/d-05-01/ will collect all .RData files of the simulation above (d=5, gamma=1, and alpha=0.505). 
3. Go to the subfolder of each design and run **gen_Graph_Table_dopar_simple.R**. This code will generate graphs of the coverage rates, CI lengths, and the computation time as well as a table of the summary statistics. 

## How to replicate the results of logistic regression
1. The codes can be run similarly as above.
2. The codes in the **logistic_single_parameter** directory compute only **Random Scale** method. So, please ignore summary statistics for other methods in the table. 


