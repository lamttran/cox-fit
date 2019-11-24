# coxfit

This is the R package for Homework 4 of Biostat 625, Fall 2019. The Cox proportional hazards model is a commonly used regression model in the context of survival data (design matrix x, death or censoring time y, and censoring vector d). If the proportional hazards assumption is met, that is the explanatory variables act solely on the baseline hazard function and not the survival/failure times, the Cox model makes no assumptions on the underlying distribution of these times. As such, the Cox model is a semi-parametric model.

Briefly, this package fits a Cox proportional hazards model to an input of survival data (in this case, the design matrix and a censoring vector), assuming the data is sorted in order of increasing survival time/time to death. 

## Getting Started
A simple example for using coxfit with simulated survival data is described below:
1. Install the package from GitHub and load/attach it to R: 
````
install.packages("devtools")
library(devtools)
install_github("lamttran/coxfit", build_vignettes = T)
library(coxfit)
````
2. Simulate some survival data:
````
set.seed(2019)
p = 5 
n = 100 
c0 = 0.5 
beta = rnorm(p, 0, 1) 
x = matrix(rnorm(n * length(beta)), ncol = length(beta)) 
yvar = rnorm(n, 0, 1) 
times = rexp(n, c0 * exp(rowSums(t(t(x) * beta)) + yvar)) 
time.censor = rexp(n, c0 * exp(beta[1] * x[, 1] + yvar)) 
censorv = ifelse(times < time.censor, 1, 0) 
y = ifelse(times < time.censor, times, time.censor) 
ord = order(y, decreasing =F) 
dataset = list(x = x[ord,], y = y[ord], censor = censorv[ord]) 
````
3. Run fit_cox and compare to true beta values:
````
fit_cox(dataset$x, dataset$censor)
beta
````

