# Death Rate Predictors Analysis
Relationship investigation between demographic variables and death rate using R.


## Dataset

The dataset was assembled in 1960.  Each observation is on a metropolitan area in the U.S.  The outcome variable of interest is the death rate per 100,000 population.   There are several kinds of predictor variables -- weather, demographic, economic, and air-pollution related.

popDens -- population per square mile

nonWh -- percent of population that is non-White

educ -- average number of years of education of adults

death -- age-adjusted deaths per 100,000 population

[Dataset Link](http://homepage.divms.uiowa.edu/~kcowles/Datasets/deathRate.txt)


## Analysis

This analysis investigates the relationship between three demographic variables and the death rate. Only above mentioned 4 variables in the dataset section
have been used for this analysis.


JAGS library has been used to fit several Bayesian models --  simple  models using one predictor variable at a time, and models with 2 predictors and then all 3 predictors.  

Performance of these models were compared using Bayesian criterion (Deviance information criterion and Posterior predictive checking) to choose the best one.

Detailed analysis can be found at [analysis](https://github.com/muneeb706/Deathrate_Predictors_Analysis/blob/main/analysis.pdf).

R-script used for analysis can be found at [r-script](https://github.com/muneeb706/Deathrate_Predictors_Analysis/blob/main/deathrate-predictors-analysis.R).
