library(fitdistrplus)
library(logspline)

deathRateData <- read.table("http://homepage.divms.uiowa.edu/~kcowles/Datasets/deathRate.txt", header=TRUE)

# determining the distribution of the variable of interest (death rate per 100,000 population)

hist(x = deathRateData$death, freq = FALSE, main="Histogram of the death rate data", xlab = "death rate", ylab = "Density")
lines(x = density(deathRateData$death), col="red")

# from above plot distribution looks like normal.

library(rjags)
load.module("glm")
load.module("dic")

# dic

## model with one predictor variable

deathRateCodeOne <- "
    data {
    
        x_mean <- mean(x)
        
        for (i in 1:n) {
            xC[i] <- x[i] - x_mean
        }
    }
    
    model {
        
        for (i in 1:n) {
            y[i] ~ dnorm(mu[i], tausq)
            ypred[i] ~ dnorm(mu[i], tausq) # predicted values
            aresid[i] <- y[i] - mu[i] # actual residuals
            arespred[i] <- ypred[i] - mu[i] # predicted residuals
            mu[i] <- beta0 + beta1*xC[i]
            
        }
        
        # calculate lag 1 autocorrelation in residuals from real data
        
        mean1 <- mean(aresid[1:(n-1)])
        mean2 <- mean(aresid[2:n])
        
        for (i in 1:(n-1)) {
            summand[i] <- (aresid[i] - mean1) * (aresid[i+1] - mean2)
        }
        
        lag1auto <- sum(summand[]) / ((n-1)*sd(aresid[1:(n-1)])*sd(aresid[2:n]))
        
        # calculate lag 1 autocorrelation in residuals from replicated/predicted data
        
        mean1pred <- mean(arespred[1:(n-1)])
        mean2pred <- mean(arespred[2:n])
        
        for (i in 1:(n-1)) {
            summandpred[i] <- (arespred[i] - mean1pred) * (arespred[i+1] - mean2pred)
        }
        
        lag1autopred <- sum(summandpred[]) / ((n-1)*sd(arespred[1:(n-1)])*sd(arespred[2:n]))
        
        # is lag1auto > lag1autopred?
        pppval <- step(lag1auto -lag1autopred)
        
        beta0 ~ dnorm(0,0.000001)
        beta1 ~ dnorm(0,0.000001)
        tausq ~ dgamma(0.001,0.001)
    
    }

"

## poulation density as a predictor variable

deathRateModelData <- list( x = deathRateData$educ,
                            y = deathRateData$death, n = length(deathRateData$death))

## fit frequentist models to get ideas for initial values for JAGS

respmat <- cbind(deathRateModelData$y)
xcent <- deathRateModelData$x - mean(deathRateModelData$x)

summary(glm( respmat ~ xcent, family = gaussian(link="identity")))

## inits for logit model
deathRateModelInits <- list(list(beta0 = 912.67, beta1 = -70.8, ypred=rnorm(60), .RNG.name="base::Wichmann-Hill", .RNG.seed=350),
                            list(beta0 = 968.03, beta1 = -4.4, ypred=rnorm(60),.RNG.name="base::Wichmann-Hill", .RNG.seed=32542),
                            list(beta0 = 995.71, beta1 = 28.8, ypred=rnorm(60), .RNG.name="base::Wichmann-Hill", .RNG.seed=5280))

deathRateModel <- jags.model(textConnection(deathRateCodeOne), data = deathRateModelData,
                             inits = deathRateModelInits, n.chains = 3, n.adapt = 0)

deathRateModelResults <- coda.samples(deathRateModel, 
                                      c("beta0", "beta1", "tausq", "lag1auto", "lag1autopred", "pppval"),
                                      n.iter = 10000)

par(mfrow=c(3,2))
traceplot(deathRateModelResults)
gelman.plot(deathRateModelResults, auto.layout = FALSE)
gelman.diag(deathRateModelResults, autoburnin=TRUE,transform=TRUE, multivariate = FALSE)

summary(window(deathRateModelResults, start=5001))

logit.dic.pd <- dic.samples(deathRateModel,n.iter = 3000,type = "pD")
logit.dic.popt <- dic.samples(deathRateModel,n.iter = 3000,type = "popt")
logit.dic.pd
logit.dic.popt



## model with two predictor variables

deathRateCodeTwo <- "
    data {
        x1_mean <- mean(x1)
        x2_mean <- mean(x2)
        for (i in 1:n) {
            x1C[i] <- x1[i] - x1_mean
            x2C[i] <- x2[i] - x2_mean
        }
    }
    model {
        for (i in 1:n) {
            y[i] ~ dnorm(mu[i], tausq)
            ypred[i] ~ dnorm(mu[i], tausq) # predicted values
            aresid[i] <- y[i] - mu[i] # actual residuals
            arespred[i] <- ypred[i] - mu[i] # predicted residuals
            mu[i] <- beta0 + beta1*x1C[i] + beta2*x2C[i]
        }

        # calculate lag 1 autocorrelation in residuals from real data
        
        mean1 <- mean(aresid[1:(n-1)])
        mean2 <- mean(aresid[2:n])
        
        for (i in 1:(n-1)) {
            summand[i] <- (aresid[i] - mean1) * (aresid[i+1] - mean2)
        }
        
        lag1auto <- sum(summand[]) / ((n-1)*sd(aresid[1:(n-1)])*sd(aresid[2:n]))
        
        # calculate lag 1 autocorrelation in residuals from replicated/predicted data
        
        mean1pred <- mean(arespred[1:(n-1)])
        mean2pred <- mean(arespred[2:n])
        
        for (i in 1:(n-1)) {
            summandpred[i] <- (arespred[i] - mean1pred) * (arespred[i+1] - mean2pred)
        }
        
        lag1autopred <- sum(summandpred[]) / ((n-1)*sd(arespred[1:(n-1)])*sd(arespred[2:n]))
        
        # is lag1auto > lag1autopred?
        
        pppval <- step(lag1auto -lag1autopred)
        beta0 ~ dnorm(0,0.000001)
        beta1 ~ dnorm(0,0.000001)
        beta2 ~ dnorm(0,0.000001)
        tausq ~ dgamma(0.001,0.001)
    }
"

## poulation density and percentage of non white population as a predictor variables

deathRateModelData <- list( x1 = deathRateData$nonWh, x2 = deathRateData$educ,
                            y = deathRateData$death, n = length(deathRateData$death))

## fit frequentist models to get ideas for initial values for JAGS

respmat <- cbind(deathRateModelData$y)
x1cent <- deathRateModelData$x1 - mean(deathRateModelData$x1)
x2cent <- deathRateModelData$x2 - mean(deathRateModelData$x2)

summary(glm( respmat ~ x1cent + x2cent, family = gaussian(link="identity")))

## inits for logit model

deathRateModelInits <- list(list(beta0 = 918.8, beta1 = 1.42, beta2 = -55.33, ypred=rnorm(60), .RNG.name="base::Wichmann-Hill", .RNG.seed=350),
                            list(beta0 = 962.0, beta1 = 6.41, beta2 = -2.62, ypred=rnorm(60),.RNG.name="base::Wichmann-Hill", .RNG.seed=32542),
                            list(beta0 = 983.6, beta1 = 8.9, beta2 = 23.7, ypred=rnorm(60), .RNG.name="base::Wichmann-Hill", .RNG.seed=5280))

deathRateModel <- jags.model(textConnection(deathRateCodeTwo), data = deathRateModelData,
                             inits = deathRateModelInits, n.chains = 3, n.adapt = 0)

deathRateModelResults <- coda.samples(deathRateModel, 
                                      c("beta0", "beta1", "beta2", "tausq", "lag1auto", "lag1autopred", "pppval"),
                                      n.iter = 10000)

par(mfrow=c(3,2))
traceplot(deathRateModelResults)
gelman.plot(deathRateModelResults, auto.layout = FALSE)
gelman.diag(deathRateModelResults, autoburnin=TRUE,transform=TRUE, multivariate = FALSE)

summary(window(deathRateModelResults, start=5001))

logit.dic.pd <- dic.samples(deathRateModel,n.iter = 3000,type = "pD")
logit.dic.popt <- dic.samples(deathRateModel,n.iter = 3000,type = "popt")
logit.dic.pd
logit.dic.popt





## model with three predictor variables

deathRateCodeThree <- "
    data {
        x1_mean <- mean(x1)
        x2_mean <- mean(x2)
        x3_mean <- mean(x3)
        for (i in 1:n) {
            x1C[i] <- x1[i] - x1_mean
            x2C[i] <- x2[i] - x2_mean
            x3C[i] <- x3[i] - x3_mean
        }
    }
    model {
        for (i in 1:n) {
            y[i] ~ dnorm(mu[i], tausq)
            ypred[i] ~ dnorm(mu[i], tausq) # predicted values
            aresid[i] <- y[i] - mu[i] # actual residuals
            arespred[i] <- ypred[i] - mu[i] # predicted residuals
            mu[i] <- beta0 + beta1*x1C[i] + beta2*x2C[i] + beta3*x3C[i]
        }

        # calculate lag 1 autocorrelation in residuals from real data
        
        mean1 <- mean(aresid[1:(n-1)])
        mean2 <- mean(aresid[2:n])
        
        for (i in 1:(n-1)) {
            summand[i] <- (aresid[i] - mean1) * (aresid[i+1] - mean2)
        }
        
        lag1auto <- sum(summand[]) / ((n-1)*sd(aresid[1:(n-1)])*sd(aresid[2:n]))
        
        # calculate lag 1 autocorrelation in residuals from replicated/predicted data
        
        mean1pred <- mean(arespred[1:(n-1)])
        mean2pred <- mean(arespred[2:n])
        
        for (i in 1:(n-1)) {
            summandpred[i] <- (arespred[i] - mean1pred) * (arespred[i+1] - mean2pred)
        }
        
        lag1autopred <- sum(summandpred[]) / ((n-1)*sd(arespred[1:(n-1)])*sd(arespred[2:n]))
        
        # is lag1auto > lag1autopred?
        
        pppval <- step(lag1auto -lag1autopred)
        beta0 ~ dnorm(0,0.000001)
        beta1 ~ dnorm(0,0.000001)
        beta2 ~ dnorm(0,0.000001)
        beta3 ~ dnorm(0,0.000001)
        tausq ~ dgamma(0.001,0.001)
    }
"


deathRateModelData <- list( x1 = deathRateData$popDens, x2 = deathRateData$nonWh, x3 = deathRateData$educ, 
                            y = deathRateData$death, n = length(deathRateData$death))

## fit frequentist models to get ideas for initial values for JAGS

respmat <- cbind(deathRateModelData$y)
x1cent <- deathRateModelData$x1 - mean(deathRateModelData$x1)
x2cent <- deathRateModelData$x2 - mean(deathRateModelData$x2)
x3cent <- deathRateModelData$x3 - mean(deathRateModelData$x3)

summary(glm( respmat ~ x1cent + x2cent + x3cent, family = gaussian(link="identity")))

## inits for logit model

deathRateModelInits <- list(list(beta0 = 919.4, beta1 = -0.0072, beta2 = 1.59, beta3 = -51.4,
                                 ypred=rnorm(60), .RNG.name="base::Wichmann-Hill", .RNG.seed=350),
                            list(beta0 = 961.32, beta1 = 0.0228, beta2 = 6.39, beta3 = 1.4,
                                 ypred=rnorm(60),.RNG.name="base::Wichmann-Hill", .RNG.seed=32542),
                            list(beta0 = 982.28, beta1 = 0.0378, beta2 = 8.79, beta3 = 27.8,
                                 ypred=rnorm(60), .RNG.name="base::Wichmann-Hill", .RNG.seed=5280))

deathRateModel <- jags.model(textConnection(deathRateCodeThree), data = deathRateModelData,
                             inits = deathRateModelInits, n.chains = 3, n.adapt = 0)

deathRateModelResults <- coda.samples(deathRateModel, 
                                      c("beta0", "beta1", "beta2", "beta3", "tausq", "lag1auto", "lag1autopred", "pppval"),
                                      n.iter = 10000)

par(mfrow=c(3,2))
traceplot(deathRateModelResults)
gelman.plot(deathRateModelResults, auto.layout = FALSE)
gelman.diag(deathRateModelResults, autoburnin=TRUE,transform=TRUE, multivariate = FALSE)

summary(window(deathRateModelResults, start=5001))

logit.dic.pd <- dic.samples(deathRateModel,n.iter = 3000,type = "pD")
logit.dic.popt <- dic.samples(deathRateModel,n.iter = 3000,type = "popt")
logit.dic.pd
logit.dic.popt