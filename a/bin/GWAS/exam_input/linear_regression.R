# /usr/bin/rscript

## Description:
## This script takes a csv file as input and output a plot to visualize 
## input statistics.
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 22 Feb 2016
## 
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <date> <psf> <infile1> <infile2> <outfile1> 


## Linear Regression model: Test DGAT1 genotype effects on FY, MY, PY ##
get_lm_fit <- function(x, y) {
    mean.pheno <- mean(y)
    fit <- lm(y ~ x-1)
    summary(fit) # show results
    coefficients(fit) # model coefficients
    confint(fit, level=0.95) # CIs for model parameters 
    fitted(fit) # predicted values
    residuals(fit) # residuals
    anova(fit) # anova table 
    vcov(fit) # covariance matrix for model parameters 
    influence(fit) # regression diagnostics 
    return(summary(fit))
}

sfy <- get_lm_fit(data$DGAT1_geno, data$FY)
smy <- get_lm_fit(data$DGAT1_geno, data$MY)
spy <- get_lm_fit(data$DGAT1_geno, data$PY)

capture.output(sfy, file = outfile2, append = FALSE)
capture.output(smy, file = outfile2, append = TRUE)
capture.output(spy, file = outfile2, append = TRUE)






