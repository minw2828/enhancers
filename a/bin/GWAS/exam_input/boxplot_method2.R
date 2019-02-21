# /usr/bin/rscript

## Description:
## This script plots an boxplot + error bar in DGAT1 effects on phenotype.
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 07 March 2015
## 
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <date> <psf> <gender> <pheno> <infile> 


## Generate boxplot using plotting method 2: base plot ##
png(filename = outfile, width = 1200, height = 1500)
par(mfrow=c(3,1))
boxplot(formula = FY~DGAT1_geno, data = data, main = "DGAT1 genotype effects on FY", xlab = "DGAT1 genotypes", ylab = "Fat Yield")
boxplot(formula = MY~DGAT1_geno, data = data, main = "DGAT1 genotype effects on MY", xlab = "DGAT1 genotypes", ylab = "Milk Yield")
boxplot(formula = PY~DGAT1_geno, data = data, main = "DGAT1 genotype effects on PY", xlab = "DGAT1 genotypes", ylab = "Protein Yield")
dev.off()



