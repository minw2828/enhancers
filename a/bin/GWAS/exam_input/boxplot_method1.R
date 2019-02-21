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


## Generate boxplot using plotting method 1: ggplot2 ##
aim <- 'Boxplot: DGAT1 genotype effects on'
Xlab  <- 'DGAT1 genotypes'

Title <- paste(aim, 'FY')
Ylab  <- 'Fat Yield'
g1 <- ggplot(data = data) +
      geom_boxplot(mapping = aes(factor(DGAT1_geno), FY, fill = factor(DGAT1_geno))) + 
      labs(title = Title, x = Xlab, y = Ylab) +
      theme(plot.title = element_text(face = "bold", size = 25),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.text.x = element_text(size = 15),
            axis.text.y = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.title = element_text(size = 15)) 

Title <- paste(aim, 'MY')
Ylab  <- 'Milk Yield'
g2 <- ggplot(data = data) +
     geom_boxplot(mapping = aes(factor(DGAT1_geno), MY, fill = factor(DGAT1_geno))) +
     labs(title = Title, x = Xlab, y = Ylab) +
     theme(plot.title = element_text(face = "bold", size = 25),
           axis.title.x = element_text(size = 20),
           axis.title.y = element_text(size = 20),
           axis.text.x = element_text(size = 15),
           axis.text.y = element_text(size = 15),
           legend.text = element_text(size = 15),
           legend.title = element_text(size = 15))

Title <- paste(aim, 'PY')
Ylab  <- 'Protein Yield'
g3 <- ggplot(data = data) +
     geom_boxplot(mapping = aes(factor(DGAT1_geno), PY, fill = factor(DGAT1_geno))) +
     labs(title = Title, x = Xlab, y = Ylab) +
     theme(plot.title = element_text(face = "bold", size = 25),
           axis.title.x = element_text(size = 20),
           axis.title.y = element_text(size = 20),
           axis.text.x = element_text(size = 15),
           axis.text.y = element_text(size = 15),
           legend.text = element_text(size = 15),
           legend.title = element_text(size = 15))


