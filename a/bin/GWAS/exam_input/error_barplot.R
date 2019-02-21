# /usr/bin/rscript

## Description:
## This script plots an error barplot of mean and sd of phenotype.
## 
## Author:
## Min Wang (min.wang@ecodev.vic.gov.au)
##
## Date Created:
## 31 March 2015
## 
## Date modified and reason: 
## 
## Execution: 
## Rscript <module_name> <date> <psf> <gender> <pheno> <infile> 


## Generate error barplot for mean and standard deviation ##
dt_FY <- as.data.table(ddply(data, .(DGAT1_geno), plyr::summarize, mean = mean(FY), sd = sd(FY)))
dt_MY <- as.data.table(ddply(data, .(DGAT1_geno), plyr::summarize, mean = mean(MY), sd = sd(MY)))
dt_PY <- as.data.table(ddply(data, .(DGAT1_geno), plyr::summarize, mean = mean(PY), sd = sd(PY)))

aim <- 'Barplot: DGAT1 genotype effects on'
SubTitle <- 'Error bar: (mean + sd, mean - sd)'
Xlab  <- 'DGAT1 genotypes'
limits <- aes(ymax = mean + sd, ymin = mean - sd)
dodge <- position_dodge(width = 0.9)

Title <- paste(aim, 'mean(FY)')
Ylab  <- 'mean(Fat Yield)'
g11 <- ggplot(dt_FY, aes(fill = factor(DGAT1_geno), y = mean, x = DGAT1_geno)) + 
       geom_bar(position = dodge, stat = "identity") + 
       geom_errorbar(limits, position = dodge, width = 0.25) + 
       ggtitle(bquote(atop(.(Title), atop(.(SubTitle))))) +
       labs(x = Xlab, y = Ylab) +
       theme(plot.title = element_text(face = "bold", size = 25),
             axis.title.x = element_text(size = 20),
             axis.title.y = element_text(size = 20),
             axis.text.x = element_text(size = 15),
             axis.text.y = element_text(size = 15),
             legend.text = element_text(size = 15),
             legend.title = element_text(size = 15)) 

Title <- paste(aim, 'mean(MY)')
Ylab  <- 'mean(Milk Yield)'
g21 <- ggplot(dt_MY, aes(fill = factor(DGAT1_geno), y = mean, x = DGAT1_geno)) +
       geom_bar(position = dodge, stat = "identity") +
       geom_errorbar(limits, position = dodge, width = 0.25) +
       ggtitle(bquote(atop(.(Title), atop(.(SubTitle))))) +
       labs(x = Xlab, y = Ylab) +
       theme(plot.title = element_text(face = "bold", size = 25),
             axis.title.x = element_text(size = 20),
             axis.title.y = element_text(size = 20),
             axis.text.x = element_text(size = 15),
             axis.text.y = element_text(size = 15),
             legend.text = element_text(size = 15),
             legend.title = element_text(size = 15))

Title <- paste(aim, 'mean(PY)')
Ylab  <- 'mean(Protein Yield)'
g31 <- ggplot(dt_PY, aes(fill = factor(DGAT1_geno), y = mean, x = DGAT1_geno)) +
       geom_bar(position = dodge, stat = "identity") +
       geom_errorbar(limits, position = dodge, width = 0.25) +
       ggtitle(bquote(atop(.(Title), atop(.(SubTitle))))) +
       labs(x = Xlab, y = Ylab) +
       theme(plot.title = element_text(face = "bold", size = 25),
             axis.title.x = element_text(size = 20),
             axis.title.y = element_text(size = 20),
             axis.text.x = element_text(size = 15),
             axis.text.y = element_text(size = 15),
             legend.text = element_text(size = 15),
             legend.title = element_text(size = 15))



