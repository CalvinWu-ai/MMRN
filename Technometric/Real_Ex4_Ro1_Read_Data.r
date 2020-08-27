### Real Data Analysis - Example 4        ###
### Runze Li, Wei Zhong & Liping Zhu      ###
### March 2011                            ###

# The Cardiomyopathy microarray dataset was once analyzed by Segal, Dahlquist
# and Conklin (2003) and Hall and Miller (2009). The goal is to identify the most influential
# genes for overexpression of a G protein-coupled receptor (Ro1) in mice. The response Y is
# the Ro1 expression level, and the predictors Xk¡¯s are other gene expression levels. Compared
# with the sample size n = 30 in this dataset, the dimension p = 6319 is very large.

############## Read in data  ##############
X             <- t(read.table("D:\\My Research-Wei Zhong\\Distance-Correlation-Screening\\Codes\\DC-SIS-LZZ-Supplements\\Data\\Ro131.csv", header=TRUE, sep=","))
Ro131_names	  <- read.table("D:\\My Research-Wei Zhong\\Distance-Correlation-Screening\\Codes\\DC-SIS-LZZ-Supplements\\Data\\Genenum.csv", header=FALSE, sep=",") # List of gene names as used in paper
Ro131_names	  <- Ro131_names[,1][c(1:6077,6079:6320)] # Need to remove name of

# Standardise
colmeans    <- apply(X,2,mean)
colsd       <- apply(X,2,sd)
X           <- t((t(X) - colmeans)/colsd)
# Remove predictor gene, which is tightly related to Ro131 as done in Segal et al
X           <- X[,c(1:6077,6079:6320)]

# Read in the response variable, entered manually
Y           <- c(143,84,98,83,153,141,191,130,744,381,1047,806,621,849,475,966,708,487,1447,1693,1731,1025,376,126,102,149,91,153,235,68)

n           <- dim(X)[1]
p           <- dim(X)[2]
d           <- floor(n/log(n))

