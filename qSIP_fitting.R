# function for qSIP
# performs calculation of average density of OTUs by fitting a nonlinear function
# y = e^-((x-a)^2)/(2*b^2))
# a is the average density of the OTU
# b is the width of the density distribution
# x are densities of the individual fractions of the SIP density gradient
density.qSIP.example <- data.frame(sample1=c(1.7600,1.7527,1.7452,1.7345,1.7270,1.7196,1.7121,1.7035,1.6961,1.6886,1.6798,1.6718,1.6638,1.6558))
# y are relative abundances based on gene copy numbers of the respective OTU in the individual fractions
abundance.qSIP.example <- data.frame(OTU1=c(0.11734270,0.08877948,0.04472409,0.03846455,0.02364541,0.04509810,0.19691188,1.00000000,0.89248736,0.93851400,0.11722801,0.08686277,0.04280901,0.02749255), OTU2=c(0.04085279,0.04390676,0.02153932,0.03129609,0.02480411,0.06086885,0.18491933,0.75019344,1.00000000,0.28620960,0.06442150,0.06428784,0.04863304,0.03725738), OTU3=c(0.03742077,0.02782354,0.02054674,0.02373387,0.02439002,0.05797520,0.20807243,0.79773130,1.00000000,0.12569966,0.04107789,0.05147142,0.02921525,0.01841257))
# function
profile.qSIP.example <- function(sample,OTU) {
  x <- density.qSIP.example[,sample]
  y <- abundance.qSIP.example[,OTU]
  nonlin_mod=nls(y~exp(-((x-a)^2)/(2*b^2)),start=list(a=x[which(y==max(y))],b=0.005), trace=FALSE, control=nls.control(maxiter=50, tol=1e-05, minFactor=1/1024,printEval=TRUE,warnOnly=TRUE))
  nonlin_mod
  plot(x,y,main=paste("OTU ",OTU,"; sample ",sample,"; R²=",round(cor(y,predict(nonlin_mod))^2,4),"; error=",round(sqrt(mean((y-predict(nonlin_mod))^2)),4)))
  points(x,predict(nonlin_mod),col="#FF0000", pch=3)
  summary(nonlin_mod)
  return(list(OTU,sample,cor(y,predict(nonlin_mod))^2, sqrt(mean((y-predict(nonlin_mod))^2)), coef(nonlin_mod)[1],coef(nonlin_mod)[2]))
}

# run for one OTU in one sample
sapply("OTU1",FUN=profile.qSIP.example,sample="sample1")
# run for all OTUs in one sample
# when running large numbers of OTUs, it is recommended to remove or comment the plot() and points() lines in the function
sample1.qSIP.example <- sapply(1:3,FUN=profile.qSIP.example,sample="sample1")
# sample1.qSIP.example contains the following information: 1 = number of OTU, 2 = sample name, 3 = R² of the fit, 4 = error of the fit, 5 = a (density), 6 = b (width of density distribution)




