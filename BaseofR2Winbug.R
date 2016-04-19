mainDir <- "E:/R code_wd/Poisson-Norm/Lkeq/PoisMd"
mainDir2 <-"E:/R output/Poisson-Norm/Lkeq/PoisMd"
BugDir<-"E:/winbugs14"
dataDir <-"E:/Simulated dataset/Poisson-Norm/Lkeq"

setwd(mainDir)
install.packages("R2WinBUGS")
library(R2WinBUGS)

#first dataset before simulation#


subDir <- "1"
mixed=read.csv(file=paste(file.path(dataDir, subDir),".csv",sep=""), header=T)
#set working directory#


dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))

#set output director#
dir.create(file.path(mainDir2, subDir), showWarnings = TRUE)
outputdirc=file.path(mainDir2, subDir)

#WINBUG CODE#
poisunit<-function() {
  #link function#
  for (i in 1:16){
      eta[i]<-trteff[Trt[i]]+blkeff[Block[i]]+unit[i]
      mu[i]<-exp(eta[i])
      y[i]~dpois(mu[i])
  }
  #random effect#
  for (i in 1:16){
    unit[i]~dnorm(0,unitprecisn)
  }
  unitsigma<-1/unitprecisn
  for (j in 1:8){
    blkeff[j]~dnorm(0,blkprecisn)
    }
  blksigma<-1/blkprecisn
  #priorS#
  for (k in 1:2){
    trteff[k]~dnorm(0,precisn)
  }
  unitprecisn~dgamma(aunit,bunit)
  blkprecisn~dgamma(ablk,bblk)
  
  #######Compute fit Statistics#######
  for (i in 1:16){
    BICobs[i]<-2*exp(eta[i])+2*loggam(y[i]+1)-2*y[i]*eta[i]
  }
  BIC<-sum(BICobs[1:16])+4*log(16)  
  
  }


filename <- file.path(outputdirc,"poisutmd.txt")
write.model(poisunit, con = filename, digits = 5)  


attach(mixed)
y = as.vector(count_ij)
Trt = as.vector(treatment)
Block = as.vector(Block)
aunit=0.001;bunit=0.001;ablk=0.001;bblk=0.001;precisn=0.001;

data=list("y","Trt","Block","aunit","bunit","ablk","bblk","precisn")
parameters = c("trteff", "blksigma", "unitsigma","BIC","DIC")
inits = list(list( trteff = c(2,    1) , blkprecisn=1, unitprecisn=0.5) , 
             list( trteff = c(1.8, 1), blkprecisn=0.5,unitprecisn=0.25))
fitpois <- bugs(data, inits, parameters, model.file = filename, 
                      n.chains=2, n.iter = 50000, n.burnin = 10000, n.thin=1, DIC=TRUE, codaPkg = TRUE , debug=FALSE,
                      program="WinBUGS", bugs.directory=BugDir,
                      working.directory=outputdirc)



mcmc.out = file.path(outputdirc,"coda1.txt")
index = file.path(outputdirc,"codaIndex.txt")
mcmc = read.coda(output.file=mcmc.out, index.file=index)
head(mcmc)

effectiveSize(mcmc)[,c(1,5,6,8,9)])
raftery.diag(mcmc, q=0.025, r=0.005, s=0.95, converge.eps=0.001)  # One chain only
geweke.diag(mcmc, frac1=0.1, frac2=0.5)	# One chain only

str(summary(mcmc))
summary(mcmc)
HPDinterval(mcmc, prob=0.95)
plot(mcmc[,c(2,3,4)],trace = TRUE, density = TRUE, smooth=FALSE, auto.layout = TRUE)



