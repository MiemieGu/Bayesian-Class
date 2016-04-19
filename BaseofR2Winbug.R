setwd("I:/R code_wd/Poisson-Norm")
install.packages("R2WinBUGS")
library(R2WinBUGS)

#first dataset before simulation#
mixed=read.csv(file="I:/Simulated R dataset/Poisson-Norm/lkeq_1.csv", header=T)
#set working directory#

mainDir <- "I:/R code_wd/Poisson-Norm"
subDir <- "lkeq01"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
setwd(file.path(mainDir, subDir))

#set output director#
mainDir2 <- "I:/R output/Poisson-Norm"
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
  
  #Compute fit Statistics#  
  
  
 
  }


filename <- file.path(outputdirc,"poisutmd.txt")
write.model(poisunit, con = filename, digits = 5)  


attach(mixed)
y = as.vector(count_ij)
Trt = as.vector(treatment)
Block = as.vector(Block)
aunit=0.001;bunit=0.001;ablk=0.001;bblk=0.001;precisn=0.001;

data=list("y","Trt","Block","aunit","bunit","ablk","bblk","precisn")
parameters = c("trteff", "blksigma", "unitsigma")
inits = list(list( trteff = c(2,    1) , blkprecisn=1, unitprecisn=0.5) , 
             list( trteff = c(1.8, 1), blkprecisn=0.5,unitprecisn=0.25))
fitpoi1 <- bugs(data, inits, parameters, model.file = filename, 
                      n.chains=2, n.iter = 50000, n.burnin = 10000, n.thin=1, DIC=TRUE, codaPkg = TRUE , debug=TRUE,
                      program="WinBUGS", bugs.directory="I:/winbugs14",
                      working.directory=outputdirc)





