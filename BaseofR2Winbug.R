#set winbugs working environment
mainDir <- "E:/R code_wd/Poisson-Norm/Lkeq/PoisMd"
mainDir2 <-"E:/R output/Poisson-Norm/Lkeq/PoisMd"
BugDir<-"E:/winbugs14"
dataDir <-"E:/Project Dataset/Poisson-Norm/Lkeq"

setwd(mainDir)
#install.packages("R2WinBUGS")
library(R2WinBUGS)

#define simulation result array
elpdls<-1:300
ePostls<-1:300
pDICls<-1:300
DICls<-1:300
BICdimxls<-1:300
BICparals<-1:300


#start of simulation run;the number of simulation iteration was defined by simudt

for (simudt in 1:300){
#first dataset before simulation#
#simudt=1
subDir =toString(simudt)
mixed=read.csv(file=paste(file.path(dataDir, subDir),".csv",sep=""), header=T)

#set working directory#
#dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
#setwd(file.path(mainDir, subDir))

#set output director#
dir.create(file.path(mainDir2, "tmpoutpt"), showWarnings = FALSE)
outputdirc=file.path(mainDir2, "tmpoutpt")


#WINBUG CODE#
poisunit<-function() {
  pi<-3.1416
  #link function#
  for (i in 1:32){
      eta[i]<-trteff[Trt[i]]+blkeff[Block[i]]+unit[i]
      mu[i]<-exp(eta[i])
      y[i]~dpois(mu[i])
  }
  #random effect#
  for (i in 1:32){
    unit[i]~dnorm(0,unitprecisn)
  }
  unitsigma<-1/unitprecisn
  for (j in 1:8){
    blkeff[j]~dnorm(0,blkprecisn)
    }
  blksigma<-1/blkprecisn
  #priorS#
  for (k in 1:4){
    trteff[k]~dnorm(0,precisn)
  }
  unitprecisn~dgamma(aunit,bunit)
  blkprecisn~dgamma(ablk,bblk)
   
  #######Compute fit Statistics#######
  for (i in 1:32){
    #calculate teh loglikelihood of y|theta 
    epost1obs[i]<-y[i]*eta[i]-mu[i]-loggam(y[i]+1)#+0.5*log(unitprecisn)
    #-0.5*(unit[i]*unit[i]*unitprecisn+log(2*pi))
  }
  
  ePOST<-sum(epost1obs[1:32])
  }

##set R2Winbug environment
filename <- file.path(outputdirc,"poisutmd.txt")
write.model(poisunit, con = filename, digits = 5)  

#read dataset into winbugs
attach(mixed)
y = as.vector(count_ij)
Trt = as.vector(treatment)
Block = as.vector(Block)

#set prior parameters
aunit=0.001;bunit=0.001;ablk=0.001;bblk=0.001;precisn=0.001;
data=list("y","Trt","Block","aunit","bunit","ablk","bblk","precisn")

#set parameter to monitor
parameters = c("ePOST","mu")
#parameters = c("blkprecisn", "unitprecisn","ePOST","mu","blkeff","trteff","eta")

#initials
inits = list(list( trteff = c(2,2,1,1) , blkprecisn=1, unitprecisn=0.5) , 
             list( trteff = c(1.8, 2,1,1), blkprecisn=0.5,unitprecisn=0.25))


#run mcmc in winbug#
fitpois <- bugs(data, inits, parameters, model.file = filename, 
                     n.chains=2, n.iter = 35000, n.burnin = 15000, n.thin=1, DIC=TRUE, codaPkg = TRUE , 
                debug=FALSE,program="WinBUGS", bugs.directory=BugDir,working.directory=outputdirc)


##read MCMC result
mcmc.out = file.path(outputdirc,"coda1.txt")
index = file.path(outputdirc,"codaIndex.txt")
mcmc = read.coda(output.file=mcmc.out, index.file=index)
head(mcmc)

## Find out the Posterior Mean of mu
postMeans=apply(mcmc,2,mean)
mu=postMeans[grep("mu",names(postMeans))]

#calculated the log-likelihood of y give posterior mean of mu
elpd1ob<-1:32
pi=3.1416
for (i in 1:32){
  elpd1ob[i]<-y[i]*log(mu[i])-mu[i]-lfactorial(y[i])#+0.5*log(unitprecisn)
  #-0.5*(unit[i]*unit[i]*unitprecisn+log(2*pi))
}
elpd=sum(elpd1ob)
elpdls[simudt]=elpd
#calculated the mean of likelihood on each iteration
ePost=postMeans["ePOST"] 
ePostls[simudt]=ePost
#Fit Statistics
pDIC=2*(elpd-ePost)
pDICls[simudt]=pDIC

DICls[simudt]=-2*(elpd-pDIC)
BICdimxls[simudt]=-2*(elpd)+32*log(32)
BICparals[simudt]=-2*(elpd)+pDIC*log(32)

}

#output result
result=cbind(elpdls,ePostls,pDICls,DICls,BICdimxls,BICparals)
write.csv(result,"result.csv")

