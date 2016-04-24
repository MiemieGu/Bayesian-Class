mainDir <- "E:/R code_wd/Poisson-Norm/Lkeq/NBMd/"
dir.create(file.path(mainDir), showWarnings = FALSE)
mainDir2 <-"E:/R output/Poisson-Norm/Lkeq/NBMd"
dir.create(file.path(mainDir2), showWarnings = FALSE)
BugDir<-"E:/winbugs14"
dataDir <-"E:/Project Dataset/Poisson-Norm/Lkeq"

setwd(mainDir)
install.packages("R2WinBUGS")
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
  
  
  #set output director#
  dir.create(file.path(mainDir2, "tmpoutpt"), showWarnings = FALSE)
  outputdirc=file.path(mainDir2, "tmpoutpt")
dir.create(file.path(outputdirc), showWarnings = FALSE)
  
  #WINBUGs CODE#
  negbin<-function() {
    #model 
      for (i in 1:32){
      eta[i]<-trteff[Trt[i]]+blkeff[Block[i]]
      lambda[i]<-exp(eta[i])
      p[i]<-(r/(r+lambda[i]))
      y[i]~dnegbin(p[i],r)
    }
    
    #random effect#
    for (j in 1:8){
      blkeff[j]~dnorm(0,blkprecisn)
    }
      
    
    #priors#
      #blocks#
      for (k in 1:4){
      trteff[k]~dnorm(0,precisn)
    }
      #phi#
    phi~dgamma(aphi,bphi)
      r<-1/phi
      #block#
    blkprecisn~dgamma(ablk,bblk)
    
    
    #######Compute fit Statistics#######
    for (i in 1:32){
      #calculate teh loglikelihood of y|theta 
      epost1obs[i]<-y[i]*log(1-p[i])+r*log(p[i])+loggam(y[i]+r)-loggam(y[i]+1)-loggam(r)
    }    
    ePOST<-sum(epost1obs[1:32])
  }
  #end of winbugs code#
  
  ##set R2Winbug environment
  filename <- file.path(outputdirc,"negbinmd.txt")
  write.model(negbin, con = filename, digits = 5)  
  
  #read dataset into winbugs
  attach(mixed)
  y=as.vector(count_ij)
  Trt = as.vector(treatment)
  Block = as.vector(Block)
  
  #set prior parameters
  aphi=10;bphi=10;ablk=0.01;bblk=0.01;precisn=0.001;
  data=list("y","Trt","Block","aphi","bphi","ablk","bblk","precisn")
  
  #set parameter to monitor
  parameters = c("ePOST","p","r")
  #parameters = c("blkprecisn", "unitprecisn","ePOST","mu","blkeff","trteff","eta")
  
  #initials
  inits = list(list( trteff = c(2,2,1,1) , blkprecisn=1, phi=0.5) , 
               list( trteff = c(1.8, 2,1,1), blkprecisn=0.5,phi=0.25))
  
  
  #run mcmc in winbug#
  fitnegbin<- bugs(data, inits, parameters, model.file = filename, 
                   n.chains=2, n.iter = 35000, n.burnin = 15000, n.thin=1, DIC=TRUE, codaPkg = TRUE , 
                   debug=TRUE,program="WinBUGS", bugs.directory=BugDir,working.directory=outputdirc)
  
  
  ##read MCMC result
  mcmc.out = file.path(outputdirc,"coda1.txt")
  index = file.path(outputdirc,"codaIndex.txt")
  mcmc = read.coda(output.file=mcmc.out, index.file=index)
  head(mcmc)
  
  ## Find out the Posterior Mean of mu
  postMeans=apply(mcmc,2,mean)
  #postSd=apply(mcmc,2,sd)
  p=postMeans[grep("p",names(postMeans))]
  ePost=postMeans["ePOST"] 
  r=postMeans["r"]
  
  #calculated the log-likelihood of y give posterior mean of mu
  elpd1ob<-1:32
  for (i in 1:32){
    elpd1ob[i]<-y[i]*log(1-p[i])+r*log(p[i])+lfactorial(y[i]+r-1)-lfactorial(y[i])-lfactorial(r-1)
  }
  elpd=sum(elpd1ob)
  elpdls[simudt]=elpd
  #calculated the mean of likelihood on each iteration
  
  ePostls[simudt]=ePost
  #Fit Statistics
  pDIC=2*(elpd-ePost)
  pDICls[simudt]=pDIC
  
  DICls[simudt]=-2*(elpd-pDIC)
  BICdimxls[simudt]=-2*(elpd)+13*log(32)
  BICparals[simudt]=-2*(elpd)+pDIC*log(32)
}

result=cbind(elpdls,ePostls,pDICls,DICls,BICdimxls,BICparals)
write.csv(result,"result.csv")

