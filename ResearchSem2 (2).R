install.packages("gsl")
library("gsl")
install.packages("f21hyper")
library("f21hyper")

install.packages("appell")

library("devtools")

install_github("larryleihua/CopulaOne")

##genhypergeo

library("devtools")
library("rootSolve")
library("hypergeo")
library("parallel")
library("doParallel")
library("foreach")
library("pryr")
library("CopulaOne")
library("appell")

hypergeo_gosper

a <- seq(0.1,1,0.1)
b <- seq(0.1,1,0.1)
a <- 0.5
b <- 0.4
x <- seq(0.1, 2000, 0.1)
y <- seq(0.1, 2000, 0.1)
OneVarF <- function(a,b,x) {
  hold <- tryCatch(1 - a/(a+b)*Re(hypergeo::hypergeo(1,b,a+b+1,1-x)), error = function(err) FALSE, warning = function(err) FALSE)
  if (!is.logical(hold) && is.finite(hold)) {
    return(hold)
  }
  return(NA)
}
OneVarFPDF <- function(a,b,x) {
  hold <- tryCatch((a*b)/(a+b)/(a+b+1)*Re(hypergeo::hypergeo(2,b+1,a+b+2,1-x,tol = 1e-06, maxiter = 10000)), error = function(err) FALSE, warning = function(err) FALSE)
  if (!is.logical(hold) && is.finite(hold)) { 
    return(hold)
  }
  return(NA)
}
TwoVarF <- function(a,b,x1,x2) {
  hold <- tryCatch(OneVarF(a,b,x1)+OneVarF(a,b,x2)-1+((a*(a+1)/((a+b)*(a+b+1)))*Re(appellf1(b,1,1,a+b+2,1-x1,1-x2)$val)), error = function(err) FALSE, warning = function(err) FALSE)
  if (!is.logical(hold) && is.finite(hold)) {
    return(hold)
  }
  return(NA)
}
TwoVarFPDF <- function(a,b,x1,x2) {
  hold <- tryCatch((((a*b*(a+1)*(b+1))/((a+b)*(a+b+1)*(a+b+2)*(a+b+3)))*Re(appellf1(b+2,2,2,a+b+4,1-x1,1-x2)$val)), error = function(err) FALSE, warning = function(err) FALSE)
  if (!is.logical(hold) && is.finite(hold)) 
    return(hold)
  return(NA)
}

no_cores <- detectCores()
registerDoParallel(no_cores)
c1 <- makeCluster(no_cores)
clusterEvalQ(c1,library(pryr))
clusterExport(c1,varlist=c("hypergeo","OneVarF","OneVarFPDF","TwoVarF","TwoVarFPDF","Copula","CopulaPDF","rootFinder","appellf1","a","b","x","y","jpGGEE","jdGGEE","dGGEE_COP"))

stopCluster(c1)

##uniroot,code in copula one,newton raphson method in r
##Find Copula?

rootFinder <- function(a,b,v)
{
  if (v == 0) {
    out <- 0
  }
  else {
    CDF <- -0.01
    DEN <- 1
    i <- 0
    t <- 0
    lower <- -1e20
    upper <- 1e20
    
    while ((i < 1000) && ((abs(v - CDF) > 1e-6))) {
      i <- i + 1
      t <- t - (CDF - v)/DEN
      if (t < lower || t > upper) {
        t <- 0.5 * (lower + upper)
      }
      DEN <- OneVarFPDF(a,b,exp(t)) * exp(t)
      CDF <- OneVarF(a,b,exp(t))
      if (CDF < v) {
        lower <- t
      }
      else {
        upper <- t
      }
      if(is.na(DEN)) {
        return(NA)
      }
    }
    out <- exp(t)
  }
  return(out)
}

##empiricalCDF - DONE
##Copula - DONE
##CopulaPDF - DONE
##MLE
##Cores And Multithreading

empiricalCDF <- function(v) {
  return(rank(v)/(length(v)+1))
}

Copula <- function(a,b,v1,v2) {
  x1 <- rootFinder(a,b,v1)
  x2 <- rootFinder(a,b,v2)
  return(TwoVarF(a,b,x1,x2))
}

CopulaPDF <- function(a,b,v1,v2) {
  x1 <- rootFinder(a,b,v1)
  x2 <- rootFinder(a,b,v2)
  c1 <- OneVarFPDF(a,b,x1)
  c2 <- OneVarFPDF(a,b,x2)
  return(TwoVarFPDF(a,b,x1,x2)/c1/c2)
}

MLE <- function(a,b,m) {
  hold <- parRapply(c1,m,function(x) {-log(CopulaPDF(a,b,x[1],x[2]))})
  hold <- hold[!is.na(hold)]
  hold <- hold[is.finite(hold)]
  hold <- hold[!is.nan(hold)]
  return(sum(hold))
}

solveMLE <- function(m) {
  hold <-  m
  corr <- cor(x=hold[,1],y=hold[,2],method="kendall")
  if(corr < 0) {
    hold[,2] <- 1-hold[,2]
  }
  matrix1 <<-  hold
  temp <- function(a) {
    return(MLE(a[1],a[2],matrix1))
  }
  return(optim(c(0.5,0.5),temp,method="L-BFGS-B", lower = 0.1, upper = 1-1e-6))
}

##CopulaOne Empirical cdf - DONE
##Time TwoVarPDF and TwoVarF against CopulaOne - DONE
##Remove loop when listing data in MLE -LEAVE
##Test against CopulaOne data
##Test against different methods in optim

a <- seq(0.1,0.9,0.1)
b <- seq(0.1,0.9,0.1)
x <- seq(0.1,200,0.1)
y <- 1

system.time(parSapply(c1,x,function(x) {mapply(function(a,b) TwoVarF(a,b,x,1),a,b)}))
system.time(parSapply(c1,x,function(x) {mapply(function(a,b) jpGGEE(a,b,x,1),a,b)}))

system.time(parSapply(c1,x,function(x) {mapply(function(a,b) TwoVarFPDF(a,b,x,1),a,b)}))
system.time(parSapply(c1,x,function(x) {mapply(function(a,b) jdGGEE(a,b,x,1),a,b)}))

system.time(parSapply(c1,x,function(x) {mapply(function(a,b) CopulaPDF(a,b,x,1),a,b)}))
system.time(parSapply(c1,x,function(x) {mapply(function(a,b) dGGEE_COP(a,b,x,1),a,b)}))

testData <- read.table("CopulaOneData.txt")

system.time(empiricalCDF(testData$V3))

system.time(uscore(testData$V3))

MLE(0.2,0.8,testData$V3,testData$V4)

dGGEE_COP(0.2,0.2,0.8,0.8)
CopulaPDF(0.8,0.8,0.2,0.2)

u1 <- empiricalCDF(testData$V3)
u2 <- empiricalCDF(testData$V4)

hold <- rbind(u1,u2)
hold

testData$V3

empiricalCDF(testData$V4)

par0 <- c(0.3,0.7)
patternpar = seq(1,2)

uscore(testData$v3)

fitCopulaOne(c(0.5,0.5),dat = hold,copula_family = "GGEE")
fitCopulaOne(par0 = par0, whichpar=seq(1,length(par0)), patternpar=seq(1,2), dat = hold, flag=1, integration=F, opt="L-BFGS-B", se=F, lower=rep(0.1, length(unique(patternpar))), upper=rep(5, length(unique(patternpar))), trace=0, factr=1e9, printlevel=0, copula_family="GGEE")

solveMLE(testData$V3,testData$V4)

##Find why our TwoVarPDF and CDF are faster
##Get Copula One working in seperate R file
##Compare to solve MLE
data("euro0306")
dat <- matrix(c(testData$V2,testData$V3),ncol=2)
dat <- uscore(dat)

fitCopulaOne(c(0.5,0.5),dat = dat,copula_family = "GGEE",lower = 0.1,upper = 1-1e-6)
solveMLE(testData$V2,testData$V3)

data("euro0306")
dat <- uscore(euro0306[,c(2,3)])[1:100,]
dat <- matrix(c(uscore(testData$V2),uscore(testData$V3)),ncol = 2)
par0 <- c(0.5, 0.5)

system.time(cat(solveMLE(dat)$par))
system.time(cat(fitCopulaOne(par0, dat=dat, copula_family="GGEE",lower = 0.1, upper = 1-1e-6)$par))

testData$V5

MLE(0.511,0.65,dat[,1],dat[,2])

dGGEE_COP(0.2,0.4,1.2,0.2)
CopulaPDF(1.2,0.2,0.2,0.4)


## Implent COpula Removal of NA values
## Get parallel computing working better
## Compare ranking methods for NA values

##Computational Speed
##Overleaf

rm(list = ls())

##Interpelation