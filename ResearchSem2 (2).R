##Loading in nessecary packages
library("devtools")
library("rootSolve")
library("hypergeo")
library("parallel")
library("doParallel")
library("foreach")
library("pryr")
library("CopulaOne")
library("appell")
##Marginal CDF
OneVarF <- function(a,b,x) {
  hold <- tryCatch(1 - a/(a+b)*Re(hypergeo::hypergeo(1,b,a+b+1,1-x)), error = function(err) FALSE, warning = function(err) FALSE)
  if (!is.logical(hold) && is.finite(hold)) {
    return(hold)
  }
  return(NA)
}
##Marginal PDF
OneVarFPDF <- function(a,b,x) {
  hold <- tryCatch((a*b)/(a+b)/(a+b+1)*Re(hypergeo::hypergeo(2,b+1,a+b+2,1-x,tol = 1e-06, maxiter = 10000)), error = function(err) FALSE, warning = function(err) FALSE)
  if (!is.logical(hold) && is.finite(hold)) { 
    return(hold)
  }
  return(NA)
}
OneVarFPDFVector <- Vectorize(OneVarFPDF,"x")
##Joint CDF
TwoVarF <- function(a,b,x1,x2) {
  hold <- tryCatch(OneVarF(a,b,x1)+OneVarF(a,b,x2)-1+((a*(a+1)/((a+b)*(a+b+1)))*Re(appellf1(b,1,1,a+b+2,1-x1,1-x2)$val)), error = function(err) FALSE, warning = function(err) FALSE)
  if (!is.logical(hold) && is.finite(hold)) {
    return(hold)
  }
  return(NA)
}
##Joint PDF
TwoVarFPDF <- function(a,b,x1,x2) {
  tem1 <- a*b*(a+1)*(b+1)/(a+b)/(a+b+1)/(a+b+2)/(a+b+3)
  tem2 <- tryCatch(Re(appell::appellf1(b+2,2,2,a+b+4,1-x1,1-x2,userflag=1)$val),error=function(err) FALSE,warning=function(err) FALSE)
  if (!is.logical(tem2) && is.finite(tem2)) { 
    return(tem1*tem2)
  }
  return(NA)
}
TwoVarFPDFVector <- Vectorize(TwoVarFPDF,c("x1","x2"))
##Newtons Method / Midpoint Method Root Finder
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
      if(is.na(DEN)) {
        return(NA)
      }
      if(is.na(CDF)) {
        return(NA)
      }
      if (CDF < v) {
        lower <- t
      }
      else {
        upper <- t
      }
    }
    out <- exp(t)
  }
  return(out)
}
rootFinderVector <- Vectorize(rootFinder,"v")
##Cubic Spline Interpolation NOT FULLY IMPLEMENTED YET
rootFinderInterp <- function(a,b,v) {
  minval <- min(v)
  maxval <- max(v)
  testvals <- seq(minval,maxval,(maxval-minval)/199)
  basevals <- unlist(parLapply(c1,1:7,MLEPartial,a=a,b=b,v=testvals))
  return(spline(x = testvals,y = basevals,xout = v,method = "hyman")$y)
}
##Empirical CDF (uniroot appears to be better)
empiricalCDF <- function(v) {
  return(rank(v)/(length(v)+1))
}
##Copula CDF
Copula <- function(a,b,v1,v2) {
  x1 <- rootFinder(a,b,v1)
  x2 <- rootFinder(a,b,v2)
  return(TwoVarF(a,b,x1,x2))
}
##Copula PDF
CopulaPDF <- function(a,b,v1,v2) {
  q1 <- tryCatch(rootFinderInterp(a,b,v1), error = function(err) FALSE, warning = function(err) FALSE)
  q2 <- sort(q1)[rank(v2)]
##  q2 <- tryCatch(rootFinderInterp(a,b,v2), error = function(err) FALSE, warning = function(err) FALSE)
  if (is.logical(q1)||is.logical(q2))
  {
    return(NA)
  } 
  else {
    if (is.finite(q1) && is.finite(q2))
    {
      tem1 <- TwoVarFPDFVector(a,b,q1,q2)
      tem2 <- OneVarFPDFVector(a,b,q1)
      tem3 <- OneVarFPDFVector(a,b,q2)
      if (is.finite(tem1) && is.finite(tem2) && is.finite(tem3))
      {
        return(tem1/tem2/tem3)
      } else
      {
        return(NA)
      }
    } else
    {
      return(NA)
    }
  }
}
##Allows Copula PDF to take a vector as input
CopulaPDFVector <- Vectorize(CopulaPDF,c("v1","v2"))
##Portion of the MLE run on a specific core
MLEPartial <- function(a,b,v,i) {
  x <- floor(length(v)/7)
  if(i < 7) {
    ourdat <- v[(((i-1)*x+1)):(i*x)]
  }
  else {
    ourdat <- v[(((i-1)*x)+1):(length(v))]
  }
  hold <- rootFinderVector(a,b,ourdat)
  return(hold)
}
##Portion of the MLE that assigns portions of the computations to the cores
MLEFinal <- function(a,b,m) {
  ##hold <- unlist(parLapply(c1,1:7,MLEPartial,a=a,b=b,m=m))
  hold <- (CopulaPDF(a,b,m[,1],m[,2]))
  hold <- -log(hold)
  hold <- hold[!is.na(hold)]
  hold <- hold[is.finite(hold)]
  hold <- hold[!is.nan(hold)]
  return(sum(hold))
}
##Coputes the actual MLE using the optim function
solveMLE <- function(a,b,m) {
  hold <-  m
  corr <- cor(x=hold[,1],y=hold[,2],method="kendall")
  if(corr < 0) {
    hold[,2] <- 1-hold[,2]
  }
  matrix1 <-  hold
  temp <- function(c) {
    return(MLEFinal(c[1],c[2],matrix1))
  }
  no_cores <- detectCores() - 1
  registerDoParallel(no_cores)
  c1 <<- makeCluster(no_cores)
  clusterEvalQ(c1,library(pryr))
  clusterExport(c1,varlist=c("hypergeo","OneVarF","OneVarFPDF","OneVarFPDFVector","TwoVarF","TwoVarFPDF","TwoVarFPDFVector","Copula","CopulaPDF","CopulaPDFVector","rootFinder","rootFinderInterp","appellf1","jpGGEE","jdGGEE","dGGEE_COP","qGGEE","dGGEE","dat","rootFinderVector"))
  answer <- (optim(c(a,b),temp,method="L-BFGS-B",lower = 0.1, upper = 5))
  stopCluster(c1)
  return(answer)
}
##Coputes the actual MLE using the NLM function
solveMLENLM <- function(a,b,m) {
  hold <-  m
  corr <- cor(x=hold[,1],y=hold[,2],method="kendall")
  if(corr < 0) {
    hold[,2] <- 1-hold[,2]
  }
  matrix1 <-  hold
  temp <- function(c) {
    return(MLEFinal(exp(c[1])+0.01,exp(c[2])+0.01,matrix1))
  }
  no_cores <- detectCores() - 1
  registerDoParallel(no_cores)
  c1 <<- makeCluster(no_cores)
  clusterEvalQ(c1,library(pryr))
  clusterExport(c1,varlist=c("hypergeo","OneVarF","OneVarFPDF","OneVarFPDFVector","TwoVarF","TwoVarFPDF","TwoVarFPDFVector","Copula","CopulaPDF","CopulaPDFVector","rootFinder","rootFinderInterp","appellf1","jpGGEE","jdGGEE","dGGEE_COP","qGGEE","dGGEE","dat","rootFinderVector"))
  answer <- (nlm(temp,c(log(a-0.01),log(b-0.01))))
  stopCluster(c1)
  answer$estimate <- exp(answer$estimate)+0.01
  return(answer)
}

################################
## Work in Progress / Testing ##
################################

testData <- read.table("CopulaOneData.txt")

data("euro0306")
dat <- uscore(euro0306[,c(2,3)])[1:100,]
dat <- matrix(c(uscore(testData$V2),uscore(testData$V3),uscore(testData$V4),uscore(testData$V5)),ncol = 2)
par0 <- c(0.5, 0.5)

system.time(cat(solveMLENLM(0.5,0.5,dat)$estimate))

system.time(cat(fitCopulaOne(par0, dat=dat, copula_family="GGEE",lower = 0.1, upper = 5,opt="nlm")$estimate))

fitCopulaOne(par0, dat=dat, copula_family="GGEE",lower = 0.1, upper = 5,opt="optim")

solveMLE(dat)


rootFinderInterp(0.5,0.5,dat[,1])
MLEFinal(0.5,0.5,dat)


stopCluster(c1)

rm(list = ls())

##Interpelation

##Simulation

a <- 0.3
b <- 0.4

set.seed(1000)
r1 <- rgamma(1000, shape = a, scale = 1)

##set.seed(3000)
r2 <- rgamma(1000, shape = b, scale = 1)

##set.seed(5000)
e11 <- rexp(1000,1)

##set.seed(7000)
e12 <- rexp(1000,1)

##set.seed(9000)
e21 <- rexp(1000,1)

##set.seed(11000)
e22 <- rexp(1000,1)

x1 <- (r1/r2)*(e11/e12)
x2 <- (r1/r2)*(e21/e22)

dat <- matrix(c(uscore(x1),uscore(x2)),ncol = 2)
system.time(solveMLENLM(dat,0.3,0.3))
range(x1,x2)

c(0.1,0.4,0.6)[rank(c(0.4,0.1,0.6))]

##Try different starting functions to see if we get better estimation on optim and nlm
##Find out which is faster, optim or nlm
##Clean the code, add comments, push to github
##Try different models - DONE
##Cubic Spline Interpolation - DONE


##Implement Cubic Spline into main file

##Compute all data points but only for q1
##Try cubic spline with transformation or make it smoother
##Try different numbers of interpolation points (difference in alpha and beta)
##AIC