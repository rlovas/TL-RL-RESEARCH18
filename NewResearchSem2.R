install.packages("gsl")
library("gsl")
install.packages("f21hyper")
library("f21hyper")

install.packages("devtools")

library("devtools")

##genhypergeo

install_github("larryleihua/CopulaOne")

library("rootSolve")
library("hypergeo")
library("parallel")
library("doParallel")
library("foreach")
library("pryr")
library("CopulaOne")

hypergeo_gosper

a <- seq(0.1,1,0.1)
b <- seq(0.1,1,0.1)
a <- 0.5
b <- 0.4
x <- seq(0.1, 2000, 0.1)
y <- seq(0.1, 2000, 0.1)
OneVarF <- function(a,b,x) {
  out <- tryCatch(1 - a/(a+b)*Re(hypergeo::hypergeo(1,b,a+b+1,1-x,tol = 1e-06, maxiter = 10000)), error = function(err) FALSE, warning = function(err) FALSE)
  if (!is.logical(out) && is.finite(out)) {
    return(out)
  }
  return(NA)
}
OneVarFPDF <- function(a,b,x) {
  out <- tryCatch((a*b)/(a+b)/(a+b+1)*Re(hypergeo::hypergeo(2,b+1,a+b+2,1-x,tol = 1e-06, maxiter = 10000)), error = function(err) FALSE, warning = function(err) FALSE)
  if (!is.logical(out) && is.finite(out)) 
    return(out)
  return(NA)
}
OneVarFGosper <- function(a,b,x) {
  return(1-((a/(a+b))*Re(hypergeo_gosper(1,b,a+b+1,1-x))))
}
OneVarFBuhring <- function(a,b,x) {
  return(1-((a/(a+b))*Re(hypergeo_buhring(1,b,a+b+1,1-x))))
}
OneVarFShanks <- function(a,b,x) {
  return(1-((a/(a+b))*Re(hypergeo_shanks(1,b,a+b+1,1-x))))
}
TwoVarF <- function(a,b,x1,x2) {
  return(OneVarF(a,b,x1)+OneVarF(a,b,x2)-1+((a*(a+1)/((a+b)*(a+b+1)))*Re(hypergeo(b,1,1,a+b+2,1-x1,1-x2))))
}
TwoVarPDF <- function(a,b,x1,x2) {
  return(1-(((a*(a+1)*b)/((a+b)*(a+b+2)*(a+b+3)*(a+b+1)))*Re(hypergeo(b+2,2,2,a+b+4,1-x1,1-x2))))
}


no_cores <- detectCores() - 1
registerDoParallel(no_cores)
c1 <- makeCluster(no_cores)
clusterEvalQ(c1,library(pryr))
clusterExport(c1,varlist=c("a","b","x","y","hypergeo","hypergeo_gosper","hypergeo_buhring","hypergeo_shanks","OneVarF","OneVarF","OneVarFGosper","OneVarFBuhring","OneVarFShanks","TwoVarF"))

stopCluster(c1)

TimeMethod1Hyper1 <- function() {
  return(parSapply(c1,x,function(x) {mapply(function(a,b) OneVarF(a,b,x),a,b)}))
}
TimeMethod2Hyper1 <- function() {
  return(parSapply(c1,x,function(x) {mapply(function(a,b) OneVarFGosper(a,b,x),a,b)}))
}
TimeMethod3Hyper1 <- function() {
  return(parSapply(c1,x,function(x) {mapply(function(a,b) OneVarFBuhring(a,b,x),a,b)}))
}
TimeMethod4Hyper1 <- function() {
  return(parSapply(c1,x,function(x) {mapply(function(a,b) OneVarFShanks(a,b,x),a,b)}))
}
tempVar1 <- TimeMethod1Hyper1()
tempVar2 <- TimeMethod2Hyper1()
tempVar3 <- TimeMethod3Hyper1()
tempVar4 <- TimeMethod4Hyper1()

TimeMethod1PDF <- function() {
  return(parSapply(c1,x,function(x) {mapply(function(a,b) OneVarFShanks(a,b,x),a,b)}))
}

tempVar1[is.na(tempVar1)]

plot(1:20000,tempVar2)
plot(1:19983,tempVar3noInf)
tempVar3noInf <- tempVar3[abs(tempVar3) != Inf]

tempVar3noInf <- tempVar3noInf[!is.na(tempVar3noInf)]

tempVar3noInf <- tempVar3noInf[tempVar3noInf < 1000]

tempVar2noInf[tempVar2noInf == Inf]

tempi <- 0

for(i in 1:10000) {
  if(tempVar3noInf[i] < tempVar3noInf[i+1]) {
    tempi = i
    break
  }
}

cleanVector <- function(v) {
  v <- parSapply(c1,v,function(x) {if(abs(x) == Inf) {0} else {x}})
  v <- parSapply(c1,v,function(x) {if(x < 1) {0} else {x}})
  return(v)
}

tempVar4cleaned <- cleanVector(tempVar4)

plot(1:20000,tempVar1)
plot(1:20000,tempVar2)
plot(20:20000,tempVar3[20:20000])
plot(1:20,tempVar4[1:20])


TimeMethod1Hyper2 <- function() {
  system.time(parSapply(c1,x,function(x) {mapply(function(a,b,y) TwoVarF(a,b,x,y),a,b,y)}))
}

newtonsMethod <- function(a,b,v,x) {
  x0 <- x
  x1 <- x0-((OneVarF(a,b,exp(x0))-v)/((OneVarFPDF(a,b,exp(x0))-v))*exp(x0))
  i <- 0
  while(abs(x1-x0) > 1e-6 && i < 1000) {
    x0 <- x1
    x1 <- x0-((OneVarF(a,b,exp(x0))-v)/((OneVarFPDF(a,b,exp(x0))-v)*exp(x0)))
    i <- i+1
  }
  return(exp(x1))
}

alternateRootFinder <- function(a,b,v) {
  iterations <- 0
  lower <- -1
  upper <- 1
  while ((OneVarF(a,b,exp(upper))-v) < 0) {
    lower <- upper
    upper <- upper*2
  }
  while ((OneVarF(a,b,exp(lower))-v) < 0) {
    upper <- lower
    lower <- lower*2
  }
  while(upper-lower > 1e-6 && iterations < 1000) {
    hold <- (upper+lower)/2
    if(OneVarF(a,b,exp(hold))-v < 0) {
      lower <- hold
    }
    else {
      upper <- hold
    }
    iterations <- iterations + 1
  }
  return(exp((upper+lower)/2))
}

alternateRootFinder(0.4,0.5,0.5)

newtonsMethod(0.4,0.5,0.9,57)

alternateRootFinder(0.4,0.5,0.5)

qGGEE(0.5,0.4,0.5)

system.time(sapply(seq(0.501,0.6,0.001), function(x) {newtonsMethod(x,0.4,0.5,1)}))

system.time(sapply(seq(0.501,0.6,0.001), function(x) {alternateRootFinder(0.4,0.5,x)}))

rootFinder(0.4,0.5,0.8)

alternateRootFinder(0.4,0.5,0.5)

qGGEE(0.9,0.4,0.5)

newtonsMethod(0.4,0.5,0.9,57)

system.time(sapply(seq(0.01,0.99,0.01), function(x) {qGGEE(0.4,0.5,x)}))

system.time(sapply(seq(0.01,0.99,0.01), function(x) {tempFunction(x,0.4,0.5)}))

##uniroot,code in copula one,newton raphson method in r
##Find Copula?


rootOneVar <- function(a,b,v) {
  return(uniroot(function(x) {1-((a/(a+b))*Re(hypergeo(1,b,a+b+1,1-x)))-v},lower=0,upper=100)$root)
}

qGGEE(0.8,0.4,0.1)
system.time(qGGEE(0.8,0.4,0.5))
rootOneVar(0.4,0.5,0.8)

system.time(rootOneVar(0.4,0.5,0.8))

Copula <- function(a,b,v1,v2) {
  x1 <- rootFinder(a,b,v1)
  x2 <- rootFinder(a,b,v2)
  return(TwoVarF(a,b,x1,x2))
}

Copula(0.9,0.5,0.2,0.7)

a = 0.4
b = 0.1
v = 0.8
x = 10000000000
1-((a/(a+b))*Re(hypergeo(1,b,a+b+1,1-x)))-v

plot(1:200,sapply(log(1:200,base = exp(1)),function(x) {1-((a/(a+b))*Re(hypergeo(1,b,a+b+1,1-x)))}))
plot(1:200,sapply(1:200,function(x) {log(1-((a/(a+b))*Re(hypergeo(1,b,a+b+1,1-x))),base = exp(1))}))
plot(1:200,sapply(1:200,function(x) {exp(1-((a/(a+b))*Re(hypergeo(1,b,a+b+1,1-x))))}))

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
    
    while ((i < 1000) && (abs(v - CDF) > 1e-6)) {
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
    }
    out <- exp(t)
  }
  return(out)
}

rootFinder(0.4,0.5,0.1)
qGGEE(0.1,0.4,0.5)
system.time(sapply(seq(0.001,0.001,0.1), function(x) {rootFinder(0.4,0.1,x)}))
system.time(sapply(seq(0.001,0.001,0.1), function(x) {rootFinderBase(0.4,0.1,x)}))

str(1)

rootFinderBase <- function(a,b,v)
{
  if (v == 0) {
    out <- 0
  }
  else {
    CDF <- -0.01
    DEN <- 1
    i <- 0
    t <- 0
    lower <- 0
    upper <- 0.5
    
    while ((i < 1000) && (abs(v - CDF) > 1e-6)) {
      i <- i + 1
      t <- t - (CDF - v)/DEN
      if (t < lower || t > upper) {
        t <- 0.5 * (lower + upper)
      }
      DEN <- OneVarFPDF(a,b,t)
      CDF <- OneVarF(a,b,t)
      if (CDF < v) {
        lower <- t
      }
      else {
        upper <- t
      }
    }
    out <- t
  }
  return(out)
}

rootFinder(0.4,0.1,0.05)
rootFinderBase(0.4,0.1,0.05)


##empiricalCDF
##Copula
##CopulaPDF
##MLE
##Cores And Multithreading
##First Three Function In Copula One understood

