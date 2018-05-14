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
  return(1-((a/(a+b))*Re(hypergeo(1,b,a+b+1,1-x))))
}
OneVarFPDF <- function(a,b,x) {
  return(1-(((a*b)/((a+b)*(1+a+b)))*Re(hypergeo(2,b+1,a+b+2,1-x))))
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

newtonsMethod <- function(v,x0,a,b) {
  x1 <- x0-((OneVarF(a,b,x0)-v)/(OneVarFPDF(a,b,x0)-v))
  if(abs(x1-x0) > 0.0001) {
    newtonsMethod(v,x1,a,b)    
  }
  else {
    return(x0)
  }
}
newtonsMethod(0.625,1,0.4,0.5)

rootOneVar(0.4,0.1,0.625)

qGGEE(0.8,0.4,0.1)

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
  x1 <- rootOneVar(a,b,v1)
  x2 <- rootOneVar(a,b,v2)
  return(TwoVarF(a,b,x1,x2))
}
a = 0.4
b = 0.1
v = 0.8
x = 10000000000
1-((a/(a+b))*Re(hypergeo(1,b,a+b+1,1-x)))-v

plot(1:20000,sapply(1:20000,function(x) {1-((a/(a+b))*Re(hypergeo(1,b,a+b+1,1-x)))}))

Copula(0.4,0.5,0.6,0.3)

lines(c(0,25,40,60,110),c(5.3,7.3,8.5,9.7,11.2))

"This is a test"