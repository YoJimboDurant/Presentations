library(plyr)
library(NADA)

smallKM = function(obs,cen){
  # this function takes NADA format data and returns KM mean
  #   method is based off of D. Helsel 2012 Book. I wrote it for bootstrapping
  #   Function does not require the overhead of the cenfit and runs approximately 10
  #   times faster. 
  
  #   Args:
  #     obs: observed value or censoring limit
  #     cen: logical value if obs values is censored or not
  # 
  #   Returns:
  #     mean value using KM statistics
  #     
  num.risk <- cumsum(table(obs))
  num.det <- table(obs * !cen)
  num.det <-num.det[rownames(num.det)!=0]
  det.values <-as.numeric(rownames(num.det))
  min.stime=min(0,min(det.values))
  num.risk <- num.risk[rownames(num.det)]
  S <- cumprod(((num.risk -num.det)/num.risk)[length(num.risk):1])
  
  sum(c(diff(c(min.stime, det.values[length(det.values):1])),0) *
        c(1, S))
}

samplesize <- rep(100,3)
detectionlimits <- c(0.5,1,3)
proportions <- c(0.1,0.5,0.4)
logMu <- rep(0,3)
logSD <- rep(1,3)

myrun <- data.frame(samplesize,detectionlimits, proportions, logMu, logSD)

makeCenData <- function(samplesize, detectionlimits, proportions, logMu, logSD){
  n <- round(samplesize * proportions, 0)
  values <- rlnorm(n, logMu, logSD)
  cen <- values < detectionlimits
  obs <- values
  obs[cen] <- detectionlimits
  data.frame(n,values,obs,cen)
}

mydata <- mdply(myrun, makeCenData)
  

KMuncen <- cenfit(mydata$value, rep(FALSE, length(mydata$value)), conf.int=0)
KMcen <- cenfit(mydata$obs, mydata$cen)

plot(KMuncen)

KMVals <- KMcen@survfit$time
KMProbs <- KMcen@survfit$surv

.x <- seq(from=0, to=max(KMVals), length=1000)

# do pdf of data
plot(.x, dlnorm(.x, logMu, logSD), type="l", xlab="Value", col="blue", ylab= "Density", ylim=c(0,1.01))
legend("topleft", "Probability Density Function", cex=0.8, col="blue", 
       lty=1, lwd=2, bty="n")

# show a probability of a value
plot(.x, dlnorm(.x, logMu, logSD), type="l", xlab="Value", col="blue", ylab= "Density", ylim=c(0,1.01))

cutpoint <- qlnorm(0.5,0,1)
polygon(c(.x[.x <= cutpoint], rev(.x[.x <= cutpoint])),        
        c(dlnorm(.x[.x <= cutpoint], 0,1), rep(0, sum(.x <= cutpoint))),        
        col = grey(.5))

legend("topleft", "Probability Density Function", cex=0.8, col="blue", 
       lty=1, lwd=2, bty="n")

text(cutpoint, dlnorm(cutpoint, logMu, logSD), paste("p(x \u2266 ", cutpoint,") = ", 
                        plnorm(cutpoint, logMu, logSD)), adj=c(-0.5,0))


plot(.x, plnorm(.x, logMu, logSD), type="l", xlab="Value", col="red", ylab= "p(x \u2266 X)", log="x", ylim=c(0,1.01))
legend("topleft", "Cumulative Distribution Function", cex=0.8, col="red", 
       lty=1, lwd=2, bty="n")
abline(h=1)

text(cutpoint, plnorm(cutpoint, logMu, logSD), paste("p(x \u2266 ", cutpoint,") = ", 
                                                     plnorm(cutpoint, logMu, logSD)), adj=c(1.1,0))

points(cutpoint, plnorm(cutpoint, logMu[1],logSD[1]), col="blue", pch=19)


#Survival Function
plot(.x, 1-plnorm(.x, logMu, logSD), type="l", xlab="Value", col="red", ylab= "p(x > X)", log="x", ylim=c(0,1.01))
abline(h=1)
legend("topright", "Survival Function", cex=0.8, col="red", 
       lty=1, lwd=2, bty="y", bg="white")


text(cutpoint, 1-plnorm(cutpoint, logMu, logSD), paste("p(x > ", cutpoint,") = ", 
                                                     plnorm(cutpoint, logMu, logSD)), adj=c(1.1,0))

points(cutpoint, plnorm(cutpoint, logMu[1],logSD[1]), col="blue", pch=19)

polygon(c(.x, rev(.x)),        
        c(1-plnorm(.x, logMu,logSD), rep(0, length(.x))),        
        col = "skyblue")

#Survival Function
plot(.x, 1-plnorm(.x, logMu, logSD), type="l", xlab="Value", col="red", ylab= "p(x > X)", log="x", ylim=c(0,1.01))
abline(h=1)
legend("topright", "Survival Function", cex=0.8, col="red", 
       lty=1, lwd=2, bty="y", bg="white")


polygon(c(.x, rev(.x)),        
        c(1-plnorm(.x, logMu,logSD), rep(0, length(.x))),        
        col = "skyblue")


plot(KMuncen, ylab= "p(x \u2266 X)")
legend("topleft", "Emperical Distribution Function (sample)", cex=0.8, col="black", 
       lty=1, lwd=2, bty="n")


plot(KMuncen, ylab= "p(x \u2266 X)")
lines(.x, plnorm(.x, logMu, logSD), type="l", col="red")
legend("topleft", c("Cumulative Distribution Function (population)", "Emperical Cumulative Distribution Function (sample)"), cex=c(0.8,0.8), col=c("red","black"), 
       lty=1, lwd=2, bty="n")

plot(KMuncen, ylab= "p(x \u2266 X)")
lines(.x, plnorm(.x, logMu, logSD), type="l", col="red")
legend("topleft", c("Cummulative Distribution Function (population)", "Emperical Cumulative Distribution Function (complete sample)"), cex=c(0.8,0.8), col=c("red","black"), 
       lty=1, lwd=2, bty="n")

plot(KMcen, ylab= "p(x \u2266 X)")
lines(.x, plnorm(.x, logMu, logSD), type="l", col="red")
legend("topleft", c("Cumulative Distribution Function (population)", "Emperical Cumulative Distribution Function (complete sample)"), cex=c(0.8,0.8), col=c("red","black"), 
       lty=1, lwd=2, bty="n")

#Get the ECDF area.

#calculate area by multiplying cumulative proportions by distance
#between knots, then summing
k=c(knots(ecdf(mydata$value)))
ecdf.area=vector("numeric",(length(mydata$values)-1))
for(i in 1:(length(mydata$values)-1)){
  ecdf.area[i]=(k[i+1]-k[i])*(sum(mydata$values<=k[i])/length(mydata$values))
}
ecdf.area=sum(ecdf.area)

kmMeanNoCen <- max(mydata$values)-sum(ecdf.area)



meancdf = integrate(function(...) 1-plnorm(...), 0, Inf, logMu, logSD)

# try with cenlimit

mydata$obs1 <- mydata$values
mydata$cen1 <- mydata$obs1 < 1

km1cen <- cenfit(mydata$obs1, mydata$cen1)
plot(km1cen)
lines(.x, plnorm(.x, logMu, logSD), col="red")
