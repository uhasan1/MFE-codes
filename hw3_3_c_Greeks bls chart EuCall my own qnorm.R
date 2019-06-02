CumPNorm <- function(q){
  d1=0.0498673470
  d2=0.0211410061
  d3=0.0032776263
  d4=0.0000380036
  d5=0.0000488906
  d6=0.0000053830
  if ((q>0)|(q==0))  N=1-0.5*(1+d1*q+d2*q^2+d3*q^3+d4*q^4+d5*q^5+d6*q^6)^(-16)
  else N=1-CumPNorm(-q)
  N
}

BlsCalldelta <- function(S0, K, r, T, sigma) {
  d1 <- (log(S0/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  BlsCalldelta= pnorm(d1)        #*exp(-delta*T)
  BlsCalldelta
}

BlsCallgamma <- function(S0, K, r, T, sigma) {
  d1 <- (log(S0/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  BlsCallgamma <- dnorm(d1)/(S0*sigma*sqrt(T))      #*exp(-delta*T) in Soorat
  BlsCallgamma
}

BlsCallvega <- function(S0, K, r, T, sigma) {
  d1 <- (log(S0/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T))
  d2 <- d1 - sigma * sqrt(T) 
  BlsCallvega = S0*dnorm(d1)*sqrt(T) #*exp(-delta*T) in Soorat
  BlsCallvega
}
  
BlsCalltheta <- function(S0, K, r, T, sigma) {
  d1 <- (log(S0/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  BlsCalltheta = (-r)*K*exp(-r*T)*pnorm(d2)- (K*exp(-r*T)*dnorm(d2)*sigma)/(2*sqrt(T)) 
  BlsCalltheta
}

BlsCallrho <- function(S0, K, r, T, sigma) {
  d1 <- (log(S0/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T))
  d2 <- d1 - sigma * sqrt(T) 
  BlsCallrho = T*K*exp(-r*T)*pnorm(d2)
  BlsCallrho
}
  
S0=seq(from=15, to=25, by=1)  

delta = BlsCalldelta(S0 = S0 ,K = 20,r = 0.04,T =0.5 ,sigma = 0.25)
plot(S0,delta,type = "b")

gamma = BlsCallgamma(S0 = S0 ,K = 20,r = 0.04,T =0.5 ,sigma = 0.25)
plot(S0,gamma,type = "b")

vega = BlsCallvega(S0 = S0 ,K = 20,r = 0.04,T =0.5 ,sigma = 0.25)
plot(S0,vega,type = "b")

theta = BlsCalltheta(S0 = S0 ,K = 20,r = 0.04,T =0.5 ,sigma = 0.25)
plot(S0,theta,type = "b")

rho = BlsCallrho(S0 = S0 ,K = 20,r = 0.04,T =0.5 ,sigma = 0.25)
plot(S0,rho,type = "b")
