
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

BlackscholesCall <- function(S0, K, r, T, sigma) {

  d1 <- (log(S0/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  
  Eucall <- S0*CumPNorm(d1) - K*exp(-r*T)*CumPNorm(d2)
  # Euput  <- K*exp(-r*T) * CumPNorm(-d2) - S0*CumPNorm(-d1)
  Eucall
}

BlackscholesCall(S0 = 88,K = 100,r = 0.04,T = 5,sigma = 0.2)
