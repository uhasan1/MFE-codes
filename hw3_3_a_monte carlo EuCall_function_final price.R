

MC_BlackscholesCall <- function(S0, K, r, T, sigma) {

n=1000000  # number of simulations
WT=sqrt(T)*rnorm(n)
ST1 = S0*(exp((T*(r-0.5*sigma^2))+(sigma*WT)))
ST2 = S0*(exp((T*(r-0.5*sigma^2))-(sigma*WT)))

ST = c(ST1,ST2)

output=ST-K
zerovec = rep(0,n)
output_positive=(zerovec<output)*output
Eucall = exp(-r*T)*output_positive
avg=mean(Eucall)
avg
# s= sd(Eucall)
# avg-1.96*s/sqrt(n)
# avg+1.96*s/sqrt(n)
}
MC_BlackscholesCall(S0 = 88,K = 100,r = 0.04,T = 5,sigma = 0.2)

