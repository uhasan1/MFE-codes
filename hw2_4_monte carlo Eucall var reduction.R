
m=2^31-1
a=7^5
b=0

n= 10000
# #first n uniform numbers
# U1= rep(NA,n)
# U1[1]=5
# for (i in 1:(length(U1)-1)){
#   U1[i+1]=(a*U1[i]+b)%%m
# }
# U1=U1/m
# # second n uniform numbers
# U2 = rep(NA,n)
# U2[1]=6
# for (i in 1:(length(U2)-1)){
#   U2[i+1]=(a*U2[i]+b)%%m
# }
# U2=U2/m
# 
# Narray1_BM = rep(NA,n)
# Narray2_BM = rep(NA,n)
# 
# Narray1_BM = sqrt(-2*log(U1))*cos(2*pi*U2)
# Narray2_BM = sqrt(-2*log(U1))*sin(2*pi*U2)

Z = rnorm(n)

S0 = 88; 
K = 100; 
sigma = 0.20; 
r = 0.04;
T = 5;

WT=sqrt(T)*Z
ST = S0*(exp((T*(r-0.5*sigma^2))+(sigma*WT)))

output=ST-K
zerovec = rep(0,n)
output_positive=(zerovec<output)*output
Eucall = exp(-r*T)*output_positive
avg=mean(Eucall)
avg
s= sd(Eucall)
s
avg+(1.96*s/sqrt(length(ST)))-(avg-(1.96*s/sqrt(length(ST))))

library(RQuantLib)
EuropeanOption("call",underlying = S0,strike = K,dividendYield = 0,riskFreeRate = r,maturity = T,volatility = sigma)

### variance reduction
WT=sqrt(T)*Z
ST1 = S0*(exp((T*(r-0.5*sigma^2))+(sigma*WT)))
ST2 = S0*(exp((T*(r-0.5*sigma^2))-(sigma*WT)))

ST = c(ST1,ST2)

output=ST-K
zerovec = rep(0,n)
output_positive=(zerovec<output)*output
Eucall = exp(-r*T)*output_positive
avg=mean(Eucall)
avg
s= sd(Eucall)
s
avg+1.96*s/sqrt(length(ST))-(avg-1.96*s/sqrt(length(ST)))
#variance reduced 30%!