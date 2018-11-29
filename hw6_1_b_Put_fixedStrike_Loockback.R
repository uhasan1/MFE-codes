S0 = 98;  
sigma = 0.12; 
r = 0.03;
K= 100;

T = 1;
h =500; 
Dt = T/h; 
n = 10000; # number of simulations
sigma = seq(from= 0.12, to = 0.48, by = 0.04)
result = matrix(nrow = 2, ncol = length(sigma))

for (i in 1:length(sigma)){
  ##### generate stock path #####
  onescol = matrix(data = 1,nrow = n, ncol = 1)
  ranmat = matrix(data = rnorm(n*h),nrow = n,ncol = h)  
  expranmat= exp((r-0.5*sigma[i]^2)*Dt+sigma[i]*sqrt(Dt)*ranmat)
  St = S0*cbind(onescol,t(apply(X = expranmat,MARGIN = 1, FUN = cumprod)))
  
  ##### payoff fixed strike loockback Put #####
  min.each.row = apply(St,1, min)
  Put.fixedK.lookback = exp(-r*T)*(pmax((K-min.each.row),0))
  Put.fixedK.lookback.price = mean(Put.fixedK.lookback)
  # s= sd(Put.fixedK.lookback)
  # Put.fixedK.lookback.price.upper = Put.fixedK.lookback.price + 1.96*s/sqrt(n)
  # Put.fixedK.lookback.price.lower = Put.fixedK.lookback.price - 1.96*s/sqrt(n)
  result[1,i] = sigma[i]
  result[2,i] = Put.fixedK.lookback.price
}

plot(x = result[1,],y = result[2,],type = "p",xlab = "sigma", ylab = "Put Lookback Price")
