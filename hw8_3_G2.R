# I have calculated all of the functions based on a face value of 1 dollar for future consistency, 
#in this case, K should be put as the percentage of Spot. 
#Then easily can multiply the final answer, either for bond or option, in the face value. 

##### zcbprice_G2_MC #####
zcbprice_G2_MC = function(r0,x0,y0,a,b,T,sigma,etta,rho,n){
  
  dt = 1/365 # time step is assumed to be a day
  h = round(T/dt) # number of intervals
  
  x = matrix(data = NA ,n,h+1)  
  y = matrix(data = NA ,n,h+1)
  r = matrix(data = NA ,n,h+1)
  
  x[,1] = x0
  y[,1] = y0
  r[,1] = r0
  
  
  W1 = matrix(data = rnorm(n*(h+1)),nrow = n,ncol = h+1) 
  W2 = matrix(data = rnorm(n*(h+1)),nrow = n,ncol = h+1) 
  W_corr =  rho * W1 + sqrt(1-rho^2) * W2

  
  # the interest rate path
  for(j in 1:h){
    dx = -a*x[,j]*dt + sigma*sqrt(dt)*W1[,j]
    dy = -b*y[,j]*dt + etta*sqrt(dt)*W_corr[,j]
    dr = dx+dy
    x[,j+1] = x[,j] + dx
    y[,j+1] = y[,j] + dy
    r[,j+1] = r[,j] + dr
  }
  
  P.0.T = mean(exp(-dt*apply(X = r,MARGIN = 1, FUN = sum)))
  
  return((P.0.T))
}


#########  put_G2_MC_MC ###########
put_G2_MC_MC = function(T, S, K,r0,x0,y0,a,b,sigma,etta,rho,n){
  if (S < T){
    cat ("bond maturity cannot be less than option maturity")
    return()
  }
  
  dt = 1/365 # time step is assumed to be a day
  h = round(T/dt) # number of intervals
  
  x = matrix(data = NA ,n,h+1)  
  y = matrix(data = NA ,n,h+1)
  r = matrix(data = NA ,n,h+1)
  
  x[,1] = x0
  y[,1] = y0
  r[,1] = r0
  
  
  W1 = matrix(data = rnorm(n*(h+1)),nrow = n,ncol = h+1) 
  W2 = matrix(data = rnorm(n*(h+1)),nrow = n,ncol = h+1) 
  W_corr =  rho * W1 + sqrt(1-rho^2) * W2
  
  
  # the interest rate path
  for(j in 1:h){
    dx = -a*x[,j]*dt + sigma*sqrt(dt)*W1[,j]
    dy = -b*y[,j]*dt + etta*sqrt(dt)*W_corr[,j]
    dr = dx+dy
    x[,j+1] = x[,j] + dx
    y[,j+1] = y[,j] + dy
    r[,j+1] = r[,j] + dr
  }
  
  zcbprice_G2_MCs = matrix(nrow=n,ncol = 1)
  
  for (i in 1:n){
    zcbprice_G2_MCs[i,1] = zcbprice_G2_MC(r0 = r[i,ncol(r)],T = S-T,x0 = x0,y0 = y0,a = a,b = b,sigma = sigma,etta = etta,rho=rho,n = n) 
  }
  
  output = K - zcbprice_G2_MCs
  output_positive=(0<output)*output
  put.discounted = mean(exp(-dt*apply(X = r,MARGIN = 1, FUN = sum))*output_positive)
  return(put.discounted)
}

put_G2_explicit_explicit = function(T, S, K,r0,x0,y0,a,b,sigma,etta,rho,phi){

  V_T = (sigma^2/a^2)*(T + (2/a)*exp(-a*T) - (1/(2*a))*exp(-2*a*T) - 3/(2*a)) + (etta^2/b^2)*(T+ (2/b)*exp(-b*T) - 1/(2*b)*exp(-2*b*T) - 3/(2*b)) + (2*rho*sigma*etta)/(a*b)*(T + 1/a*(exp(-a*T)-1) + (1/b)*(exp(-b*T)-1) - (1/(a+b))*(exp(-(a+b)*T)-1)) 
  P_T = exp(-phi*T - (1/a)*(1-exp(-a*T))*x0 - (1/b)*(1-exp(-b*T))*y0 + 1/2*V_T)
  
  V_S = (sigma^2/a^2)*(S + (2/a)*exp(-a*S) - (1/(2*a))*exp(-2*a*S) - 3/(2*a)) + (etta^2/b^2)*(S+ (2/b)*exp(-b*S) - 1/(2*b)*exp(-2*b*S) - 3/(2*b)) + (2*rho*sigma*etta)/(a*b)*(S + 1/a*(exp(-a*S)-1) + (1/b)*(exp(-b*S)-1) - (1/(a+b))*(exp(-(a+b)*S)-1)) 
  P_S = exp(-phi*S - (1/a)*(1-exp(-a*S))*x0 - (1/b)*(1-exp(-b*S))*y0 + 1/2*V_S)
           
           
  SIGMA_SQUARED =( sigma^2/(2*a^3))*(1-exp(-a*(S-T)))^2*(1-exp(-2*a*T)) + etta^2/(2*b^3)*(1-exp(-b*(S-T)))^2*(1-exp(-2*b*T)) + 2*rho*sigma*etta/(a*b*(a+b))*(1-exp(-a*(S-T)))*(1-exp(-b*(S-T)))*(1-exp(-(a+b)*T))

  d1 = log((K*P_T)/(P_S))/sqrt(SIGMA_SQUARED) - 1/2* sqrt(SIGMA_SQUARED)
  d2 = log((K*P_T)/(P_S))/sqrt(SIGMA_SQUARED) + 1/2* sqrt(SIGMA_SQUARED)

  put = - P_S * pnorm(d1) + K * P_T * pnorm(d2)
  return(put)
}



put_G2_MC_MC(T = 0.5,S = 1,K = 0.95,r0 = 0.03,x0 = 0,y0 = 0,a = 0.1,b = 0.3,sigma = 0.03,etta = 0.08,rho = 0.7,n = 100)
put_G2_explicit_explicit(T = 0.5,S = 1,K = 0.95,r0 = 0.03,phi=0.03,x0 = 0,y0 = 0,a = 0.1,b = 0.3,sigma = 0.03,etta = 0.08,rho = 0.7)

# the answers are very close, the monte carlo method has a large standard deviation though and increasing
# the number of simulations will tremendously increase the time of computation since it is a 2 step monte carlo
