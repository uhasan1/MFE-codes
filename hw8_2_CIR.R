# I have calculated all of the functions based on a face value of 1 dollar for future consistency, 
#in this case, K should be put as the percentage of Spot. 
#Then easily can multiply the final answer, either for bond or option, in the face value. 

r0 = 0.05
sigma = 0.12
rbar = 0.055
kappa = 0.92

n = 100    # number of simulations
T = 0.5
S = 1
K = 0.98
dt = 1/252
######## zcbprice_CIR_MC ########
zcbprice_CIR_MC = function(r0,T,sigma,kappa,rbar,n){
  
  dt = 1/252 # time step is assumed to be a day
  h = round(T/dt) # number of intervals
  
  r = matrix(data = NA ,n,h+1)  # matrix to hold short rate paths
  r[,1] = r0
  
  W = matrix(data = rnorm(n*(h+1)),nrow = n,ncol = h+1) 
  # the interest rate path
  for(j in 1:h){
    dr = kappa*(rbar-r[,j])*dt + sigma*sqrt(r[,j])*sqrt(dt)*W[,j]
    r[,j+1] = r[,j] + dr
  }
  
  P.0.T = mean(exp(-dt*apply(X = r,MARGIN = 1, FUN = sum)))
  
  return((P.0.T))
}


zcbprice_CIR_explicit = function(r0,T,sigma,kappa,rbar){
  h1= sqrt(kappa^2+2*sigma^2)
  h2= (kappa+h1)/2
  h3= (2*kappa*rbar)/sigma^2
  
  A= ((h1*exp(h2*T))/(h2*(exp(h1*T)-1)+h1))^h3
  B= (exp(h1*T)-1)/(h2*(exp(h1*T)-1)+h1)
  
  P.0.T = A*exp(-B*r0)
  return(P.0.T)
}

zcbprice_CIR_MC(r0 = r0,T = 0.5,sigma = sigma,kappa = kappa,rbar = rbar,n = 1000)
zcbprice_CIR_explicit(r0 = r0,T = 0.5,sigma = sigma,kappa = kappa,rbar = rbar)
  


##### call_CIR_MC_MC #####
call_CIR_MC_MC = function(T, S, K, r0, sigma, kappa, rbar, n){
  if (S < T){
    cat ("bond maturity cannot be less than option maturity")
    return()
  }
  
  dt = 1/252 # time step is assumed to be a day
  h = round(T/dt) # number of intervals
  
  r1 = matrix(data = NA ,nrow = n,ncol = h+1)  # matrix to hold short rate paths
  r1[,1] = r0
  
  W = matrix(data = rnorm(n*(h+1)),nrow = n,ncol = h+1) 
  # the interest rate path
  for(j in 1:h){
    dr = kappa*(rbar-r1[,j])*dt + sigma*sqrt(r1[,j])*sqrt(dt)*W[,j]
    r1[,j+1] = r1[,j] + dr
  }
  zcbprice_CIR_MCs = matrix(nrow=n,ncol = 1)
  for (i in 1:n){
    zcbprice_CIR_MCs[i,1] = zcbprice_CIR_MC(T = S-T, r0 = r1[i,ncol(r1)], sigma = sigma,kappa = kappa,rbar = rbar,n = n)
  }
  output = zcbprice_CIR_MCs - K
  output_positive=(0<output)*output
  call.discounted = mean(exp(-dt*apply(X = r1,MARGIN = 1, FUN = sum))*output_positive)
  return(mean(call.discounted))
}


####### call_CIR_explicit #######
call_CIR_explicit = function(T ,S ,K = 0.98,r0 = r0,sigma = sigma,kappa = kappa,rbar = rbar){
theta = sqrt(kappa^2 + 2*sigma^2)
phi = 2*theta/(sigma^2*(exp(theta*T)-1))
gamma = (theta+kappa)/(sigma^2)

h1 = sqrt(kappa^2+2*sigma^2)
h2 = (kappa+h1)/2
h3 = (2*kappa*rbar)/sigma^2


A = ((h1*exp(h2*T))/(h2*(exp(h1*T)-1)+h1))^h3
B = (exp(h1*T)-1)/(h2*(exp(h1*T)-1)+h1)
r_star = (log(A/K))/B

x_p= 2*r_star*(phi+gamma+B) 
df_p= 4*kappa*rbar/(sigma^2)
ncp_p= 2*phi^2*r0*exp(theta*T)/(phi+gamma+B) 

x_k= 2*r_star*(phi+gamma)
df_k= (4*kappa*rbar)/(sigma^2)
ncp_k= (2*(phi^2)*r0*exp(theta*T))/(phi+gamma)

P_T = zcbprice_CIR_explicit(kappa = kappa, rbar = rbar, r0 = r0, sigma = sigma,T = T) 
P_S = zcbprice_CIR_explicit(kappa = kappa, rbar = rbar, r0 = r0, sigma = sigma,T = S) 

call_CIR_explicit =  P_S * pchisq(q = x_p,df = df_p,ncp = ncp_p) - K * P_T * pchisq(q = x_k,df = df_k,ncp = ncp_k) 
return(call_CIR_explicit)
}

###### finite difference method #####
r_lower= 0
dr = 0.001
dt =1/365
n_r = (r0 - r_lower)/dr
r_upper = r0+(n_r*dr)
size = 2 * n_r + 1
r_series = seq(from = r_lower, to = r_upper, by = dr)

bond = matrix(data = 0, nrow = size, ncol = 1)
call = matrix(data = 0, nrow = size, ncol = 1)


for (i in (1:size)){
  bond[i,1]= zcbprice_CIR_explicit(r0 = r_series[i], kappa = kappa, rbar =  rbar, sigma =  sigma, T= S-T)
  call[i,1]= pmax(bond[i]-K, 0)
}

A= matrix(data = 0, nrow = size, ncol = size)

for (j in  1:(n_r*2-1)){
  A[j+1,j+2] = (-dt)/(2*dr) * (kappa*(rbar - j*dr) + j* sigma^2)
  A[j+1,j+1] = 1+ j* (dt/dr) * sigma^2 + j*dt*dr
  A[j+1,j] = dt/(2*dr) * (kappa*(rbar - j*dr) - (j*sigma^2))
}

A[size,size-1]= 1
A[size,size]= -1
A[1,1] = 1
A[1,2] = -1


X = call

for (i in n:1){
  B= matrix(data = 0, nrow = size , ncol = 1)
  B[1,1] = bond[1,1]- bond[2,1] 
  B[2:n_r*2,1] = X[2:n_r*2,1]
  B[n_r*2+1,1] = 0
  X = solve(A,B)
}

#middle part of the matrix
# my answer is not consistent with the two others.

X[n_r-1,1]

##################

call_CIR_explicit(T = 0.5,S = 1,K = 0.98,r0 = r0,sigma = sigma,kappa = kappa,rbar = rbar)
call_CIR_MC_MC(T = 0.5,S = 1,K = 0.98,r0 = r0,sigma = sigma,kappa = kappa,rbar = rbar,n = 1000)

# the monte carlo and explicit method yield the same result, though I couldn't get the same answer in finite difference method.

