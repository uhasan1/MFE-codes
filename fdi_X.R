fdi_X = function(S0=10, strike=10, r=0.04, T=0.5, sigma,m = 100 ,Dt = 0.002,type = c("call", "put"), style = c("european", "american")) {
  # m determines the coputation rate, don't set it too high!
  
  h = T / Dt
  m = 2*trunc((m+1)/2)                         
  X0 = log(S0)
  Xmin = X0 - 3*sigma*sqrt(T)
  Xmax = X0 + 3*sigma*sqrt(T)
  dX = (Xmax - Xmin)/m
  
  ## Set stock price limits to +/- 3 standard deviations.
  
  X.col = seq(from = Xmin, to = Xmax, length.out = m+1)
  f = matrix(NA, nrow=h+1,ncol =  m+1)
  
  if (type == "call"){
    f[h+1,] = pmax(exp(X.col)-strike,0)
    f[,m+1] = exp(X.col[m+1])-strike
    f[,1] = 0
  }
  else if (type == "put"){
    f[h+1,] = pmax(strike-exp(X.col),0)
    f[,m+1] = 0
    f[,1] = strike-exp(X.col[1])
  }
  else {
    cat("wrong input for option type! only 'call' or 'put' please! ")
  }
  
  
  pu = Dt/(2*dX)*(r - 0.5*sigma^2) - Dt/(2*dX^2)*sigma^2
  pm = 1 + Dt/dX^2*sigma^2 + r*Dt
  pd = -Dt/(2*dX)*(r - 0.5*sigma^2) - Dt/(2*dX^2)*sigma^2
  
  for (i in (h:1)) {             # Iterate from end to beginning.
    A = matrix(0, nrow=m-1, ncol=m-1)
    B = f[i+1,(2:m)]
    for (j in 1:(m-1)) {
      if (j == 1) {                
        A[j,1:2] = c(pm, pd)
        B[1] = B[1] - pu*f[i,1]
      }
      else if (j < m-1)
        A[j,c(j-1,j,j+1)] = c(pu, pm, pd)
      else if (j == m-1){                          
        A[j,c(j-1,j)] = c(pu, pm)
        B[m-1] = B[m-1] - pd*f[i,m+1]
      }
    }
    
    f[i,(2:m)] = solve(A, B)
    
    if (type == 'put' && style == 'american')
      f[i,] = pmax(f[i,], strike - exp(X.col))
  }
  
  return(f[1, m/2+1])
}
