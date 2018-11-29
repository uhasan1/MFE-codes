fdi_S = function(S0, strike, r, T, sigma, m = 100 ,h = 1000,type = c("call", "put"), style = c("european", "american")) {
  # m determines the convergance rate, don't set it too much!
  
  Dt = T / h
  m = 2*trunc((m+1)/2)                         
  
  S0.col = seq(from = 0, to = 2*S0, length.out = m+1)     
  
  f = matrix(data = NA, nrow=h+1,ncol = m+1)
  if (type == "call"){
    f[h+1,] = pmax(S0.col - strike, 0)
    f[,m+1] = S0.col[m+1] - strike
    f[,1] = 0
  }
  else if (type == "put"){
    f[h+1,] = pmax(strike - S0.col, 0)
    f[,m+1] = 0
    f[,1] = strike
  }
  else {
    cat("wrong input for option type! only 'call' or 'put' please! ")
  }
  
  
  for (i in h:1) {             # Iterate from last col to first col
    A = matrix(0, nrow=m-1, ncol=m-1)
    B = f[i+1,2:m]
    for (j in 1:(m-1)) {
      pu = 0.5*Dt*(r*j - sigma^2*j^2) 
      pm = 1 + Dt*(sigma^2*j^2 + r)
      pd = -0.5*Dt*(r*j + sigma^2*j^2)
      if (j == 1) {                  
        A[j,c(1,2)] = c(pm, pd)
        B[1] = B[1] - pu*f[i,1]
      }
      else if (j < m-1){
        A[j,c(j-1,j,j+1)] = c(pu, pm, pd)
      }
      else if (j == m-1){                           
        A[j,c(j-1,j)] = c(pu, pm)
        B[m-1] = B[m-1] - pd*f[i,m+1]
      }
    }
    f[i,2:m] = solve(A, B)
    
    if (type == "put" && style == "american")
      f[i,] = pmax(f[i,], strike - S0.col)
  }
  
  return(f[1, m/2+1])
}