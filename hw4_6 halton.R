

number2digits=function(n,base){ 
  #first digit in output is the least significant 
  digit=n%%base 
  if (n<base) 
    digit 
  else 
    c(digit,number2digits((n-digit)/base,base)) 
} 

digits2number=function(digits,base){ 
  #first digit in input should be the most significant 
  result=0 
  for (digit in digits) 
    result=(base*result)+digit 
  result 
} 

halton=function(n,base=2){ 
  #generate n first digits of the van der Corput sequence 
  output=NA*1:n 
  for(i in 1:n){ 
    digits=number2digits(i,base) 
    output[i]=digits2number(digits,base)/base^length(digits) 
  } 
  output 
} 

halton=function(n,base=2){ 
  #generate n first digits of the van der Corput sequence 
  output=NA*1:n 
  for(i in 1:n){ 
    digits=number2digits(i,base) 
    output[i]=digits2number(digits,base)/base^length(digits) 
  } 
  output 
} 

Halton_Call <- function(S0, K, r, T, sigma,N,b1,b2) {
  
  n=N  # number of simulations
  
  H1=halton(n = N,base=b1)
  H2=halton(n = N,base=b2)
  
  Z1 = rep(NA,n)
  Z2 = rep(NA,n)
  
  Z1 = sqrt(-2*log(H1))*cos(2*pi*H2)
  Z2 = sqrt(-2*log(H1))*sin(2*pi*H2)
  
  Z=c(Z1,Z2)
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
}

Halton_Call(S0 = 88,K = 100,r = 0.04,T = 5,sigma = 0.2,N = 10000,b1 = 2,b2 = 5)

