# in this question, again, just take care of cube root in R
X0 = 1;  

T = 3;
h =100; 
Dt = T/h; 
n = 10000; # number of simulations

# W= matrix(data = Narray1_BM,nrow = n, ncol = h)
# Z= matrix(data = Narray2_BM,nrow = n, ncol = h)


W = matrix(rnorm(n*h),nrow = n,ncol = h)
Z = matrix(rnorm(n*h),nrow = n,ncol = h)

X = matrix(nrow=n, ncol = h+1)

Y=exp(-0.08*T+(1/3)*sqrt(T)*W[,h]+(3/4)*sqrt(T)*Z[,h])
expectedY = mean(sign(Y) * (abs(Y)^(1/3)))

X[,1]=X0

  for (j in 1:h){
    X[,j+1]=X[,j]+(0.25)*X[,j]*Dt + (1/3)*X[,j]*sqrt(Dt)*W[,j] - (3/4)*X[,j]*sqrt(Dt)*Z[,j] 
  }

expectedX = mean(sign((X[,h+1])) * abs((X[,h+1]))^(1/3))

cat("E1 is ", expectedX, " and E2 is ", expectedY,sep = "")
