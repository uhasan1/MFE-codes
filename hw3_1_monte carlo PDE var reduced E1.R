X0 = 1;  

T = 2;
h =1000; # how many steps
Dt = T/h; 
n = 1000; # how many simulations

W= matrix(data = rnorm(n*h),nrow = n, ncol = h)

X_1 = matrix(nrow=n, ncol = h+1)
X_2 = matrix(nrow=n, ncol = h+1)

X_1[,1]=X0
X_2[,1]=X0


  for (j in 1:h){
    X_1[,j+1]=X_1[,j]+(0.2-0.5*X_1[,j])*Dt + (2/3)*sqrt(Dt)*W[,j]
    X_2[,j+1]=X_2[,j]+(0.2-0.5*X_2[,j])*Dt - (2/3)*sqrt(Dt)*W[,j]
  }

X=rbind(X_1,X_2)
# take care of cube root in R
# take care that E(x^1/3)!= E(x)^1/3
X2cubesquared = sign((X[,h+1])) * abs((X[,h+1]))^(1/3)
expectedX2cubesquared = mean(X2cubesquared)
sdX= sd(X2cubesquared)
# note that nrow(x)= 2*nrow(x_1), this way variance is reduced.
Lconfint=expectedX2cubesquared-1.96*sdX/sqrt(nrow(X))
Uconfint=expectedX2cubesquared+1.96*sdX/sqrt(nrow(X))
cat("E(X2^(1/3)) is ", expectedX2cubesquared, " within the confidence interval [", Lconfint," , ",Uconfint, "]",sep = "")
