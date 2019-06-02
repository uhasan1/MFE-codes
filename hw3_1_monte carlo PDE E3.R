X0 = 1;  
Y0 = 0.75; 

# to increase precision change n
T = 2;
h =1000; # number of steps
Dt = T/h; 
n = 10000; # number of simulations
t=0;

# W= matrix(data = Narray1_BM,nrow = n, ncol = h)
# Z= matrix(data = Narray2_BM,nrow = n, ncol = h)


W = matrix(rnorm(n*h),nrow = n,ncol = h)
Z = matrix(rnorm(n*h),nrow = n,ncol = h)

X = matrix(nrow=n, ncol = h+1)
Y = matrix(nrow=n, ncol = h+1)

X[,1]=X0
Y[,1]=Y0


  for (j in 1:h){
    X[,j+1]=X[,j]+(0.2-0.5*X[,j])*Dt + (2/3)*sqrt(Dt)*W[,j]
    Y[,j+1]=Y[,j]+((2/(1+t))*Y[,j]+(1+t^3)/3)*Dt+((1+t^3)/3)*sqrt(Dt)*Z[,j]
    t=t+Dt
  }


Result= X[,h+1]*Y[,h+1]*(X[,h+1]>1)
expectedResult=mean(Result)
expectedResult
sdResult= sd(Result)
Res_confint = c(expectedResult-1.96*sdResult/sqrt(nrow(X)),expectedResult+1.96*sdResult/sqrt(nrow(X)))
Res_confint
