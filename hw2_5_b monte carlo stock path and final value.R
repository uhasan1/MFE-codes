n= 6000

# m=2^31-1
# a=7^5
# b=0
# 
# 
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

S0 = 88;  
sigma = 0.18; 
r = 0.04;

T = 10;
h =1000; 
Dt = T/h; 
n = 6; # number of simulations

onescol = matrix(data = 1,nrow = n, ncol = 1)
ranmat = matrix(data = rnorm(n*h),nrow = n,ncol = h)  
expranmat= exp((r-0.5*sigma^2)*Dt+sigma*sqrt(Dt)*ranmat)
St = S0*cbind(onescol,t(apply(X = expranmat,MARGIN = 1, FUN = cumprod)))

# another way! turn * to +  same answer
zeroscol = matrix(data = 0,nrow = n, ncol = 1)
expranmat2 = (r-0.5*sigma^2)*Dt+sigma*sqrt(Dt)*ranmat
St2 = exp(log(S0)+cbind(zeroscol,t(apply(X = expranmat2,MARGIN = 1, FUN = cumsum))))

matplot(t(St),type = "l",xlab = "total number of time intervals", ylab = "stock price", main="stock path")

#----------- UNCOMMENT TO HAVE BOTH GRAPHS ON SAME PLOT ---------------#

T= seq(100,1000,by=100)
points(T,rowMeans(ST),main = "E(ST) vs T ",ylab = "S(T)",col="red",pch=15)

# repeat with sigma=0.35, stock price path will have greater variation and a higher range for final value

