# this time generating my own normal

m=2^31-1
a=7^5
b=0

n= 1000000 # equal to n*h later
#first n uniform numbers
U1= rep(NA,n)
U1[1]=5
for (i in 1:(length(U1)-1)){
  U1[i+1]=(a*U1[i]+b)%%m
}
U1=U1/m
# second n uniform numbers
U2 = rep(NA,n)
U2[1]=6
for (i in 1:(length(U2)-1)){
  U2[i+1]=(a*U2[i]+b)%%m
}
U2=U2/m

Narray1_BM = rep(NA,n)
Narray2_BM = rep(NA,n)

Narray1_BM = sqrt(-2*log(U1))*cos(2*pi*U2)
Narray2_BM = sqrt(-2*log(U1))*sin(2*pi*U2)



X0 = 1;  
Y0 = 0.75; 

T = 3;
h =1000; 
Dt = T/h; 
n = 1000; # number of simulations
t=0;

W= matrix(data = Narray1_BM,nrow = n, ncol = h)
Z= matrix(data = Narray2_BM,nrow = n, ncol = h)


# W = matrix(rnorm(n*h),nrow = n,ncol = h)
#Z = matrix(rnorm(n*h),nrow = n,ncol = h)

Y_1 = matrix(nrow=n, ncol = h+1)
Y_2 = matrix(nrow=n, ncol = h+1)

Y_1[,1]=Y0
Y_2[,1]=Y0


  for (j in 1:h){
    Y_1[,j+1]=Y_1[,j]+((2/(1+t))*Y_1[,j]+(1+t^3)/3)*Dt+((1+t^3)/3)*sqrt(Dt)*Z[,j]
    Y_2[,j+1]=Y_2[,j]+((2/(1+t))*Y_2[,j]+(1+t^3)/3)*Dt-((1+t^3)/3)*sqrt(Dt)*Z[,j]
    t=t+Dt
  }


Y=rbind(Y_1,Y_2)
expectedY = mean(Y[,h+1])
sdY= sd(Y[,h+1])
YLconfint = expectedY-1.96*sdY/sqrt(nrow(Y))
YUconfint = expectedY+1.96*sdY/sqrt(nrow(Y))
cat("E(Y3) is ", expectedY, " within the confidence interval [", YLconfint," , ",YUconfint, "]",sep = "")

