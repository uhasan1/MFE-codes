m=2^31-1
a=7^5
b=0

n= 1000
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

S0 = 88;  
sigma = 0.18; 
r = 0.04;

ST=matrix(nrow = 10,ncol = n)

for (T in 1:10){
  WT=sqrt(T)*Narray2_BM
  ST[T,] = S0*(exp((T*(r-0.5*sigma^2))+(sigma*WT)))
}
# Now, there are 10 rows in ST, each contain a series of possible final value of ST at time T #
T= seq(1,10,by=1)
plot(T,rowMeans(ST, dims = 1),main = "E(ST) vs T ",ylab = "S(T)")

# repeat with sigma=0.35, there will be no change since it is final expected value
