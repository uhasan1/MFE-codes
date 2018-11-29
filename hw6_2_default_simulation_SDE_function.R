mydefault <- function (landa1 = 0.2, landa2 = 0.4, T = 5) 
{
# the value of jumps is selected from a double exponential distribution; search jump in the code to change it
library(smoothmest) # in order to make double exponential
library(LSMonteCarlo)
V0 = 20000
L0 = 22000
mu = -0.1
sigma = 0.2
gamma = -0.4
landa1 = landa1

r0 = 0.02
delta = 0.25
landa2 = landa2
alpha = 0.7
epsilon = 0.95

T =T; 
h =100;# steps
Dt = T/h; 
nsim = 10000; # number of simulations


R= r0 + delta*landa2

r= R/12
n =12*T
PMT = (L0*r)/(1-(1+r)^(-n))
a = PMT/r
b = a - L0;   # also: PMT/r - L0
              # also: PMT/(r*(1+r)^n)
c = 1+r

t = seq(from = 0, to = T, by = Dt)
Lt = matrix(data = a-b*c^(12*t),nrow = nsim,ncol = h+1,byrow = T)


beta  = (epsilon - alpha)/T
qt = alpha + beta*t

# we are using the other way this time in order to turn * to +  

# onescol = matrix(data = 1,nrow = nsim, ncol = 1)
# ranmat = matrix(data = rnorm(nsim*h),nrow = nsim,ncol = h)  
# expranmat= exp((mu-0.5*sigma^2)*Dt+sigma*sqrt(Dt)*ranmat)
# St = S0*cbind(onescol,t(apply(X = expranmat,MARGIN = 1, FUN = cumprod)))

# another way! turn * to +  same answer
zeroscol = matrix(data = 0,nrow = nsim, ncol = 1)
Nt1 = matrix(data = rpois(n = nsim*h,lambda = landa1*Dt),nrow = nsim , ncol = h)  # number of jumps which occur in each interval
Nt2 = cbind(zeroscol,matrix(data = rpois(n = nsim*h,lambda = landa2*Dt),nrow = nsim , ncol = h)) # number of times that shit happens!
LnYt = matrix(data = rdoublex(nsim*h, mu = gamma,lambda = 0.1),nrow = nsim , ncol = h) # jump sizes

ranmat = matrix(data = rnorm(nsim*h) , nrow = nsim , ncol = h)
expranmat2 = (mu-0.5*sigma^2)*Dt + sigma*sqrt(Dt)*ranmat + Nt1*(LnYt) #
Vt = exp(log(V0)+cbind(zeroscol,t(apply(X = expranmat2,MARGIN = 1, FUN = cumsum)))) # simulationg the path for collateral

path =  (Vt<qt*Lt)  # where (Vt<qt*Lt)
pathnew = LSMonteCarlo::firstValueRow(path) # first time where (Vt<qt*Lt)
pathrecovery = pathnew * pmax((Lt - epsilon*Vt),0) # how much we recover when (Vt<qt*Lt)

Nt2new= LSMonteCarlo::firstValueRow(Nt2) # first time when shit happens!
Nt2recovery = Nt2new * abs(Lt - epsilon*Vt) # how much we recover when shit happens

#remove(expranmat2,LnYt,Lt,Nt1,ranmat,path,Nt2,Vt,zeroscol) # remove unwanyed variables for memory space

default = Nt2new|pathnew # when default 
default = LSMonteCarlo::firstValueRow(default) #when default for the first time
default_prob = sum(default)/nsim #number of defaults divided by total number of simulations

default_time = matrix(data = seq(from=0, to=T,by=Dt),nrow = nsim, ncol=h+1,byrow = T)
default_time = default * default_time # keep only times of default
expected_default_time = sum(default_time)/sum(default) # only if a default happened ehen is that

recovery = pmax(pathrecovery,Nt2recovery)  # recovery is maximum of the two
recovery= default*recovery*exp(-default_time*r0) # how much we recover when first time default
expected_default_recovery = mean(recovery)/default_prob

mydefault <- list("D" = expected_default_recovery, "Prob" = default_prob, "Et" = expected_default_time)

return(mydefault)

}
