# this program needs revision. it does not give correct answers


# what I have done is to calculate risk neutral probability of price at t=0.2
# I set the strike price equal to this expected value
# calculate the price of american and european option with T-t to maturity
# discount the price back by DF= exp(-t*r)

# constant factors in question
library(OptionPricing)
library(LSMonteCarlo)
sigma = 0.2
r = 0.06
delta = 0
S0 = 65
# strike = 60  what was this?! I assumed: strike=S(t) if not uncomment this line
t = 0.2
T =1

h = 50
n = 1000

K = S0*exp(t*r) #risk neutral price growth
# we would get the same result if this price came from the simulation

euput = OptionPricing::BS_EP(T= T-t,K=K, sigma = sigma,r = r,S0 = K)
fs_Eu_put = exp(-r*t)*euput[1]

# I am getting a wrong answer. I will look into it later if possible.
# amput = my_Put_LSMC(S0 = K, sigma = sigma, strike = K, r = r,delta = 0,T = T-t, n = n, h = h, terms = 3,method = "monomials")$price
amput = AmerPutLSM(Spot = K,sigma = sigma,n = 1000,m = 50,Strike = K,r = r,dr = 0,mT = T-t)
fs_Am_put = exp(-r*t)*amput
