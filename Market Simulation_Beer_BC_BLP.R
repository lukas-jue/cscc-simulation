##################################################################################################
#####Simulation & Estimation - Nash with Budget###################################################
#####BLP-type of indirect utility#################################################################
##################################################################################################

rm(list=ls())

# LOAD LIBRARIES
library(xtable)
library(devtools)
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(bayesm)
library(ggplot2)
library(tikzDevice)
library(plyr)
library(latex2exp)
library(FixedPoint)
library(dplyr)


###Increase memory capacities
memory.limit(size=100000000)

load("Estimation_Data_Beer_20170423.Rdata")
products = c("Amstel Extra Lata 37,5 cl","Amstel Extra Lata 33 cl","Amstel Lata 37,5 cl","Amstel Lata 33 cl","Amstel Cl?sica Lata 33 cl",
             "Cruzcampo Lata 33 cl","Estrella Damm Lata 33 cl","Estrella Galicia Lata 33 cl","Heineken Lata 33 cl","Mahou 5 Estrellas Lata 33 cl",
             "Mahou Cl?sica Lata 33 cl","San Miguel Lata 33 cl","Voll Damm Lata 33 cl","Steinburg (Marca Blanca Mercadona) Lata 33 cl",
             "Marca Blanca Carrefour Lata 33 cl")

#reorder that price is first in the design matrix
for (i in 1:length(E_Data$lgtdata)){
  E_Data$lgtdata[[i]]$X <- E_Data$lgtdata[[i]]$X[,c(16,1:15)]
}

#Number of players
nplayers = 15

##########################################
#######Run BLP Budget model###############
##########################################
###Load complete sampler now...
Rcpp::sourceCpp("rhierMnlRwMixture_rcpp_loop_Illus_BLP_type.cpp",showOutput = FALSE)
source('rhierMnlRwMixture_main_untuned_BC.R')

#number of constrained coefficients (budget & price)
nvar_c = 2
#position of price coefficient in design matrix
pr=1

###Prior setting
Amu = diag(1/10, nrow = nvar_c, ncol = nvar_c)
mustarbarc = matrix(rep(0, nvar_c), nrow = nvar_c)
nu = 15 + nvar_c
V = nu * diag(nvar_c)*0.5

Prior = list(ncomp=1, Amu = Amu, mustarbarc = mustarbarc, nu = nu, V = V)
Mcmc = list(R=30000, keep=3)#, s=1.6)
#,s=c(0.1,0.5,0.5,0.5)
out_BC = rhierMnlRwMixture_SR(Data=E_Data,Prior=Prior,Mcmc=Mcmc,nvar_c=nvar_c,pr=pr,starting_budget = log(0.74))

betastar_HB_BC = out_BC$betadraw
compdraw_HB = out_BC$nmix$compdraw
probdraw_HB = out_BC$nmix$probdraw
rejection = out_BC$rejection
loglike_BC = out_BC$loglike

###Compute rejection rate of sampler 
rej_rate_indi = apply(rejection,2,mean)
summary(rej_rate_indi)
rej_rate_agg = mean(rej_rate_indi)
rej_rate_agg


########################
###Get rid of burnin####
########################

# check visually how much burn-in is required
plot(out_BC$loglike, type="l")

burnin = 2000
R = dim(betastar_HB_BC)[3]

betastar_HB_BC = betastar_HB_BC[,,(burnin+1):R]
compdraw_HB = compdraw_HB[(burnin+1):R]
probdraw_HB = probdraw_HB[(burnin+1):R]
rejection = rejection[(burnin+1):R,]
loglike_BC = loglike_BC[(burnin+1):R]
plot(loglike_BC, type="l")

R = dim(betastar_HB_BC)[3]
N = dim(betastar_HB_BC)[1]

###Heterogeneity distribution lower level model non-smoothed  
l=10
betastar_BC_LLMns <- array(0,dim=c(R*l,dim(betastar_HB_BC)[2]))
index_r = rep(rank(runif(R)),l)
index_n = rep(rank(runif(N)),round((R*l)/N)+1)

#data generation
for(i in 1:(R*l)){
  betastar_BC_LLMns[i,] = betastar_HB_BC[index_n[i],,index_r[i]]
}
#transform betastardraws to betadraws
beta_BC_LLMns = betastar_BC_LLMns
beta_BC_LLMns[,1] = exp(betastar_BC_LLMns[,1])
beta_BC_LLMns[,2] = exp(betastar_BC_LLMns[,2])

summary(beta_BC_LLMns)

#######################  
###Nash optimization
#######################

### Load functions for Fixed-point approach
source('Morrow_Skerlos_Implementations_BC_BLP_MarkupEquations_Final_Efficient.R')
Rcpp::sourceCpp("Speed++_MS_BC_BLP_Efficient.cpp",showOutput = FALSE)

###########################################################################
# Nash Equilibrium Prices and Shares for all 15 brands
###########################################################################

# ingredients
nplayers = 15 + 1 ### 15 inside + outside
min_p = 0.22
prices_init = c(rep(min_p,nplayers - 1),0)
designBase = rbind(diag(nplayers-1),rep(0,nplayers-1))
Xdata = cbind(prices_init,designBase); colnames(Xdata)[1] = "PRICE" 
Ownership = array(0,dim=c(nplayers,nplayers))
Ownership[1:(nplayers-1),1:(nplayers-1)] = diag((nplayers-1))

# define which brands compete on price, i.e. are owned by different companies
# without any specification, all brands compete (owned by different companies)
#Ownership[1,2] <- 1
#Ownership[2,1] <- 1

inside_row=which(rowSums(Ownership) != 0, arr.ind = TRUE)
p0=Xdata[inside_row,"PRICE"]
costBase = as.vector(prices_init*0.9)
MC = costBase[inside_row]

### Run Fixed-Point algorithm with Xi-markup equation (reliable and fast: See Table 3 in paper for comparison)
p_Markup_Xi_FixedPoint_BC_BLP = FixedPoint(Function = function(price_vec) FixedPoint_BLP_Xi(price_vec,MC=MC,
                                           ownership=Ownership,Xdata=Xdata,beta_draws=beta_BC_LLMns,pr=1), 
                                           Inputs = p0, MaxIter = 10000, ConvergenceMetricThreshold = 1e-10, Method = "Anderson")

p_Markup_Xi_FixedPoint_BC_BLP$FixedPoint

# Overview optimal prices & shares
Optimal_prices <- array(0,dim=c(nplayers,2))
rownames(Optimal_prices) = c(products[1:(nplayers-1)],"Outside")
colnames(Optimal_prices) = c("Equilibrium Price","Equilibrium Shares")
# Save equi-prices
Optimal_prices[,"Equilibrium Price"] <-c(p_Markup_Xi_FixedPoint_BC_BLP$FixedPoint,0)
# Check for all of the 15 brands -> change one, keep all the others fixed
computeShares_BC_BLP <- function(prices, beta, design, pr = 1) {
  fullDesign <- cbind(prices,design) ###put prices to the last position here
  probabilities_BC_BLP_log_cpp(beta,fullDesign,pr)
}
Optimal_prices[,"Equilibrium Shares"] <- as.vector(computeShares_BC_BLP(Optimal_prices[,"Equilibrium Price"],
                                                                        beta_BC_LLMns,designBase,pr=1))
round(Optimal_prices,2)

# store results in data frame
res_matrix <- data.frame(
  product = c(products[1:(nplayers-1)],"Outside"),
  equi_price = as.numeric(c(p_Markup_Xi_FixedPoint_BC_BLP$FixedPoint,0)),
  equi_share = as.vector(computeShares_BC_BLP(Optimal_prices[,"Equilibrium Price"],
                                              beta_BC_LLMns,designBase,pr=1)))

# create list to store different competitive scenarios
scenarios <- list("full_comp" = res_matrix)

###########################################################################
# Ownership according to brand names (ie all Amstel go to the same company)
###########################################################################

# ingredients
nplayers = 15 + 1 ### 15 inside + outside
min_p = 0.22
prices_init = c(rep(min_p,nplayers - 1),0)
designBase = rbind(diag(nplayers-1),rep(0,nplayers-1))
Xdata = cbind(prices_init,designBase); colnames(Xdata)[1] = "PRICE" 
Ownership = array(0,dim=c(nplayers,nplayers))
Ownership[1:(nplayers-1),1:(nplayers-1)] = diag((nplayers-1))

# define which brands compete on price, i.e. are owned by different companies
# Amstel (five brands)
for (i in 1:5){
  for (k in 1:5){
    Ownership[i,k] <- 1
  }
}
# Estrella (two brands)
for (i in 7:8){
  for (k in 7:8){
    Ownership[i,k] <- 1
  }
}
# Mahou (two brands)
for (i in 10:11){
  for (k in 10:11){
    Ownership[i,k] <- 1
  }
}

inside_row=which(rowSums(Ownership) != 0, arr.ind = TRUE)
p0=Xdata[inside_row,"PRICE"]
costBase = as.vector(prices_init*0.9)
MC = costBase[inside_row]

### Run Fixed-Point algorithm with Xi-markup equation (reliable and fast: See Table 3 in paper for comparison)
p_Markup_Xi_FixedPoint_BC_BLP = FixedPoint(Function = function(price_vec) FixedPoint_BLP_Xi(price_vec,MC=MC,
                                           ownership=Ownership,Xdata=Xdata,beta_draws=beta_BC_LLMns,pr=1), 
                                           Inputs = p0, MaxIter = 10000, ConvergenceMetricThreshold = 1e-10, Method = "Anderson")

p_Markup_Xi_FixedPoint_BC_BLP$FixedPoint

Optimal_prices <- array(0,dim=c(nplayers,2))
rownames(Optimal_prices) = c(products[1:(nplayers-1)],"Outside")
colnames(Optimal_prices) = c("Equilibrium Price","Equilibrium Shares")
# Save equi-prices
Optimal_prices[,"Equilibrium Price"] <-c(p_Markup_Xi_FixedPoint_BC_BLP$FixedPoint,0)
computeShares_BC_BLP <- function(prices, beta, design, pr = 1) {
  fullDesign <- cbind(prices,design) ###put prices to the last position here
  probabilities_BC_BLP_log_cpp(beta,fullDesign,pr)
}
Optimal_prices[,"Equilibrium Shares"] <- as.vector(computeShares_BC_BLP(Optimal_prices[,"Equilibrium Price"],
                                                                        beta_BC_LLMns,designBase,pr=1))

round(Optimal_prices,2)

# store results in data frame
res_matrix <- data.frame(
  product = c(products[1:(nplayers-1)],"Outside"),
  equi_price = as.numeric(c(p_Markup_Xi_FixedPoint_BC_BLP$FixedPoint,0)),
  equi_share = as.vector(computeShares_BC_BLP(Optimal_prices[,"Equilibrium Price"],
                                       beta_BC_LLMns,designBase,pr=1)))


# add to results list
scenarios[["brand_comp"]] <- res_matrix

scenarios[["brand_comp"]] %>% 
  ggplot(aes(x = product, y = equi_price)) +
    geom_bar(stat = "identity")+
    #geom_bar(aes(y = equi_share), stat = "identity")+
    theme(axis.text.x = element_text(angle = 90))

# display results of full competition vs. 
cbind(scenarios[["full_comp"]], scenarios[["brand_comp"]])
