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
library(ggpmisc)
library(tikzDevice)
library(plyr)
library(latex2exp)
library(FixedPoint)
library(dplyr)
library(reshape2)
library(tidyr)
library(ggcorrplot)

load("100k_MCMC_corrected.RData")

###Increase memory capacities
memory.limit(size=100000000)

load("Estimation_Data_Beer_20170423.Rdata")
products = c("Amstel Extra Lata 37,5 cl","Amstel Extra Lata 33 cl","Amstel Lata 37,5 cl","Amstel Lata 33 cl","Amstel Clasica Lata 33 cl",
             "Cruzcampo Lata 33 cl","Estrella Damm Lata 33 cl","Estrella Galicia Lata 33 cl","Heineken Lata 33 cl","Mahou 5 Estrellas Lata 33 cl",
             "Mahou Clasica Lata 33 cl","San Miguel Lata 33 cl","Voll Damm Lata 33 cl","Steinburg (Marca Blanca Mercadona) Lata 33 cl",
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

data.frame(loglikelihood = out_BC$loglike) %>% 
  mutate(index = row_number()) %>% 
    ggplot(aes(x = index, y = loglikelihood)) +
      geom_line(alpha = 0.7) +
      geom_smooth(method = "lm", se = FALSE) +
      stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~~~")), 
                   parse=TRUE,label.x.npc = "right", label.y.npc = 0.1,
                   output.type = "expression") +
  theme_bw()
      

burnin = 3100
R = dim(betastar_HB_BC)[3]

betastar_HB_BC = betastar_HB_BC[,,(burnin+1):R]
compdraw_HB = compdraw_HB[(burnin+1):R]
probdraw_HB = probdraw_HB[(burnin+1):R]
rejection = rejection[(burnin+1):R,]
loglike_BC = loglike_BC[(burnin+1):R]
plot(loglike_BC, type="l")

data.frame(loglikelihood = loglike_BC) %>% 
  mutate(index = row_number()) %>% 
  ggplot(aes(x = index, y = loglikelihood)) +
  geom_line(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~~~")), 
               parse=TRUE,label.x.npc = "right", label.y.npc = 0.01,
               output.type = "expression") +
  theme_bw()


lmdata <- data.frame(loglikelihood = loglike_BC) %>% 
  mutate(index = row_number())
lm(loglike_BC ~ index, data = lmdata) %>% 
          summary()

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


###########################################################################
# Only two Amstel, two Estrella and one Heineken (Ownership according to brands)
###########################################################################

# ingredients
nplayers = 5 + 1 
min_p = 0.22
prices_init = c(rep(min_p,nplayers - 1),0)
designBase = rbind(diag(nplayers-1),rep(0,nplayers-1))
Xdata = cbind(prices_init,designBase); colnames(Xdata)[1] = "PRICE" 
Ownership = array(0,dim=c(nplayers,nplayers))
Ownership[1:(nplayers-1),1:(nplayers-1)] = diag((nplayers-1))

# define which brands compete on price, i.e. are owned by different companies
# Amstel (two brands)
for (i in 1:2){
  for (k in 1:2){
    Ownership[i,k] <- 1
  }
}
# Estrella (two brands)
for (i in 3:4){
  for (k in 3:4){
    Ownership[i,k] <- 1
  }
}


inside_row=which(rowSums(Ownership) != 0, arr.ind = TRUE)
p0=Xdata[inside_row,"PRICE"]
costBase = as.vector(prices_init*0.9)
MC = costBase[inside_row]

### Run Fixed-Point algorithm with Xi-markup equation (reliable and fast: See Table 3 in paper for comparison)
p_Markup_Xi_FixedPoint_BC_BLP = FixedPoint(Function = function(price_vec) FixedPoint_BLP_Xi(price_vec,MC=MC,
                                                                                            ownership=Ownership,Xdata=Xdata,beta_draws=beta_BC_LLMns[,(c(1,2,4,7,9,10,11))],pr=1), 
                                           Inputs = p0, MaxIter = 10000, ConvergenceMetricThreshold = 1e-10, Method = "Anderson")

p_Markup_Xi_FixedPoint_BC_BLP$FixedPoint

Optimal_prices <- array(0,dim=c(nplayers,2))
rownames(Optimal_prices) = c(products[c(2,5,7,8,9)],"Outside")
colnames(Optimal_prices) = c("Equilibrium Price","Equilibrium Shares")
# Save equi-prices
Optimal_prices[,"Equilibrium Price"] <-c(p_Markup_Xi_FixedPoint_BC_BLP$FixedPoint,0)
computeShares_BC_BLP <- function(prices, beta, design, pr = 1) {
  fullDesign <- cbind(prices,design) ###put prices to the last position here
  probabilities_BC_BLP_log_cpp(beta,fullDesign,pr)
}
Optimal_prices[,"Equilibrium Shares"] <- as.vector(computeShares_BC_BLP(Optimal_prices[,"Equilibrium Price"],
                                                                        beta_BC_LLMns[,c(1,2,4,7,9,10,11)],designBase,pr=1))

round(Optimal_prices,2)

# store results in data frame
res_matrix <- data.frame(
  product = rownames(Optimal_prices),
  equi_price = as.numeric(c(p_Markup_Xi_FixedPoint_BC_BLP$FixedPoint,0)),
  equi_share = as.vector(computeShares_BC_BLP(Optimal_prices[,"Equilibrium Price"],
                                              beta_BC_LLMns[,c(1,2,4,7,9,10,11)],designBase,pr=1)))


# add to results list
scenarios[["five_full_comp"]] <- res_matrix

###########################################################################
# Only two Amstel, two Estrella and one Heineken (Merger Amstel and Heineken, bc highly correlated betas)
###########################################################################

# ingredients
nplayers = 5 + 1 
min_p = 0.22
prices_init = c(rep(min_p,nplayers - 1),0)
designBase = rbind(diag(nplayers-1),rep(0,nplayers-1))
Xdata = cbind(prices_init,designBase); colnames(Xdata)[1] = "PRICE" 
Ownership = array(0,dim=c(nplayers,nplayers))
Ownership[1:(nplayers-1),1:(nplayers-1)] = diag((nplayers-1))

# define which brands compete on price, i.e. are owned by different companies
# Amstel (two brands)
for (i in c(1:2,5)){
  for (k in c(1:2,5)){
    Ownership[i,k] <- 1
  }
}
# Estrella (two brands)
for (i in 3:4){
  for (k in 3:4){
    Ownership[i,k] <- 1
  }
}

inside_row=which(rowSums(Ownership) != 0, arr.ind = TRUE)
p0=Xdata[inside_row,"PRICE"]
costBase = as.vector(prices_init*0.9)
MC = costBase[inside_row]

### Run Fixed-Point algorithm with Xi-markup equation (reliable and fast: See Table 3 in paper for comparison)
p_Markup_Xi_FixedPoint_BC_BLP = FixedPoint(Function = function(price_vec) FixedPoint_BLP_Xi(price_vec,MC=MC,
                                                                                            ownership=Ownership,Xdata=Xdata,beta_draws=beta_BC_LLMns[,(c(1,2,4,7,9,10,11))],pr=1), 
                                           Inputs = p0, MaxIter = 10000, ConvergenceMetricThreshold = 1e-10, Method = "Anderson")

p_Markup_Xi_FixedPoint_BC_BLP$FixedPoint

Optimal_prices <- array(0,dim=c(nplayers,2))
rownames(Optimal_prices) = c(products[c(2,5,7,8,9)],"Outside")
colnames(Optimal_prices) = c("Equilibrium Price","Equilibrium Shares")
# Save equi-prices
Optimal_prices[,"Equilibrium Price"] <-c(p_Markup_Xi_FixedPoint_BC_BLP$FixedPoint,0)
computeShares_BC_BLP <- function(prices, beta, design, pr = 1) {
  fullDesign <- cbind(prices,design) ###put prices to the last position here
  probabilities_BC_BLP_log_cpp(beta,fullDesign,pr)
}
Optimal_prices[,"Equilibrium Shares"] <- as.vector(computeShares_BC_BLP(Optimal_prices[,"Equilibrium Price"],
                                                                        beta_BC_LLMns[,c(1,2,4,7,9,10,11)],designBase,pr=1))

round(Optimal_prices,2)

# store results in data frame
res_matrix <- data.frame(
  product = rownames(Optimal_prices),
  equi_price = as.numeric(c(p_Markup_Xi_FixedPoint_BC_BLP$FixedPoint,0)),
  equi_share = as.vector(computeShares_BC_BLP(Optimal_prices[,"Equilibrium Price"],
                                              beta_BC_LLMns[,c(1,2,4,7,9,10,11)],designBase,pr=1)))


# add to results list
scenarios[["five_merge_comp"]] <- res_matrix

###########################################################################
# Heineken and Estrella Merger
###########################################################################

# ingredients
nplayers = 5 + 1 
min_p = 0.22
prices_init = c(rep(min_p,nplayers - 1),0)
designBase = rbind(diag(nplayers-1),rep(0,nplayers-1))
Xdata = cbind(prices_init,designBase); colnames(Xdata)[1] = "PRICE" 
Ownership = array(0,dim=c(nplayers,nplayers))
Ownership[1:(nplayers-1),1:(nplayers-1)] = diag((nplayers-1))

# define which brands compete on price, i.e. are owned by different companies
# Amstel (two brands)
for (i in c(1:2)){
  for (k in c(1:2)){
    Ownership[i,k] <- 1
  }
}
# Estrella (two brands)
for (i in c(3:4,5)){
  for (k in c(3:4,5)){
    Ownership[i,k] <- 1
  }
}

inside_row=which(rowSums(Ownership) != 0, arr.ind = TRUE)
p0=Xdata[inside_row,"PRICE"]
costBase = as.vector(prices_init*0.9)
MC = costBase[inside_row]

### Run Fixed-Point algorithm with Xi-markup equation (reliable and fast: See Table 3 in paper for comparison)
p_Markup_Xi_FixedPoint_BC_BLP = FixedPoint(Function = function(price_vec) FixedPoint_BLP_Xi(price_vec,MC=MC,
                                                                                            ownership=Ownership,Xdata=Xdata,beta_draws=beta_BC_LLMns[,(c(1,2,4,7,9,10,11))],pr=1), 
                                           Inputs = p0, MaxIter = 10000, ConvergenceMetricThreshold = 1e-10, Method = "Anderson")

p_Markup_Xi_FixedPoint_BC_BLP$FixedPoint

Optimal_prices <- array(0,dim=c(nplayers,2))
rownames(Optimal_prices) = c(products[c(2,5,7,8,9)],"Outside")
colnames(Optimal_prices) = c("Equilibrium Price","Equilibrium Shares")
# Save equi-prices
Optimal_prices[,"Equilibrium Price"] <-c(p_Markup_Xi_FixedPoint_BC_BLP$FixedPoint,0)
computeShares_BC_BLP <- function(prices, beta, design, pr = 1) {
  fullDesign <- cbind(prices,design) ###put prices to the last position here
  probabilities_BC_BLP_log_cpp(beta,fullDesign,pr)
}
Optimal_prices[,"Equilibrium Shares"] <- as.vector(computeShares_BC_BLP(Optimal_prices[,"Equilibrium Price"],
                                                                        beta_BC_LLMns[,c(1,2,4,7,9,10,11)],designBase,pr=1))

round(Optimal_prices,2)

# store results in data frame
res_matrix <- data.frame(
  product = rownames(Optimal_prices),
  equi_price = as.numeric(c(p_Markup_Xi_FixedPoint_BC_BLP$FixedPoint,0)),
  equi_share = as.vector(computeShares_BC_BLP(Optimal_prices[,"Equilibrium Price"],
                                              beta_BC_LLMns[,c(1,2,4,7,9,10,11)],designBase,pr=1)))


# add to results list
scenarios[["five_merge_comp_estr"]] <- res_matrix

###########################################################################
# Only two Amstel, two Estrella and one Heineken (Full Competition)
###########################################################################

# ingredients
nplayers = 5 + 1 
min_p = 0.22
prices_init = c(rep(min_p,nplayers - 1),0)
designBase = rbind(diag(nplayers-1),rep(0,nplayers-1))
Xdata = cbind(prices_init,designBase); colnames(Xdata)[1] = "PRICE" 
Ownership = array(0,dim=c(nplayers,nplayers))
Ownership[1:(nplayers-1),1:(nplayers-1)] = diag((nplayers-1))

inside_row=which(rowSums(Ownership) != 0, arr.ind = TRUE)
p0=Xdata[inside_row,"PRICE"]
costBase = as.vector(prices_init*0.9)
MC = costBase[inside_row]

### Run Fixed-Point algorithm with Xi-markup equation (reliable and fast: See Table 3 in paper for comparison)
p_Markup_Xi_FixedPoint_BC_BLP = FixedPoint(Function = function(price_vec) FixedPoint_BLP_Xi(price_vec,MC=MC,
                                                                                            ownership=Ownership,Xdata=Xdata,beta_draws=beta_BC_LLMns[,(c(1,2,4,7,9,10,11))],pr=1), 
                                           Inputs = p0, MaxIter = 10000, ConvergenceMetricThreshold = 1e-10, Method = "Anderson")

p_Markup_Xi_FixedPoint_BC_BLP$FixedPoint

Optimal_prices <- array(0,dim=c(nplayers,2))
rownames(Optimal_prices) = c(products[c(2,5,7,8,9)],"Outside")
colnames(Optimal_prices) = c("Equilibrium Price","Equilibrium Shares")
# Save equi-prices
Optimal_prices[,"Equilibrium Price"] <-c(p_Markup_Xi_FixedPoint_BC_BLP$FixedPoint,0)
computeShares_BC_BLP <- function(prices, beta, design, pr = 1) {
  fullDesign <- cbind(prices,design) ###put prices to the last position here
  probabilities_BC_BLP_log_cpp(beta,fullDesign,pr)
}
Optimal_prices[,"Equilibrium Shares"] <- as.vector(computeShares_BC_BLP(Optimal_prices[,"Equilibrium Price"],
                                                                        beta_BC_LLMns[,c(1,2,4,7,9,10,11)],designBase,pr=1))

round(Optimal_prices,2)

# store results in data frame
res_matrix <- data.frame(
  product = rownames(Optimal_prices),
  equi_price = as.numeric(c(p_Markup_Xi_FixedPoint_BC_BLP$FixedPoint,0)),
  equi_share = as.vector(computeShares_BC_BLP(Optimal_prices[,"Equilibrium Price"],
                                              beta_BC_LLMns[,c(1,2,4,7,9,10,11)],designBase,pr=1)))


# add to results list
scenarios[["five_full_comp"]] <- res_matrix

###########################################################################
# Only two Amstel, two Estrella and one Heineken (Monopoly)
###########################################################################

# ingredients
nplayers = 5 + 1 
min_p = 0.22
prices_init = c(rep(min_p,nplayers - 1),0)
designBase = rbind(diag(nplayers-1),rep(0,nplayers-1))
Xdata = cbind(prices_init,designBase); colnames(Xdata)[1] = "PRICE" 
Ownership = array(0,dim=c(nplayers,nplayers))
Ownership[1:(nplayers-1),1:(nplayers-1)] = diag((nplayers-1))

# define which brands compete on price, i.e. are owned by different companies
# All five belong to the same owner
for (i in 1:5){
  for (k in 1:5){
    Ownership[i,k] <- 1
  }
}

inside_row=which(rowSums(Ownership) != 0, arr.ind = TRUE)
p0=Xdata[inside_row,"PRICE"]
costBase = as.vector(prices_init*0.9)
MC = costBase[inside_row]

### Run Fixed-Point algorithm with Xi-markup equation (reliable and fast: See Table 3 in paper for comparison)
p_Markup_Xi_FixedPoint_BC_BLP = FixedPoint(Function = function(price_vec) FixedPoint_BLP_Xi(price_vec,MC=MC,
                                                                                            ownership=Ownership,Xdata=Xdata,beta_draws=beta_BC_LLMns[,(c(1,2,4,7,9,10,11))],pr=1), 
                                           Inputs = p0, MaxIter = 10000, ConvergenceMetricThreshold = 1e-10, Method = "Anderson")

p_Markup_Xi_FixedPoint_BC_BLP$FixedPoint

Optimal_prices <- array(0,dim=c(nplayers,2))
rownames(Optimal_prices) = c(products[c(2,5,7,8,9)],"Outside")
colnames(Optimal_prices) = c("Equilibrium Price","Equilibrium Shares")
# Save equi-prices
Optimal_prices[,"Equilibrium Price"] <-c(p_Markup_Xi_FixedPoint_BC_BLP$FixedPoint,0)
computeShares_BC_BLP <- function(prices, beta, design, pr = 1) {
  fullDesign <- cbind(prices,design) ###put prices to the last position here
  probabilities_BC_BLP_log_cpp(beta,fullDesign,pr)
}
Optimal_prices[,"Equilibrium Shares"] <- as.vector(computeShares_BC_BLP(Optimal_prices[,"Equilibrium Price"],
                                                                        beta_BC_LLMns[,c(1,2,4,7,9,10,11)],designBase,pr=1))

round(Optimal_prices,2)

# store results in data frame
res_matrix <- data.frame(
  product = rownames(Optimal_prices),
  equi_price = as.numeric(c(p_Markup_Xi_FixedPoint_BC_BLP$FixedPoint,0)),
  equi_share = as.vector(computeShares_BC_BLP(Optimal_prices[,"Equilibrium Price"],
                                              beta_BC_LLMns[,c(1,2,4,7,9,10,11)],designBase,pr=1)))


# add to results list
scenarios[["five_monopoly"]] <- res_matrix

########################################################
# Compute Average Prices in the Market, Weighted by Market Share
########################################################

w_a_p_full_comp <- t(scenarios[["full_comp"]]$equi_price) %*% scenarios[["full_comp"]]$equi_share
w_a_p_brand_comp <- t(scenarios[["brand_comp"]]$equi_price) %*% scenarios[["brand_comp"]]$equi_share

w_a_p_five_full_comp <- t(scenarios[["five_full_comp"]]$equi_price) %*% scenarios[["five_full_comp"]]$equi_share
w_a_p_five_brand_comp <- t(scenarios[["five_brand_comp"]]$equi_price) %*% scenarios[["five_brand_comp"]]$equi_share
w_a_p_five_monopoly <- t(scenarios[["five_monopoly"]]$equi_price) %*% scenarios[["five_monopoly"]]$equi_share
w_a_p_five_merge <- t(scenarios[["five_merge_comp"]]$equi_price) %*% scenarios[["five_merge_comp"]]$equi_share
w_a_p_five_merge_estr <- t(scenarios[["five_merge_comp_estr"]]$equi_price) %*% scenarios[["five_merge_comp_estr"]]$equi_share

w_a_p_brand_comp / w_a_p_full_comp
w_a_p_five_brand_comp / w_a_p_five_full_comp

# compute average annual welfare loss for typical German consumer with the merger
welf_loss_avg <- 102/(1/3)*(w_a_p_five_merge - w_a_p_five_brand_comp)
welf_loss_total <- welf_loss_avg * 82000000


########################################################
# Compute contribution margin * market share for all scenarios
########################################################
MC <- rep(0.198, 15)
for (i in 1:length(scenarios)) { #
  # position of outside good
  out_pos <- length(scenarios[[i]]$equi_price)
  
  # compute contribubtion margin times market share for every brand in every scenario in the list
  scenarios[[i]]$CMxMS <- c((scenarios[[i]][-out_pos,"equi_price"] - MC[1:out_pos-1]) * scenarios[[i]][-out_pos,"equi_share"], NA)
  
}

# compute CMxMS for both Amstel and Heineken before and after merger
CMxMS_full_comp <- sum(scenarios$five_full_comp[c(1,2,5),"CMxMS"])
CMxMS_before <- sum(scenarios$five_brand_comp[c(1,2,5),"CMxMS"])
CMxMS_after <- sum(scenarios$five_merge_comp[c(1,2,5),"CMxMS"])
CMxMS_monop <- sum(scenarios$five_monopoly[c(1,2,5),"CMxMS"])

CMxMS_before_estr <- sum(scenarios$five_brand_comp[c(3:5),"CMxMS"])
CMxMS_after_estr <- sum(scenarios$five_merge_comp_estr[c(3:5),"CMxMS"])


CMxMS_after / CMxMS_before
CMxMS_monop / CMxMS_full_comp
CMxMS_monop / CMxMS_after

CMxMS_after_estr / CMxMS_before_estr

########################################################
# Plotting
########################################################

# corrplot of betas for all 15 brands
colnames(beta_BC_LLMns) <- c("Budget", "Price", products)
windows()
cor(beta_BC_LLMns) %>% 
  ggcorrplot(type = "lower", lab = TRUE)

cor(beta_BC_LLMns[,c(1,2,4,7,9,10,11)]) %>% 
  ggcorrplot(type = "lower", lab = TRUE)


# histogram of all betas
beta_BC_LLMns %>% 
    data.frame() %>% 
      gather(key = "beer_brand") %>% 
        ggplot(aes(value)) +
          geom_histogram() +
          facet_wrap(~beer_brand, scales = "free") +
          theme_bw()

# density plot of five brands (all in one)
beta_BC_LLMns[,c(4,7,9,10,11)] %>% 
  data.frame() %>% 
  gather(key = "beer_brand") %>% 
  ggplot(aes(value)) +
  geom_density(aes(fill = beer_brand), position="identity", alpha = 0.3) +
  xlim(c(-10, 20)) +
  ylim(c(0, 0.095)) +
  scale_fill_manual(values = c("#FD0505", "#FF9A9A", "#000CD3", "#00C1EA", "green")) +
  theme_bw() +
  theme(legend.position = c(0.2, 0.8),
        legend.direction = "vertical")

# density plot of five brands (separate)
beta_BC_LLMns[,c(2,4,7,9,10,11)] %>% 
  data.frame() %>% 
  gather(key = "beer_brand") %>% 
  ggplot(aes(value)) +
  geom_density(aes(fill = beer_brand), position="identity") +
  xlim(c(-10, 25)) +
  ylim(c(0, 0.13)) +
  scale_fill_manual(values = c("#FD0505", "#FF9A9A", "#000CD3", "#00C1EA", "green", "grey"), guide = FALSE) +
  facet_wrap(~beer_brand) +
  theme_bw()

# basic plot of one scenario
scenarios[["brand_comp"]] %>% 
  melt() %>% 
    ggplot(aes(x = reorder(product, value) , y = value)) +
      geom_bar(stat = "identity", position = "dodge") +
      coord_flip() +
      theme(axis.text.x = element_text(angle = 90)) +
      facet_wrap(~variable, scales = "free_x") +
      theme_bw()

# uses all competitive situations in the scenarios list and plots them side by side
windows()
scenarios[c("brand_comp", "full_comp")] %>% 
  melt() %>% 
  rename(comp_scenario = L1) %>% 
  ggplot(aes(x = reorder(product, value) , y = value, fill = comp_scenario)) +
  geom_bar(stat = "identity", colour="black", position = "dodge") +
  coord_flip() +
  labs(x = "Beer Brand",
       y = NULL) +
  facet_wrap(~variable, scales = "free_x") +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw()

 # only five brands #c("five_full_comp", "five_brand_comp", "five_monopoly", "five_merge_comp")
scenarios[c("five_full_comp","five_brand_comp", "five_merge_comp", "five_monopoly")] %>% 
  melt() %>% 
  rename(comp_scenario = L1) %>% 
  ggplot(aes(x = reorder(product, value) , y = value, fill = comp_scenario)) +
  geom_bar(stat = "identity", colour="black", position = "dodge") +
  coord_flip() +
  labs(x = "Beer Brand",
       y = NULL) +
  facet_wrap(~variable, scales = "free_x") +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw()

# price-share scatterplot, w/o outside option (either for single scenario or all)
#scenarios[["brand_comp"]] %>% 
scenarios[c("five_full_comp","five_brand_comp", "five_merge_comp", "five_monopoly")] %>% 
bind_rows() %>% 
  filter(product != "Outside") %>% 
    ggplot(aes(x = equi_price, y = equi_share)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    stat_poly_eq(formula = y ~ x, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~~~")), 
                 parse=TRUE,label.x.npc = "right",
                 output.type = "expression") +
    theme_bw()

cor(scenarios[["brand_comp"]]$equi_price, scenarios[["brand_comp"]]$equi_share)

# bivariate denisity plot
hist(beta_BC_LLMns[,c(4,7)])
beta_BC_LLMns[,c(4,7)] %>% 
  data.frame() %>% 
  gather(key = "beer_brand") %>% 
  ggplot(aes(value)) +
  geom_density(aes(fill = beer_brand), position="identity")

data.frame(beta_BC_LLMns[,c(4,7)])%>% 
  ggplot(aes_string(x = "Amstel.Extra.Lata.33.cl", y = "Amstel.Clasica.Lata.33.cl")) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  scale_fill_continuous(low="lavenderblush", high="blue")+
  geom_abline(slope = 1) +
  labs(fill = "density") +
  #theme(legend.title = element_text("density")) +
  theme_bw()
