################################################################
### Functions for Computing Nash-equilibrium in the budget model
### As in Morrow & Skerlos 2011#################################
### BLP-type of indirect utility function#######################
### More efficient ways#########################################
################################################################

###requires "Speed++_MS_BC_BLP_Efficient.cpp" file 

### Write BLP-markup equation first
FOC_BLP_Markup = function(price_vec,MC,ownership,Xdata,beta_draws,pr) {
  
  # Note: "price_vec" should have same length as "MC". I.e., only the prices that are optimized
  
  X=Xdata

  position=pr

  FOC_result= logit_FOC_NE_BLP_markup(beta_draws,X,price_vec,MC,position,ownership)

  #Results:
  return(FOC_result)
}


logit_FOC_NE_BLP_markup = function(beta_draws,X,price_vec,MC,position,ownership) {
  
  # Indicate inside goods
  inside_row = which(rowSums(ownership) != 0, arr.ind = TRUE)
  # update prices, only inside goods
  X[,position][inside_row] = price_vec
  
  # get ingredients from loop over MCMC draws 
  out_ingred = ingredients_BC_BLP_draws_log_cpp(beta_draws, X, position)
  # Compute average choice probs over draws
  logit_probas.mean = apply(out_ingred$Logitdraws,2,mean)
  # compute jacobian using Proposition 2.3, p.12 (Morrow & Skerlos book) 
  Lambda = apply(out_ingred$Lambdabig, c(1,2), mean) 
  Gamma_tilde = apply(out_ingred$Gammabig, c(1,2), mean)*ownership
  logit_jacobian.mean = Lambda - Gamma_tilde
  
  # strip outside good from matrices
  strip_col = which(colSums(ownership) == 0, arr.ind = TRUE)
  strip_row = which(rowSums(ownership) == 0, arr.ind = TRUE)
  ownership_strip = ownership[-strip_row,-strip_col]

  # check whether Jacobian is non-singular at prices
  if(f(logit_jacobian.mean[-strip_row,-strip_col])){
    FOC = price_vec-MC-(-solve(ownership_strip*logit_jacobian.mean[-strip_row,-strip_col])%*%logit_probas.mean[-strip_row]) ### Similar as Eq. (6) in Morrow and Skerlos, 2011
  }else{
    FOC = -1*1e50
  }
  #Results:
  logit_FOC_NE_BLP_markup=FOC
  return(logit_FOC_NE_BLP_markup)
  
}

### Write another Fixed-Point Equation (more reliable and fast equilibrium computations), see Morrow and Skerlos, 2011 part 2.5
FOC_BLP_Xi = function(price_vec,MC,ownership,Xdata,beta_draws,pr) {
  
  # Note: "price_vec" should have same length as "MC". I.e., only the prices that are optimized
  
  X=Xdata
  position=pr
  
  FOC_result=logit_FOC_NE_Xi_markup(beta_draws,X,price_vec,MC,position,ownership)
  
  #Results:
  return(FOC_result)
}


logit_FOC_NE_Xi_markup = function(beta_draws,X,price_vec,MC,position,ownership) {
  
  #Indicate inside goods
  inside_row=which(rowSums(ownership) != 0, arr.ind = TRUE)
  #update prices, only inside goods
  X[,position][inside_row]=price_vec
  
  # get ingredients from loop over MCMC draws 
  out_ingred = ingredients_BC_BLP_draws_log_cpp(beta_draws, X, position)
  # Compute average choice probs over draws
  logit_probas.mean = apply(out_ingred$Logitdraws,2,mean)
  # Compute Xi as proposed in Equ. (8), Morrow and Skerlos, 2011
  Lambda = apply(out_ingred$Lambdabig, c(1,2), mean) 
  Gamma_tilde = apply(out_ingred$Gammabig, c(1,2), mean)*ownership
  strip_col = which(colSums(ownership) == 0, arr.ind = TRUE)
  strip_row = which(rowSums(ownership) == 0, arr.ind = TRUE)
  Lambda_strip = Lambda[-strip_row,-strip_col]
  Gamma_tilde_strip = Gamma_tilde[-strip_row,-strip_col]
  # check whether lambda is non-singular
  if( f(Lambda_strip) ){
    Xi = solve(Lambda_strip)%*%t(Gamma_tilde_strip)%*%(price_vec-MC)-solve(Lambda_strip)%*%logit_probas.mean[-strip_row]
    # compute FOC
    FOC = price_vec-MC-Xi ### Similar as Eq. (8) in Morrow and Skerlos, 2011
  }else{
    FOC = -1*1e50
  }
  #Results:
  logit_FOC_NE_Xi_markup=FOC
  return(logit_FOC_NE_Xi_markup)
}

###################################
###Re-write Xi Fixed-Point equation
###################################
FixedPoint_BLP_Xi = function(price_vec,MC,ownership,Xdata,beta_draws,pr) {
  
  # Note: "price_vec" must have same length as "MC". I.e., only the prices that are optimized
  
  X=Xdata
  position=pr
  
  Output=logit_FixedPoint_Xi_markup(beta_draws,X,price_vec,MC,position,ownership)
  
  #Results:
  return(Output)
}

logit_FixedPoint_Xi_markup = function(beta_draws,X,price_vec,MC,position,ownership) {
  
  #Indicate inside goods
  inside_row=which(rowSums(ownership) != 0, arr.ind = TRUE)
  #update prices, only inside goods
  X[,position][inside_row]=price_vec
  
  # get ingredients from loop over MCMC draws 
  out_ingred = ingredients_BC_BLP_draws_log_cpp(beta_draws, X, position)
  # Compute average choice probs over draws
  logit_probas.mean = apply(out_ingred$Logitdraws,2,mean)
  # Compute Xi as proposed in Equ. (8), Morrow and Skerlos, 2011
  Lambda = apply(out_ingred$Lambdabig, c(1,2), mean) 
  Gamma_tilde = apply(out_ingred$Gammabig, c(1,2), mean)*ownership
  strip_col = which(colSums(ownership) == 0, arr.ind = TRUE)
  strip_row = which(rowSums(ownership) == 0, arr.ind = TRUE)
  Lambda_strip = Lambda[-strip_row,-strip_col]
  Gamma_tilde_strip = Gamma_tilde[-strip_row,-strip_col]
  # check whether lambda is non-singular
  if( f(Lambda_strip) ){
    Xi = solve(Lambda_strip)%*%t(Gamma_tilde_strip)%*%(price_vec-MC)-solve(Lambda_strip)%*%logit_probas.mean[-strip_row]
    # compute FOC
    y_out = MC+Xi ### Similar as Eq. (8) in Morrow and Skerlos, 2011
  }else{
    y_out = -1*1e50
  }
  return(y_out)
}

#############################################
###INFER MARGINAL COSTS BASED ON BLP-MARKUP##
#############################################

MC_BLP_logit_NE = function(beta_draws,X,position,ownership) {
  
  #update price
  price_vec = X[,position]
  
  # get ingredients from loop over MCMC draws 
  out_ingred = ingredients_BC_BLP_draws_log_cpp(beta_draws, X, position)
  # Compute average choice probs over draws
  logit_probas.mean = apply(out_ingred$Logitdraws,2,mean)
  # compute jacobian using Proposition 2.3, p.12 (Morrow & Skerlos book) 
  Lambda = apply(out_ingred$Lambdabig, c(1,2), mean) 
  Gamma_tilde = apply(out_ingred$Gammabig, c(1,2), mean)*ownership
  logit_jacobian.mean = Lambda - Gamma_tilde
  
  # strip outside good from matrices
  strip_col = which(colSums(ownership) == 0, arr.ind = TRUE)
  strip_row = which(rowSums(ownership) == 0, arr.ind = TRUE)
  ownership_strip = ownership[-strip_row,-strip_col]
  logit_jacobian.mean = logit_jacobian.mean[-strip_row,-strip_col]
  logit_probas.mean =  logit_probas.mean[-strip_row]
  price_vec = price_vec[-strip_row]
  
  # check whether Jacobian is non-singular at prices
  if(f(logit_jacobian.mean)){
    # Infer MC as in Morrow & Skerlos 2011, Equ. (6)
    MC = price_vec + t(solve(ownership_strip*logit_jacobian.mean))%*%logit_probas.mean ### Need transpose here???
  }else{
    MC = -1*1e50
  }
  
  Margin = price_vec-MC 
  #Results:
  MC_logit_NE <- NULL
  MC_logit_NE$MC = MC
  MC_logit_NE$Margin = Margin

  return(MC_logit_NE)
}


### function to check whether matrix is non-singular
f <- function(m) class(try(solve(m),silent=T))=="matrix"


