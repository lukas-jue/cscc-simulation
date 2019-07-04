// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace arma; // use the Armadillo library for matrix computations
using namespace Rcpp;

//Market Share Computations (simple, pa=1) with BC using BLP-type of indirect utility function
//FIRST entry in constrained vector corresponds to BC parameter
//pr-th row of design contains prices
//[[Rcpp::export]]
mat probabilities_BC_BLP_log_cpp(mat const& beta, mat const& design_joint, int const& pr){
  int pr_cpp = pr - 1; //c++ indexing
  int draws = beta.n_rows;
  int nplay = design_joint.n_rows;
  int nvar = beta.n_cols;
  vec Xbeta = zeros(nplay);
  mat X_trans = zeros(nplay,nvar);
  vec budget_copies = zeros(nplay);
  vec Xbeta_stab = zeros(nplay);
  vec max_temp = zeros(nplay);
  vec numerator = zeros(nplay);
  vec denominator = zeros(nplay);
  mat market_share_draws = zeros(draws,nplay);
  vec zero_share = zeros(nplay);
  zero_share(nplay-1) = 1; // everything goes to the outside good now  

  //Compute stabilized choice probs now...
  for(int r = 0; r<draws; r++){
    double budget = beta(r,0); //update budget for r-th draw
    budget_copies = rep(budget,nplay);
    X_trans = design_joint;
    X_trans(span::all,pr_cpp) = vectorise(log(budget_copies-design_joint(span::all,pr_cpp))); 
    //kick out alternatives outside budget that will be NaN elements in Xbeta
    for(int n=0; n<nplay; n++){
      if (design_joint(n,pr_cpp) <= budget) {
        Xbeta(n) = sum(X_trans(n,span::all) % beta(r,span(1,nvar-1))); // ??check this against *-operator??
      } else {
        Xbeta(n) = -1*datum::inf; //set Xbeta to -Inf if price is larger than budget
      }
    }
    if(max(Xbeta)==-1*datum::inf){
      market_share_draws(r,span()) = trans(zero_share);
    } else {
      max_temp.fill(max(Xbeta)); // ensure that max is not -Inf
      Xbeta_stab = Xbeta - max_temp;
      numerator = Xbeta_stab;
      denominator.fill(log(sum(exp(numerator))));
      market_share_draws(r,span()) = trans(exp((numerator - denominator)));
    }
  }
  vec exp_ms = arma::conv_to< vec >::from(mean(market_share_draws)); //computes mean for each column
  return exp_ms;
}

//Market Share Computations (simple, pa=1) with BC using BLP-type of indirect utility function
//FIRST entry in constrained vector corresponds to BC parameter
//pr-th row of design contains prices
//[[Rcpp::export]]
mat probabilities_BC_BLP_draws_log_cpp(mat const& beta, mat const& design_joint, int const& pr){
  int pr_cpp = pr - 1; //c++ indexing
  int draws = beta.n_rows;
  int nplay = design_joint.n_rows;
  int nvar = beta.n_cols;
  vec Xbeta = zeros(nplay);
  mat X_trans = zeros(nplay,nvar);
  vec budget_copies = zeros(nplay);
  vec Xbeta_stab = zeros(nplay);
  vec max_temp = zeros(nplay);
  vec numerator = zeros(nplay);
  vec denominator = zeros(nplay);
  mat market_share_draws = zeros(draws,nplay);
  vec zero_share = zeros(nplay);
  zero_share(nplay-1) = 1; // everything goes to the outside good now  
  
  //Compute stabilized choice probs now...
  for(int r = 0; r<draws; r++){
    double budget = beta(r,0); //update budget for r-th draw
    budget_copies = rep(budget,nplay);
    X_trans = design_joint;
    X_trans(span::all,pr_cpp) = vectorise(log(budget_copies-design_joint(span::all,pr_cpp))); 
    //kick out alternatives outside budget that will be NaN elements in Xbeta
    for(int n=0; n<nplay; n++){
      if (design_joint(n,pr_cpp) <= budget) {
        Xbeta(n) = sum(X_trans(n,span::all) % beta(r,span(1,nvar-1))); // ??check this against *-operator??
      } else {
        Xbeta(n) = -1*datum::inf; //set Xbeta to -Inf if price is larger than budget
      }
    }
    if(max(Xbeta)==-1*datum::inf){
      market_share_draws(r,span()) = trans(zero_share);
    } else {
      max_temp.fill(max(Xbeta)); // ensure that max is not -Inf
      Xbeta_stab = Xbeta - max_temp;
      numerator = Xbeta_stab;
      denominator.fill(log(sum(exp(numerator))));
      market_share_draws(r,span()) = trans(exp((numerator - denominator)));
    }
  }
  return market_share_draws;
}



//Compute ingredients used for FOC over draws 
//[[Rcpp::export]]
List ingredients_BC_BLP_draws_log_cpp(mat const& beta, mat const& design_joint, int const& pr){
  
  int pr_cpp = pr - 1; //c++ indexing
  int draws = beta.n_rows;
  int nplay = design_joint.n_rows;
  int nvar = beta.n_cols;
  vec Xbeta = zeros(nplay);
  mat X_trans = zeros(nplay,nvar);
  vec budget_copies = zeros(nplay);
  vec Xbeta_stab = zeros(nplay);
  vec max_temp = zeros(nplay);
  vec numerator = zeros(nplay);
  vec denominator = zeros(nplay);
  vec probs_at_draw = zeros(nplay);
  mat market_share_draws = zeros(draws,nplay);
  vec zero_share = zeros(nplay);
  zero_share(nplay-1) = 1; // everything goes to the outside good now  
  vec Dwj_pj_at_draw = zeros(nplay);
  cube Lambda_at_draw = zeros(nplay, nplay, draws);
  cube Gamma_at_draw = zeros(nplay, nplay, draws);
  mat Help = zeros(nplay, nplay);
  mat Help_Gamma = zeros(nplay, nplay);
  
  //Compute stabilized choice probs now...
  for(int r = 0; r<draws; r++){
    double budget = beta(r,0); //update budget for r-th draw
    budget_copies = rep(budget,nplay);
    double price_coeff = beta(r,(pr_cpp+1));
    X_trans = design_joint;
    X_trans(span::all,pr_cpp) = vectorise(log(budget_copies-design_joint(span::all,pr_cpp))); 
    //kick out alternatives outside budget that will be NaN elements in Xbeta
    for(int n=0; n<nplay; n++){
      if (design_joint(n,pr_cpp) < budget) {
        Xbeta(n) = sum(X_trans(n,span::all) % beta(r,span(1,nvar-1))); // ??check this against *-operator??
        Dwj_pj_at_draw(n) = -1*(price_coeff/(budget-design_joint(n,pr_cpp)));
      } else if (design_joint(n,pr_cpp) == budget){
        Xbeta(n) = sum(X_trans(n,span::all) % beta(r,span(1,nvar-1))); // ??check this against *-operator??
        Dwj_pj_at_draw(n) = 0;
      } else {
        Xbeta(n) = -1*datum::inf; //set Xbeta to -Inf if price is larger than budget
        Dwj_pj_at_draw(n) = 0;
      }
    }
    if(max(Xbeta)==-1*datum::inf){
      probs_at_draw = zero_share;
      market_share_draws(r,span()) = trans(probs_at_draw);
    } else {
      max_temp.fill(max(Xbeta)); // ensure that max is not -Inf
      Xbeta_stab = Xbeta - max_temp;
      numerator = Xbeta_stab;
      denominator.fill(log(sum(exp(numerator))));
      probs_at_draw = exp((numerator - denominator));
      market_share_draws(r,span()) = trans(probs_at_draw);
    }
    // fill Lambda & Gamma based on computed choice probs
    for(int n=0; n<nplay; n++){
      Lambda_at_draw(n,n,r) = Dwj_pj_at_draw(n)*probs_at_draw(n);
      Help(n,span::all) = trans(Dwj_pj_at_draw);
    }
    Help_Gamma = probs_at_draw*trans(probs_at_draw);
    Gamma_at_draw.slice(r) = Help_Gamma % Help;
  }
  
  List out = List::create(Named("Logitdraws") = market_share_draws,
                          Named("Lambdabig") = Lambda_at_draw,
                          Named("Gammabig") = Gamma_at_draw);
  return out;
}















//Market Share Computations (simple, pa=1) with BC
//FIRST entry in constrained vector corresponds to BC parameter
//pr-th row of design contains prices
//[[Rcpp::export]]
mat probabilities_BC_log_CORR_cpp(mat const& beta, mat const& design_joint, int const& pr){
  int pr_cpp = pr - 1; //c++ indexing
  int draws = beta.n_rows;
  int nplay = design_joint.n_rows;
  int nvar = beta.n_cols;
  mat Xbeta = zeros(draws,nplay);
  vec zero_share = zeros(nplay);
  zero_share(nplay-1) = 1; // everything goes to the outside good now  
  //Xbeta
  mat beta_red = beta(span::all,span(1,nvar-1)); //Compute Xbeta
  Xbeta = beta_red * trans(design_joint);
  mat market_share_draws = zeros(draws,nplay);
  //Compute stabilized choice probs now...
  for(int r = 0; r<draws; r++){
    double budget = beta(r,0); //update budget for r-th draw
    vec Xbeta_temp = zeros(nplay); 
    Xbeta_temp = arma::conv_to< vec >::from(Xbeta.row(r)); //Initialize here to kick out alts not considered
    if (sum(exp(Xbeta_temp(span(0,nplay-2))))==0) { // if xbeta of inside brands --> -infty
      market_share_draws(r,span()) = trans(zero_share);
    } else {
      vec constraint_ind = ones(nplay);
      for(int n=0; n<nplay; n++){
        if (design_joint(n,pr_cpp) > budget) {
          constraint_ind(n) = 0; //set indicator to zero if price is larger than budget
          Xbeta_temp(n) = -1*datum::inf;
        }
      }
      if (sum(constraint_ind)==0) {
        market_share_draws(r,span()) = trans(zero_share);
      } else {
        vec Xbeta_temp_stab = zeros(nplay);
        vec max_temp = zeros(nplay);
        max_temp.fill(max(Xbeta_temp));
        Xbeta_temp_stab = Xbeta_temp - max_temp;
        vec numerator = log(exp(Xbeta_temp_stab)%constraint_ind);
        vec denominator = zeros(nplay);
        denominator.fill(log(sum(exp(numerator))));
        market_share_draws(r,span()) = trans(exp((numerator - denominator)));
      }
    }
  }
  vec exp_ms = arma::conv_to< vec >::from(mean(market_share_draws)); //computes mean for each column
  return exp_ms;
}


//Market Share Computations (simple, pa=1) with BC
//FIRST entry in constrained vector corresponds to BC parameter
//pr-th row of design contains prices
//[[Rcpp::export]]
mat probabilities_BC_log_NoSTAB_cpp(mat const& beta, mat const& design_joint, int const& pr){
  int pr_cpp = pr - 1; //c++ indexing
  int draws = beta.n_rows;
  int nplay = design_joint.n_rows;
  int nvar = beta.n_cols;
  mat Xbeta = zeros(draws,nplay);
  //Xbeta
  mat beta_red = beta(span::all,span(1,nvar-1)); //Compute Xbeta
  Xbeta = beta_red * trans(design_joint);
  mat market_share_draws = zeros(draws,nplay);
  //Compute stabilized choice probs now...
  for(int r = 0; r<draws; r++){
    double budget = beta(r,0); //update budget for r-th draw
    vec constraint_ind = ones(nplay);
    for(int n=0; n<nplay; n++){
      if (design_joint(n,pr_cpp) > budget) {
        constraint_ind(n) = 0; //set indicator to zero if price is larger than budget
      }
    }
    // vec Xbeta_temp = zeros(nplay); 
    // Xbeta_temp = arma::conv_to< vec >::from(Xbeta.row(r))%constraint_ind; //Do NOT choose maximum xbeta not affordable
    // vec Xbeta_temp_stab = zeros(nplay);
    // vec max_temp = zeros(nplay);
    // max_temp.fill(max(Xbeta_temp));
    // Xbeta_temp_stab = Xbeta_temp - max_temp;
    
    vec Xbeta_temp = zeros(nplay);
    Xbeta_temp = arma::conv_to< vec >::from(Xbeta.row(r));
    vec numerator = log(exp(Xbeta_temp)%constraint_ind);
    vec denominator = zeros(nplay);
    denominator.fill(log(sum(exp(numerator))));
    market_share_draws(r,span()) = trans(exp((numerator - denominator)));
  }
  vec exp_ms = arma::conv_to< vec >::from(mean(market_share_draws)); //computes mean for each column
  return exp_ms;
}


