// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <mvnorm.h>
#include <wishart.h>

#define M_PIl          3.141592653589793238462643383279502884L /* pi */

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

// multinomial
IntegerVector oneMultinomCalt(NumericVector probs) {
  int k = probs.size();
  IntegerVector ans(k);
  rmultinom(1, probs.begin(), k, ans.begin());
  return(ans);
}

// which function
SEXP mywhich(SEXP x, SEXP y) {
  
  //For each supported type, turn it into the 'real' type and
  //perform the operation. We can use TYPEOF to check the type.
  switch(TYPEOF(x)){
  case REALSXP: { 
    Environment base("package:base"); 
    Function f("which");
    NumericVector answer = f(as<NumericVector>(y)(0) == as<NumericVector>(x));
    //                                           ^^^
    return wrap(answer);
  }
  case INTSXP: { 
    Environment base("package:base"); 
    Function f("which");
    IntegerVector answer = f(as<IntegerVector>(y)(0) == as<IntegerVector>(x));
    //                                           ^^^
    return wrap(answer);
  }
  default: {
    stop("Only integer and numeric vectors are supported");
  }
  }}

// vec to mat
arma::mat vec2mat(arma::vec x, int nrow, int ncol) {
  arma::mat y(x);
  y.reshape(nrow, ncol);
  return y;
}

Function dat_by_clus("dat_by_cluster");

IntegerVector one(1,1);

// [[Rcpp::export]]
List htdp_update(const NumericMatrix dat, IntegerVector& indic, List& mu, List& V, 
                  const double nu, const arma::mat diag_a, const double rho, 
                  const arma::vec xi, const NumericVector inact_prop) {
  
  const int n = dat.nrow();
  const int p = dat.ncol();
  
  for (int i = 0; i < n; i++) { 
    
    IntegerVector cluster = sort_unique(indic);
    int k = cluster.size();
    IntegerVector nj = table(indic);
    
    const IntegerVector int_curr_indic(1, indic(i));
    const String curr_indic(indic(i));
    
    arma::vec curr_x = dat(i,_);
    
    const List dbc = dat_by_clus(Named("dat")=dat, Named("indic")=indic, Named("clus")=cluster);
    
    NumericVector log_likeli(k);
    NumericVector likeli(k);
    
    const arma::mat xi_mat = vec2mat(xi, p, 1);
    const arma::mat curr_x_mat = vec2mat(curr_x, p ,1);
    
    // calculate log likelihood
    for (int j=0; j<k; j++) {
      
      NumericMatrix curr_dat = dbc(j);
      
      if (curr_dat.ncol() != 1) { 
        
        arma::vec xi_star = (as<arma::vec>(colSums(curr_dat)) + rho*xi) / (rho + nj(j));
        const arma::mat xi_star_mat = vec2mat(xi_star, p, 1);
        const arma::mat curr_dat_mat = as<arma::mat>(curr_dat);
        arma::mat A_star = 2*nu*diag_a + rho*(xi_mat*xi_mat.t()) + curr_dat_mat.t()*curr_dat_mat -
          (xi_star_mat * xi_star_mat.t()) * (rho + nj(j));
        const arma::mat curr_centered = curr_x_mat - xi_star_mat;
        const arma::mat x_xi_star = curr_centered * curr_centered.t();
        
        const double lmvg_1 = lgamma((nu+p+nj(j))/2);
        const double lmvg_2 = lgamma((nu+nj(j))/2);
        log_likeli(j) = p/2*log((rho+nj(j))/(rho+nj(j)+1)) - p/2*log(M_PIl) + 
          lmvg_1 - lmvg_2 + (nu+p-1+nj(j))/2 * log(arma::det(A_star)) -
          (nu+p+nj(j))/2*log(arma::det(A_star + (rho+nj(j))/(rho+nj(j)+1)*x_xi_star));
        
      } else {
        
        arma::vec xi_star = (as<arma::vec>(transpose(curr_dat)) + rho*xi) / (rho + nj(j));
        const arma::mat xi_star_mat = vec2mat(xi_star, p, 1);
        const arma::mat curr_dat_mat = as<arma::mat>(curr_dat);
        arma::mat A_star = 2*nu*diag_a + rho*(xi_mat*xi_mat.t()) + curr_dat_mat*curr_dat_mat.t() -
          (xi_star_mat * xi_star_mat.t()) * (rho + nj(j));
        const arma::mat curr_centered = curr_x_mat - xi_star_mat;
        const arma::mat x_xi_star = curr_centered * curr_centered.t();
        
        const double lmvg_1 = lgamma((nu+p+nj(j))/2);
        const double lmvg_2 = lgamma((nu+nj(j))/2);
        log_likeli(j) = p/2*log((rho+nj(j))/(rho+nj(j)+1)) - p/2*log(M_PIl) + 
          lmvg_1 - lmvg_2 + (nu+p-1+nj(j))/2 * log(arma::det(A_star)) -
          (nu+p+nj(j))/2*log(arma::det(A_star + (rho+nj(j))/(rho+nj(j)+1)*x_xi_star));
        
      }
    }
    
    // calculate likelihood
    likeli = exp(log_likeli);
    
    
    // nj adjusting
    const IntegerVector seq_cluster = seq_len(k);
    const IntegerVector vec_idx = mywhich(cluster, int_curr_indic);
    const int idx = vec_idx(0);
    
    const arma::vec temp_mu = as<arma::vec>(mu[idx-1]);
    
    if (nj[curr_indic]==1) { 
      
      likeli = likeli[seq_cluster != idx];
      V = V[seq_cluster != idx];
      mu = mu[seq_cluster != idx];
      nj = nj[seq_cluster != idx];
      
      k = k - 1;
      
    } else {
      nj[curr_indic] = nj[curr_indic] - 1;
    }
    const CharacterVector active_names = nj.names();
    
    
    // calculate active and inactive proportions
    NumericVector prop(k + 1);
    for (int l = 0; l < k; l++) prop(l) = likeli(l)*nj(l);
    
    prop(k) = inact_prop(i);
    prop = prop/sum(prop);
    
    // update indicator
    const IntegerVector rand_mult = oneMultinomCalt(prop);
    const IntegerVector vec_position = mywhich(rand_mult, one);
    const int int_position = vec_position(0);
    
    
    if (int_position <= k) {
      
      const std::string char_new_cluster = as<std::string>(active_names[int_position-1]);
      const int int_new_cluster = std::stoi(char_new_cluster);
      indic(i) = int_new_cluster;
      
    } else {
      
      indic(i) = max(cluster) + 1;
      
      arma::vec part1 = curr_x-temp_mu;
      arma::vec part2 = temp_mu-xi;
      arma::mat part1_mat = vec2mat(part1, p, 1);
      arma::mat part2_mat = vec2mat(part2, p, 1);
      
      arma::mat part3 = part1_mat * part1_mat.t() + rho*(part2_mat * part2_mat.t()) + 2*nu*diag_a;
      
      arma::mat new_V = riwish(nu+p+1, part3);
      //arma::mat sigma = new_V / (1+rho);
      //NumericMatrix sigma = as<NumericMatrix>(new_V)/(1+rho);
      
      arma::mat new_mu = rmvnorm(1, (curr_x + rho*xi)/(rho+1), new_V/(1+rho));
      
      mu.push_back(new_mu);
      V.push_back(new_V);
      
    }
    
    
  }
  
  
  return List::create(indic, mu, V);
  //return List::create(nj);
  
}




