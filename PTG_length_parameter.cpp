// Author: Minoru Kusaba
// Affiliation: School of Multidisciplinary Sciences Department of Statistical Science SOKENDAI
// Contact: kusaba@ism.ac.jp
// File name: PTG_length_parameter.cpp
// Task: define functions for sampling length-scale prameters in GTM_LDLV algorithm
// Last update: 2019/08/22

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;  // use the Armadillo library for matrix computations
using namespace Rcpp;

// create matrix with r vectors
// [[Rcpp::export]]
mat get_m(vec r,int K){
  mat m(K,K,fill::zeros);
  for(int i=0; i<K; ++i){
    for(int j=0; j<K; ++j){
      m(i,j) = r(i);
    }
  }
  return m;
}

// give Ch matrix(covariance matrix for h)
// [[Rcpp::export]]
mat give_Ch(mat m,mat d,double epsilon,double L){
  mat tm = m.t();
  mat Ch = pow((2*exp(m+tm)/(exp(2*m)+exp(2*tm))),L/2)%exp(-d/(exp(2*m)+exp(2*tm)));
  Ch.diag() = Ch.diag() * (1 + epsilon);
  return Ch;
}

// conbine above functions(instance function for creating Ch)
// [[Rcpp::export]]
mat get_Ch(vec r,mat d,double epsilon,double L){
  int K = r.n_elem;
  mat m = get_m(r,K);
  return give_Ch(m,d,epsilon,L);
}

// give log density of posterior distribution of r
// [[Rcpp::export]]
double give_logDens(vec r,mat H,mat inCst,mat d,double epsilon,double L){
  double D = H.n_rows;
  mat Ch = get_Ch(r,d,epsilon,L);
  mat inCh = pinv(Ch); //genelized inverse matrix
  double logdet = real(log_det(Ch));
  double Tird = sum((inCst * r) % r);
  double logDens = (-D/2)*logdet-(0.5)*trace(H * inCh * H.t())-(0.5)*Tird;
  return logDens;
}


// core function to derive gradient of posterior distribution of r
// [[Rcpp::export]]
mat give_dCh(mat m,mat d,double L){
  mat tm = m.t();
  mat dCh = (L/2) * pow((2*exp(m+tm)/(exp(2*m)+exp(2*tm))),L/2) % exp(-d/(exp(2*m)+exp(2*tm))) % (exp(4*tm)-exp(4*m) + (1/L)*4*d%exp(2*m)) / pow((exp(2*m)+exp(2*tm)),2);
  return dCh;
}

// core function to derive hessian of posterior distribution of r
// [[Rcpp::export]]
mat give_ddCh(mat m,mat d,double L){
  mat tm = m.t();
  mat ddCh = ((L/2) * pow((2*exp(m+tm)/(exp(2*m)+exp(2*tm))),L/2) % exp(-d/(exp(2*m)+exp(2*tm))) / pow((exp(2*m)+exp(2*tm)),3)) % ( (exp(4*tm)-exp(4*m) + (1/L)*4*d%exp(2*m))%(2*d%exp(2*tm)/(exp(2*m)+exp(2*tm)) + (L/2)*(exp(2*m)-exp(2*tm))) + 4*exp(2*m + 2*tm)%(exp(2*m)+exp(2*tm)-(1/L)*4*d) );                           ;
  return ddCh;
}

// for each elements(ddCh)
// [[Rcpp::export]]
mat give_ddChij(mat ddCh,int i,int j,int K){
  mat m(K,K,fill::zeros);
  m(i,j) = ddCh(i,j);
  mat x = m + m.t();
  return x;
}

// for each elements(dCh)
// [[Rcpp::export]]
mat give_dChi(mat dCh,int i,int K){
  mat m(K,K,fill::zeros);
  m.row(i) = dCh.row(i);
  mat x = m + m.t();
  return x;
}

// give gradient of posterior distribution of r
// [[Rcpp::export]]
vec give_dr(vec r,mat d,mat H,mat inCst,double epsilon,double L){
  int K = r.n_elem;
  double D = H.n_rows;
  vec predr(K,fill::zeros);
  mat m = get_m(r,K);
  mat Ch = give_Ch(m,d,epsilon,L);
  mat iCh = pinv(Ch); //genelized inverse matrix
  mat dCh = give_dCh(m,d,L);
  for(int i=0; i<K; ++i){
    mat dChi = give_dChi(dCh,i,K);
    double X = (-D/2)*trace(iCh * dChi) + 0.5*trace(H * iCh * dChi * iCh * H.t());
    predr(i) = X;
  }
  vec dr = predr - inCst * r;
  return dr;
}

// give hessian of posterior distribution of r
// [[Rcpp::export]]
mat give_ddr(vec r,mat d,mat H,mat inCst,double epsilon,double L){
  int K = r.n_elem;
  double D = H.n_rows;
  mat x(K,K,fill::zeros);
  mat m = get_m(r,K);
  mat Ch = give_Ch(m,d,epsilon,L);
  mat iCh = pinv(Ch); //genelized inverse matrix
  mat ddCh = give_ddCh(m,d,L);
  mat dCh = give_dCh(m,d,L);
  for(int i=0; i<K; ++i){
    for(int j=0; j<(i+1); ++j){
      mat ddChij = give_ddChij(ddCh,i,j,K);
      mat dChi = give_dChi(dCh,i,K);
      mat dChj = give_dChi(dCh,j,K);
      x(i,j) = (-D/2)*trace(iCh * (ddChij-dChj*iCh*dChi)) + 0.5*trace(H*iCh*(ddChij-dChj*iCh*dChi-dChi*iCh*dChj)*iCh*H.t()) - inCst(i,j);
    }
  }
  mat A = x + x.t();
  A.diag() = x.diag();
  return A;
}