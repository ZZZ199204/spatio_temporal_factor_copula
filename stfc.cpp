#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <vector>
#include <complex>
#include <string>
#include <numeric>
#include <algorithm>

// [[Rcpp::depends(RcppGSL)]]

// [[Rcpp::depends(RcppArmadillo)]]

using std::pow; using std::log;
using std::exp; using std::vector;
using std::complex; using std::size_t;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;

// [[Rcpp::export]]
double nllk_p_C(NumericVector par, NumericMatrix dat, NumericMatrix centers, int Kcen, double g){
    NumericVector tt = dat(Rcpp::_, 0);
    NumericVector y = dat(Rcpp::_, 1);
    NumericVector LON = dat(Rcpp::_, 2);
    NumericVector LAT = dat(Rcpp::_, 3);
    
    const int M = 12;  // number of time periods (monthly data) per year
    const double w = 2*PI/M;
    double tmp = 0.0;
    
    vector<double> weig1;
    
    int i,j;
  
    for(j=0;j<Kcen+1;++j){
        weig1.push_back(par[j]); // the last one is for intercept for hx
    }
    
    double a0 = par[1+Kcen];
    double a1 = par[2+Kcen];
    double a2 = par[3+Kcen];
    double a3 = par[4+Kcen];
    
    int dd = dat.nrow();
    double sumj = 0.0;
    double pvec = 0.0;
    double nllk = 0.0;
    double tmpx = 0.0;
    
    for(i=0;i<dd;++i){
        sumj = weig1[Kcen];
        for(j=0;j<Kcen;++j){
            sumj +=  weig1[j]*exp(-g*(pow((LON[i] - centers(j,0)),2) + pow((LAT[i] - centers(j,1)),2)));
            //tmpx = pow((LON[i] - centers(j,0)),2) + pow((LAT[i] - centers(j,1)),2);
            //sumj +=  weig1[j]*(tmpx*log(pow(tmpx, 0.5)));
        }
        tmp = exp(sumj + a0 + a1*tt[i] + a2*sin(w*tt[i]) + a3*cos(w*tt[i]));
        
        if(y[i] == 1){
            pvec = tmp / (1 + tmp);
        }else{
            pvec = 1 / (1 + tmp);
        }
         // Rcpp::Rcout << pvec << std::endl;
         // Rcpp::Rcout << " i = " << i << " centers(j,0) " << j << " " << centers(j,0) << std::endl;
        if(!R_finite(pvec)){
            pvec = 0.0001;
        }else if(pvec > 0.9999999){
            pvec = 0.9999999;
        }else if(pvec < 0.0000001){
            pvec = 0.0000001;
        }
        nllk -= log(pvec);
    }
    return nllk;
}

