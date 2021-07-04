#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]  

using namespace Rcpp;
// [[Rcpp::export]]
List clustercpp (arma::mat rawcl,
                 arma::mat X,
                 arma::mat Res,
                 arma::mat invX,
                 double q
                ){
    arma::colvec unique_cl;
    arma::uvec index;
    arma::mat subX;
    unique_cl = unique(rawcl);
    int i;
    int N0=unique_cl.n_rows;
    arma::mat hmeat;
    arma::mat hmeat2;
    arma::mat meat;
    arma::mat vcov;

    for(i=0;i<=N0;i++){
        index = find(rawcl == unique_cl[i]);
        subX = X.rows(index);
        hmeat = subX.t()*Res.rows(index);
        hmeat2 = hmeat*hmeat.t();
        if(i==0){
            meat = hmeat2;
        }else{
            meat = meat + hmeat2;
        }
    }
    vcov = q*invX*meat*invX;
    List output;
    output["unique_cl"] = unique_cl;
    output["meat"] = meat;
    output["vcov"] = vcov;
    return(output);
}

// [[Rcpp::export]]
List ivclustercpp (arma::mat rawcl,
                 arma::mat X,
                 arma::mat Res,
                 arma::mat invxPzx,
                 double q,
                 arma::mat invzzx
                ){
    arma::colvec unique_cl;
    arma::uvec index;
    arma::mat subX;
    unique_cl = unique(rawcl);
    int i;
    int N0=unique_cl.n_rows;
    arma::mat hmeat;
    arma::mat hmeat2;
    arma::mat meat;
    arma::mat vcov;

    for(i=0;i<=N0;i++){
        index = find(rawcl == unique_cl[i]);
        subX = X.rows(index);
        hmeat = subX.t()*Res.rows(index);
        hmeat2 = hmeat*hmeat.t();
        if(i==0){
            meat = hmeat2;
        }else{
            meat = meat + hmeat2;
        }
    }

    vcov = q*invxPzx*(invzzx.t()*meat*invzzx)*invxPzx;
    List output;
    output["unique_cl"] = unique_cl;
    output["meat"] = meat;
    output["vcov"] = vcov;
    return(output);
}