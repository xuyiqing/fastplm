#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]  

using namespace Rcpp;

// [[Rcpp::export]]
List solveiv (arma::mat Y,
              arma::mat X,
              arma::mat Z
              ){

    arma::mat invzz ;
    arma::mat invxx;
    arma::mat XZ;
    arma::mat ZY;
    arma::mat XPzX;
    arma::mat invXPzX;
    arma::mat invzzx;
    arma::mat ivcoef; 
    arma::mat Y_hat; 
    arma::mat residual;
    
    invxx = inv_sympd(X.t()*X);
    invzz = inv_sympd(Z.t()*Z);
    XZ = X.t()*Z;
    XPzX = XZ*invzz*XZ.t();
    invXPzX = inv_sympd(XPzX);
    ZY = Z.t()*Y;
    invzzx = invzz*XZ.t();
    ivcoef = invXPzX * invzzx.t() * ZY;
    Y_hat = X * ivcoef;
    residual = Y - Y_hat;
    
    List output;
    
    output["invxx"] = invxx; 
    output["invzz"] = invzz;
    output["invXPzX"] = invXPzX;
    output["ZY"] = ZY;
    output["XZ"] = XZ;
    output["coef"] = ivcoef;
    output["residual"] = residual;
    output["invzzx"] = invzzx;
    output["Y_hat"] = Y_hat;
    return(output);
}

// [[Rcpp::export]]
arma::mat solvecpp (arma::mat X){
    arma::mat invxx;
    invxx = inv_sympd(X.t()*X);
    return(invxx);
}

// [[Rcpp::export]]
List solvegmm (arma::mat Y,
                    arma::mat X,
                    arma::mat Z, 
                    arma::mat u_hat
                    ){
    int L = Z.n_cols ;                  
    arma::mat gmm_coef;
    arma::mat XZ;
    arma::mat ZY;
    arma::mat invmeat;
    arma::mat meat;
    arma::mat meat_sub;
    arma::mat solve1;
    arma::mat solve2;
    arma::mat u_hat_gmm;
    arma::mat hansen_sub;
    arma::mat hansen;

    XZ = X.t()*Z;
    ZY = Z.t()*Y;
    meat_sub = Z % repmat(u_hat,1,L);
    meat = meat_sub.t()*meat_sub;
    invmeat = inv(meat);
    solve1 = XZ*invmeat*XZ.t();
    solve2 = XZ*invmeat*ZY;
    gmm_coef = solve(solve1, solve2);
    u_hat_gmm = Y - X*gmm_coef;
    hansen_sub = u_hat_gmm.t()*Z;
    hansen = hansen_sub * invmeat * hansen_sub.t();

    List output;
    output["invmeat"] = invmeat;
    output["coef"] = gmm_coef;
    output["meat"] = meat;
    output["u_hat"] = u_hat_gmm;
    output["hansen"] = hansen;
    return(output);
}

// [[Rcpp::export]]
List solvegmm_meat (arma::mat Y,
                    arma::mat X,
                    arma::mat Z, 
                    arma::mat meat
                    ){
    arma::mat gmm_coef;
    arma::mat XZ;
    arma::mat ZY;
    arma::mat invmeat;
    arma::mat solve1;
    arma::mat solve2;
    arma::mat u_hat_gmm;
    arma::mat hansen_sub;
    arma::mat hansen;

    XZ = X.t()*Z;
    ZY = Z.t()*Y;
    invmeat = inv(meat);
    solve1 = XZ*invmeat*XZ.t();
    solve2 = XZ*invmeat*ZY;
    gmm_coef = solve(solve1, solve2);
    u_hat_gmm = Y - X*gmm_coef;
    hansen_sub = u_hat_gmm.t()*Z;
    hansen = hansen_sub * invmeat * hansen_sub.t();

    List output;
    output["invmeat"] = invmeat;
    output["coef"] = gmm_coef;
    output["meat"] = meat;
    output["u_hat"] = u_hat_gmm;
    output["hansen"] = hansen;
    return(output);
}


// [[Rcpp::export]]
List simpleols (arma::mat Y,
                arma::mat X
                ){
    arma::mat coef;
    arma::mat u_hat;
    arma::mat Y_hat;

    coef=solve(X,Y);
    Y_hat=X*coef;
    u_hat= Y - Y_hat;
    List output;
    output["coef"] = coef;
    output["u_hat"] = u_hat;
    output["Y_hat"] = Y_hat;
    return(output);
}

