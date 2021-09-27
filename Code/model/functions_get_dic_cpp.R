c.get_loglik =  cxxfunction(
  signature(Yy="numeric", ALPHA="numeric", DELTA="numeric", Xitm1="numeric",
            DIMS="integer"),
  body='
/*Dims is c(n,p,TT)
  X pxTxn; Y pxpxT; 
  */
  Rcpp::IntegerVector dims(DIMS);
  Rcpp::NumericVector YY(Yy);
  arma::cube Y(YY.begin(),dims[0],dims[0],dims[2]); // n X n X TT
  Rcpp::NumericVector X_EST(Xitm1);
  arma::cube X_est(X_EST.begin(),dims[1],dims[2],dims[0]); //p X TT X n
  double Alpha = Rcpp::as<double>(ALPHA);
  double Delta = Rcpp::as<double>(DELTA);
  double loglik = 0;
  double d_tij = 0;
  double eta_tij = 0;
  
  // calculate log(y1|x1)
  for(int i=0; i<(dims[0]-1); i++){// n
    for(int j=(i+1); j<dims[0]; j++){
      d_tij = arma::norm(X_est.slice(i).col(0)-X_est.slice(j).col(0),2);
      eta_tij = Alpha - d_tij;
      loglik = loglik + Y.slice(0)(i,j)*eta_tij - log(1+exp(eta_tij));
    }
  }
  // calculate log(yt|xt), t>=2
  for(int tt=1; tt < dims[2]; tt++){ //TT
    for(int i=0; i<(dims[0]-1); i++){// n
      for(int j=(i+1); j<dims[0]; j++){
        d_tij = arma::norm(X_est.slice(i).col(tt)-X_est.slice(j).col(tt),2);
        eta_tij = Alpha - d_tij + Delta*Y.slice(tt-1)(i,j);
        loglik = loglik + Y.slice(tt)(i,j)*eta_tij - log(1+exp(eta_tij));
      }
    }
  }
  return Rcpp::List::create(loglik);
  ', plugin="RcppArmadillo")

