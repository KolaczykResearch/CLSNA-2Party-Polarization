# This code is developed upon the source code for paper: Sewell&Chen, Latent Space Models for Dynamic Networks

# ------------ Define function to draw samples for X, alpha, delta-------------------------- #
c.update1 <-  cxxfunction(
  signature(Xitm1="numeric",DIMS="integer",TUNEXs="numeric",Yy="numeric",
            ALPHA="numeric",DELTA="numeric",TUNEA="numeric", TUNED="numeric",
            WW="numeric",t2X="numeric",s2X="numeric",
            xiALPHA="numeric",xiDELTA="numeric",nuALPHA="numeric",
            nuDELTA="numeric",CAUCHY="integer",
            RNORMS="numeric",RNORMSAD="numeric", PIX="numeric", IT='integer', SUBWINDOW='integer'),
  body='
/*Dims is c(n,p,TT,K)
  Z nxT; X pxTxn; Y pxpxT; 
  */
  Rcpp::IntegerVector dims(DIMS);
  int Cauchy = Rcpp::as<int>(CAUCHY);
  
  Rcpp::NumericVector XOLD(Xitm1);
  arma::cube Xold(XOLD.begin(),dims[1],dims[2],dims[0]); //p X TT X n
  arma::cube Xnew = arma::zeros(dims[1],dims[2],dims[0]);

  for(int i=0;i<dims(0);i++){ // for each agent i
    Xnew.slice(i) = Xold.slice(i);
  }

  arma::mat Xconc = arma::zeros(dims(0)*dims(2),dims(1)); //(nXTT) X p
  arma::colvec centerings = arma::zeros(dims(1),1); //column vector of length p
  
  Rcpp::NumericVector rnormsVec(RNORMS); 
  arma::cube rnorms(rnormsVec.begin(),dims(1),dims(2),dims(0));
  arma::colvec rnormsAD = Rcpp::as<arma::colvec>(RNORMSAD);
  arma::colvec pi = Rcpp::as<arma::colvec>(PIX);

  Rcpp::NumericVector YY(Yy);
  arma::cube Y(YY.begin(),dims[0],dims[0],dims[2]); // n X n X TT
  
  double Alpha = Rcpp::as<double>(ALPHA);
  double Delta = Rcpp::as<double>(DELTA);
  double xiAlpha = Rcpp::as<double>(xiALPHA);
  double xiDelta = Rcpp::as<double>(xiDELTA);
  double nuAlpha = Rcpp::as<double>(nuALPHA);
  double nuDelta = Rcpp::as<double>(nuDELTA);
  double AlphaNew =0, DeltaNew =0;
  Rcpp::NumericVector TUNEXsVec(TUNEXs); 
  arma::mat tuneXs(TUNEXsVec.begin(), dims(0), dims(2)); //n X TT
  //double tunex = Rcpp::as<double>(TUNEX);
  //double tuneAD = Rcpp::as<double>(TUNEAD);
  double tuneA = Rcpp::as<double>(TUNEA);
  double tuneD = Rcpp::as<double>(TUNED);

  double t2 = Rcpp::as<double>(t2X);
  double s2 = Rcpp::as<double>(s2X);
  
  arma::colvec ww = Rcpp::as<arma::colvec>(WW);
  
  double AccProb =0;
  double dz=0, dx=0, uu=0, count1=0, count2=0;
  arma::mat insides = arma::zeros(1,1);
  arma::colvec muit = arma::zeros(dims[1],1); // p 
  arma::colvec avg1 = arma::zeros(dims[1],1);
  arma::colvec avg2 = arma::zeros(dims[1],1);
  arma::colvec avg1_old = arma::zeros(dims[1],1);
  arma::colvec avg2_old = arma::zeros(dims[1],1);
  arma::colvec Ai_1 = arma::zeros(dims[1],1);
  arma::colvec Ai_2 = arma::zeros(dims[1],1);
  arma::colvec Ai_1_old = arma::zeros(dims[1],1);
  arma::colvec Ai_2_old = arma::zeros(dims[1],1);
  arma::colvec AccRate = arma::zeros(dims(0)*dims(2)+4,1); // n*TT + 4
  int it = Rcpp::as<int>(IT); // current iteration
  int subwindow = Rcpp::as<int>(SUBWINDOW); // current window index
  arma::vec vector = {1, 1}; // for update of tuning params
  /*-------------------- Latent Positions-------------------*/
  
  for(int tt=0; tt < dims[2]; tt++){ // TT
    for(int i=0; i<dims[0]; i++){// n

      AccProb=0; // X_ti

      // propose Xnew.slice(i).col(tt) by adding normal perturbation
      if(Cauchy<0.5){
        /*---Normal Random Walk---*/
        Xnew.slice(i).col(tt) = Xold.slice(i).col(tt) + exp(tuneXs(i,tt))*rnorms.slice(i).col(tt);
      }else{
        /*---Cauchy Random Walk---*/
        for(int ell=0; ell<dims(1); ell++){
          uu = arma::randu();
          Xnew.slice(i)(ell,tt) = Xold.slice(i)(ell,tt) + exp(tuneXs(i,tt))*tan(PI*(uu-0.5));
        }
      }

      // calculate the AccProb for Xnew.slice(i).col(tt)
      for(int j=0; j<dims[0]; j++){ // n, log(pt,ij)+log(pt, ji) term
        if(j != i){
          dx = arma::norm(Xold.slice(i).col(tt)-Xnew.slice(j).col(tt),2); // s(z_ti^old, z_tj^new) 
          dz = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(j).col(tt),2); // s(z_ti^new, z_tj^new), updated distance between i and j
          if(tt==0){ //no delta term
            AccProb += (dx-dz)*(Y.slice(tt)(j,i)+Y.slice(tt)(i,j));
            AccProb += 2*log(1+exp(Alpha-dx));
            AccProb -= 2*log(1+exp(Alpha-dz));
          }
          if(tt>0){
            AccProb += (dx-dz)*(Y.slice(tt)(j,i)+Y.slice(tt)(i,j));
            AccProb += 2*log(1+exp(Alpha+Delta*Y.slice(tt-1)(i,j)-dx));
            AccProb -= 2*log(1+exp(Alpha+Delta*Y.slice(tt-1)(i,j)-dz));
          }
        }
      }
      
      // normal term
      if(tt == 0){
        // add term regarding N(Z_t,i| 0, tau2)
        insides = trans(Xnew.slice(i).col(tt))*(Xnew.slice(i).col(tt))/t2;
        AccProb += -0.5*insides(0,0); //from 1X1 matrix to a number
        insides = trans(Xold.slice(i).col(tt))*(Xold.slice(i).col(tt))/t2;
        AccProb -= -0.5*insides(0,0);
      }

      if(tt > 0){
        // add term regarding N(Z_t,i| Z_t-1,i + gamma_1*A_i_1 + gamma_2*A_i_2, sigma2)
        // find neighbors of i, calculate attractors
        avg1 = arma::zeros(dims[1],1);
        avg2 = arma::zeros(dims[1],1);
        count1 = 0;
        count2 = 0;
        for(int j=0; j<dims[0]; j++){
          if(Y.slice(tt-1)(i, j)==1){
            if(pi(i) == pi(j)){
              avg1 += Xnew.slice(j).col(tt-1);
              count1 += 1;
            }else{
              avg2 += Xnew.slice(j).col(tt-1);
              count2 += 1;
            }
          } 
        }
        if(count1 != 0){
          Ai_1 = avg1/count1 - Xnew.slice(i).col(tt-1);
        }else{
          Ai_1 = arma::zeros(dims[1],1); //zeros
        }
        if(count2 != 0){
          Ai_2 = avg2/count2 - Xnew.slice(i).col(tt-1);
        }else{
          Ai_2 = arma::zeros(dims[1],1); //zeros
        }
        if(pi(i)==1){
          muit = Xnew.slice(i).col(tt-1)+ww(0)*Ai_1+ww(2)*Ai_2;
        }else{
          muit = Xnew.slice(i).col(tt-1)+ww(1)*Ai_1+ww(2)*Ai_2;
        }
        insides = trans(Xnew.slice(i).col(tt)-muit)*(Xnew.slice(i).col(tt)-muit)/s2;
        AccProb += -0.5*insides(0,0);
    
        insides = trans(Xold.slice(i).col(tt)-muit)*(Xold.slice(i).col(tt)-muit)/s2;
        AccProb -= -0.5*insides(0,0);  
      }

      if(tt < dims[2]-1){
        // add term regarding N(Z_t+1,i| Z_t,i + gamma_i*A_i + gamma_i_2*A_i_2, sigma2)
        avg1 = arma::zeros(dims[1],1);
        avg2 = arma::zeros(dims[1],1);
        count1 = 0;
        count2 = 0;
        for(int j=0; j<dims[0]; j++){
          if(Y.slice(tt)(i, j)==1){
            if(pi(i) == pi(j)){
              avg1 += Xnew.slice(j).col(tt);
              count1 += 1;
            }else{
              avg2 += Xnew.slice(j).col(tt);
              count2 += 1;
            }
          } 
        }
        if(count1 != 0){
          Ai_1 = avg1/count1 - Xnew.slice(i).col(tt);
          Ai_1_old = avg1/count1 - Xold.slice(i).col(tt);
        }else{
          Ai_1 = arma::zeros(dims[1],1); //zeros
          Ai_1_old = arma::zeros(dims[1],1);
        }
        if(count2 != 0){
          Ai_2 = avg2/count2 - Xnew.slice(i).col(tt);
          Ai_2_old = avg2/count2 - Xold.slice(i).col(tt);
        }else{
          Ai_2 = arma::zeros(dims[1],1); //zeros
          Ai_2_old = arma::zeros(dims[1],1); //zeros
        }
        if(pi(i)==1){
          muit = Xnew.slice(i).col(tt)+ww(0)*Ai_1+ww(2)*Ai_2;
        }else{
          muit = Xnew.slice(i).col(tt)+ww(1)*Ai_1+ww(2)*Ai_2;
        }
        insides = trans(Xnew.slice(i).col(tt+1)-muit)*(Xnew.slice(i).col(tt+1)-muit)/s2;
        AccProb += -0.5*insides(0,0); 
        
        if(pi(i)==1){
          muit =  Xold.slice(i).col(tt)+ww(0)*Ai_1_old+ww(2)*Ai_2_old;
        }else{
          muit =  Xold.slice(i).col(tt)+ww(1)*Ai_1_old+ww(2)*Ai_2_old;
        }
        insides = trans(Xnew.slice(i).col(tt+1)-muit)*(Xnew.slice(i).col(tt+1)-muit)/s2;
        AccProb -= -0.5*insides(0,0);  
      }
      
      if(tt < dims[2]-1){
        // add term regarding sum_{j in N_i} log N(Z_t+1,j| Z_t,j + gamma_j*A_j + gamma_j_2*A_j_2, sigma2)
        for(int j=0; j<dims[0]; j++){
          if(Y.slice(tt)(i, j)==1){ //j is a neighbor of i at time tt
            // add term N(Z_t+1,j| Z_t,j + gamma_j*A_j + gamma_j_2*A_j_2, sigma2)
            avg1 = arma::zeros(dims[1],1);
            avg1_old = arma::zeros(dims[1],1);
            avg2 = arma::zeros(dims[1],1);
            avg2_old = arma::zeros(dims[1],1);
            count1 = 0;
            count2 = 0;
            for(int k=0; k<dims[0]; k++){ //for node j, find average of neighbors by memberships
              if(Y.slice(tt)(j, k)==1){
                if(pi(j) == pi(k)){
                  if(k==i){
                    avg1 += Xnew.slice(k).col(tt);
                    avg1_old += Xold.slice(k).col(tt);
                    count1 += 1;
                  }else{
                    avg1 += Xnew.slice(k).col(tt);
                    avg1_old += Xnew.slice(k).col(tt);
                    count1 += 1;
                  }
                }else{
                  if(k==i){
                    avg2 += Xnew.slice(k).col(tt);
                    avg2_old += Xold.slice(k).col(tt);
                    count2 += 1;
                  }else{
                    avg2 += Xnew.slice(k).col(tt);
                    avg2_old += Xnew.slice(k).col(tt);
                    count2 += 1;
                  }
                }
              } 
            }
            if(count1 != 0){
              Ai_1 = avg1/count1 - Xnew.slice(j).col(tt);
              Ai_1_old = avg1_old/count1 - Xnew.slice(j).col(tt);
            }else{
              Ai_1 = arma::zeros(dims[1],1); //zeros
              Ai_1_old = arma::zeros(dims[1],1);
            }
            if(count2 != 0){
              Ai_2 = avg2/count2 - Xnew.slice(j).col(tt);
              Ai_2_old = avg2_old/count2 - Xnew.slice(j).col(tt);
            }else{
              Ai_2 = arma::zeros(dims[1],1); //zeros
              Ai_2_old = arma::zeros(dims[1],1); //zeros
            }
            muit = Xnew.slice(j).col(tt)+ww(pi(j)-1)*Ai_1+ww(2)*Ai_2;
            insides = trans(Xnew.slice(j).col(tt+1)-muit)*(Xnew.slice(j).col(tt+1)-muit)/s2;
            AccProb += -0.5*insides(0,0); 
            
            muit = Xnew.slice(j).col(tt)+ww(pi(j)-1)*Ai_1_old+ww(2)*Ai_2_old;
            insides = trans(Xnew.slice(j).col(tt+1)-muit)*(Xnew.slice(j).col(tt+1)-muit)/s2;
            AccProb -= -0.5*insides(0,0);  
          } 
        }
      }
      
      uu= arma::randu();
      if(uu<exp(AccProb)){
        AccRate(4+tt*dims(0)+i) = 1; //4+tt*n + i
      }else{
        Xnew.slice(i).col(tt) = Xold.slice(i).col(tt);
      }
    
      if(subwindow!=1 && tt==0){
        Xnew.slice(i).col(tt) = Xold.slice(i).col(tt);
      }
      // update tuneXs.row(i).col(tt)
      vector = {1, exp(AccProb)};
      tuneXs(i, tt) = tuneXs(i, tt) + 1/pow(it,0.8)*(arma::min(vector)-0.234);
    }
  }
  
  /*---Centering---*/
  // concatenate all latent positions
  for(int i=0; i<dims(0); i++){ // n
    for(int tt=0; tt<dims(2); tt++){ // TT
      Xconc.row(i*dims(2)+tt) = trans(Xnew.slice(i).col(tt)); // Xconc is of dimension (n*TT) X 2
    }
  }
  // find centers for each dimensions
  for(int ell=0; ell<dims(1); ell++){ //p
    centerings(ell) = sum(Xconc.col(ell))/(dims(0)*dims(2)); //n*TT, mean of each dimension of latent positions
  }
  for(int i=0;i<dims(0);i++){ //n 
    for(int tt=0;tt<dims(2);tt++){ // TT
      for(int ell=0;ell<dims(1);ell++){ //p
        Xnew.slice(i)(ell,tt) = Xnew.slice(i)(ell,tt) - centerings(ell);
      }
    }
  }
 
  /*-------------------- Alpha ------------------------------*/

  AccProb=0;

  // propose AlphaNew (AlphaNew) by adding normal perturbation
  if(Cauchy<0.5){
    // AlphaNew = Alpha + tuneAD*arma::randn();
    // AlphaNew = Alpha + tuneAD*rnormsAD(0);
    // AlphaNew = Alpha + tuneA*rnormsAD(0);
    AlphaNew = Alpha + exp(tuneA)*rnormsAD(0);
  }else{
    uu = arma::randu();
    AlphaNew = Alpha + tuneA*tan(PI*(uu-0.5));
  }
  
  for(int tt=0;tt<dims(2);tt++){
    for(int i = 0; i < dims(0); i++){
      for(int j = 0; j < dims(0); j++){
        if(i != j){
          dx = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(j).col(tt), 2);
          if(tt==0){
            if(subwindow==1){
              AccProb += Y.slice(tt)(i,j)*(AlphaNew-Alpha);
              AccProb += log(1+exp(Alpha-dx));
              AccProb -= log(1+exp(AlphaNew-dx));
            }
          }
          if(tt>0){
            AccProb += Y.slice(tt)(i,j)*(AlphaNew-Alpha);
            AccProb += log(1+exp(Alpha+Delta*Y.slice(tt-1)(i,j)-dx));
            AccProb -= log(1+exp(AlphaNew+Delta*Y.slice(tt-1)(i,j)-dx));
          }
        }
      }
    }
  }
  
  AccProb += -0.5*(AlphaNew-nuAlpha)*(AlphaNew-nuAlpha)/xiAlpha;
  AccProb -= -0.5*(Alpha-nuAlpha)*(Alpha-nuAlpha)/xiAlpha;
  
  uu= arma::randu();
  if(uu<exp(AccProb)){
    AccRate(0) = 1;
  }else{
    AlphaNew = Alpha;
  }

  // update tuneA 
  vector = {1, exp(AccProb)};
  tuneA = tuneA + 1/pow(it,0.8)*(arma::min(vector)-0.234);

  /*-------------------- Delta ------------------------------*/

  AccProb=0;

  // propose DeltaNew (DeltaNew) by adding normal perturbation
  if(Cauchy<0.5){
    // DeltaNew = Delta + tuneAD*arma::randn();
    // DeltaNew = Delta + tuneAD*rnormsAD(1);
    //DeltaNew = Delta + tuneD*rnormsAD(1);
    DeltaNew = Delta + exp(tuneD)*rnormsAD(1);
  }else{
    uu = arma::randu();
    DeltaNew = Delta + tuneD*tan(PI*(uu-0.5));
  }
  
  for(int tt=1; tt<dims(2); tt++){
    for(int i = 0; i < dims(0); i++){
      for(int j = 0; j < dims(0); j++){
        if(i != j){
          dx = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(j).col(tt),2);
          AccProb += Y.slice(tt)(i,j)*(DeltaNew-Delta)*Y.slice(tt-1)(i,j);
          AccProb += log(1+exp(AlphaNew+Delta*Y.slice(tt-1)(i,j)-dx));
          AccProb -= log(1+exp(AlphaNew+DeltaNew*Y.slice(tt-1)(i,j)-dx));
        }
      }
    }
  }
  
  AccProb += -0.5*(DeltaNew-nuDelta)*(DeltaNew-nuDelta)/xiDelta;
  AccProb -= -0.5*(Delta-nuDelta)*(Delta-nuDelta)/xiDelta;
  
  uu= arma::randu();
  if(uu<exp(AccProb)){
    AccRate(1) = 1;
  }else{
    DeltaNew = Delta;
  }

  // update tuneD
  vector = {1, exp(AccProb)};
  tuneD = tuneD + 1/pow(it,0.8)*(arma::min(vector)-0.234);

  return Rcpp::List::create(Xnew, AlphaNew, DeltaNew, AccRate, tuneA, tuneD, tuneXs);
  
  ', plugin="RcppArmadillo")


c.t2s2Parms <-  cxxfunction(
  signature(DATA="numeric",DIMS="integer",THETAT="numeric",
            THETAS="numeric",PHIT="numeric",PHIS="numeric", Yy="numeric", 
            WW="numeric", PIX="numeric"),
  body='
// Dims is c(n,p,TT,K)
  Rcpp::IntegerVector dims(DIMS);
  Rcpp::NumericVector Xvec(DATA);
  arma::cube X(Xvec.begin(),dims[1],dims[2],dims[0]);
  double thetaT=Rcpp::as<double>(THETAT);
  double phiT=Rcpp::as<double>(PHIT);
  double thetaS=Rcpp::as<double>(THETAS);
  double phiS=Rcpp::as<double>(PHIS);
  double shapeT=0, scaleT=0, shapeS=0, scaleS=0;
  
  Rcpp::NumericVector YY(Yy);
  arma::cube Y(YY.begin(),dims[0],dims[0],dims[2]); // n X n X TT
  arma::mat insides = arma::zeros(1,1);
  arma::colvec ww = Rcpp::as<arma::colvec>(WW);
  double count1=0;
  double count2=0;
  arma::colvec avg1 = arma::zeros(dims[1],1);
  arma::colvec avg2 = arma::zeros(dims[1],1);
  arma::colvec Ai_1 = arma::zeros(dims[1],1);
  arma::colvec Ai_2 = arma::zeros(dims[1],1);
  arma::colvec pi = Rcpp::as<arma::colvec>(PIX);

  Rcpp::NumericVector A1vec(dims[1]*dims[2]*dims[0]);
  Rcpp::NumericVector A2vec(dims[1]*dims[2]*dims[0]);
  arma::cube Ai_1s(A1vec.begin(),dims[1],dims[2],dims[0]); // p x (TT) x n
  arma::cube Ai_2s(A2vec.begin(),dims[1],dims[2],dims[0]);

  // ---------------------------------------
  
  shapeT = thetaT + 0.5*dims(0)*dims(1);
  scaleT = phiT;
  shapeS = thetaS + 0.5*dims(0)*dims(1)*(dims(2)-1);
  scaleS = phiS;
  for(int i =0; i<dims(0); i++){
    insides = 0.5*trans(X.slice(i).col(0))*(X.slice(i).col(0));
    scaleT += insides(0,0);
    for(int tt=1;tt<dims(2);tt++){
      avg1 = arma::zeros(dims[1],1);
      avg2 = arma::zeros(dims[1],1);
      count1 = 0;
      count2 = 0;
      for(int j=0; j<dims[0]; j++){
        if(Y.slice(tt-1)(i, j)==1){
          if(pi(i) == pi(j)){
            avg1 += X.slice(j).col(tt-1);
            count1 += 1;
          }else{
            avg2 += X.slice(j).col(tt-1);
            count2 += 1;
          }
        } 
      }
      if(count1>0){
        Ai_1 = avg1/count1 - X.slice(i).col(tt-1);
      }else{
        Ai_1 = arma::zeros(dims[1],1);
      }
      if(count2>0){
        Ai_2 = avg2/count2 - X.slice(i).col(tt-1);
      }else{
        Ai_2 = arma::zeros(dims[1],1);
      }
      Ai_1s.slice(i).col(tt) = Ai_1;
      Ai_2s.slice(i).col(tt) = Ai_2;
      insides = 0.5*trans(X.slice(i).col(tt)-X.slice(i).col(tt-1)-ww(pi(i)-1)*Ai_1-ww(2)*Ai_2)*
      (X.slice(i).col(tt)-X.slice(i).col(tt-1)-ww(pi(i)-1)*Ai_1-ww(2)*Ai_2);
      scaleS += insides(0,0);
    }
  }
  
  return Rcpp::List::create(shapeT, scaleT, shapeS, scaleS, Ai_1s, Ai_2s);
  
  ', plugin="RcppArmadillo")

c.update2Gamma_norm <-  cxxfunction(
  signature(Data="numeric",DIMS="integer",Yy="numeric", TUNEG1="numeric", TUNEG2="numeric",
            WWOld="numeric", s2X="numeric", NUGamma1='numeric', XIGamma1='numeric', 
            NUGamma2='numeric', XIGamma2='numeric', IT='integer', RNGX='numeric', 
            PIX="numeric", A1='numeric', A2='numeric'),
  body='
  /*Dims is c(n,p,TT,K)
  Z nxT; X pxTxn; Y pxpxT; 
  */
  Rcpp::IntegerVector dims(DIMS);
  
  Rcpp::NumericVector Xvec(Data);
  arma::cube X(Xvec.begin(),dims[1],dims[2],dims[0]);
  
  Rcpp::NumericVector YY(Yy);
  arma::cube Y(YY.begin(),dims[0],dims[0],dims[2]);
  
  arma::colvec wwOld = Rcpp::as<arma::colvec>(WWOld);
  // arma::colvec wwNew = Rcpp::as<arma::colvec>(WWNew);
  arma::colvec wwNew = arma::zeros(dims[0], 1); // n  
  // arma::colvec wwNew = arma::zeros(2, 1) // for gamma1, gamma2

  double tuneG1 = Rcpp::as<double>(TUNEG1);
  double tuneG2 = Rcpp::as<double>(TUNEG2);
  double nuGamma1 = Rcpp::as<double>(NUGamma1);
  double xiGamma1 = Rcpp::as<double>(XIGamma1);
  double nuGamma2 = Rcpp::as<double>(NUGamma2);
  double xiGamma2 = Rcpp::as<double>(XIGamma2);
  int it = Rcpp::as<int>(IT); // current iteration
  arma::colvec RNGs = Rcpp::as<arma::colvec>(RNGX);

  double AccProb =0,dx=0, uu=0;
  arma::mat insides = arma::zeros(1,1);
  arma::colvec muit = arma::zeros(dims[1],1);
  arma::colvec AccRate = arma::zeros(2, 1);
  arma::colvec Ai_1 = arma::zeros(dims[1],1);
  arma::colvec Ai_2 = arma::zeros(dims[1],1);
  arma::colvec pi = Rcpp::as<arma::colvec>(PIX);
  double s2 = Rcpp::as<double>(s2X);

  Rcpp::NumericVector A1vec(A1);
  Rcpp::NumericVector A2vec(A2);
  arma::cube Ai_1s(A1vec.begin(),dims[1],dims[2],dims[0]);
  arma::cube Ai_2s(A2vec.begin(),dims[1],dims[2],dims[0]);

  /*------------Update gamma1---------------------------*/

  // propose wwNew(0), i.e. new gamma1 by adding normal perturbation
  wwNew(0) = wwOld(0) + exp(tuneG1)*RNGs(0);

  // calculate acceptance prob for wwNew(0) with wwOld(0)
  for(int tt=1;tt<dims(2);tt++){
    for(int i = 0; i < dims(0); i++){
      Ai_1 = Ai_1s.slice(i).col(tt);
      Ai_2 = Ai_2s.slice(i).col(tt);
      muit = X.slice(i).col(tt-1)+wwNew(0)*Ai_1+wwOld(1)*Ai_2;
      insides = trans(X.slice(i).col(tt)-muit)*(X.slice(i).col(tt)-muit)/s2;
      AccProb += -0.5*insides(0,0);
      muit = X.slice(i).col(tt-1)+wwOld(0)*Ai_1 + wwOld(1)*Ai_2;
      insides = trans(X.slice(i).col(tt)-muit)*(X.slice(i).col(tt)-muit)/s2;
      AccProb -= -0.5*insides(0,0);
    }
  }
  AccProb += -0.5*(wwNew(0)-nuGamma1)*(wwNew(0)-nuGamma1)/xiGamma1;
  AccProb -= -0.5*(wwOld(0)-nuGamma1)*(wwOld(0)-nuGamma1)/xiGamma1;
  
  uu= arma::randu();
  if(uu<exp(AccProb)){
    AccRate(0) = 1;
  }else{
    wwNew(0) = wwOld(0);
  }
  // update tuneG1
  arma::vec vector = {1, exp(AccProb)};
  tuneG1 = tuneG1 + 1/pow(it,0.8)*(arma::min(vector)-0.234);

  /*------------Update gamma2---------------------------*/
  // propose wwNew(1), i.e. new gamma2 by adding normal perturbation
  wwNew(1) = wwOld(1) + exp(tuneG2)*RNGs(1);
  
  // calculate acceptance prob for wwNew(1) with wwOld(1)
  for(int tt=1;tt<dims(2);tt++){
    for(int i = 0; i < dims(0); i++){
      Ai_1 = Ai_1s.slice(i).col(tt);
      Ai_2 = Ai_2s.slice(i).col(tt);
      muit = X.slice(i).col(tt-1)+wwNew(0)*Ai_1+wwNew(1)*Ai_2;
      insides = trans(X.slice(i).col(tt)-muit)*(X.slice(i).col(tt)-muit)/s2;
      AccProb += -0.5*insides(0,0);
      muit = X.slice(i).col(tt-1)+wwNew(0)*Ai_1 + wwOld(1)*Ai_2;
      insides = trans(X.slice(i).col(tt)-muit)*(X.slice(i).col(tt)-muit)/s2;
      AccProb -= -0.5*insides(0,0);
    }
  }
  AccProb += -0.5*(wwNew(1)-nuGamma2)*(wwNew(1)-nuGamma2)/xiGamma2;
  AccProb -= -0.5*(wwOld(1)-nuGamma2)*(wwOld(1)-nuGamma2)/xiGamma2;
  
  uu= arma::randu();
  if(uu<exp(AccProb)){
    AccRate(1) = 1;
  }else{
    wwNew(1) = wwOld(1);
  }

  // update tuneG2 
  vector = {1, exp(AccProb)};
  tuneG2 = tuneG2 + 1/pow(it,0.8)*(arma::min(vector)-0.234);
  
  return Rcpp::List::create(wwNew, AccRate, tuneG1, tuneG2);
  
  ', plugin="RcppArmadillo")

