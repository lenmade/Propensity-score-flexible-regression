###############################################################
##
## Supporting function for PS unified framework
##
## WIN 10, R 3.2.4
## April 22, 2016, Huzhahng Ma0
##
###############################################################
library(mvtnorm);  # For function rmvnorm;
library(boot);  # For function inv.logit
library(MatchIt);  # matchit

calc.trunc.basis <- function(xi, Knots, dgr) {
  # calculate truncated power basis
  #
  # Args
  #   xi: data point at which basis is to be calculated
  #   Knots: spline knots
  #   dgr: spline degree
  # Return
  #   A vector of spline basis
  
  tmp1 <- sapply( c(0 : dgr), function(x) xi^x );
  tmp2 <- xi-Knots;
  tmp2 <- tmp2^dgr*(tmp2 >= 0);
  ans <- c( tmp1, tmp2 );
  
  ans;
}


set.basis.name <- function (nbasis, prefix) {
  # Set basis names 
  #
  # Args:
  # n.basis: number of basis used
  #
  # Returns:
  # A character vector of basis names with defined prefix
  
  ans <- vector(mode = "character", length = 0)
  for (i in 1 : nbasis) {
    ans <- c(ans, paste(prefix, i, sep = ""))
  }
  ans
}

calc.unif.sp <- function( dat, weight, dgr, nknots, Knots0=NULL, Knots1=NULL, 
                          pelty0, pelty1, ps.form, pscoef.name, Fcoef0.name, Fcoef1.name,
                          Rcoef0.name, Rcoef1.name, Fname, R0name, R1name, O.name, M.name, 
                          X.name, var.verbose=T ) {
  # main function to calculate causal effect within unified framework using penalized spline
  #
  # Args
  #   dat: data used
  #   weight: PS weight
  #   dgr: spline degree
  #   nknots: number of knots, must be specified when Knots is NULL
  #   Knots0, Knots1: user provided spline knots for control and treated group, respectively, by default is NULL
  #   pelty0, pelty1: penalty parameters for treatment groups
  #   Q0.mat, Q1.mat: matrice associated with  
  #   ps.form: PS model
  #   pscoef.name: names of PS coefficients
  #   Fcoef0.name, Fcoef1.name: coef names fixed effect for control and treated groups, respectively
  #   Rcoef0.name, Rcoef1.name: coef names random effect for control and treated groups, respectively
  #   Fname: Fxied covariate name
  #   R0name, R1name: covariate names in control and treated groups in the S matrix
  #   O.name: outcome name
  #   M.name: weighted average treatment effect
  #   X.name: covariate names in PS model, including the intercept
  #   var.verbose: calcuate variance. TRUE by default
  
  Q1.mat <- Q0.mat <- cbind( matrix( 0, ncol=length(Fname), nrow=(length(Fname) + length(R1name)) ),
                             rbind( matrix( 0, ncol=length(R1name), nrow=length(Fname) ), diag( length(R1name)) ) );
  
  # estimated PS
  out.ps <- ps.model( dat=dat, as.formula(form.ps) );
  dat$ps <- out.ps$ps.hat;  # estimated ps
  dat$ltps1 <- logit( dat$ps );
  omega <- sapply( dat$ps, calc.omega, weight = weight ); # generic omega function
  beta.hat <- coef( out.ps$fm );  # PS model coefficient
  names( beta.hat ) <- pscoef.name;
  
  # control group
  Z <- 0;
  index0 <- which( dat[, Tname] == Z );  # index for the control group
  if ( is.null(Knots0) ) {
    Knots0 <- quantile( dat[ index0, ]$ltps1, probs=seq( 0, 1, length=nknots+2 ) );
    Knots0 <- Knots0[ -c( 1, length(Knots0) ) ];    
  }
  S0 <- as.data.frame( t( sapply( dat$ltps1, calc.trunc.basis, Knots=Knots0, dgr=dgr ) ) );
  colnames( S0 ) <- c( Fname, R0name );
  X0 <- as.matrix( S0[ index0, ] );
  coef0 <- solve( t(X0) %*% X0 + pelty0*Q0.mat ) %*% t(X0) %*% dat[ index0, O.name ];
  names( coef0 ) <- c( Fcoef0.name, Rcoef0.name );
  alpha0.hat <- as.matrix( S0 ) %*% coef0;
  
  # treated group
  Z <- 1;
  index1 <- which( dat[, Tname] == Z );  # index for the treated group
  if ( is.null(Knots1) ) {
    Knots1 <- quantile( dat[ index1, ]$ltps1, probs=seq( 0, 1, length=nknots+2 ) );
    Knots1 <- Knots1[ -c( 1, length(Knots1) ) ];    
  }
  S1 <- as.data.frame( t( sapply( dat$ltps1, calc.trunc.basis, Knots=Knots1, dgr=dgr ) ) );
  colnames( S1 ) <- c( Fname, R1name );
  X1 <- as.matrix( S1[ index1, ] );
  coef1 <- solve( t(X1) %*% X1 + pelty1*Q1.mat ) %*% t(X1) %*% dat[ index1, O.name ];
  names( coef1 ) <- c( Fcoef1.name, Rcoef1.name );
  alpha1.hat <- as.matrix( S1 ) %*% coef1;
  
  est.indv <- alpha1.hat - alpha0.hat;
  est.sp <- sum( omega * est.indv ) / sum( omega );
  names( est.sp ) <- M.name;
  
  if ( var.verbose ) {
    param.hat <- c( est.sp, coef1, coef0, beta.hat );  # estimated beta
    Q.hat <- matrix(0, nrow=1 + 2*(1 + dgr + nknots) + length(beta.hat), ncol=1 + 2*(1 + dgr + nknots) + length(beta.hat) );
    diag(Q.hat) <- c( 0, pelty1*c( rep(0, 1 + dgr), rep(1, nknots) ), 
                      pelty0*c( rep(0, 1 + dgr), rep(1, nknots) ),
                      rep( 0, length(beta.hat) ) );
    
    #print( Knots1 );
    #print( Knots0 );
    
    var.coef <- calc.var( dat=dat, param=param.hat, weight=weight, Knots0=Knots0, Knots1=Knots1, 
                        dgr=dgr, Q.mat=Q.hat, O.name=O.name, X.name=X.name, Tname=Tname, 
                        mu.name=M.name, Fcoef0.name=Fcoef0.name, Fcoef1.name=Fcoef1.name,
                        Rcoef0.name=Rcoef0.name, Rcoef1.name=Rcoef1.name, ps.coef.name=pscoef.name );
    
    std.sp <- sqrt( var.coef[1, 1] );
    
    return( list(est.sp=est.sp, std.sp=std.sp, coef1=coef1, coef0=coef0, var.coef=var.coef) );
  }
  return( list(est.sp=est.sp, coef1=coef1, coef0=coef0) );
} 

calc.var <- function( dat, param, weight, Knots0, Knots1, dgr, Q.mat, O.name, X.name, Tname, 
                      mu.name, Fcoef0.name, Fcoef1.name, Rcoef0.name, Rcoef1.name, ps.coef.name) {
  # calcuate the var-covar adjusting for PS estimation
  #
  # Args
  #   dat: data 
  #   param: a named vector of all parameters
  #   weight: name of propensity score weight
  #   Knots0, Knots1: knots for contorl and treated group, respectively
  #   dgr: degree of truncated power basis
  #   Q.mat: Q matrix in estimating equation
  #   O.name: outcome name
  #   X.name: covariate name in propensity score model
  #   Tname: treatment indicator name
  #   mu.name: name for estimated treatment effect
  #   Fcoef0.name, Fcoef1.name: name of fixed effect coefficient for control and treated groups, respectively
  #   Rcoef0.name, Rcoef1.name: name of random effect coefficient for control and treated groups, respectively
  #   ps.coef.name: name of propensity score model coefficient
  #
  # Returns
  #   Var-covar matrix of all parameters
  
  Amat <- Bmat <- 0;
  n <- nrow(dat);
  for ( i in 1 : n ) {
    this.dat <- dat[ i, ];
    Yi <- as.numeric( this.dat[ , O.name ] );
    Zi <- as.numeric( this.dat[ , Tname] );
    Xi <- as.numeric( this.dat[ , X.name] );
    
    phi <- calc.phi( param=param, Xi=Xi, Zi=Zi, Yi=Yi, dgr=dgr, weight=weight, Knots1=Knots1, Knots0=Knots0,
                     Q.mat=Q.mat, ndata=n, mu.name=mu.name, Fcoef1.name=Fcoef1.name, Fcoef0.name=Fcoef0.name, 
                     Rcoef1.name=Rcoef1.name, Rcoef0.name=Rcoef0.name, ps.coef.name=ps.coef.name );
    Bmat <- Bmat + outer( phi, phi );
    
    phi.deriv <- jacobian( calc.phi, param, Xi=Xi, Zi=Zi, Yi=Yi, dgr=dgr, weight=weight, Knots1=Knots1, Knots0=Knots0,
                           Q.mat=Q.mat, ndata=n, mu.name=mu.name, Fcoef1.name=Fcoef1.name, Fcoef0.name=Fcoef0.name, 
                           Rcoef1.name=Rcoef1.name, Rcoef0.name=Rcoef0.name, ps.coef.name=ps.coef.name );
    Amat <- Amat + phi.deriv;
  }
  Amat <- Amat/n;
  Bmat <- Bmat/n;
  Amat.inv <- solve(Amat);
  var.mat <- ( Amat.inv %*% Bmat %*% t(Amat.inv)) / n;
  colnames( var.mat ) <- rownames( var.mat ) <- c(mu.name, Fcoef1.name, Rcoef1.name, Fcoef0.name, Rcoef0.name, ps.coef.name );
  
  return( var.mat );
}


calc.phi <- function( param, Xi, Zi, Yi, dgr, weight, Knots1, Knots0, Q.mat, ndata, mu.name, 
                      Fcoef1.name, Fcoef0.name, Rcoef1.name, Rcoef0.name, ps.coef.name ) {
  
  # calculate the phi function in the estimating equation for balance diagnosis
  #
  # Args
  #   param: a named vector of parameters
  #   Xi: covariate in PS model
  #   Zi: treatment indicator
  #   Yi: the outcome
  #   dgr: degree of spline
  #   weight: PS weight type
  #   Knots0, Knots1: spline knots for treated and control groups, respectively
  #   Q.mat: the combinded penalyt matrix for all parameter
  #   ndata: number of subject, adjust for penalty term in estimating equation
  #   mu.name: name for estimated treatment effect
  #   Fcoef0.name, Fcoef1.name: name of fixed effect coefficient for control and treated groups, respectively
  #   Rcoef0.name, Rcoef1.name: name of random effect coefficient for control and treated groups, respectively
  #   ps.coef.name: name of propensity score model coefficient
  # 
  # Retruns
  #   A vector of estimating equation for each subject
  
  mu <- param[ mu.name ];
  w1 <- param[ c(Fcoef1.name, Rcoef1.name) ];
  w0 <- param[ c(Fcoef0.name, Rcoef0.name) ];
  beta <- param[ ps.coef.name ];
  
  ei <- calc.ps.Xbeta( Xmat = Xi, beta = beta );
  omegai <- calc.omega( ps = ei, weight = weight );
  ltei <- logit( ei );
  
  G0i <- calc.trunc.basis(xi = ltei, Knots = Knots0, dgr=dgr);
  G1i <- calc.trunc.basis(xi = ltei, Knots = Knots1, dgr=dgr);

  alpha1i <- G1i %*% w1;
  alpha0i <- G0i %*% w0;
  est <- alpha1i -  alpha0i;
  
  ans1 <- omegai*( est - mu );
  ans2 <- -2*Zi*(Yi - t(G1i) %*% w1) * G1i;
  ans3 <- -2*(1-Zi)*(Yi - t(G0i) %*% w0) * G0i;
  ans4 <- (Zi-ei)*Xi;
  ans.penalty <- as.numeric( 2*Q.mat %*% param / ndata );
  ans <- c( ans1, ans2, ans3, ans4) + as.numeric( ans.penalty );
  
  return( ans );
}

ps.model <- function(dat, form) {
  # ps.model to calculate propensity score
  # 
  # Args:
  #   dat: the data from which estimated ps to be calculated
  #   form: the formula used in propensity score model
  #
  # Return:
  #   Estimated propensity score and glm fitting
  
  fm <- glm(form, data = dat, family = binomial(link = "logit"));
  ps.hat <- as.numeric(predict(fm, newdata = dat, type = "response"));
  ps.hat <- pmin(pmax(0.000001, ps.hat), 0.999999);  # ps.hat cannot be exactly 0 or 1
  return(list(ps.hat = ps.hat, fm = fm));
}

calc.omega <- function(ps, weight, delta=0.002, K=4) {
  # Calculate omega for each weighting method
  #
  # Args:
  #   ps: estimated proopensity score,
  #   weight: weighting method.
  #   K: wighting coefficient for trapezoidal weighting
  #
  # Return:
  #   A vector of length(ps)
  
  ans <- 0;
  if (weight == "IPW") {
    ans <- 1;
  } else if (weight == "MW") {
    ans <- calc.omega.MW(ps, delta);
  } else if (weight == "ATT") {
    ans <- ps;
  } else if (weight == "ATC") {
    ans <- 1-ps;
  } else if (weight == "OVERLAP") {
    ans <- 4*ps*(1-ps);
  } else if (weight == "TRAPEZOIDAL") {
    ans <- calc.omega.trapzd(ps=ps, delta=delta, K=K);
  } else {
    stop("Error in calc.omega: weight method does not exist!");
  }
  
  ans;
}


calc.omega.MW <- function (ps, delta) {
  # Calculate omega value for MW method
  #
  # Args:
  # ps: propensity score, scalar
  # delta: closeness to non-differential point
  #
  # Return:
  # Value of omega, scalar
  
  ans <- 0;
  if (ps <= 0.5 - delta) {
    ans <- 2*ps;
  } else if (ps >= 0.5 + delta) {
    ans <- 2*(1 - ps);
  } else {
    ans <- approx.omega.MW(ps, delta);
  }
  
  ans;
}


approx.omega.MW <- function (ps, delta) {
  # Approximate omega function of MW at non-differential point, eta(ps)
  #
  # Args:
  # ps: propensity score, scalar
  # delta: closeness to non-differential point
  #
  # Return:
  # Approximated omega at non-differential point, scalar
  
  A <- solve.A.MW(delta);
  ans <- rowVec(A) %*% c(1, ps, ps^2, ps^3);
  
  ans;
}

solve.A.MW <- function(delta) {
  # Get coefficients for approximated cubic polynomial for omega (eta(e)) function of MW
  #
  # Args:
  #   delta: pre-defined closeness to midpiece of ps.
  #
  # Return
  #   A vecotor of solved coefficients
  
  if ( delta < 0.00001 ) { 
    stop("*** ERROR in solve.a: delta too small ***");
  }
  
  tmp1 <- 0.5 - delta;
  tmp2 <- 0.5 + delta;
  
  D <- matrix(c(1, tmp1, tmp1^2, tmp1^3,
                0, 1, 2*tmp1, 3*tmp1^2,
                1, tmp2, tmp2^2, tmp2^3,
                0, 1, 2*tmp2, 3*tmp2^2),
              ncol = 4, nrow = 4, byrow = TRUE);
  C <- 2*c(tmp1, 1, tmp1, -1);
  A <- solve(D) %*% C;  # coefficients of cubic polynomial
  
  A;
}

calc.omega.trapzd <- function (ps, delta, K) {
  # Calculate omega value for trapezoidal weighting method
  #
  # Args:
  # ps: propensity score, scalar
  # delta: closeness to non-differential point
  # K: trapezoidal weighting coffecient
  #
  # Return:
  # Value of omega, scalar
  
  ans <- 0;
  if ( (0 < ps) & (ps <= 1/K - delta) ) {
    ans <- K*ps;
  } else if ( (1/K + delta <= ps) & (ps <= 1 - 1/K - delta) ) {
    ans <- 1;
  } else if ( (1 - 1/K + delta <= ps) & (ps < 1) ) {
    ans <- K*(1 - ps);
  } else {
    ans <- approx.omega.trapzd(ps, delta, K);
  }
  
  ans;
}

approx.omega.trapzd <- function (ps, delta, K) {
  # Approximate omega function of trapezoidal weight at non-differential point, eta(ps)
  #
  # Args:
  # ps: propensity score, scalar
  # delta: closeness to non-differential point
  # K: trapezoidal weight coefficient
  #
  # Return:
  # Approximated omega at non-differential point, scalar
  
  A <- 0;
  if ( (1/K - delta < ps) & (ps < 1/K + delta) ) {
    A <- solve.A.trapzd1st(delta, K);
  } else {
    A <- solve.A.trapzd2nd(delta, K);
  }  
  ans <- rowVec(A) %*% c(1, ps, ps^2, ps^3);
  
  ans;
}

solve.A.trapzd1st <- function (delta = delta, K = K) {
  # Get coefficients for approximated cubic polynomial for omega (eta(e)) function of 
  # trapezoidal weighting at the first non-differential pivot
  #
  # Args:
  #   delta: pre-defined closeness to midpiece of ps.
  #   K: coefficient of trapezoidal weight
  #
  # Return
  #   A vecotor solved coefficients
  
  if ( delta < 0.00001 ) { 
    stop("*** ERROR in solve.a: delta too small ***");
  }
  
  tmp1 <- 1/K - delta;
  tmp2 <- 1/K + delta;
  
  D <- matrix(c(1, tmp1, tmp1^2, tmp1^3,
                0, 1, 2*tmp1, 3*tmp1^2,
                1, tmp2, tmp2^2, tmp2^3,
                0, 1, 2*tmp2, 3*tmp2^2),
              ncol = 4, nrow = 4, byrow = TRUE);
  C <- 2*c(K*tmp1, K, 1, 0);
  A <- solve(D) %*% C;  # coefficients of cubic polynomial
  
  A;
}

solve.A.trapzd2nd <- function (delta = delta, K = K) {
  # Get coefficients for approximated cubic polynomial for omega (eta(e)) function of 
  # trapezoidal weighting at the second non-differential pivot
  #
  # Args:
  #   delta: pre-defined closeness to midpiece of ps.
  #   K: coefficient of trapezoidal weight
  #
  # Return
  #   A vecotor solved coefficients
  
  if ( delta < 0.00001 ) { 
    stop("*** ERROR in solve.a: delta too small ***");
  }
  
  tmp1 <- 1 - 1/K - delta;
  tmp2 <- 1 - 1/K + delta;
  
  D <- matrix(c(1, tmp1, tmp1^2, tmp1^3,
                0, 1, 2*tmp1, 3*tmp1^2,
                1, tmp2, tmp2^2, tmp2^3,
                0, 1, 2*tmp2, 3*tmp2^2),
              ncol = 4, nrow = 4, byrow = TRUE);
  C <- 2*c(1, 0, K*(1/K - delta), -K);
  A <- solve(D) %*% C;  # coefficients of cubic polynomial
  
  A;
}

calc.ps.Xbeta <- function(Xmat, beta) {
  # Calculate the propensity score, e(X, beta)
  #
  # Args:
  #   X: a vector or a matrix with each column being X_i
  #   beta: the coefficient, of the same length as the nrow(X)
  #
  # Return:
  #   A scalar of matching weight
  
  Xmat <- as.matrix(Xmat);
  tmp <- as.numeric(rowVec(beta) %*% Xmat);
  tmp <- exp(tmp);
  names(tmp) <- NULL;
  
  return( tmp/(1 + tmp) );
}

rowVec <- function(x) {
  t(x);
}

colVec <- function(x) {
  t(t(x));
}
