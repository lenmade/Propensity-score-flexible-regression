#######################################################################################################################################
##
## Sample R code to implement the method Flexible Regression Method in Propensity Score Analysis 
## Note: due to the version of associated R packages, warning message may come up, but it's been tested not impact analysis result.
#######################################################################################################################################
rm(list=ls(all=T));
setwd( '/Users/huzhangmao/Dropbox/PS_Reg/psReg_simulation_4Li copy/nonparam' );

source( 'func_psReg.R' );
###########################################################################################################
## sample data
###########################################################################################################
data <- read.csv( "sample_data.csv" );

#############################################################################################################
## Following setting is for the sample data with 
##    1. Z is the treatment indicator, Z=0 or 1;
##    2. Y is the continuous outcome;
##    3. Variates: X0=1, X1 and X2 are contious variable, X3 and X4 are dichotomous variables 
#############################################################################################################
setup <- list( weight = 'MW',  # Specify propensity score weighting method
               sp.degree = 1,  # Specify the degree of spline function
               nknots = 15 );  # Specify number of knots in spline function

form.ps <- as.formula( 'Z ~ X1 + X2 + X3 + X4' );
Xname <- c( 'X0', 'X1', 'X2', 'X3', 'X4' );
Tname <- 'Z';
Oname <- 'Y';
Fname <- c( 'X0', sapply( c(1 : setup$sp.degree), function(x) paste( 'ltps', x, sep='') ) );
R0name <- set.basis.name( nbasis = setup$nknots, prefix = 'S0_' );
R1name <- set.basis.name( nbasis = setup$nknots, prefix = 'S1_' );
Mname <- 'mu';
Fcoef0.name <- sapply( c(0 : setup$sp.degree), function(x) paste( 'alpha0', x, sep='' ) );
Fcoef1.name <- sapply( c(0 : setup$sp.degree), function(x) paste( 'alpha1', x, sep='' ) );
Rcoef0.name <- sapply( c(1 : setup$nknots), function(x) paste( 'u0', x, sep='_' ) );
Rcoef1.name <- sapply( c(1 : setup$nknots), function(x) paste( 'u1', x, sep='_' ) );
pscoef.name <- sapply( c(1 : length(Xname)), function(x) paste( 'beta', x-1, sep='' ) );

#############################################################################################################
## Prepare data to calculate penalty parameter for penalized regression
#############################################################################################################
ps.hat <- ps.model(dat=data, form=form.ps )$ps.hat; 
data$ltps1 <- logit( ps.hat );
# control group
Z <- 0;
index0 <- which( data[, Tname] == Z );  # index for the control group
Knots0 <- quantile( data[ index0, ]$ltps1, probs=seq( 0, 1, length=setup$nknots+2 ) );
Knots0 <- Knots0[ -c( 1, length(Knots0) ) ];
S0 <- as.data.frame( t( sapply( data[ index0, ]$ltps1, calc.trunc.basis, Knots=Knots0, dgr=setup$sp.degree ) ) );
colnames( S0 ) <- c( Fname, set.basis.name( nbasis = setup$nknots, prefix = 'S_' ) );
dat0 <- cbind( data[ index0, c( Oname, Tname) ], S0 );
#write.csv( dat0, file=paste( "data_Z", Z, ".csv", sep="" ), row.names=FALSE );  # data to be used for SAS PROC GLIMMIX

# treated group
Z <- 1;
index1 <- which( data[, Tname] == Z );  # index for the treated group
Knots1 <- quantile( data[ index1, ]$ltps1, probs=seq( 0, 1, length=setup$nknots+2 ) );
Knots1 <- Knots1[ -c( 1, length(Knots1) ) ];    
S1 <- as.data.frame( t( sapply( data[ index1, ]$ltps1, calc.trunc.basis, Knots=Knots1, dgr=setup$sp.degree ) ) );
colnames( S1 ) <- c( Fname, set.basis.name( nbasis = setup$nknots, prefix = 'S_' ) );
dat1 <- cbind( data[ index1, c( Oname, Tname) ], S1 );
#write.csv( dat1, file=paste( "data_Z", Z, ".csv", sep="" ), row.names=FALSE );  # data to be used for SAS PROC GLIMMIX

#############################################################################################################
## Implementation of flexible regression method
#############################################################################################################
pelty0 <- 30.9;  # penalty parameter for penalized spline, it is data specific.
pelty1 <- 29.2;  # penalty parameter for penalized spline, it is data specific.
var.verbose <- T;  # calculate variance?
res.sp <- calc.unif.sp( dat=data, weight=setup$weight, dgr=setup$sp.degree, nknots=setup$nknots, Knots0=NULL, Knots1=NULL, 
                        pelty0=pelty0, pelty1=pelty1, ps.form=form.ps, pscoef.name=pscoef.name, Fcoef0.name=Fcoef0.name,
                        Fcoef1.name=Fcoef1.name, Rcoef0.name=Rcoef0.name, Rcoef1.name=Rcoef1.name, Fname=Fname, R0name=R0name, 
                        R1name=R1name, O.name=Oname, M.name=Mname, X.name=Xname, var.verbose=var.verbose );
est.sp <- res.sp$est.sp;
std.sp <- res.sp$std.sp;
