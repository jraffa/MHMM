## Version 0.004
## This code is from Raffa JD and Dubin JA (2015) "Multivariate Longitudinal Data Analysis with Mixed Effects Hidden Markov Models."
## It contains two main functions to simulate multivariate response longitudinal data arises from a hidden Markov Model, 
## and implementation of a MCMC approach to estimate the parameters from such models.  For further details see:
## the original paper, supplementary material and example.R / README files.
## Jesse Raffa (jraffa@uw.edu): http://github.com/jraffa/MHMM
##
#######
#	genSim generates simulated MHMM dataset
#	code written for clarity, maybe not efficient
#	arguments:
# 	N: Number of study subjects (integer)
#	ni: Number of observation per subject (fixed for all subjects,, integer)
# 	hs: Number of hidden states (integer)
#	tau1: means for each hidden state for Poisson response (vector of length hs)
#	ti: number of 'days' for the Poisson response (fixed for all observation times, integer)
#   tau2: means for each  hidden state for Normally distributed response (vector of length hs)
#	residinvvar: inverse of the residual variance (numeric)
#	reSigma: random effects covariance matrix (2x2 matrix)
#	P0: hidden states transition probability matrix (hs x hs matrix) for placebo or baseline group
#	Pi0: initial probability vector (vector, length hs) for placebo or baseline group
#	w: random effects design matrix -- uses default as in paper (separate but correlated REs), requires R list of length 2, with each element containing a 2 x 1 row vector
#	rx: binary vector of treatment indicators of length N; if not specified assumes not treatment effect!
#	P1: hidden state tpm for treatment group (same as P0, requires rx to be specified); if left null, assumes no difference
#	Pi1: hidden state initial probability vector for treatment group (same as Pi0, requires rx to be specified); if left null, assumes no difference
#	fitRx: a logical vector of length two: fitRx[1] is whether to fit initial treatment probabilities separately for Rx and Control; fitRx[2] fits separate tpms
#
# 	Returns as a list everything you need to start computation: 
#	y1: generated response 1 (Poisson) as a matrix, 1 row per subject
#	y2: generated response 2 (Normal) as as matrix, 1 row per subject
#	mi: N x ni matrix of number of days in each week
#	N: Number of subjects, as specified
#	ni: number of observation times per subject
#	Z: generated hidden states (N x ni matrix) -- you don't have these normally
#	randomEffects: subject specific random effects (N x 2 matrix -- you don't have these normally either)
#	rx: treatment dummy variables
#	pars: parameter vector used to simulate data (useful for thinking specifying possible starting values)
#######
mhmmversion <- 0.004;

genSimpois <- function(N=354,ni=6,hs=3,tau1=c(-3,1,2.5),ti=7,tau2=c(1.2,1.2,1.1),
residinvvar=10,reSigma=matrix(c(1.5,0.2,0.2,0.1),nr=2),P0=matrix(c(0.85,0.6,0.1,0.1,0.2,0.1,0.05,0.2,0.8),nr=3),Pi0=c(0.6,0.25,0.15),
P1=matrix(c(0.85,0.6,0.1,0.1,0.2,0.1,0.05,0.2,0.8),nr=3),Pi1=c(0.6,0.25,0.15),
rx.=NULL,fitRx=c(FALSE,FALSE),w=list(matrix(c(1,0),nr=1),matrix(c(0,1),nr=1)))  {
	Res <- matrix(NA,nr=N,nc=dim(reSigma)[1]);
	reSigmachol <- chol(reSigma);
	Z <- matrix(NA,nr=N,nc=ni);
	
	y1 <- matrix(NA,nr=N,nc=ni);
	y2 <- matrix(NA,nr=N,nc=ni);
	for(i in 1:N) {
		if(fitRx[1] & !is.null(rx.) & !is.null(Pi1)) {
			if(rx.[i]==1) {
				Z[i,1] <- sample(1:hs,1,replace=TRUE,Pi1); #intial hidden states
			} else {
				Z[i,1] <- sample(1:hs,1,replace=TRUE,Pi0);
			}
		} else {
			Z[i,1] <- sample(1:hs,1,replace=TRUE,Pi0); #intial hidden states
		}
		Res[i,] <- t(reSigmachol)%*%rnorm(dim(reSigma)[1])
		tmpneta1 <- sum(tau1[1:Z[i,1]]) + w[[1]]%*%Res[i,] + log(ti);
		tmpneta2 <- sum(tau2[1:Z[i,1]]) + w[[2]]%*%Res[i,];
		y1[i,1] <- rpois(1,exp(tmpneta1))
		y2[i,1] <- rnorm(1,tmpneta2,sqrt(1/residinvvar));
		for(j in 2:ni) {
			if((fitRx[2] & !is.null(rx.) & !is.null(P1))) {
				if(rx.[i]==1) {
					Z[i,j] <- sample(1:hs,1,replace=TRUE,P1[Z[i,j-1],]);
				} else {
					Z[i,j] <- sample(1:hs,1,replace=TRUE,P0[Z[i,j-1],]);
				}
			} else {
				Z[i,j] <- sample(1:hs,1,replace=TRUE,P0[Z[i,j-1],]);
			}
			tmpneta1 <- sum(tau1[1:Z[i,j]]) + w[[1]]%*%Res[i,] + log(ti);
			tmpneta2 <- sum(tau2[1:Z[i,j]]) + w[[2]]%*%Res[i,];
			y1[i,j] <- rpois(1,exp(tmpneta1));
			y2[i,j] <- rnorm(1,tmpneta2,sqrt(1/residinvvar));
		}
	}
	
	if(fitRx[1] & fitRx[2]) { #Both Rx Effects
		mpars <- c(Pi0[1:(hs-1)],as.numeric(t(P0))[-seq(hs^2,1,-(hs+1))],tau1,tau2,as.numeric(reSigma),residinvvar,as.numeric(t(P1))[-seq(hs^2,1,-(hs+1))],Pi1[1:(hs-1)]);
	}
	if(!fitRx[2] & fitRx[1]) { #Only Initial Prob Rx Effect
		mpars <- c(Pi0[1:(hs-1)],as.numeric(t(P0))[-seq(hs^2,1,-(hs+1))],tau1,tau2,as.numeric(reSigma),residinvvar,Pi1[1:(hs-1)]);
	}
	if(fitRx[2] & !fitRx[1]) { #Only TPM Rx Effect
		mpars <- c(Pi0[1:(hs-1)],as.numeric(t(P0))[-seq(hs^2,1,-(hs+1))],tau1,tau2,as.numeric(reSigma),residinvvar,as.numeric(t(P1))[-seq(hs^2,1,-(hs+1))]);
	}
	if(!fitRx[1] & !fitRx[2])  { #No Rx Effect (Default)
		mpars <- c(Pi0[1:(hs-1)],as.numeric(t(P0))[-seq(hs^2,1,-(hs+1))],tau1,tau2,as.numeric(reSigma),residinvvar);
	}
	return(list(y1=y1,y2=y2,ti=matrix(rep(ti,N*ni),nr=N,nc=ni),N=N,ni=rep(ni,N),Z=Z,randomEffects=Res,rx=rx.,hs=hs,pars=mpars,fitRx=fitRx))
}

library(Rcpp); library(RcppArmadillo); library(inline); library(MCMCpack); library(mvtnorm);

## The following are C++ code that will be compiled on your system.  
## Note running this code requires R-tools in windows, or a compiler with other operating systems.
## It will spawn a compiler process in order to compile the C++ code.


## The following two blocks of code implement the sample procedure for the hidden states, Z_i

inccode.pois <- 'List  HMMlablpois(NumericVector y1, NumericVector y2, NumericVector states, NumericVector lmb, NumericMatrix p, NumericVector pi, NumericVector ni, NumericVector lti) {
//Code adapted from Zucchini and MacDonald (2009)
NumericVector x1(y1); NumericVector x2(y2);
NumericVector ltii(lti);
NumericVector nt(ni); NumericVector lambda(lmb); NumericMatrix gamma(p);
NumericVector mt(states); NumericVector delta(pi);
int n = x1.size();
int m = (lambda.size()-1)/2;
NumericMatrix lalpha(m,n);
NumericMatrix lbeta(m,n);
NumericMatrix allprobs(n,m);
arma::mat gma = Rcpp::as<arma::mat>(gamma);
for(int q=0; q<n; q++) {
	for(int r=0; r<m; r++) {
		allprobs(q,r) = ::exp(::Rf_dpois(x1[q],lambda[r]*::exp(ltii[q]),true) + ::Rf_dnorm4(x2[q],lambda[r+m],lambda[2*m],true));
		//Rcpp::Rcout << allprobs(q,r) << std::endl;
	}
	
}
NumericMatrix foo(n,m);
NumericVector lscale(n);
foo.row(0) = delta*allprobs.row(0);
NumericVector footmp(foo.row(0));
double sumfoo = std::accumulate(footmp.begin(),footmp.end(), 0.0);
//if(sumfoo==0 && n==1) {
if(sumfoo==0 || std::isnan(sumfoo) || std::isinf(sumfoo)) {
	Rcpp::Rcout << "underflow error (usually poor initial values)" << std::endl;
	sumfoo = ::pow(10,-100); //Not usually necessary when using this code in MCMC, but if you have problems, may be worth checking.
	footmp[m-1] = ::pow(10,-100);
}
	
foo.row(0) = footmp/sumfoo;
lscale[0] = ::log(sumfoo);
NumericVector logfootmp(foo.row(0));
std::transform(logfootmp.begin(), logfootmp.end(), logfootmp.begin(), ::log);
lalpha.column(0) = logfootmp+lscale[0];
if(n>1) {
for(int i=1; i<n; i++) {
	NumericVector foa(foo.row(i-1));
	arma::colvec fooa = Rcpp::as<arma::colvec>(foa);
	NumericVector ttt(m);
	ttt = arma::trans(fooa)*gma;
	for(int j=0; j<m; j++) {
		foo(i,j) = ttt[j]*allprobs(i,j);
	}
	NumericVector footmp(foo.row(i));
	double sumfoo = std::accumulate(footmp.begin(),footmp.end(), 0.0);
	lscale[i] = lscale[i-1] + ::log(sumfoo);
	foo.row(i) = footmp/sumfoo;
	NumericVector logfootmp(foo.row(i));
	std::transform(logfootmp.begin(), logfootmp.end(), logfootmp.begin(), ::log);
	lalpha.column(i) = logfootmp+lscale[i];
}
}
	
List ret; ret["lalpha"] = lalpha; ret["lbeta"] = 0;
return ret;
}
';




code5.pois <-'
NumericMatrix y1(y1m); NumericMatrix y2(y2m);
NumericVector beta(betaa); IntegerVector ni(nis); NumericMatrix con(ccc);
NumericMatrix ltii(lti);

LogicalVector fitRx(R_fitRx);
NumericVector rx(rx1);
Rcpp::List pidx(paridx);

arma::mat co = Rcpp::as<arma::mat>(con);
IntegerVector states(m);
int me = states[0];
arma::mat pf(arma::zeros(me,me));

int ctt = 0;
IntegerVector ppidx = pidx["P0"];

for(int i=0; i<me; i++) {
	for(int j=0; j<me; j++) {
		if(i!=j) {
			pf(j,i) = beta[ppidx[ctt]-1];
			ctt++;
		}
}
}

pf = pf.t();
for(int i=0; i<me; i++) {
		pf(i,i) = 1- arma::sum(pf.row(i));
}

NumericMatrix p(as<Rcpp::NumericMatrix>(wrap(pf)));
NumericMatrix pt(me,me);
if(!fitRx[1]) {
	pt = p;
} else {

	arma::mat pft(arma::zeros(me,me));
	int ctt = 0;
	IntegerVector ppidx1 = pidx["P1"];
	for(int i=0; i<me; i++) {
		for(int j=0; j<me; j++) {
			if(i!=j) {
				pft(j,i) = beta[ppidx1[ctt]-1];
				ctt++;
			}
	}
	}
	pft = pft.t();
	for(int i=0; i<me; i++) {
		pft(i,i) = 1- arma::sum(pft.row(i));
	}
	pt = (as<Rcpp::NumericMatrix>(wrap(pft)));
}
NumericVector pint(me);
IntegerVector pintidx = pidx["pi0"];
double pitotal = 0;
for(int i = 0; i<me-1; i++) {
	pint[i] = beta[pintidx[i]-1];
	pitotal = pitotal + pint[i];
}
pint[me-1] = 1-pitotal;
NumericVector pintt(me);
if(!fitRx[0]) {
	pintt = pint;
} else {
	IntegerVector pinttidx = pidx["pi1"];
	double pitotal = 0;
	for(int i = 0; i<me-1; i++) {
		pintt[i] = beta[pinttidx[i]-1];
		pitotal = pitotal + pintt[i];
	}
	pintt[me-1] = 1-pitotal;
}
IntegerVector invsigidx = pidx["invsigma"];
NumericVector sdd = NumericVector::create(::pow(beta[invsigidx[0]-1],-0.5)); // use sd instead of inverse variance
NumericMatrix ptmp(y1.nrow(),y1.ncol());
NumericMatrix ptmp2(y1.nrow(),y1.ncol());
IntegerVector tau1idx = pidx["tau1"];
IntegerVector tau2idx = pidx["tau2"];
NumericVector poisb = beta[tau1idx-1];
arma::colvec pb = Rcpp::as<arma::colvec>(poisb);
NumericVector normb = beta[tau2idx-1];
arma::colvec nb = Rcpp::as<arma::colvec>(normb);
NumericVector basepois1(me);
NumericVector basenorm(me);
basepois1 = co*pb;
basenorm = co*nb;
NumericVector lmb(2*me+1);
for(int i=0; i<me; i++) {
	lmb[i] = basepois1[i];
	lmb[me+i] = basenorm[i];
}
lmb[2*me] = sdd[0];
IntegerMatrix Z(y1.nrow(),y1.ncol());

for(int i=0; i<y1.nrow(); i++) {
IntegerVector hsseq =	seq_len(me);
	NumericVector tmpy1(ni[i]);
	NumericVector tmpy2(ni[i]);
	for(int j=0; j<ni[i]; j++) {
		tmpy1[j] = y1(i,j);
		tmpy2[j] = y2(i,j);
	}
	NumericVector lmbtmp(lmb);
	IntegerVector re1pidx = pidx["re1"];
	IntegerVector re2pidx = pidx["re2"];
	for(int cf=0; cf<me; cf++) {
		lmbtmp[cf] = ::exp(basepois1[cf] +beta[re1pidx[i]-1]);
		lmbtmp[cf+me] = basenorm[cf] + beta[re2pidx[i]-1];
	}
	NumericMatrix pii(me,me);
	NumericVector piii(me);
	if(rx[i]==1 && (fitRx[0] || fitRx[1])) {
		pii = pt;
		piii = pintt;
	} else {
		pii=p;
		piii=pint;
	}
	//Rcpp::Rcout << i << std::endl;
	List out = HMMlablpois(tmpy1,tmpy2,me,lmbtmp,pii,piii,ni[i],ltii.row(i));
	NumericMatrix tmpla = out["lalpha"];
	NumericVector tmplc(tmpla.column((ni[i]-1)));
	double tmpc = *std::max_element(tmplc.begin(),tmplc.end());
	tmplc = tmplc-tmpc;
	std::transform(tmplc.begin(),tmplc.end(),tmplc.begin(),::exp);
	//double ss = ::log(std::accumulate(tmplc.begin(),tmplc.end(), 0.0));
	//double tmpllk = tmpc+ss;
	IntegerVector tmpZZ(y1.ncol());
	NumericVector tmpprob(tmpla.column((ni[i]-1)));
	std::transform(tmpprob.begin(),tmpprob.end(),tmpprob.begin(),::exp);
	double tmpnorm = std::accumulate(tmpprob.begin(),tmpprob.end(), 0.0);
	tmpprob = tmpprob/tmpnorm;
	if(is_nan(tmpprob)[0] || is_infinite(tmpprob)[0]) {
			Rcpp::Rcout << "underflow error 2 (usually poor initial values) simulating fake hidden states.  if error continues beyond the first few iterations, there is likely something very wrong." << std::endl;
			for(int hq=0; hq<me; hq++) {
				tmpprob[hq] = ::pow(me,-1);
			}
		}
	RNGScope scope;
	IntegerVector t1 = Rcpp::RcppArmadillo::sample(hsseq, 1, false, tmpprob);
	tmpZZ[ni[i]-1] = t1[0];
	for(int qz = ni[i]-2; qz>=0; qz--) {
		NumericVector tmpprobi(tmpla.column(qz));
		std::transform(tmpprobi.begin(),tmpprobi.end(),tmpprobi.begin(),::exp);
		double tmpnormi1 = std::accumulate(tmpprobi.begin(),tmpprobi.end(), 0.0);
		tmpprobi = tmpprobi/tmpnormi1;
		NumericVector mtp(pii.column(tmpZZ[qz+1]-1));
		for(int qy = 0; qy<me; qy++) {
			tmpprobi[qy] = tmpprobi[qy]*mtp[qy];
		}
		double tmpnormi = std::accumulate(tmpprobi.begin(),tmpprobi.end(), 0.0);
		tmpprobi = tmpprobi/tmpnormi;
		if(is_nan(tmpprobi)[0] || is_infinite(tmpprobi)[0]) {
			Rcpp::Rcout << "underflow error 3 (usually poor initial values) simulating fake hidden states.  if error continues beyond the first few iterations, there is likely something very wrong" << std::endl;
			for(int hq=0; hq<me; hq++) {
				tmpprobi[hq] = ::pow(me,-1);
			}
		}
		IntegerVector tx = Rcpp::RcppArmadillo::sample(hsseq, 1, false, tmpprobi);
		tmpZZ[qz] = tx[0];
	}
	Z.row(i) = tmpZZ;
	}
	return Z;  // returns the matrix of simulated Zi
';




pc.com.pois <- cxxfunction(signature( y1m = "matrix", y2m = "matrix", betaa="numeric",m="integer",nis="integer",ccc="matrix",lti="matrix",rx1="integer",paridx="List",R_fitRx="logical"),code5.pois,plugin = "RcppArmadillo",includes=c('#include <RcppArmadilloExtensions/sample.h>','#include <cmath>',inccode.pois));


#This code is involved in sampling the random effects b_i

code.pois <- 'Environment mvtnorm("package:mvtnorm");
	Function dmvnorm = mvtnorm["dmvnorm"];
	NumericMatrix y1(y1m);
	NumericMatrix y2(y2m);
	NumericMatrix zz(zzz);
	NumericMatrix ZZ(ZZZ);
	NumericMatrix ltii(lti);
	IntegerVector nii(nis);
	int maxni = y1.ncol();
	arma::mat ZZa = Rcpp::as<arma::mat>(ZZ);
	Rcpp::List pidx(paridx);
	int N = y1.nrow();
	int ni = y1.ncol();
	NumericMatrix nrs(nre);
	NumericMatrix ors(ore);
	NumericVector obeta(b1);
	NumericVector nbeta(b2);
	NumericMatrix sig(2,2);
	IntegerVector reSigidx = pidx["Sigma"];
	sig(0,0) = nbeta[reSigidx[0]-1];
	sig(0,1) = nbeta[reSigidx[1]-1];
	sig(1,0) = nbeta[reSigidx[2]-1];
	sig(1,1) = nbeta[reSigidx[3]-1];
	IntegerVector errorSigidx = pidx["invsigma"];
	NumericVector sd = NumericVector::create(::pow(obeta[errorSigidx[0]-1],-0.5));
	NumericVector zeros(2);
	NumericVector npp(N);
	NumericVector opp(N);
	IntegerVector states(m);
	int me = states[0];
	NumericVector  tmpbetac(me);
	NumericVector tmpbetan(me);
	IntegerVector tau1idx = pidx["tau1"];
	IntegerVector tau2idx = pidx["tau2"];
	for(int i=0; i<me; i++) {		
		tmpbetac[i] = obeta[tau1idx[i]-1];
		tmpbetan[i] = obeta[tau2idx[i]-1];
	}
	arma::colvec tmpbetacc = Rcpp::as<arma::colvec> (tmpbetac);
	arma::colvec tmpbetann = Rcpp::as<arma::colvec> (tmpbetan);
	NumericVector fixedc(N*ni);
	fixedc = ZZa*tmpbetacc;
	NumericVector fixedn(N*ni);
	fixedn = ZZa*tmpbetann;
	NumericMatrix fixednn(N,ni);
	NumericMatrix fixedpp(N,ni);
	int counter1 = 0;
	for(int q=0; q<maxni; q++) {
		for(int r=0; r<N; r++) {
			if(nii[r]>q) {
			fixednn(r,q) = fixedn[counter1];
			fixedpp(r,q) = fixedc[counter1];
			counter1++;
			} else {  // Do not think necessary for this code.
				fixednn(r,q) = -10000;
				fixedpp(r,q) = -20000;
			}
			
		}
	}
	NumericVector orep = dmvnorm(ors,zeros,sig,true);
	NumericVector nrep = dmvnorm(nrs,zeros,sig,true);
	for (int i=0; i<N; i++) {
		NumericVector nnb(nii[i]);
		NumericVector npb(nii[i]);
		NumericVector onb(nii[i]);
		NumericVector opb(nii[i]);
		NumericVector eopb(nii[i]);
		NumericVector enpb(nii[i]);
		NumericVector ocpb(nii[i]);
		NumericVector ncpb(nii[i]);
		NumericVector onpb(nii[i]);
		NumericVector nnpb(nii[i]);
		for (int j=0; j<nii[i]; j++) {
				ocpb(j) = ::Rf_dpois(y1(i,j),::exp(ors(i,0) + fixedpp(i,j) + ltii(i,j)),true);
				ncpb(j) = ::Rf_dpois(y1(i,j),::exp(nrs(i,0) + fixedpp(i,j) + ltii(i,j)),true);
				onpb(j) = ::Rf_dnorm4(y2(i,j),ors(i,1) + fixednn(i,j),sd[0],true);
				nnpb(j) = ::Rf_dnorm4(y2(i,j),nrs(i,1) + fixednn(i,j),sd[0],true);
			}
		npp(i) = std::accumulate(ncpb.begin(),ncpb.end(), 0.0) + std::accumulate(nnpb.begin(),nnpb.end(), 0.0) + nrep(i);
		opp(i) = std::accumulate(ocpb.begin(),ocpb.end(), 0.0) + std::accumulate(onpb.begin(),onpb.end(), 0.0) + orep(i);
	}
	List ret; ret["npp"] = npp; ret["opp"] = opp;
	return ret;
';
	
	
library(mvtnorm);
re.n.pois <- cxxfunction(signature( y1m = "matrix", y2m = "matrix", zzz= 'matrix', ZZZ= 'matrix', b1 = "numeric",b2 = "numeric",nre = "matrix",ore="matrix",nis="numeric",lti="matrix",paridx="list",m="integer"),code.pois,plugin = "RcppArmadillo");

#Optional, but may increase speed a little
library(compiler)
enableJIT(3);

###
# simulated: Simulate from the posterior of the MVHMM through Gibbs Sampling
# arguments:
# y1.: a vector of length N * ni of Poisson responses (include NAs for missing data after dropout; this software will not currently handle intermittently missing data)
# y2.: a vector of length N * ni of normal responses (include NAs for missing data after dropout;  this software will not currently handle intermittently missing data)
# ti.: a N x ni matrix of days in a week for Poisson counts
# inits.: initial values to start the simulation.  Can include/exclude the random effects, but the other parameter values must be specified
# nsim.: number of MCMC iterations to run
# report1: how often to report progress of MCMC
# ksamp.: thinning parameter: keep every ksamp. samples
# sig.y.1.: covariance matrix for MV normal proposal for MH sampler for Poisson coefficients
# sig.re.:  covariance matrix for MV normal proposal for MH sampler for random effects
# N.: sample size (# of subjects)
# ni.: vector of length N of number of observation times for each subject
# hs:  Number of hidden states in the HMM
# rx.: treatment/covariate vector of length N
# fitRx: logical vector of length two for initial probability vector, and tpm respectively.
# scales.: scaling constant for sig.y.1. and sig.re., respectively
# id: vector of length N * ni of ids for all observations; probably just rep(1:N.,max(ni))
# hyperpar: specify hyper-parameters for inverse-Wishart (indices 1-3) and gamma (4-5)
# run.: for keeping tracking chains run in parallel; not supported in this code.
#
# Returns:
# betaa: Matrix of Posterior Samples
# Z: Matrix of hidden state samples	
# accp: # of accepted moves for MH sampler for Poisson coefficients
# pidx: List of Parameter indices
#
###


simulated.pois <- function(y1.,y2.,ti.,inits.,nsim.,report1.=1000,
ksamp.=1,sig.y.1.=diag(hs),sig.re.=diag(2),N.,ni.,
hs=3,rx.=NULL,fitRx=c(FALSE,FALSE),scales.=c(0.01,0.01,0.01),id=rep(1:N.,6),hyperpar=c(3,1,1,.001,.0002),run.=1) {
	y1m <- matrix(y1.,nr=N.,nc=max(ni.));
	y2m <- matrix(y2.,nr=N.,nc=max(ni.));
	
	if(floor(nsim./ksamp.)!=ceiling(nsim./ksamp.)) {
		stop("nsim is not a multiple of ksamp");
	}
	if(!all.equal(!is.na(y1m),!is.na(y2m))) {
		stop("responses not matching with respect to missing data");
	}
	inc <- !is.na(as.numeric(y1m));
	if(sum(fitRx)>0 & is.null(rx.)) {
		stop("No treatment vector passed to fit to");
	}
	#Parameter index
	pi0idx <- 1:(hs-1);
	P0idx <- (max(pi0idx)+1):(max(pi0idx)+hs*(hs-1));
	tau1idx <- (max(P0idx)+1):(max(P0idx)+hs);
	tau2idx <- (max(tau1idx)+1):(max(tau1idx)+hs);
	Sigmaidx <- (max(tau2idx+1)):(max(tau2idx)+4);
	invsigmaidx <- (max(Sigmaidx)+1);

	if(fitRx[1] & fitRx[2]) {
		P1idx<- (invsigmaidx+1):(invsigmaidx+hs*(hs-1));
		pi1idx <- (max(P1idx)+1):(max(P1idx)+hs-1);
		re1idx <- (max(pi1idx)+1):(max(pi1idx)+N.);
		re2idx <- (max(re1idx)+1):(max(re1idx)+N.);
		pidx <- list(pi0=pi0idx,P0=P0idx,tau1=tau1idx,tau2=tau2idx,Sigma=Sigmaidx,invsigma=invsigmaidx,P1=P1idx,pi1=pi1idx,re1=re1idx,re2=re2idx);
	} else if(!fitRx[1] & fitRx[2]) {
		P1idx<- (invsigmaidx+1):(invsigmaidx+hs*(hs-1));
		re1idx <- (max(P1idx)+1):(max(P1idx)+N.);
		re2idx <- (max(re1idx)+1):(max(re1idx)+N.);
		pidx <- list(pi0=pi0idx,P0=P0idx,tau1=tau1idx,tau2=tau2idx,Sigma=Sigmaidx,invsigma=invsigmaidx,P1=P1idx,re1=re1idx,re2=re2idx);
	} else if(fitRx[1] & !fitRx[2]) {
		pi1idx <- (max(invsigmaidx)+1):(max(invsigmaidx)+hs-1);
		re1idx <- (max(pi1idx)+1):(max(pi1idx)+N.);
		re2idx <- (max(re1idx)+1):(max(re1idx)+N.);
		pidx <- list(pi0=pi0idx,P0=P0idx,tau1=tau1idx,tau2=tau2idx,Sigma=Sigmaidx,invsigma=invsigmaidx,pi1=pi1idx,re1=re1idx,re2=re2idx);
	} else {
		re1idx <- (invsigmaidx+1):(invsigmaidx+N.);
		re2idx <- (max(re1idx)+1):(max(re1idx)+N.);
		pidx <- list(pi0=pi0idx,P0=P0idx,tau1=tau1idx,tau2=tau2idx,Sigma=Sigmaidx,invsigma=invsigmaidx,re1=re1idx,re2=re2idx);
	}
	#matrix to hold parameter/REs from MCMC ouput
	betaa <- matrix(NA,nr=nsim./ksamp.,nc=max(pidx$re2));
	if(!fitRx[1] & !fitRx[2]) {
		basepars <- max(pidx$invsigma);
		rx. <- rep(0,N.);
	} 
	if(fitRx[1] & fitRx[2]) {
		basepars <- max(pidx$pi1)
	}
	if(fitRx[1] & !fitRx[2]) {
		basepars <- max(pidx$pi1)
	}
	if(!fitRx[1] & fitRx[2]) {
		basepars <- max(pidx$P1)
	}
	betaa[1,1:basepars] <- inits.[1:basepars];
	if(length(inits.)==(basepars+N.*2)) { #If init values for random effects are passed, use them, otherwise start at zero
	
	betaa[1,c(pidx$re1,pidx$re2)] <- inits.[c(pidx$re1,pidx$re2)]
	} else {
		betaa[1,(basepars+1):(basepars+2*N.)] <-0;
	}

	o.betaa <- betaa[1,]; #Start 
	n.betaa <- rep(NA,length(betaa[1,]));
	a.param2 <- hyperpar[4];
	b.param2 <- hyperpar[5];
	scale1 <- scales.[1]; # For RE MH Simulation
	scale2 <- scales.[2]; # For Poisson \taus
	U.re <- chol(scale1*sig.re.);
	acc <- 0;
	accp <- 0;
	acc2 <- 0;
	y1.na <- y1.[inc];
	y2.na <- y2.[inc];
	id.na <- id[inc];
	ti.na <- as.numeric(ti.[inc]);
	lti.na <- log(ti.na);
	lti <- log(ti.);
	ccc <- matrix(NA,nr=hs,nc=hs);
	for(q in 1:hs) {
		ccc[q,] <- c(rep(1,q),rep(0,hs-q));
	}
	Zmat <- matrix(NA,nr=nsim./ksamp.,nc=max(ni.)*N.)
for(i in 2:nsim.) {
	#Simulate Sigma
	#print(i);
	ress <- t(cbind(o.betaa[re1idx],o.betaa[re2idx]))
	cvvv <- ress%*%t(ress);
	n.betaa[Sigmaidx] <- as.numeric(riwish(N.+hyperpar[1],diag(hyperpar[2:3]) + cvvv));
	#simulated hidden states (Z)
	#print(y1m);
	#print(y2m);
	#print(o.betaa)
	#print(lti);
	#print(ni.);
	zz2 <- pc.com.pois(y1m,y2m,o.betaa,hs,ni.,ccc,lti,rx.,pidx,fitRx)
	zz2[zz2==0] <- NA;
	zz.na <- na.omit(as.numeric(zz2));
	#Simulate P and Pi
	if(fitRx[2]) {
		zz2.rx <- zz2[rx.==1,];
		zzt.rx <- table(factor(zz2.rx[,-ncol(zz2.rx)],levels=as.character(c(1:hs))),factor(zz2.rx[,-1],levels=as.character(c(1:hs))));
		zz2.px <- zz2[rx.==0,];
		zzt.px <- table(factor(zz2.px[,-ncol(zz2.px)],levels=as.character(c(1:hs))),factor(zz2.px[,-1],levels=as.character(c(1:hs))));
		ctt <- 1;
		for(ps in 1:hs) {
			n.betaa[P0idx[ctt:(ctt+hs-2)]] <- rdirichlet(1,zzt.px[ps,]+rep(1,hs))[-ps];
			n.betaa[P1idx[ctt:(ctt+hs-2)]] <- rdirichlet(1,zzt.rx[ps,]+rep(1,hs))[-ps];
			ctt <- ctt + hs-1;
		}
	} else {
		zzt <- table(factor(zz2[,-ncol(zz2)],levels=as.character(c(1:hs))),factor(zz2[,-1],levels=as.character(c(1:hs))));
		ctt <- 1;
		for(ps in 1:hs) {
			n.betaa[P0idx[ctt:(ctt+hs-2)]] <- rdirichlet(1,zzt[ps,]+rep(1,hs))[-ps];
			ctt <- ctt + hs-1;
		}
	}
	if(fitRx[1]) {
	    zz2.rx <- zz2[rx.==1,];
		zz2.px <- zz2[rx.==0,];
		n.betaa[pi0idx] <- rdirichlet(1,as.numeric(table(factor(zz2.px[,1],levels=as.character(c(1:hs))))+1))[1:(hs-1)];
		n.betaa[pi1idx] <- rdirichlet(1,as.numeric(table(factor(zz2.rx[,1],levels=as.character(c(1:hs))))+1))[1:(hs-1)];
	} else {
		n.betaa[pi0idx] <- rdirichlet(1,as.numeric(table(factor(zz2[,1],levels=as.character(c(1:hs))))+1))[1:(hs-1)];
	}
	ZZZ <- matrix(NA,nr=length(zz.na),nc=hs);
	ZZZ[,1] <- 1;
	for(s in 2:hs) {
		ZZZ[,s] <- as.numeric(1*zz.na>(s-1));
	}
	new.res <- matrix(NA,nr=N.,nc=2);
	tmp.re <- matrix(rnorm(2*N.,0,1),nr=N.,nc=2);
	for(qqq in 1:N.) {
		new.res[qqq,] <- cbind(o.betaa[re1idx[qqq]],o.betaa[re2idx[qqq]]) + as.vector(t(U.re)%*%tmp.re[qqq,]);
	}
	old.res <- matrix(c(o.betaa[re1idx],o.betaa[re2idx]),nr=N.);
	mh.s <- re.n.pois(y1m,y2m,zz2,ZZZ,o.betaa[1:basepars],n.betaa[1:basepars],new.res,old.res,ni.,lti,pidx,hs);
	npp <- mh.s$npp;
	opp <- mh.s$opp;
		
	alp <- exp(mh.s$npp-mh.s$opp);
	
	alp <- ifelse(alp>1,1,alp);
	u <- runif(N.);
	app.new <- as.numeric(u<=alp);
	n.betaa[re1idx] <- app.new*new.res[,1] + (1-app.new)*o.betaa[re1idx];
	n.betaa[re2idx] <- app.new*new.res[,2] + (1-app.new)*o.betaa[re2idx];
	
	#Simulate tau_1
	new.betaa <- rmvnorm(1,o.betaa[tau1idx],scale2*sig.y.1.);
	npp <- sum(dpois(y1.na,exp(n.betaa[re1idx[id.na]] +ZZZ%*%t(new.betaa) + lti.na),TRUE));
	if(is.na(npp)) {
		message(paste("Warning cannot sample from tau_1", o.betaa,sep=""));
		stop("");
	}
	 opp <- sum(dpois(y1.na,exp(n.betaa[re1idx[id.na]] +ZZZ%*%o.betaa[tau1idx] + lti.na),TRUE));
	 u <- runif(1);
	 if(u<=min(c(1,exp(npp-opp)))) {
			 n.betaa[tau1idx] <- new.betaa;
			 accp <- accp +1;
		 } else {
			 n.betaa[tau1idx] <- o.betaa[tau1idx];
		 }
		 rm(npp);
		 rm(opp); rm(new.betaa);
	#Simualte tau_2
	bbcov <- solve(crossprod(ZZZ))
	bbmean <- bbcov%*%t(ZZZ)%*%(y2.na-n.betaa[re2idx[id.na]]);
	bbcov <- 1/o.betaa[invsigmaidx]*bbcov;
	n.betaa[tau2idx] <- rmvnorm(1,bbmean,bbcov); 
	#Simulate 1/sigma^2_e
	n.betaa[invsigmaidx] <- rgamma(1,(sum(ni.)/2 + a.param2),b.param2 + sum((y2.na-n.betaa[re2idx[id.na]] -ZZZ%*%n.betaa[tau2idx])^2)/2); #tau
	if(floor(i/ksamp.)==ceiling(i/ksamp.)) {
		betaa[(i-1)/ksamp.+1,] <- n.betaa;
		Zmat[(i-1)/ksamp.+1,] <- as.vector(zz2);
	}
	o.betaa <- n.betaa;
	rm(n.betaa);
	n.betaa <- rep(NA,length(o.betaa));
	if(ceiling(i/report1.)==floor(i/report1.)) {
		message("iteration: ", i, " ", date(), " Run: ", run.);
	}
	
	
	
}
	return(list(betaa=betaa,Z=Zmat,accp=accp,pidx=pidx));
}


