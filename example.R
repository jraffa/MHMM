## This is a brief example of how one may use the functions contained in ff-final.R to first simulate 
## multivariate response longitudinal data from an HMM with genSim() and subsequently provide parameter
## estimates via a MCMC-based approach with simulated()
## For further details see the original paper and supplementary material.
##
## For specific information about the functions, see the comments in ff-final.R.
## Jesse Raffa (jraffa@uw.edu): http://github.com/jraffa/MHMM


#Load R Functions and compile C++ code
file.name <- 'sim-hs3-1'; #What you want R to save the posterior sample files to. (it will add the .Rdata extension)
source("ff-final.R");
# Set seed
set.seed(25562);

#Generate Simulated Data under default settings
daf <- genSim();

#Run Preliminary MCMC
x <- simulated(y1.=as.numeric(daf$y1),y2.=as.numeric(daf$y2),mi.=daf$mi,inits.=c(daf$pars),nsim.=40000,ksamp.=40,
N.=daf$N,ni=rep(6,daf$N),rx.=daf$rx,fitRx=daf$fitRx,report1.=500,id=rep(1:daf$N,6))


#Get updated proposals
fracc <- 2
prelim.samps <- dim(x$betaa)[1];
cov.beta.y1 <- cov(x$betaa[ceiling(prelim.samps/fracc+1):prelim.samps, x$pidx$tau1])


no.pars <- max(c(x$pidx$invsigma,x$pidx$P1,x$pidx$pi1))
cov.re <- list();
for(i in 1:daf$N) {
	cov.re[[i]] <- cov(cbind(x$betaa[(ceiling(prelim.samps/fracc)+1):(prelim.samps),no.pars+i],x$betaa[(ceiling(prelim.samps/fracc)+1):(prelim.samps/fracc),no.pars+i+daf$N]))
}

sig.re <- Reduce('+',cov.re)/daf$N

#Run MCMC for a longer period of time to sample from Posterior
samp <- simulated(y1.=as.numeric(daf$y1),y2.=as.numeric(daf$y2),mi.=daf$mi,inits.=colMeans(x$betaa[(ceiling(prelim.samps/fracc)+1):(prelim.samps),]),nsim.=3*10^6,ksamp.=60,
N.=daf$N,ni=rep(6,daf$N),rx.=daf$rx,fitRx=daf$fitRx,report1.=500,id=rep(1:daf$N,6), sig.y.1.=cov.beta.y1,sig.re.=sig.re,scales=c(0.1,0.1))


#Save the Posterior Samples to Disk
save.image(paste(file.name,".rdata",sep=""));

#Output the posterior means of the samples after some burnin period (note: \Sigma_{12} is reported twice)
colMeans(samp$betaa[30001:50000,1:no.pars]);


#See samp$pidx for the indices these correspond with.  For example, this reports the non-diagonal entries for the transition probability matrix, P.
colMeans(samp$betaa[30001:50000,samp$pidx$P0]);



