parameter.est.allobs = function(z.targ,z.pert,theta.pert,theta.mean.prior=NULL,theta.covm.prior=NULL,augment.cov=TRUE,pnames=NULL) {
#### ESTIMATE PARAMETERS OF THE MODEL
#### Y = L * THETA + MU + NOISE
####
#### INPUT:
#		Z.TARG[TMAX.TARG,NVAR]: 		THE TARGET TIME SERIES
#		Z.PERT[TMAX.PERT,NVAR,NEXP]: 	PERTURATION RUNS
#		THETA.PERT[NEXP,NPARM]: 		PARAMETERS IN THE PERTURBATION RUNS
#		THETA.MEAN.PRIOR[NPARM]: 		MEAN OF THE PRIOR DIST ON THETA
#		THETA.COVM.PRIOR[NPARM,NPARM]: 	COV. MATRIX OF THE PRIOR DIST ON THETA
#		TIME.MEAN: 						PROCESS TIME MEANS (T) OR PROCESS INDIVIDUAL TIME STEPS (F)
#		PNAMES[NPARM]: 					NAMES OF THE PARAMETERS (OPTIONAL)


tmax.targ = dim(z.targ)[1]
tmax.pert = dim(z.pert)[1]
nvar      = dim(z.pert)[2]
nexp      = dim(z.pert)[3]
nparm     = dim(theta.pert)[2]
if (nvar != dim(z.targ)[2]) stop('z.targ and z.pert have inconsistent number of variables')
if (nexp != dim(theta.pert)[1]) stop('inconsistent number of experiments in theta.pert and z.pert')

if (!is.null(pnames)) {
	if (nparm != length(pnames)) stop('number of parameters inconsistent between theta.pert and pnames')
} else {
	pnames = paste('P',1:nparm,sep='')
}

##### CHECK FOR ENSEMBLE COLLAPSE
l.collapse = any( sqrt(diag(cov(theta.pert)))/abs(colMeans(theta.pert)) < 1.e-5 )

if (nvar > tmax.pert * nexp - nparm - 1 | l.collapse) return(
list(theta.mean = NA, theta.mean.end = NA,
     theta.covm = NA, theta.covm.end = NA,
     theta.sd   = NA, theta.sd.end   = NA,
     intercept  = NA, lmat           = NA,
     cov.noise  = NA, l.collapse     = l.collapse,
     parm.est   = rep(NA,nparm), 
     parm.std   = rep(NA,nparm)))

if (is.null(theta.mean.prior)) theta.mean.prior = colMeans(theta.pert)
if (is.null(theta.covm.prior)) theta.covm.prior = cov(theta.pert)



##############################################################
######## ESTIMATE THETA IN CHUNKS
##############################################################
if (tmax.targ %% tmax.pert != 0) lchunk = FALSE else lchunk = TRUE
nchunk = floor(tmax.targ / tmax.pert)

theta.mean      = array(NA,dim=c(nparm,nchunk+1))
theta.covm      = array(NA,dim=c(nparm,nparm,nchunk+1))
rownames(theta.mean) = pnames
theta.mean[ ,1] = theta.mean.prior
theta.covm[,,1] = theta.covm.prior

if (lchunk) for ( nc in 1:nchunk) {
	no.pic             = 1:tmax.pert + (nc-1) * tmax.pert
	kf.lm.list         = kf.lm(z.targ[no.pic,,drop=FALSE],z.pert,theta.pert,theta.mean[,nc],theta.covm[,,nc],augment.cov=augment.cov,pnames=pnames)
	theta.mean[ ,nc+1] = kf.lm.list$theta.mean.update
	theta.covm[,,nc+1] = kf.lm.list$theta.covm.update
} 

##############################################################
######## ESTIMATE THETA IN ONE STEP
##############################################################
kf.lm.list     = kf.lm(z.targ,z.pert,theta.pert,theta.mean[,1],theta.covm[,,1],augment.cov=augment.cov,pnames=pnames)
theta.mean.end = kf.lm.list$theta.mean.update
theta.covm.end = kf.lm.list$theta.covm.update

##############################################################
######## COMPUTE STANDARD DEVIATIONS
##############################################################
theta.sd.end   = sqrt(diag(theta.covm.end))
theta.sd       = array(NA,dim=c(nparm,nchunk+1))
for (nc in 0:nchunk+1) theta.sd[,nc] = sqrt(diag(theta.covm[,,nc]))

##############################################################
######## GENERATE NEW PARAMETER VALUES FOR THE NEXT ITERATION
##############################################################
cov.eigen = eigen(theta.covm.end)
theta.new = as.numeric(theta.mean.end) + cov.eigen$vectors %*% diag(sqrt(cov.eigen$values),nrow=nparm,ncol=nparm) %*% array(rnorm(nparm*nexp),dim=c(nparm,nexp))

det.noise = sum(log(eigen(kf.lm.list$cov.noise)$values))

list(theta.mean = theta.mean, theta.mean.end = theta.mean.end,
     theta.covm = theta.covm, theta.covm.end = theta.covm.end,
     theta.sd   = theta.sd  , theta.sd.end   = theta.sd.end  ,
     intercept  = kf.lm.list$intercept , lmat           = kf.lm.list$lmat          ,
     rsqr       = kf.lm.list$rsqr      , can.cor.sqr    = kf.lm.list$can.cor.sqr   ,
     cov.noise  = kf.lm.list$cov.noise , l.collapse      = l.collapse   ,
     parm.est   = as.numeric(theta.mean.end), parm.std = as.numeric(theta.sd.end)  ,
     parm.update = theta.new, det.noise = det.noise)
     
}



