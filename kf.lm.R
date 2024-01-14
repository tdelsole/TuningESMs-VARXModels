kf.lm = function(z.targ,z.pert,theta.pert,theta.mean.prior,theta.covm.prior,lambda=0,augment.cov=TRUE,pnames=NULL) {
#### ESTIMATE PARAMETERS OF THE MODEL
#### Y = L * THETA + MU + NOISE
####
#### INPUT:
#		Z.TARG[TMAX.TARG,NVAR]: 		THE TARGET TIME SERIES
#		Z.PERT[TMAX.PERT,NVAR,NEXP]: 	PERTURATION RUNS
#		THETA.PERT[NEXP,NPARM]: 		PARAMETERS IN THE PERTURBATION RUNS
#		THETA.MEAN.PRIOR[NPARM]: 		MEAN OF THE PRIOR DIST ON THETA
#		THETA.COVM.PRIOR[NPARM,NPARM]: 	COV. MATRIX OF THE PRIOR DIST ON THETA
#		LAMBDA:							REGULARIZATION PARAMETER R = LAMBDA * I
#		AUGMENT.COV:					LOGICAL: NOISE COV ESTIMATED FROM PERTURBATION ONLY (FALSE), OR PERT+TARG (TRUE)
#		PNAMES[NPARM]: 					NAMES OF THE PARAMETERS (OPTIONAL)
#### COMMENTS
# 1) THE TARGET AND PERTURBATIONS CAN BE DIFFERENT LENGTHS
# 2) PERTURBATION RUNS CAN BE DIFFERENT LENGTHS, BUT FOR EFFICIENCY THIS ALGORITHM ASSUMES THEY ARE THE SAME
# 3) VERIFIED THAT RECURSIVE = BATCH (7/2/22), PROVIDED AUGMENT.COV = FALSE (BECAUSE COV[Y|THETA] DEPENDS ON ALL TARGET DATA)

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
list(parm.est   = NA , parm.std        = NA        ,
     intercept  = NA , lmat            = NA        ,
     cov.noise  = NA , l.collapse      = l.collapse,
     theta.mean.udpate = NA, 
     theta.covm.udpate = NA ))

##############################################################
######## ESTIMATE PARAMETERS OF P[Y|THETA]
##############################################################
y = NULL
x = NULL
j = rep(1,tmax.pert)
for (ne in 1:nexp) {
	y = rbind(y,array(z.pert[,,ne],dim(z.pert)[1:2]))
	x = rbind(x,j %*% t(theta.pert[ne,]))
}
colnames(x) = pnames
xy.lm       = lm(y~x)
intercept   = coef(xy.lm)[1,]
lmat        = t(coef(xy.lm)[-1,])
cov.noise   = (t(residuals(xy.lm)) %*% residuals(xy.lm)) / xy.lm$df.residual
# cov.theta   = cov(parm.pert)
cov.theta   = cov(theta.pert)

##############################################################
######## COMPUTE MULTIVARIATE R-SQUARE
##############################################################
gev.list     = gev(cov.noise,cov.noise + lmat %*% cov.theta %*% t(lmat))
cca.list     = gev(lmat %*% cov.theta %*% t(lmat),cov.noise + lmat %*% cov.theta %*% t(lmat))
rsqr         = 1 - prod(gev.list$lambda)


##############################################################
######## AUGMENT COVARIANCE MATRIX USING TARGET DATA
##############################################################
if (augment.cov) cov.noise = (cov.noise * xy.lm$df.residual + cov(z.targ) * (tmax.targ-1))/(xy.lm$df.residual + (tmax.targ-1))

obs               = colMeans(z.targ)
cov.ytot          = lmat %*% theta.covm.prior %*% t(lmat) + cov.noise/tmax.targ
diag(cov.ytot)    = diag(cov.ytot) + lambda
cov.ytot.inv      = chol2inv(chol(cov.ytot))
theta.mean.update = theta.mean.prior + theta.covm.prior %*% t(lmat) %*% cov.ytot.inv %*% ( obs - intercept - lmat %*% theta.mean.prior )
theta.covm.update = theta.covm.prior - theta.covm.prior %*% t(lmat) %*% cov.ytot.inv %*% lmat %*% theta.covm.prior

##############################################################
######## COMPUTE STANDARD DEVIATIONS
##############################################################
parm.est   = as.numeric(theta.mean.update)
parm.std   = as.numeric(sqrt(diag(theta.covm.update)))
names(parm.est)             = pnames
names(parm.std)             = pnames
rownames(theta.covm.update) = pnames
colnames(theta.covm.update) = pnames


list(parm.est   = parm.est  , parm.std       = parm.std      ,
     intercept  = intercept , lmat           = lmat          ,
     cov.noise  = cov.noise , l.collapse     = l.collapse    ,
     cov.ytot   = cov.ytot  , obs            = obs           ,
     rsqr       = rsqr      , can.cor.sqr    = cca.list$lambda,
     theta.mean.update = theta.mean.update, 
     theta.covm.update = theta.covm.update )

	
}