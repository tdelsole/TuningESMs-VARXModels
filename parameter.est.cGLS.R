parameter.est.cGLS = function(ts.list,pnames=NULL) {
#### ESTIMATE PARAMETERS OF THE MODEL
#### Y = L * THETA + MU + NOISE
#### USING GENERALIZED LEAST SQUARES, CONDITIONAL ON ESTIMATES OF L, MU, COV_NOISE FROM PERTURBATION RUNS
#### INPUT:
##		TS.LIST(TSERIES[NTOT,SDIM], FIRST.STEP, THETA)[NENS+1]: 
##			LIST CONTAINING TARGET AND MODEL RUN TIME SERIES AND ASSOCIATED METADATA
##			TS.LIST[[1]] IS THE TARGET.  TS.LIST[[2]],...,TS.LIST[[E+1]] ARE PERTURBATION RUNS
##			THETA: PARAMETER VALUES ASSOCIATED WITH TIME SERIES
##			NTOT IS THE LENGHT OF THE TIME SERIES (CAN DIFFER ACROSS ELEMENTS)
##			SDIM IS THE NUMBER OF VARIABLES/EOFS/LAPLACIANS (MUST BE SAME ACROSS TS.LIST[[1]],...,TS.LIST[[E+1]])

nreal       = length(ts.list)
nparm.base  = length(ts.list[[1]]$theta)

if (is.null(pnames)) {
	pnames.say = paste('P',1:nparm.base,sep='')
} else {
	pnames.say = pnames
}

############################################	
########## SET UP XMAT AND YMAT	
############################################	
xmat = NULL
ymat = NULL
for (nr in 2:nreal) {
	theta.pic = ts.list[[nr]]$theta
	yget      = as.matrix(ts.list[[nr]]$tseries)
	ndim      = dim(yget)[1]
	sdim      = dim(yget)[2]
	nparm     = length(theta.pic)
	if (nr == 2) {sdim.base = sdim}
	if (sdim.base  != sdim ) stop('inconsistent spatial dimesion')
	if (nparm.base != nparm) stop('inconsistent number of parameter')
		
	ymat = rbind(ymat,yget)
	xmat = rbind(xmat,matrix(theta.pic,nrow=ndim,ncol=nparm,byrow=TRUE))
}

colnames(ymat) = paste('L',1:sdim,sep='')
colnames(xmat) = pnames.say

############################################	
########## ESTIMATE PARAMETERS OF P[Y | THETA]	
############################################
xy.lm      = lm(ymat~xmat)
lmat.trans = coef(xy.lm)[-1,]
mu.cond    = coef(xy.lm)[1,]
covm.cond  = (t(residuals(xy.lm)) %*% residuals(xy.lm)) / xy.lm$df.residual

############################################	
########## ESTIMATE THETA0	
############################################
z.targ        = as.matrix(ts.list[[1]]$tseries)
n.targ        = dim(z.targ)[1]

covm.cond     = covm.cond /n.targ
covm.cond.inv = chol2inv(chol(covm.cond))

lt.covm       = lmat.trans %*% covm.cond.inv
lt.covm.l.inv = chol2inv(chol(lt.covm %*% t(lmat.trans)))
theta.o.mean  = lt.covm.l.inv %*% lt.covm %*% ( colMeans(z.targ) - mu.cond )
theta.o.covm  = lt.covm.l.inv

list(theta.o.mean = theta.o.mean, theta.o.covm = theta.o.covm, theta.o.sd = sqrt(diag(theta.o.covm)))

}



