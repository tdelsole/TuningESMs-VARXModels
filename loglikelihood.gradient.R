loglikelihood.gradient = function(theta,output = 'likelihood') {
#### COMPUTE LIKELIHOOD FUNCTION FOR PARAMETER ESTIMATION
#### INPUT:
#### 	THETA[NPARM]: PARAMETER VALUES
####	OUTPUT: CHARACTER STRING SPECIFYING THE OUTPUT
####			DEFAULT = 'likelihood' to produce log-likelihood and gradient
####			= 'gamma' to produce gamma matrix (noise cov. matrix)
####			= 'varx' to produce all parameters in VARX model
####			= 'information' to produce information matrix
#### ENVIRONMENT VARIABLES:
#### 	XYMAT.LIST: LIST CONTAINING OUTPUT FROM PARAMETER.EST.XYMAT	
####
	
	## FORCE XMAT.O.F TO BE MATRIX, SO THAT NUMBER OF FORCINGS = NUMBER OF COLUMNS	
	xmat.o.f  = as.matrix(xymat.list$xmat.targ.f)
	num.f     = dim(xymat.list$xmat.targ.f)[2]
	nrow.targ = dim(xymat.list$ymat.targ)[1]
	nrow.pert = dim(xymat.list$ymat.pert)[1]
	
	## CONSTRUCT XMAT USING THETA VALUE
	xmat.o.ft = NULL
	for (nf in 1:num.f) xmat.o.ft = cbind(xmat.o.ft,xymat.list$xmat.targ.f[,nf] %*% t(theta)) 
	xmat = rbind(cbind(xymat.list$xmat.targ.ar,xymat.list$xmat.targ.f,xmat.o.ft),
	             cbind(xymat.list$xmat.pert.ar,xymat.list$xmat.pert.f,xymat.list$xmat.pert.ft))
	             
	### CONSTRUCT YMAT
	ymat = rbind(xymat.list$ymat.targ,xymat.list$ymat.pert)
	
	### SOLVE FOR LS ESTIMATE OF BETA.HAT, CONDITIONED ON THETA; COMPUTE RESIDUALS
	ntot  = dim(xmat)[1]
	xy.lm = lm(ymat~xmat-1)
	gamma.hat = (t(residuals(xy.lm)) %*% residuals(xy.lm)) / ntot
	
	### COMPUTE GRADIENT
	gamma.inv    = chol2inv(chol(gamma.hat))
	ncol.ar      = dim(xymat.list$xmat.targ.ar)[2]
	ncol.f       = dim(xymat.list$xmat.targ.f)[2]
	ncol.ft      = dim(xymat.list$xmat.targ.ft)[2]
	sdim         = dim(xymat.list$ymat.targ)[2]
	
	if (is.null(ncol.ar)) {
		ncol.ar  = 0
		ar.trans = NA
	} else {
		ar.trans = coef(xy.lm)[1:ncol.ar,]
	}
	ce.trans     = coef(xy.lm)[1:ncol.f  + ncol.ar,,drop=FALSE]
	cce.trans    = coef(xy.lm)[1:ncol.ft + ncol.ar + ncol.f,,drop=FALSE]
	
	dim(cce.trans) = c(nparm,ncol.f,sdim)
	umat           = 0
	vvec           = 0

	for (j in 1:ncol.f) for (jp in 1:ncol.f) umat = umat + 
	   cce.trans[,j,] %*% gamma.inv %*% t(cce.trans[,jp,]) * sum(xymat.list$xmat.targ.f[,j] * xymat.list$xmat.targ.f[,jp])
	   
	res = xymat.list$ymat.targ - xymat.list$xmat.targ.f %*% ce.trans
	if (ncol.ar != 0) res = res - xymat.list$xmat.targ.ar %*% ar.trans	
	for (j in 1:ncol.f) vvec = vvec + cce.trans[,j,] %*% gamma.inv %*% t(res) %*% xymat.list$xmat.targ.f[,j]
	   
	# for (j in 1:ncol.f) vvec = vvec + 
	   # cce.trans[,j,] %*% gamma.inv %*% 
	   # t(xymat.list$ymat.targ - xymat.list$xmat.targ.ar %*% ar.trans - xymat.list$xmat.targ.f %*% ce.trans) %*% xymat.list$xmat.targ.f[,j] 
	   
	if (!is.null(xymat.list$prior.mean)) {
		umat = umat + xymat.list$prior.covm.inv
		vvec = vvec + xymat.list$prior.covm.inv %*% xymat.list$prior.mean
	}
	
	gradient = 2 * (umat %*% theta - vvec)
	
	if (output == 'information') {
		singular = FALSE
		if (dim(xmat)[1] < dim(xmat)[2]) singular = TRUE
		if (!singular) xtx.inv  = chol2inv(chol( t(xmat) %*% xmat))

		gmat = t(cbind(xymat.list$xmat.targ.ar,xymat.list$xmat.targ.f,xmat.o.ft)) %*% xymat.list$xmat.targ.f
		
		if (is.null(xymat.list$prior.mean)) info.mat = 0 else info.mat = xymat.list$prior.covm.inv

		for (j in 1:ncol.f) for (jp in 1:ncol.f) {
			if (j == ncol.f & jp == ncol.f) fctr1 = sum(xymat.list$xmat.targ.f[,j] * xymat.list$xmat.targ.f[,jp])
			if (j == ncol.f & jp == ncol.f) fctr2 =    as.numeric(gmat[,j] %*% xtx.inv %*% gmat[,jp])
			fctr = sum(xmat.o.f[,j] * xmat.o.f[,jp]) - as.numeric(gmat[,j] %*% xtx.inv %*% gmat[,jp])
			info.mat = info.mat + fctr * cce.trans[,j,] %*% gamma.inv %*% t(cce.trans[,jp,])
			if (is.null(xymat.list$prior.mean) & abs(fctr) < 1.e-6) singular = TRUE
		}
		if (singular) info.mat = matrix(Inf,nrow=length(theta),ncol=length(theta)) else info.mat = chol2inv(chol(info.mat))
		
	}
	
	
	if (output == 'likelihood') {
		gev.list   = gev(gamma.hat,gamma.first.guess)
		output.say = ntot * sum(log(gev.list$lambda))
		if (!is.null(xymat.list$prior.mean)) output.say = output.say + 
		   (theta - xymat.list$prior.mean) %*% xymat.list$prior.covm.inv %*% (theta - xymat.list$prior.mean)
		attr(output.say,'gradient') = gradient
	} else if (output == 'gamma') {
		output.say = gamma.hat
	} else if (output == 'varx') {
		dim(cce.trans) = c(nparm*ncol.f,sdim)
		output.say = list(gamma=gamma.hat,ar.trans=ar.trans,ce.trans=ce.trans,cce.trans=cce.trans,
		    xmat=xmat,ymat=ymat,beta.hat=coef(xy.lm),residuals=residuals(xy.lm),
		    ncol.ar = ncol.ar, ncol.f = ncol.f, ncol.ft = ncol.ft,
		    nrow.targ  = nrow.targ , nrow.pert = nrow.pert,
		    umat = umat, vvec = vvec)
	} else if (output == 'information') {
		output.say = list(info.mat=info.mat,fctr1=fctr1,fctr2=fctr2)
	}
	
	output.say
	
}

