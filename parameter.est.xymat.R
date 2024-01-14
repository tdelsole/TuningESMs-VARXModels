parameter.est.xymat = function(ts.list,p.order,nharm,npoly=0,period=12,prior.mean = NULL,prior.covm = NULL) {
#### SET UP X AND Y MATRICES FOR Y = XB + E FOR PARAMETER ESTIMATION.
#### YO(T) = A1 * YO(T-1) + ... + AP * YO(T-P) + c1 * F1(T) + ... cJ * FJ(T) + C1 * THETA.O * F1(T) + ... + CJ * THETA.O * FJ(T)
#### YE(T) = A1 * YE(T-1) + ... + AP * YE(T-P) + c1 * F1(T) + ... cJ * FJ(T) + C1 * THETA.O * F1(T) + ... + CJ * THETA.O * FJ(T)
#### INPUT:
##		TS.LIST(TSERIES[NTOT,SDIM], FIRST.STEP, THETA)[NENS+1]: 
##			LIST CONTAINING TARGET AND MODEL RUN TIME SERIES AND ASSOCIATED METADATA
##			TS.LIST[[1]] IS THE TARGET.  TS.LIST[[2]],...,TS.LIST[[E+1]] ARE MODEL RUNS
##			FIRST.STEP: INTEGER IN {1, 2, ..., PERIOD} FOR PHASE OF THE CYCLE AT T=1.  
##			THETA: PARAMETER VALUES ASSOCIATED WITH TIME SERIES
##			NTOT IS THE LENGHT OF THE TIME SERIES (CAN DIFFER ACROSS ELEMENTS)
##			SDIM IS THE NUMBER OF VARIABLES/EOFS/LAPLACIANS (MUST BE SAME ACROSS TS.LIST[[1]],...,TS.LIST[[E+1]])
##		P.ORDER: ORDER OF THE VAR MODEL
##		NHARM: NUMBER OF HARMONICS OF THE CYCLE
##		PERIOD: PERIOD OF THE CYCLES

###############################
## HARD TO SET UP POLYNOMIAL WHILE ALLOWING OBS AND ENS TO HAVE DIFFERENT TIME LENGHTS
## SEEMS LIKE I NEED TO DEFINE A SINGLE TIME SERIES, AND THEN IDENTIFY THE START TIME OF YOBS AND YENS
## REMEMBER PERTURBATION RUNS MAY START IN DIFFERENT MONTHS AND YEARS.
## INCLUDING A TIME VECTOR FOR OBS AND PERTURBATION RUNS WOULD RESOLVE ALL THESE ISSUES
## BUT THIS MEANS EVERYONE WHO RUNS IT WILL HAVE THE BURDEN OF SPECIFYING A DATE VECTOR
## THIS COULD BE AVOIDED BY DEFAULTING TO ALL TIME SERIES START ON SAME MONTH AND YEAR
###############################
nreal     = length(ts.list)

###############################
######## DEFINE POLYNOMIAL OVER THE LONGEST LENGTH OF THE TIME SERIES
######## THIS ASSUMES ALL TIME SERIES START AT THE SAME DATE!!!!!!
###############################
if (npoly > 0) {
	tmax      = 0
	for (nr in 1:nreal) tmax = max(tmax,dim(as.matrix(ts.list[[nr]]$tseries))[1])
	t.poly    = poly(1:tmax,npoly,simple=TRUE)	
}

###############################
######## SET UP XMAT AND YMAT
###############################
ymat.e    = NULL
xmat.e.ar = NULL
xmat.e.f  = NULL
xmat.e.ft = NULL
for (nr in 1:nreal) {
	first.step = ts.list[[nr]]$first.step
	theta.pic  = ts.list[[nr]]$theta
	xget       = as.matrix(ts.list[[nr]]$tseries)
	ntot       = dim(xget)[1]
	sdim       = dim(xget)[2]
	if (nr == 1) sdim1 = sdim else if (sdim != sdim1) stop('sdim inconsistent across time series')
	npic       = (1+p.order):ntot
	
	xmat.ar     = NULL
	xmat.fonly  = NULL
	xmat.ftheta = NULL 
	if (p.order > 0) for (lag in 1:p.order) xmat.ar = cbind(xmat.ar,xget[npic-lag,])

	if (nharm   > 0) for (nh  in 1:nharm  ) {
		arg         = 2*pi*(npic -1 + first.step - 1)*nh/period
		xmat.fonly  = cbind(xmat.fonly ,cos(arg),sin(arg))
		xmat.ftheta = cbind(xmat.ftheta,cos(arg) %*% t(theta.pic),sin(arg) %*% t(theta.pic))
	}
	if (nharm > 0) colnames(xmat.fonly ) = paste(c('cos','sin'),rep(1:nharm,each=2),sep='')
	ncol.fonly = dim(xmat.fonly)[2]
	
	if (npoly > 0) for (np in 1:npoly) {
		xmat.fonly  = cbind(xmat.fonly ,t.poly[npic,np])
		xmat.ftheta = cbind(xmat.ftheta,t.poly[npic,np] %*% t(theta.pic))
	}
	if (npoly > 0) colnames(xmat.fonly)[1:npoly+ncol.fonly] = paste('poly',1:npoly,sep='')
	
	intercept   = rep(1,length(npic))
	xmat.fonly  = cbind(xmat.fonly,intercept)
	xmat.ftheta = cbind(xmat.ftheta,intercept %*% t(theta.pic))
	
	if (nr == 1) {
		xmat.o.ar = xmat.ar
		xmat.o.f  = xmat.fonly
		xmat.o.ft = xmat.ftheta
		ymat.o    = xget[npic,]
	} else {
		xmat.e.ar = rbind(xmat.e.ar,xmat.ar)
		xmat.e.f  = rbind(xmat.e.f ,xmat.fonly)
		xmat.e.ft = rbind(xmat.e.ft,xmat.ftheta)
		ymat.e    = rbind(ymat.e   ,xget[npic,])
	}
}

if (dim(xmat.o.f )[2] != dim(xmat.e.f )[2]) stop('xmat.f  has inconsistent number of columns between o and e')
if (dim(xmat.o.ft)[2] != dim(xmat.e.ft)[2]) stop('xmat.ft has inconsistent number of columns between o and e')
if (p.order > 0) if (dim(xmat.o.ar)[2] != dim(xmat.e.ar)[2]) stop('xmat.ar has inconsistent number of columns between o and e')
if (dim(xmat.o.ft)[2] != nparm * dim(xmat.o.f )[2]) stop('number of columns in xmat.o.ft is inconsistent')

ymat.o = as.matrix(ymat.o)
ymat.e = as.matrix(ymat.e)

if (is.null(prior.mean)) {
	prior.covm.inv = NULL
} else {
	prior.covm.inv = chol2inv(chol(prior.covm))
}


list(ymat.targ = ymat.o, xmat.targ.ar = xmat.o.ar, xmat.targ.f = xmat.o.f, xmat.targ.ft = xmat.o.ft, 
     ymat.pert = ymat.e, xmat.pert.ar = xmat.e.ar, xmat.pert.f = xmat.e.f, xmat.pert.ft = xmat.e.ft,
     prior.mean = prior.mean, prior.covm = prior.covm, prior.covm.inv = prior.covm.inv)

}



