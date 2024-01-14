parameter.est.enkf = function(z.targ,z.pert,parm.pert,lambda.pic = 'lambda.min') {
#### ESTIMATE PARAMETERS FROM ENKF, EXPRESSED AS RIDGE REGRESSION
# INPUT:
#	Z.TARG[NTIME,NSPACE]: 		TARGET DATA
#	Z.PERT[NTIME,NSPACE,NEXP]: 	PERTURBATION RUNS FOR 'NEXP' CASES
#	PARM.PERT[NEXP,NPARM]: 		PERTUBATION VALUES FOR 'NEXP' CASES
#	LAMBDA.PIC = 'lambda.min' or 'lambda.1se' FOR SELECTING LAMBDA IN CV.GLMNET
# OUTPUT:
#	PARM.EST[NPARM]: ESTIMATE OF EACH PARAMETER
#	PARM.CI[NPARM,2]: 2.5 - 97.5% INTERVAL OF EACH PARAMETER
#	PARM.UPDATE[NPARM,NEXP]: ENSEMBLE OF UPDATED PARAMETERS, CONSISTENT WITH KF


tmax.targ = dim(z.targ)[1]
tmax.pert = dim(z.pert)[1]
nvar      = dim(z.pert)[2]
nexp      = dim(z.pert)[3]
nparm     = dim(parm.pert)[2]
if (nvar != dim(z.targ)[2]) stop('z.targ and z.pert have inconsistent number of variables')
if (nexp != dim(parm.pert)[1]) stop('inconsistent number of experiments in parm.pert and z.pert')


o.vec          = colMeans(z.targ)						# [nspace]
stdev          = sqrt(rowMeans((t(z.targ) - o.vec)^2))	# [nspace]

f.mat          = colMeans(z.pert)						# [nspace      ,nexp]
f.mat          = rbind(f.mat,t(parm.pert))				# [nspace+nparm,nexp]
mu.b           = rowMeans(f.mat)						# [nspace+nparm]
f.mat          = (f.mat - mu.b)/sqrt(nexp-1)			# [nspace+nparm,nexp]
x.mat          = f.mat[1:nvar,]/stdev					# [nspace      ,nexp]
y.vec          = (o.vec - mu.b[1:nvar])/stdev			# [nspace]

#############################
###### SPECIFIC LAMBDA BOUNDS
#############################
stdv.y         = sd(y.vec)
nsamp          = length(y.vec)
x.svd          = svd(x.mat)
lambda.bot     = x.svd$d[length(x.svd$d)-1]^2 * stdv.y / nsamp
lambda.top     = x.svd$d[1]^2 * stdv.y / nsamp
lambda.seq     = 10^seq(from=ceiling(log10(lambda.top))+1,to=floor(log10(lambda.bot)),length.out=100)


#############################
###### GLMNET
#############################
cv.glmnet.list = cv.glmnet(x.mat,y.vec,alpha=0,standardize=FALSE,intercept=FALSE,thresh=1.e-20,lambda=lambda.seq)
# cv.glmnet.list = cv.glmnet(x.mat,y.vec,alpha=0,standardize=FALSE,intercept=FALSE,thresh=1.e-20)
plot(cv.glmnet.list)
lambda         = cv.glmnet.list[[lambda.pic]]
beta.ridge     = coef(cv.glmnet.list,s=lambda.pic)[-1]
parm.est.ridge = f.mat[-(1:nvar),] %*% beta.ridge + mu.b[-(1:nvar)]

### COMPUTE BETA BY BRUTE FORCE
# https://stats.stackexchange.com/questions/129179/why-is-glmnet-ridge-regression-giving-me-a-different-answer-than-manual-calculat

y.stdv           = sqrt(mean((y.vec-mean(y.vec))^2))
l.star           = nrow(x.mat)/y.stdv * lambda
d.mat            = chol2inv(chol( t(x.mat) %*% x.mat + l.star * diag(1,ncol=ncol(x.mat),nrow=ncol(x.mat))))
beta.ridge.brute = as.numeric(d.mat %*% (t(x.mat) %*% y.vec))
cov.analysis     = l.star * f.mat[-(1:nvar),] %*% d.mat %*% t(f.mat[-(1:nvar),])
parm.stdv.ridge  = sqrt(diag(cov.analysis))
parm.ci.ridge    = as.numeric(parm.est.ridge) + parm.stdv.ridge %*% t(c(-2,2))

### COMPUTE ANALYSIS ENSEMBLE (NOT EFFICIENT-- CHOLESKY AND EIGEN ARE BOTH COMPUTED)
d.eigen          = eigen( l.star * d.mat)
d.sqrt           = d.eigen$vectors %*% (t(d.eigen$vectors) * sqrt(d.eigen$values) )
parm.update      = f.mat[-(1:nvar),] %*% d.sqrt * sqrt(nexp-1) + as.numeric(parm.est.ridge)

### CHECK THAT BRUTE FORCE RIDGE AGREES WITH GLMNET RIDGE, EXCEPT WHEN GLMNET SETS BETA = 0
# if (any(abs(beta.ridge) > 1.e-10)) if (any(abs(beta.ridge.brute-beta.ridge) > 1.e-5)) {
	# print(cbind(beta.ridge,beta.ridge.brute,beta.ridge-beta.ridge.brute))
	# print('inconsistent beta.ridge')
# }
list(parm.est=parm.est.ridge,parm.ci=parm.ci.ridge,parm.update=parm.update,parm.stdv.ridge=parm.stdv.ridge,cv.glmnet.list=cv.glmnet.list)	

}