rm(list=ls())

lplotfile     = TRUE	### DIRECT FIGURES TO FILE (TRUE) OR ON-SCREEN (FALSE)
npoly         = 0		### ORDER OF THE POLYNOMIAL TO REMOVE
max.lag       = 1		### ORDER OF THE AR MODEL (CALLED 'P' IN THE PAPER)
nharm         = 5		### NUMBER OF ANNUAL HARMONICS (CALLED 'H' IN THE PAPER)
nlap.varx     = 20		### NUMBER OF SPHERICAL HARMONICS IN THE VARX MODEL (CALLED 'S' IN THE PAPER)
alpha         = 0.05	### SIGNIFICANCE LEVEL
lpic.plot     = 1		### THE LAPLACIAN FOR PLOTTING TIME SERIES
max.time.pert =  2 * 12	### LENGTH OF PERTURBATION RUN (IN MONTHS)
max.time.targ = 50 * 12	### LENGTH OF TARGET TIME SERIES (IN MONTHS)
nens.use      = 'all'	### NUMBER OF ENSEMBLE MEMBERS TO USE (A NUMBER, OR 'ALL')
ne.skip       = 0		### LENGTH OF PERTURBATION RUN TO SKIP (IN MONTHS)
no.skip       = 0		### LENGTH OF TARGET RUN TO SKIP (IN MONTHS)

check.analyticals = FALSE	### SHOULD 'NLM' CHECK GRADIENT FUNCTION?  
###	(recommendation: set this FALSE-- otherwise it is prone to false errors near unstable solutions)

pcntrl        = c(0.28,200)				### VALUE OF THE CESM PARAMETERS IN TARGET RUN
pnames        = c('clubbgamma','dcs')	### NAME OF THE CESM PARAMETERS

ptarg         = pcntrl					### 'PTARG' STANDS FOR 'PARAMETERS IN TARGET RUN'
nparm         = length(pnames)			### NUMBER OF TUNABLE PARAMETERS IN CESM

# dir.Rlib   = '/Users/delsole/R/delsole_tools/'### DIRECTORY OF R FUNCTION
# dir.data   = '/Users/delsole/data/CESM2/'		### DIRECTORY OF CESM DATA
dir.Rlib   = NULL								### DIRECTORY OF R FUNCTION
dir.data   = './'								### DIRECTORY OF CESM DATA


source(paste(dir.Rlib,'pdf.eps.R',sep=''))						### FUNCTION FOR MAKING PLOTS
source(paste(dir.Rlib,'parameter.est.xymat.R',sep=''))
source(paste(dir.Rlib,'loglikelihood.gradient.R',sep=''))
source(paste(dir.Rlib,'parameter.est.cGLS.R',sep=''))
source(paste(dir.Rlib,'parameter.est.allobs.R',sep=''))
source(paste(dir.Rlib,'kf.lm.R',sep=''))
source(paste(dir.Rlib,'gev.R',sep=''))
source(paste(dir.Rlib,'parameter.est.enkf.R',sep=''))


library(glmnet)

pnames.short = pnames
for (np in 1:nparm) {
	if (pnames[np] == 'clubbgamma') pnames.short[np] = 'gamma'
}


###############################
######## SPECIFY PREFIX
###############################
prefix.emulate = paste(pnames[1],ptarg[1],'.',pnames[2],ptarg[2],'.TT',max.time.targ,'.TP',max.time.pert,sep='')

prefix.compare = paste('cesm.compareVARX.p',max.lag,'.S',nlap.varx,'.No',max.time.targ/12,'.Ne',max.time.pert/12,sep='')


###############################
######## GET CONTROL DATA
###############################
fcontrol   = 'lapl.GLO.v0.b.e21.B1850.f19_g17.CMIP6-piControl-2deg.001.cam.h0.SST.000101-050012.RData'
fread    = paste(dir.data,fcontrol,sep='')
load(fread)
pc.ctrl  = pc.list$pc
ntot.ctr = dim(pc.ctrl)[1]
nlap     = dim(pc.ctrl)[2]

if (max.time.targ > ntot.ctr) stop('max.time.targ is too big')

###############################
######## GET PERTURBATION PARAMETERS
####### ASSUMES FILE NAME IS OF THE STANDARD FORM 
####### lapl.GLO.v0.clubb_gamma_0.26_micro_mg_100.D-6.SST.000101-001012.RData
###############################
fnames   = list.files(path=dir.data,pattern='lapl.GLO.v0.clubb_gamma*')
nexp     = length(fnames)
prun.all = array(NA,dim=c(nexp,nparm))
colnames(prun.all) = pnames
for (n in 1:nexp) prun.all[n,1] = as.numeric(strsplit(fnames[n],split='_')[[1]][3])
for (n in 1:nexp) prun.all[n,2] = as.numeric(strsplit(strsplit(fnames[n],split='_')[[1]][6],split='[.]')[[1]][1])

###############################
######## GET PERTURBATION RUNS
###############################
for (n in 1:nexp) {
	fread = paste(dir.data,fnames[n],sep='')
	load(fread)
	if (n == 1) {
		ntot.per = dim(pc.list$pc)[1]
		pc.pert  = array(NA,dim=c(ntot.per,nlap,nexp))
		if (max.time.pert > ntot.per) stop('max.time.pert is too large')
	}
	pc.pert[,,n] = pc.list$pc
}

pvals       = unique(prun.all[,1])
col.vals    = colorRampPalette(c('royalblue','red'))(length(pvals))
col.pic     = NULL
for (n in 1:nexp) col.pic = c(col.pic,col.vals[which(prun.all[n,1] == pvals)])

if (nens.use == 'all') nens.use = nexp



###############################
######## PLOT RAW TIME SERIES
###############################
yrange = range(pc.ctrl[,lpic.plot],pc.pert[,lpic.plot,])
xrange = c(1,ntot.per+10)

par(mfcol=c(2,1),mar=c(5,5,3,1))
plot(1,1,type='n',xlim=xrange,ylim=yrange,xlab='months',ylab=paste('Laplacian',lpic.plot,sep=''))
for (n in 1:nexp) lines(pc.pert[,lpic.plot,n],col=col.pic[n],lwd=2)
for (n in 1:nexp) text(ntot.per,pc.pert[ntot.per,lpic.plot,n],prun.all[n,1],pos=4,col=col.pic[n])
lines(pc.ctrl[,lpic.plot],col=1,lwd=3)


###############################
######## DEFINE PRIOR FOR REGULARIZED METHODS
###############################
theta.mean.prior       = colMeans(prun.all)
theta.covm.prior       = cov(prun.all)

theta.mean.prior=c(0.3,400)
theta.covm.prior=diag(c(0.1^2,300^2))

theta.mean.prior = NULL
theta.covm.prior = NULL





###############################
######## STORE DATA IN A LIST (USED IN REMAINDER OF CODE)
###############################
tso     = pc.ctrl[1:max.time.targ + 12*no.skip,1:nlap.varx]
tse     = pc.pert[1:max.time.pert + 12*ne.skip,1:nlap.varx,1:nens.use]

ts.list = list(list(tseries = tso, first.step = 1, theta = rep(NA,nparm)))
for (ne in 1:nens.use) ts.list[ne+1] = list(list(tseries=tse[,,ne], first.step = 1, theta = prun.all[ne,]))


###############################
######## DEFINE GAMMA MATRIX, WHICH IS USED AS THE REFERENCE IN THE LIKELIHOOD FUNCTION.
######## THIS AVOIDS COMPUTING LOGARITHMS OF DETERMINANTS. INSTEAD, LOGS OF EIGENVALUES ARE COMPUTED
###############################
xymat.list          = parameter.est.xymat(ts.list,max.lag,nharm,npoly,period=12,prior.mean = theta.mean.prior, prior.covm = theta.covm.prior) 
cGLS                = parameter.est.cGLS(ts.list,pnames)
gamma.first.guess   = loglikelihood.gradient(cGLS$theta.o.mean,output = 'gamma')


################################################################
######## EVALUATE ALL METHODS
###############################################################
lead.all       = c( 2, 3, 4, 5,10)
ens.all        = c(10,15,20,25,30,35)
nlead          = length(lead.all)
nens           = length(ens.all)
theta.mean.all = array(NA,dim=c(nparm,5,nens,nlead))
theta.stdv.all = array(NA,dim=c(nparm,5,nens,nlead))
theta.mean.prior = c(0.3,400)
theta.covm.prior = diag(c(0.1^2,300^2))



for (nl.loop in 1:nlead) for (ne.loop in 1:nens) {
	set.seed(1)
	lead = lead.all[nl.loop]*12
	ens  =  ens.all[ne.loop]
	
	tso     = pc.ctrl[1:max.time.targ,1:nlap.varx ]
	tse     = pc.pert[1:lead         ,1:nlap.varx,1:ens]
	
	ts.list = list(list(tseries = tso, first.step = 1, theta = rep(NA,nparm)))
	for (ne in 1:ens) ts.list[ne+1] = list(list(tseries=tse[,,ne], first.step = 1, theta = prun.all[ne,]))
	
	cGLS   = parameter.est.cGLS(ts.list,pnames)
	allobs = parameter.est.allobs(tso,tse,prun.all[1:ens,],theta.mean.prior=theta.mean.prior,theta.covm.prior=theta.covm.prior,augment.cov=FALSE)
	ridge  = parameter.est.enkf(tso[1:max.time.pert,],tse[,,1:ens],prun.all[1:ens,])
	
	### THE NEXT THREE FUNCTIONS, CALLED IN SEQUENCE, EVALUATE THE FULL-MLE AND THE INFORMATION MATRIX
	xymat.list  = parameter.est.xymat(ts.list, max.lag, nharm, npoly, period=12, prior.mean = NULL, prior.covm = NULL)
	nlm.list    = nlm(loglikelihood.gradient,cGLS$theta.o.mean,steptol=1.e-10,check.analyticals=check.analyticals)		
	info.list   = loglikelihood.gradient(nlm.list$estimate,output='information')
	
	theta.mean.all[,1,ne.loop,nl.loop] = nlm.list$estimate
	theta.mean.all[,2,ne.loop,nl.loop] = allobs$theta.mean.end
	theta.mean.all[,3,ne.loop,nl.loop] = ridge$parm.est
	theta.mean.all[,4,ne.loop,nl.loop] = cGLS$theta.o.mean
	theta.mean.all[,5,ne.loop,nl.loop] = cGLS$theta.o.mean
	
	
	theta.stdv.all[,1,ne.loop,nl.loop] = sqrt(diag(info.list$info.mat))
	theta.stdv.all[,2,ne.loop,nl.loop] = allobs$theta.sd.end
	theta.stdv.all[,3,ne.loop,nl.loop] = ridge$parm.stdv.ridge
	theta.stdv.all[,4,ne.loop,nl.loop] = cGLS$theta.o.sd
	theta.stdv.all[,5,ne.loop,nl.loop] = cGLS$theta.o.sd
	print(paste(ne.loop,nl.loop,paste(signif(theta.mean.all[2,1:5,ne.loop,nl.loop],3),collapse=' '),
	      paste(signif(theta.stdv.all[2,1:5,ne.loop,nl.loop],3),collapse=' ')))

}

methods  = c('MLE','AllObs','Ridge','First Guess','cGLS')
nmethods = length(methods)

prefix.compare2 = paste('cesm.compareVARX.p',max.lag,'.S',nlap.varx,'.No',max.time.targ/12,'.H',nharm,sep='')


fout = paste(prefix.compare2,'.ESTvRunlength',sep='')
if (lplotfile) pdf.eps(fout,'pdf',width=8.5,height=7)
par(mfcol=c(nparm,1),mar=c(5,5,3,1))
for (np in 1:nparm) {
	if (np == 1) yrange = c(0.25,0.30) else yrange = c(100,600)
	xrange = range(lead.all %*% t(ens.all))
	plot(1,1,type='n',xlim=xrange,ylim=yrange,xlab='total length of perturbation runs (yrs)',ylab=pnames[np])
	abline(h=pcntrl[np],col='grey',lwd=2)
	for (nm in nmethods:1) for (nl in 1:nlead) arrows(
	   lead.all[nl] * ens.all,theta.mean.all[np,nm,,nl]+2*theta.stdv.all[np,nm,,nl],
	   lead.all[nl] * ens.all,theta.mean.all[np,nm,,nl]-2*theta.stdv.all[np,nm,,nl],col=nm,code=3,angle=90,length=0.035)
	for (nm in nmethods:1) for (nl in 1:nlead) points(lead.all[nl] * ens.all,theta.mean.all[np,nm,,nl],pch=19,col=nm,cex=0.5)
	for (nm in 1:nmethods) for (nl in 1:nlead) for (ne in 1:nens) if (theta.mean.all[np,nm,ne,nl] < yrange[1]) points(lead.all[nl]*ens.all[ne],yrange[1],pch=4,col=nm)
	for (nm in 1:nmethods) for (nl in 1:nlead) for (ne in 1:nens) if (theta.mean.all[np,nm,ne,nl] > yrange[2]) points(lead.all[nl]*ens.all[ne],yrange[2],pch=4,col=nm)
	if (np==1) legend('bottomright',legend=methods,col=1:nmethods,pch=19,lwd=2,ncol=2)
	ftitle.top = paste(pnames[np],'Estimates')
	ftitle.bot = paste('p= ',max.lag,'; S= ',nlap.varx,'; Nobs= ',max.time.targ/12,' yr; H= ',nharm,sep='')
	title(main=ftitle.top,line=2.0)
	title(main=ftitle.bot,line=0.5)		
}
if (lplotfile) dev.off()

#####################################
######## SAVE DATA
#####################################
fsave     = paste(prefix.compare2,'.RData',sep='')
list.save = list(theta.mean.all = theta.mean.all, theta.stdv.all = theta.stdv.all,
  max.lag  = max.lag, nlap.varx = nlap.varx, max.time.targ = max.time.targ,
  lead.all = lead.all, ens.all = ens.all, parm.true = pcntrl, pnames = pnames, 
  methods = methods)
save(list.save,file=fsave)


################################################
######### CHECK REPRODUCIBILITY
################################################
list.save2 = list.save
load(paste(prefix.compare2,'.check.RData',sep=''))
print(all.equal(list.save,list.save2))
