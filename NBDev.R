NBDev<-function(counts,design,log.offset,nb.disp,print.progress=TRUE){
n<-ncol(counts); 

### Function used to evaluate quasi-negative binomial likelihood for saturated model
SAT.LIKE<-function(counts,disp){
means<-counts
like<-disp*log(disp/(disp+means))
like[counts!=0]<-like[counts!=0]+counts[counts!=0]*log(means[counts!=0]/(disp+means[counts!=0]))
-sum(like)
}

### Function used to evaluate quasi-negative binomial likelihood at current parameter estimates
LIKE<-function(parms,design,counts,disp,est.offset){
means<-as.vector(exp(design%*%parms)*est.offset)
like<-disp*log(disp/(disp+means))
like[counts!=0]<-like[counts!=0]+counts[counts!=0]*log(means[counts!=0]/(disp+means[counts!=0]))
-sum(like)
}

### Function used to evaluate quasi-negative binomial likelihood gradient at current parameter estimates
GRAD<-function(parms,design,counts,disp,est.offset){
means<-as.vector(exp(design%*%parms)*est.offset)
colSums(-(counts-means*(disp+counts)/(disp+means))*design)
}

if(is.null(log.offset)) log.offset<-rep(0,ncol(counts));
est.offset<-exp(log.offset)

deviance.vector<-rep(NA,nrow(counts)); means<-matrix(NA,nrow(counts),ncol(counts));
parms<-matrix(NA,nrow(counts),ncol(design))

### For each gene and given design matrix, find model parameters (for mean structure) that optimize quasi-likelihood
for(i in 1:nrow(counts)){
### If wanted, provide running progress update (eventually once every 5000 genes) 
if(i%in%c(2,10,100,500,1000,2500,5000*(1:200))&print.progress) print(paste("Analyzing Gene #",i))

### Find initial parameter estimates using a linear model
if(ncol(design)>1)init.parms<-lm(log(counts[i,]+1)~design[,-1],offset=log(est.offset))$coefficients
if(ncol(design)==1)init.parms<-lm(log(counts[i,]+1)~1,offset=log(est.offset))$coefficients

### Find optimum parameter estimates
opt<-optim(init.parms,fn=LIKE,gr=GRAD,method="BFGS",design=design,
control=list(reltol=1e-25,maxit=1000),counts=counts[i,],disp=1/nb.disp[i],
est.offset=est.offset)

### Save optimized means (used in Pearson's dispersion estimator)
means[i,]<-as.vector(exp(design%*%opt$par)*est.offset)
parms[i,]<-opt$par

### Save deviance (used to compute LRT for comparing models and also deviance dispersion estimator)
deviance.vector[i]<-2*(opt$value-SAT.LIKE(counts[i,],1/nb.disp[i])) ##SAT.LIKE and opt$value are !negative! likelihoods
}

return(list(dev=deviance.vector,means=means,parms=parms))
}








