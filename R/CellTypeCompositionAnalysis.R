library(lme4)
library(Matrix)
library(numDeriv)

# cell count table
Y = as.matrix(read.table("Y.txt",header=T,sep="\t",as.is=T))
colnames(Y)=unlist(strsplit(readLines("Y.txt",n=1),"\t")) # rename cell types

# meta data table
metadata = read.table("metadata.txt",header=T,sep="\t",as.is=T)

# number of samples / number of cell types
nsamples = nrow(Y)
ncells = ncol(Y)

# repeating the meta data table by the number of cell types
metadataExp=cbind(metadata[rep(match(rownames(Y),as.character(metadata$Sanger)),ncells),],Celltype=rep(colnames(Y),rep(nsamples,ncells)))

# poisson regression model
res.prop=glmer(I(c(Y))~
(1|Celltype)
+(1|Sample)
+(1|Donor)
+(1|Sex)
+Age
+(1|Covid_status)
+(1|Tissue)
+(1|X10x)
+(1|Technician)
+(1|Batch)
+(1|Date_process)
+(1|Ethnicity)

+(1|Sample:Celltype)
+(1|Donor:Celltype)
+(1|Sex:Celltype)
+(Age-1|Celltype)
+(1|Tissue:Celltype)
+(1|Ethnicity:Celltype)
+(1|X10x:Celltype)
+(1|Technician:Celltype)
+(1|Batch:Celltype)
+(1|Date_process:Celltype)
+(1|Covid_status:Celltype)
,
family=poisson,data=metadataExp,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

# standard errors of standard deviations (squre root of the variance parameters)
devfun = update(res.prop, devFunOnly=T)
pars = getME(res.prop, c("theta","fixef"))
hess = hessian(devfun, unlist(pars))
sdse.prop = data.frame(sd=unlist(pars), se=sqrt(diag(solve(hess))))

# posterior means and their standard deviations
res.prop.ranef = ranef(res.prop)

# Plots directory
system("mkdir Plots")

# Forest plot
source("Forest.R")
png(width=1200,height=1000,res=300,file="Plots/forest.png")
par(mar=c(3,6,1,1),mgp=c(1.2,0.5,0),family="Liberation Sans")
Forest(sdse.prop[seq(nrow(sdse.prop))%in%grep("Celltype",rownames(sdse.prop)),],xlim=c(-0.5,1.5))
dev.off()

# Posterior mean and standard deviation
source("getCondVal.R")
postmean = cbind(
    getCondVal(res.prop.ranef,"Covid_status:Celltype",ncells)[[1]][,c(2,1,3)],
    NA,
    getCondVal(res.prop.ranef,"Tissue:Celltype",ncells)[[1]],
    NA,
    getCondVal(res.prop.ranef,"Celltype",ncells)[[1]][,1,drop=F], # effect sizes for Age
    NA,
    getCondVal(res.prop.ranef,"Ethnicity:Celltype",ncells)[[1]]
)

lfsr = cbind(
    getCondVal(res.prop.ranef,"Covid_status:Celltype",ncells)[[2]][,c(2,1,3)],
    NA,
    getCondVal(res.prop.ranef,"Tissue:Celltype",ncells)[[2]],
    NA,
    getCondVal(res.prop.ranef,"Celltype",ncells)[[2]][,1,drop=F], # effect sizes for Age
    NA,
    getCondVal(res.prop.ranef,"Ethnicity:Celltype",ncells)[[2]]
)

# dotplot
source("col.rb.R")
source("drawDendrogram.R")
source("Dotplot.R")
png(width=2000,height=1800,file="Plots/dot.png",res=300)
par(mar=c(4,8,2,12),family="Liberation Sans")
Dotplot(postmean, SORT=c(F,T),zlim=c(log(1/3),log(3)),ltsr=1-lfsr)
dev.off()



