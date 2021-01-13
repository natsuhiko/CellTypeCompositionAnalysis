library(lme4)
library(Matrix)
library(numDeriv)

# cell count table
Y = as.matrix(read.table("../Data/Y.txt",header=T,sep="\t",as.is=T))
colnames(Y)=unlist(strsplit(readLines("../Data/Y.txt",n=1),"\t")) # rename cell types


# meta data table
metadata = read.table("../Data/metadata.txt",header=T,sep="\t",as.is=T)
metadata=metadata[match(rownames(Y),metadata$Sample),] # matching Y and metadata by Sample ID

# !!! EXAMPLE DATA SPECIFIC SCALING !!!
metadata$Cell_viability[is.na(metadata$Cell_viability)]=mean(metadata$Cell_viability,na.rm=T)
metadata$Age=scale(log10(metadata$Age))
metadata$Cell_viability=scale(metadata$Cell_viability)
metadata$Pool=as.character(metadata$Pool)

# number of samples / number of cell types
nsamples = nrow(Y)
ncells = ncol(Y)


# repeating the meta data table by the number of cell types
metadataExp=cbind(metadata[rep(match(rownames(Y),as.character(metadata$Sample)),ncells),],Celltype=rep(colnames(Y),rep(nsamples,ncells)))


# poisson regression model
# !!! EXAMPLE DATA SPECIFIC MODEL FORMULA !!!
res.prop=glmer(I(c(Y))~
(1|Celltype)
+(1|Sample)
+(1|Patient)
+(1|Sex)
+Age
+(1|COVID_status)
+(1|Tissue)
+(1|X10x)
#+(1|Technician)
+(1|Pool)
+(1|Date_of_process)
+(1|Ethnicity)
+(1|Smoker)
+Cell_viability

+(1|Sample:Celltype)
+(1|Patient:Celltype)
+(1|Sex:Celltype)
+(Age-1|Celltype)
+(1|Tissue:Celltype)
+(1|Ethnicity:Celltype)
+(1|X10x:Celltype)
#+(1|Technician:Celltype)
+(1|Pool:Celltype)
+(1|Date_of_process:Celltype)
+(1|COVID_status:Celltype)
+(1|Smoker:Celltype)
+(Cell_viability-1|Celltype)
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
# !!! EXAMPLE DATA SPECIFIC POSTERIOR COMPUTATION !!!
source("getCondVal.R")
postmean = cbind(
    getCondVal(res.prop.ranef,"Tissue:Celltype",ncells,celltype=colnames(Y))[[1]],
    NA,
    getCondVal(res.prop.ranef,"COVID_status:Celltype",ncells,celltype=colnames(Y))[[1]][,c(2,1,3)],
    NA,
    getCondVal(res.prop.ranef,"Celltype",ncells,celltype=colnames(Y))[[1]][,1,drop=F], # effect sizes for Cell viability
    NA,
    getCondVal(res.prop.ranef,"Celltype",ncells,celltype=colnames(Y))[[1]][,2,drop=F] # effect size for Age
)

lfsr = cbind(
    getCondVal(res.prop.ranef,"Tissue:Celltype",ncells,celltype=colnames(Y))[[2]],
    NA,
    getCondVal(res.prop.ranef,"COVID_status:Celltype",ncells,celltype=colnames(Y))[[2]][,c(2,1,3)],
    NA,
    getCondVal(res.prop.ranef,"Celltype",ncells,celltype=colnames(Y))[[2]][,1,drop=F], # SD for Cell viability
    NA,
    getCondVal(res.prop.ranef,"Celltype",ncells,celltype=colnames(Y))[[2]][,2,drop=F] # SD for Age
)

# Dotplot
source("col.rb.R")
source("drawDendrogram.R")
source("Dotplot.R")
png(width=2000,height=2200,file="Plots/dot.png",res=300)
par(mar=c(4,8,2,12),family="Liberation Sans")
Dotplot(postmean, SORT=c(F,T),zlim=c(log(1/3),log(3)),ltsr=1-lfsr, cex=0.8)
dev.off()



