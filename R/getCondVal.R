# obtaining posterior mean and
# local false sign rate (lfsr)

getCondVal <- function(res.prop.ranef, id, ncells, nfactors=2){
	tmp = data.frame(res.prop.ranef)[data.frame(res.prop.ranef)[[1]]==id,]
	if(length(grep(":",tmp$grp))==0){
		cnam = matrix(as.character(tmp$term),ncells)[1,]	
		rnam = matrix(as.character(tmp$grp),ncells)[,1]	
	}else if(nfactors==2){
		cnam = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),2)[1,],ncells)[1,]
		rnam = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),2)[2,],ncells)[,1]
	}else if(nfactors==3){
		cnam1 = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[1,],ncells)[1,]
		cnam2 = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[2,],ncells)[1,]
		cnam = paste(cnam1,cnam2,sep=":")
		rnam = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[3,],ncells)[,1]
	}
	condval = matrix(tmp$condval,ncells)
	condsd  = matrix(tmp$condsd, ncells)
	rownames(condval)=rownames(condsd)=rnam
	colnames(condval)=colnames(condsd)=cnam
	lfsr = pnorm(condval,0,condsd)
	lfsr[lfsr>0.5]=1-lfsr[lfsr>0.5]
	list(condval=condval, lfsr=lfsr)
}

getCondVal=function(res.prop.ranef, id, ncells, nfactors=2, celltype){
        tmp = data.frame(res.prop.ranef)[data.frame(res.prop.ranef)[[1]]==id,]
        if(length(grep(":",tmp$grp))==0){
                cnam = matrix(as.character(tmp$term),ncells)[1,]
                rnam = matrix(as.character(tmp$grp),ncells)[,1]
        }else if(nfactors==2){
a=unique(matrix(unlist(strsplit(as.character(tmp$grp),":")),2)[1,])
b=unique(matrix(unlist(strsplit(as.character(tmp$grp),":")),2)[2,])
tmp=tmp[match(paste(a[rep(1:length(a),rep(length(b),length(a)))],b[rep(1:length(b),length(a))],sep=":"),tmp$grp),]
tmp$grp=paste(a[rep(1:length(a),rep(length(b),length(a)))],b[rep(1:length(b),length(a))],sep=":")
                cnam = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),2)[1,],ncells)[1,]
                rnam = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),2)[2,],ncells)[,1]
        }else if(nfactors==3){
a=unique(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[1,])
b=unique(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[2,])
C=unique(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[3,])
ab=paste(a[rep(1:length(a),rep(length(b),length(a)))],b[rep(1:length(b),length(a))],sep=":")
abc=paste(ab[rep(1:length(ab),rep(length(C),length(ab)))],C[rep(1:length(C),length(ab))],sep=":")
tmp=tmp[match(abc,tmp$grp),]
tmp$grp=abc
                cnam1 = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[1,],ncells)[1,]
                cnam2 = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[2,],ncells)[1,]
                cnam = paste(cnam1,cnam2,sep=":")
                rnam = matrix(matrix(unlist(strsplit(as.character(tmp$grp),":")),3)[3,],ncells)[,1]
        }
        condval = matrix(tmp$condval,ncells)
        condsd  = matrix(tmp$condsd, ncells)
	condval[is.na(condval)]=0
	condsd[is.na(condsd)]=1
        rownames(condval)=rownames(condsd)=rnam
	condval=condval[match(celltype,rnam),]
	condsd =condsd[match(celltype,rnam),]
        colnames(condval)=colnames(condsd)=cnam
        lfsr = pnorm(condval,0,condsd)
        lfsr[is.na(lfsr)]=1
        lfsr[lfsr>0.5]=1-lfsr[lfsr>0.5]
        list(condval=condval,  lfsr=lfsr)
}
