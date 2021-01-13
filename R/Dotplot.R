Dotplot <-
function(X, ltsr, col=col.rb, cex=1, zlim=c(-1,1), xlab="", ylab="", lab=rep("",ncol(X)), SORT=c(F,F), beta, border=F, add=F, xlim=c(0,ncol(X)), tick=F, side=2, w=0.03, ww=50, labx=NULL, laby=NULL){
	dotplot=function(x,y,z,col,add,zlim,axes,xlab,ylab,xlim=range(x),ylim=range(y),ltsr,cex=1){
		z=ceiling((z-zlim[1])/diff(range(zlim))*99)+1
		ltsr=-log10((1-ltsr)*2)*8
		par(xpd=NA)
		plot(rep(x,length(y)),rep(y,rep(length(x),length(y))),cex=sqrt(ltsr)*cex, pch=20, col=col[z],
			xlim=range(x)+c(0,0.5),axes=axes,xlab=xlab,ylab=ylab)
		par(xpd=F)
	}
    
    X[X>zlim[2]]=zlim[2]
    X[X<zlim[1]]=zlim[1]
    ltsr[ltsr>0.9999]=0.9999
    N0=N=ncol(X)
    M0=M=nrow(X)
	M1=N1=0

	if(SORT[1]){
		dendx = hclust(dist(t(X)))
		X = X[, dendx$ord]
		ltsr = ltsr[, dendx$ord]
	}
	if(SORT[2]){
		dendy = hclust(dist(X))
		X = X[dendy$ord,]
		ltsr = ltsr[dendy$ord,]
	}
	if(!is.null(labx)){
		X=rbind(X,array(NA,c(ncol(labx)+1, N)))
		rownames(X)[(M+2):(M+ncol(labx)+1)]=colnames(labx)
		M1=ncol(labx)+1
		M=M0+M1
	}
	if(!is.null(laby)){
		X=cbind(array(NA,c(M, ncol(laby)+1)),X)
		colnames(X)[(1):(ncol(laby))]=colnames(laby)
		N1=ncol(laby)+1
		N=N0+N1
	}

	if(ncol(X)>1){
		dotplot(1:N-0.5, 1:M-0.5, t(X[nrow(X):1,,drop=F]),col=col,zlim=zlim,axes=F,xlab=xlab,ylab=ylab,add=add, xlim=xlim,ltsr=t(ltsr[nrow(ltsr):1,,drop=F]), cex=cex)
	}else{
		dotplot(1:2-0.5, 1:M-0.5, t(cbind(X,NA)[nrow(X):1,,drop=F]),col=col,zlim=zlim,axes=F,xlab=xlab,ylab=ylab,add=add, xlim=xlim+c(0,1),ltsr=t(cbind(ltsr,NA)[nrow(X):1,,drop=F]), cex=cex)
	}
	if(!is.null(labx)){
		labx=cbind(NA,matrix(as.numeric(as.factor(unlist(lapply(labx,as.character)))),N0))
		if(SORT[1]){image(N1:N, 1:M1-0.5, (labx[dendx$ord,ncol(labx):1]),col=c(col24,col24)[1:max(labx,na.rm=T)], add=T)
		}else{image(N1:N, 1:M1-0.5, (labx[,ncol(labx):1]),col=c(col24,col24)[1:max(labx,na.rm=T)], add=T)}
	}
	if(!is.null(laby)){
		laby=cbind(matrix(as.numeric(as.factor(unlist(lapply(laby,as.character)))),M0),NA)
		if(SORT[2]){ image(1:N1-0.5, M1:M, t(laby[rev(dendy$ord),]),col=c(col24,col24)[1:max(laby,na.rm=T)], add=T)
		}else{image(1:N1-0.5, M1:M, t(laby[nrow(laby):1,]),col=c(col24,col24)[1:max(laby,na.rm=T)], add=T)}
	}

	if(border){
		for(i in (-1):(N)){segments(i+0.5,0.5,i+0.5,M+0.5)}
		for(i in (-1):(M)){segments(0.5,i+0.5,N+0.5,i+0.5)}
	}
	mtext(lab,1,at=1:N)
	if(tick){
		flag=!is.na(rownames(X))
		par(xpd=NA)
		segments(N,c(M:1-0.5)[flag],N+0.1,c(M:1-0.5)[flag])
		segments(N+0.1,c(M:1-0.5)[flag],N+1, c(M:1-0.5)[flag]*w + (rank((M:1-0.5)[flag])/sum(flag)*M-ww)*(1-w))
		text(N+1, c(M:1-0.5)[flag]*w + (rank((M:1-0.5)[flag])/sum(flag)*M-ww)*(1-w), rownames(X)[flag], pos=4, offs=0.1, cex=0.75)
		par(xpd=F)
	}else{
		mtext(rownames(X),at=M:1-0.5,side,las=2,line=0)
	}
	par(xpd=NA)
	text(1:ncol(X)-0.2,rep(-.5,ncol(X)),colnames(X),pos=2,srt=45,offs=0)
	if(SORT[1]) drawDendrogram(dendx, c(N1+0.5,M), c(N0-1,3))
	if(SORT[2]) drawDendrogram(dendy, c(M-0.5,N), c(-(M0-1),3), T)
	#Gauge
	ybase=-1
	rect(rep(N+1+3-0.25,100),(seq(100)-1)/100*6+10+ybase,rep(N+1+3+0.25,100),seq(100)/100*6+10+ybase,border=NA,col=col.rb)
	text(N+5.5,16.7+ybase,"Fold change",offs=0.5)
	text(rep(N+2+3,3), c(0.02,0.5,0.97)*6+10+ybase, c(paste("<1/",exp(zlim[2]),sep=""),"1",paste(">",exp(zlim[2]),sep="")),pos=4,offs=-0.5)
	#LTSR
	points(rep(N+1+3,5),1:5*1.2+ybase,cex=sqrt(-log10((1-c(0.501,0.9,0.99,0.999,0.9999))*2)*4)*cex,pch=20)
	text(rep(N+2+3,5),  1:5*1.2+ybase,c("0.5","0.9","0.99","0.999",">0.9999"),pos=4,offs=-0.3)
	text(N+3+3,7.2+ybase,"LTSR",offs=0.5)
	par(xpd=F)
}
