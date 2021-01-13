Forest <-
function(x,ord=T,labs=rownames(x),xlim=c(-1,1.5)){
	x=x[rownames(x)!="theta.Celltype.(Intercept)",]
	if(ord){labs=labs[rev(order(x[,1]))]; x=x[rev(order(x[,1])),]}
	labs=gsub("Celltype.","",gsub("X10x","10x_kit",gsub("Sample","residual",gsub(":Celltype.\\(Intercept\\)","",gsub("theta.","",rownames(x))))))
	x=rbind(x[labs!="residual",],x[labs=="residual",])
	labs=c(labs[labs!="residual"],labs[labs=="residual"])
	n=nrow(x)
	x=cbind(x,x[,1]-x[,2]*1.96,x[,1]+x[,2]*1.96);
	plot(x[,1],n:1,type="n",axes=F,xlab="",ylab="",xlim=xlim)
	points(x[,1],n:1,pch=15)
	axis(1)
	#segments(xlim[1],(n:1)[x[,3]<xlim[1]],xlim[1]+0.05,(n:1)[x[,3]<xlim[1]]+0.1)
	#segments(xlim[1],(n:1)[x[,3]<xlim[1]],xlim[1]+0.05,(n:1)[x[,3]<xlim[1]]-0.1)
	#x[x[,3]<xlim[1],3]=xlim[1]
	segments(x[,3],n:1,x[,4],n:1)
	abline(v=0,lty=2)
	par(xpd=NA)
	text(rep(xlim[1],nrow(x)),n:1,labs,pos=2)
	x
}
