sdplot = function(x, y, sdx, sdy, col, xlab="", ylab=""){
	col=c(col)
	plot(x, y, pch=20, col=col, xlim=range(c(x-sdx*2,x+sdx*2)), ylim=range(c(y-sdy*2,y+sdy*2)),axes=F, xlab=xlab, ylab=ylab)
	segments(x, y-sdy*2, x, y+sdy*2, col=col)
	segments(x-sdx*2, y, x+sdx*2, y, col=col)
	axis(2,las=2)
	axis(1)
	abline(0,1,lty=2)
	text(min(x-2*sdx),max(y+sdy),substitute(paste(italic(R),"=",x), list(x = round(cor(x, y)*100)/100)),pos=4)
}
