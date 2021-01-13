drawDendrogram <-
function(tree, a, b, transpose=F){
	m = tree$merge
	o = tree$order
	h = tree$height/max(tree$height)
	
	# number of leavs
	K = length(o)
	# number of nodes
	N = nrow(m)
	
	x = (matrix(match(-m, o), N)-1)/(K-1)
	y = matrix(rep(0, N*2), N)
	
	#plot(1,1,xlim=c(0,1),ylim=c(0,1),type="n")
	if(transpose==F){
		for(i in 1:N){
			if(is.na(x[i,1])){
				x[i,1] = mean(x[m[i,1],])
				y[i,1] = h[m[i,1]]
			}
			if(is.na(x[i,2])){
				x[i,2] = mean(x[m[i,2],])
				y[i,2] = h[m[i,2]]
			}
			lines(x[i,c(1,1,2,2)]*b[1]+a[1], c(y[i,1],h[i],h[i],y[i,2])*b[2]+a[2])
		}
	}else{
		for(i in 1:N){
			if(is.na(x[i,1])){
				x[i,1] = mean(x[m[i,1],])
				y[i,1] = h[m[i,1]]
			}
			if(is.na(x[i,2])){
				x[i,2] = mean(x[m[i,2],])
				y[i,2] = h[m[i,2]]
			}
			lines(c(y[i,1],h[i],h[i],y[i,2])*b[2]+a[2], x[i,c(1,1,2,2)]*b[1]+a[1])
		}
	}
}
