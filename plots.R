library(fitdistrplus)

plot_norm<-function(plot_dens=T,n=500,mu=100,sig=1.5,S=1.5,L=2.5){
	set.seed(260390)
	dat<-rnorm(n,mu,sig)
	h<-hist(dat,ann=F,main="",font.axis=2,col="red",ylim=c(0,50),breaks=50)
	mtext(text="Fluorescence level (AU)",side=1,line=L,cex=S,font=2)
	mtext(text="Counts",side=2,line=L,cex=S,font=2)

	if(plot_dens){
		params<-fitdist(data=dat,method="mle",dist="norm")
		xfit<-seq(min(dat)-1,max(dat)+1,length=1000)
		yfit<-dnorm(xfit,params$estimate[1],params$estimate[2])
		yfit<-yfit*diff(h$mids[1:2]*length(dat))
		lines(xfit,yfit,col="blue",lwd=3)
	}
}





plot_coin_flip<-function(n_flip=10,grid_size=1e5,p=0.7,S=2,L=1.5){
	set.seed(42)
	flips<-rep(NA,n_flip)
	for(i in 1:n_flip) flips[i]<-rbinom(n=1,size=1,prob=p)

	p_grid<-seq(from=0,to=1,length.out=grid_size)
	if(n_flip==0){
		likelihood<-rep(1,grid_size)
	}else{
		likelihood<-dbinom(x=sum(flips),size=n_flip,prob=p_grid)
	}

	plot(p_grid
	     ,likelihood
	     ,type="l"
	     ,xlab="probability of H (p)"
	     ,ylab=""
	     ,yaxt="n"
	     ,cex.lab=S
	     ,font.lab=2
	     ,font.axis=2)
	mtext(side=2,text="likelihood",cex=S,line=L,font=2)
	mtext(side=3,text=paste("H:",sum(flips)," T:",n_flip-sum(flips)),cex=S,line=L,font=2)

}
