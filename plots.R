library(fitdistrplus)
library(bayestestR)
source("martincolourscale.R")

plot_norm<-function(plot_dens=T,n=500,mu=100,sig=1.5,S=1.5,L=2.5){
	set.seed(260390)
	dat<-rnorm(n,mu,sig)
	h<-hist(dat,ann=F,main="",font.axis=2,col=martincolourscale[1],ylim=c(0,50),breaks=50,yaxt="n")
	Axis(side=2, labels=F)
	mtext(text="Fluorescence level (AU)",side=1,line=L,cex=S,font=2)

	if(plot_dens){
		params<-fitdist(data=dat,method="mle",dist="norm")
		xfit<-seq(min(dat)-1,max(dat)+1,length=1000)
		yfit<-dnorm(xfit,params$estimate[1],params$estimate[2])
		yfit<-yfit*diff(h$mids[1:2]*length(dat))
		lines(xfit,yfit,col=martincolourscale[3],lwd=3)
		mtext(text="PDF",side=2,line=L,cex=S,font=2)

	}else{
		mtext(text="Counts",side=2,line=L,cex=S,font=2)
	}
}

plot_norm_dens<-function(S=2){
	par(mar=c(1,5,1,1))
	gridd<-seq(from=1,to=100,length.out=2000)
	plot(x=gridd
		,ylim=c(0,0.45)
		,col="white"
		,xaxt="n"
	 	,xlab=""
	     	,ylab="PDF"
		,cex=S
		,font.axis=2
	    	,font.lab=2
		,cex.axis=S
		,cex.lab=S)
	sig<-c(1,5,10)
	mu<-c(75,20,50)
	for(i in 1:3){
		points(dnorm(gridd,mean=mu[i],sd=sig[i]),col=martincolourscale[i],type="l",lwd=3)
	}
	legend("topright",legend="Normal(μ,σ)",bty="n",cex=S,text.font=2)
}


plot_binom_dens<-function(S=2,S_pnt=1.5){
	par(mar=c(1,5,1,1))
	gridd<-seq(from=1,to=100,length.out=100)
	plot(x=gridd
		,ylim=c(0,0.2)
		,col="white"
		,xaxt="n"
	 	,xlab=""
	     	,ylab="PMF"
		,cex=S
		,font.axis=2
	    	,font.lab=2
		,cex.axis=S
		,cex.lab=S)
	n=c(40,40,80)
	p=c(0.5,0.7,0.8)
	for(i in 1:3){
		dens<-dbinom(gridd,size=n[i],prob=p[i])
		points(dens,col=martincolourscale[i],type="p",cex=S_pnt,pch=19)
	}
	legend("topright",legend="Binomial(n,p)",bty="n",cex=S,text.font=2)
}

plot_beta_dens<-function(S=2){
	par(mar=c(1,5,1,1))
	gridd<-seq(from=0,to=1,length.out=2000)
	plot(x=gridd
		,ylim=c(0,1.8)
		,col="white"
		,xaxt="n"
	 	,xlab=""
	     	,ylab="PDF"
		,cex=S
		,font.axis=2
	    	,font.lab=2
		,cex.axis=S
		,cex.lab=S)
	Alpha<-c(0.1,1,3)
	Beta<-c(0.1,1,2)
	for(i in 1:3){
		points(dbeta(gridd,shape1=Alpha[i],shape2=Beta[i]),col=martincolourscale[i],type="l",lwd=3)
	}
	legend("topright",legend="Beta(α,β)",bty="n",cex=S,text.font=2)
}

plot_binom<-function(S=2,S_pnt=1.5){
	par(mar=c(5,5,1,1))
	gridd<-seq(from=0,to=1,length.out=2)	
	cent<-barplot(height=dbinom(gridd,size=1,prob=0.7)
		,cex=S_pnt
		,cex.axis=S
		,cex.lab=S
		,font.axis=2
		,font.lab=2
		,xlab=""
		,ylab="PMF"
		,col=martincolourscale[1]
		,pch=19
		,ylim=c(0,1))
	legend("topleft",inset=c(-0.1,0),legend="Binomial(n = 1,p = 0.7)",bty="n",cex=S,text.font=2)
	axis(side=1,at=cent,labels=c("Tails","Heads"),line=1.5,cex.axis=S,font.axis=2)
}

plot_binom_sample<-function(trials=20,S=2,S_pnt=1.5){
	set.seed(3)
	par(mar=c(5,5,1,1))
	samples<-rbinom(n=trials,size=1,prob=0.7)
	samples<-table(factor(samples,levels=0:1))
	ylims<-max(samples)*1.2
	cent<-barplot(height=samples
		,cex=S_pnt
		,cex.axis=S
		,cex.lab=S
		,font.axis=2
		,font.lab=2
		,xlab=""
		,xaxt="n"
		,ylab="Number of Heads and Tails"
		,col=martincolourscale[1]
		,pch=19
		,ylim=c(0,ylims))
	legend("topleft",inset=c(-0.1,0),legend=paste("Binomial(n = ",trials,",p = 0.7)"),bty="n",cex=S,text.font=2)
	legend("topleft",inset=c(-0.1,0.1),legend=paste("emp p= ",round(samples[2]/trials,digits=2)),bty="n",cex=S,text.font=2)
	axis(side=1,at=cent,labels=c("Tails","Heads"),line=1.5,cex.axis=S,font.axis=2)
}



plot_coin_flip<-function(n_flip=10,grid_size=1e5,p=0.7,S=2,L=1.5){
	set.seed(42)
	flips<-rbinom(n=n_flip,size=1,prob=p)

	p_grid<-seq(from=0,to=1,length.out=grid_size)
	if(n_flip==0){
		likelihood<-rep(1,grid_size)
		MLE<-NA
	}else{
		likelihood<-dbinom(x=sum(flips),size=n_flip,prob=p_grid)
		MLE<-p_grid[which.max(likelihood)]
		MLE<-round(MLE,digits=2)
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
	mtext(side=3,text=paste("Flips:",n_flip," H:",sum(flips)," T:",n_flip-sum(flips)),cex=S,line=L,font=2)
	legend("topright",legend=paste("MLE:",MLE),bty="n",cex=S,text.font=2)

}

plot_coin_flip_prior<-function(n_flip=50,a=1,b=1,grid_size=1e5,p=0.7,S=2,L=1.5){
	set.seed(42)

	par(mfrow=c(1,3))

	flips<-rbinom(n=n_flip,size=1,prob=p)

	p_grid<-seq(from=0,to=1,length.out=grid_size)
	likelihood<-dbinom(x=sum(flips),size=n_flip,prob=p_grid)
	MLE<-p_grid[which.max(likelihood)]
	MLE<-round(MLE,digits=2)

	prior<-dbeta(x=p_grid,shape1=a,shape2=b)

	posterior<-likelihood*prior
	posterior<-posterior/sum(posterior)
	MAP<-p_grid[which.max(posterior)]
	MAP<-round(MAP,digits=2)

	plot(p_grid
	     ,prior
	     ,type="l"
	     ,xlab="probability of H (p)"
	     ,ylab=""
	     ,yaxt="n"
	     ,cex.lab=S
	     ,cex.axis=S
	     ,font.lab=2
	     ,font.axis=2)
	mtext(side=2,text="prior",cex=S,line=L,font=2)
	
	plot(p_grid
	     ,likelihood
	     ,type="l"
	     ,xlab="probability of H (p)"
	     ,ylab=""
	     ,yaxt="n"
	     ,cex.lab=S
	     ,cex.axis=S
	     ,font.lab=2
	     ,font.axis=2)
	mtext(side=2,text="likelihood",cex=S,line=L,font=2)
	mtext(side=3,text=paste("Flips:",n_flip),cex=S,line=L,font=2)
	legend("topright",legend=paste("MLE:",MLE),bty="n",cex=S,text.font=2)

	plot(p_grid
	     ,posterior
	     ,type="l"
	     ,xlab="probability of H (p)"
	     ,ylab=""
	     ,yaxt="n"
	     ,cex.lab=S
	     ,cex.axis=S
	     ,font.lab=2
	     ,font.axis=2)
	mtext(side=2,text="posterior",cex=S,line=L,font=2)
	legend("topright",legend=paste("MAP:",MAP),bty="n",cex=S,text.font=2)
}



plot_coin_flip_CI<-function(n_flip=10,grid_size=1e5,CI_bound=0.89,p=0.7,S=2,L=1.5){
	set.seed(42)
	flips<-rbinom(n=n_flip,size=1,prob=p)

	p_grid<-seq(from=0,to=1,length.out=grid_size)
	likelihood<-dbinom(x=sum(flips),size=n_flip,prob=p_grid)
	likelihood<-likelihood/sum(likelihood)

	MLE<-p_grid[which.max(likelihood)]
	MLE<-round(MLE,digits=2)

	c_likelihood<-cumsum(likelihood)/sum(likelihood)
	
	upper_CI_ind<-(which(c_likelihood>=(1-CI_bound)/2))[1]
	lower_CI_ind<-(which(c_likelihood>=1-(1-CI_bound)/2))[1]-1

	lower_CI<-p_grid[lower_CI_ind]
	upper_CI<-p_grid[upper_CI_ind]

	
	plot(p_grid
	     ,likelihood
	     ,type="l"
	     ,xlab="probability of H (p)"
	     ,ylab=""
	     ,yaxt="n"
	     ,cex.lab=S
	     ,font.lab=2
	     ,font.axis=2)

	#shade AUC corresponding to CI
	upper_lim<-likelihood[lower_CI_ind:upper_CI_ind]
	lower_lim<-rep(0,length(upper_lim))
	x_cord<-p_grid[lower_CI_ind:upper_CI_ind]

	polygon(c(x_cord,rev(x_cord)),c(upper_lim,rev(lower_lim)),border=F,col="red")

	mtext(side=2,text="posterior",cex=S,line=L,font=2)
	mtext(side=3,text=paste("Flips:",n_flip," H:",sum(flips)," T:",n_flip-sum(flips)),cex=S,line=L,font=2)
	legend("topright",legend=paste("MAP:",MLE),bty="n",cex=S,text.font=2)


}

