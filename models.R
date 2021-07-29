library(bayestestR)

model1="model {
         
		 pB ~ dbeta(1,1)
         
		 n ~ dbinom(pB,Ntot)
		 	  
         
}"

model2="model {
         	
		for(j in 1:4){
			  n[j] ~ dbinom(pB[j],Ntot[j])
		 }	  
         
		 for(i in 1:4){
			 pB[i] ~ dbeta(alpha,beta)
		 }
         
		alpha ~ dexp(0.001) 
		beta ~ dexp(0.001)
		 
         
}"

set.seed(27)
Ntot=c(10,100,100,100)
#Ntot=c(10000,10000,10000,10000)
#p<-c(0.7,0.4,0.5,0.8)
p<-c(0.7,0.7,0.7,0.7)
n<-rep(NA,4)
for(i in 1:4)	n[i]<-rbinom(n=1,size=Ntot[i],prob=p[i])

data1=list(n=sum(n),Ntot=sum(Ntot));
data2=list(n=n,Ntot=Ntot);
varnames1=c("pB")
varnames2=c("pB","alpha","beta")
burn_in=1000;
steps=100000;
thin=100;

library(rjags)
fileConn=file("model1.tmp")
writeLines(model1,fileConn);
close(fileConn)

fileConn=file("model2.tmp")
writeLines(model2,fileConn);
close(fileConn)

m1=jags.model(file="model1.tmp",data=data1);
m2=jags.model(file="model2.tmp",data=data2);
update(m1,burn_in)
update(m2,burn_in)
draw1=jags.samples(m1,steps,thin=thin,variable.names=varnames1)
draw2=jags.samples(m2,steps,thin=thin,variable.names=varnames2)

plot_post<-function(draw=draw2,S=2){
	par(mar=c(5,5,5,3))
	plot(n/Ntot
	     ,ylim=c(0,1)
	     ,pch=19
	     ,col=rgb(0,0,1,0.5)
	     ,cex=S
	     ,xaxt="n"
	     ,xlab=""
	     ,ylab="p"
	     ,cex.lab=S
	     ,cex.axis=S
	     ,font.axis=2
	     ,font.lab=2
	     ,main="n samples")	
	mtext(text=Ntot,at=1:4,side=3,cex=S,font=2)

	MAP<-apply(draw2$pB,1,map_estimate)
	ci_hdi<-apply(draw$pB,1,ci,method="HDI")
	l_ci<-sapply(ci_hdi,"[[",2)
	u_ci<-sapply(ci_hdi,"[[",3)



	points(MAP
	       ,col=rgb(1,0,0,0.5)
	       ,pch=17
	       ,cex=S)

	segments(x0=1:4,y0=l_ci,y1=u_ci,lwd=2)

	leg_text<-c("empirical","MAP")
	cols<-c(rgb(0,0,1,0.5),rgb(1,0,0,0.5))
	legend("bottomright",legend=leg_text,bty="n",fill=cols,cex=S,text.font=2)
}

plot_prior_est<-function(draw=draw2,S=2){

	Alpha<-map_estimate(draw$alpha[1,,1])
	Beta<-map_estimate(draw$beta[1,,1])
	Alpha_ci<-ci(draw$alpha[1,,1],method="HDI")
	Beta_ci<-ci(draw$beta[1,,1],methpds="HDI")

	par(mar=c(5,5,5,3))
	plot(c(Alpha,Beta)
	     ,ylim=c(0,5000)
	     ,pch=19
	     ,col=rgb(0,0,1,0.5)
	     ,cex=S
	     ,xaxt="n"
	     ,xlab=""
	     ,ylab="estimated value"
	     ,cex.lab=S
	     ,cex.axis=S
	     ,font.axis=2
	     ,font.lab=2
	     ,main="")	
	mtext(text=c("alpha","beta"),at=1:2,side=1,cex=S,font=2,line=1.5)
	segments(x0=1:2,y0=c(Alpha_ci[[2]],Beta_ci[[2]]),y1=c(Alpha_ci[[3]],Beta_ci[[3]]),lwd=2)
}

	
plot_prior<-function(draw=draw2,S=2){
	par(mar=c(5,5,5,3))
	Alpha<-map_estimate(draw$alpha[1,,1])
	Beta<-map_estimate(draw$beta[1,,1])
	p_grid<-seq(from=0,to=1,length.out=1e4)
	dens<-dbeta(p_grid,shape1=Alpha,shape2=Beta)
	MAP<-p_grid[which.max(dens)]
	plot(p_grid
	     ,dens
	     ,type="l"
	     ,yaxt="n"
             ,ylab="Prior"
	     ,xlab="p"
	     ,font.axis=2
	     ,font.lab=2
	     ,cex.lab=S
	     ,cex.axis=S)
}

plot_stuff<-function(S=2.5){
	par(mfrow=c(1,3))
	plot_post(S=S)
	plot_prior_est(S=S)
	plot_prior(S=S)
}
