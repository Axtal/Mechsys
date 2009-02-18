# ttriax.R

doplot <- function(iele=0)
{
	d <- read.table(paste("ttriax_ele",iele,".cal",sep=""),header=T)
	op <- par(mar=c(3.5,3.5,0.5,0.5),mfcol=c(2,2))
	plot(d$Ed*100,-d$q/d$p, type='l',col='red',lwd=2,xlab=expression(epsilon[d]~'[%]'), ylab='-q/p',mgp=c(1.5,0.5,0)); grid()
	plot(d$Ed*100,d$Ev*100, type='l',col='red',lwd=2,xlab=expression(epsilon[d]~'[%]'), ylab=expression(epsilon[v]~'[%]'),mgp=c(1.5,0.5,0)); grid()
	plot(-d$Ev*100,-d$q/d$p,type='l',col='red',lwd=2,xlab=expression(-epsilon[v]~'[%]'),ylab='-q/p',mgp=c(1.5,0.5,0)); grid()
	plot(-d$p,d$Ev*100,     type='l',col='red',lwd=2,ylab=expression(epsilon[v]~'[%]'), xlab='-p', mgp=c(1.5,0.5,0)); grid()
	par(op)
}
