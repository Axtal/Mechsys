# ttriax.R

doplot <- function(iele=0)
{
	c  <- read.table(paste("ttriax_ele",iele,".cal",sep=""),header=T)
	d  <- read.table("~/mechsys_cvs/data/FCH/FCH.TTP.01.dat",header=T)
	op <- par(mar=c(3.5,3.5,0.5,0.5),mfcol=c(2,2))

	plot  (c$Ed*100,-c$q/c$p, type='l',col='red',lwd=2,xlab=expression(epsilon[d]~'[%]'), ylab='-q/p',mgp=c(1.5,0.5,0)); grid()
	points(d$Ed    , d$q/d$p, type='b')

	plot  (c$Ed*100, c$Ev*100, type='l',col='red',lwd=2,xlab=expression(epsilon[d]~'[%]'), ylab=expression(epsilon[v]~'[%]'),mgp=c(1.5,0.5,0)); grid()
	points(d$Ed    ,-d$Ev    , type='b')

	plot  (-c$Ev*100,-c$q/c$p, type='l',col='red',lwd=2,xlab=expression(-epsilon[v]~'[%]'),ylab='-q/p',mgp=c(1.5,0.5,0)); grid()
	points( d$Ev    , d$q/d$p, type='b')

	plot  (-c$p, c$Ev*100, type='l',col='red',lwd=2,ylab=expression(epsilon[v]~'[%]'), xlab='-p', mgp=c(1.5,0.5,0)); grid()
	points( d$p,-d$Ev    , type='b')

	par(op)
}
