x <- seq(-4, 4, length=1000)
hx <- dnorm(x)
rb = 1
lb = qnorm(1.5-pnorm(rb))
plot(x, hx, type="l", lty=1, xlab="",
     ylab="Density", main="Standard Normal Density",axes=FALSE)
#abline(v=0,lty=2)
#abline(v=rb,lty=2,col=2)
#abline(v=lb,lty=2,col=3)
i = c(which(x>=0&x<=lb))
polygon(c(0,x[i],lb),c(0,hx[i],0),col='gray50')
i = c(which(x>=rb))
polygon(c(rb,x[i],4),c(0,hx[i],0),col='gray50')

axis(1,at = c(-4,0,lb,rb,4),labels = c("-4",'0',expression(tilde(z)),'z',"4"),pos=0)
axis(2)

