plot.funpca <-function(x,...){
    f = x$f
    y = x$y
    fi = x$fi
    di = x$di
    fd1 = x$fd1
    f.ucb = x$f.ucb
    f.lcb = x$f.lcb
    fd1.ucb = x$fd1.ucb
    fd1.lcb = x$fd1.lcb

    index = seq(0,1,length.out=length(f))
    nindex = length(index) 

    par(mfrow=c(2,2),mar=c(2,2,2,2))
    # fi and f
    plot(index,f,type="l",col=1,ylim=range(y),xlab="",ylab="",main="Functional Data"); 
    for(i in 1:ncol(fi)) lines(index,fi[,i],col=8)
    lines(index,f,lwd=4)

    # di
    plot(index,di[,1],type="l",col=0,ylim=range(di),xlab="",ylab="",main="Subject Specific Deviations"); 
    for(i in 1:ncol(di)) lines(index,di[,i],col=8)
    lines(index,rep(0,nindex),col=1,lwd=4)

    # f and Confidence bands
    plot(index,f,type="l",col=0,lwd=4,ylim=range(c(f.ucb,f.lcb)),xlab="",ylab="",main="f and 95% CB"); 
    polygon(c(index, rev(index)), c(f.lcb, rev(f.ucb)),col="gray", lty = 0)
    lines(index,f,lwd=4)

    # fd1 and Confidence bands
    plot(index,fd1,type="l",col=0,lwd=4,ylim=range(c(fd1.ucb,fd1.lcb)),xlab="",ylab="",main="1st Der of f and 95% CB"); 
    polygon(c(index, rev(index)), c(fd1.lcb, rev(fd1.ucb)),col="gray", lty = 0)
    lines(index,fd1,lwd=4)
}



