plot.funpca <-function(x,...){
    f = x$f
    y = x$y

    plot(f,ylim = range(y),col=0,xlab="",ylab="f estimation");
    apply(y,2,lines,col=8)
    lines(f,col=2,lwd=4)
}



