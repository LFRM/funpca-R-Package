funpca.default <-
function(DATA, k = NULL, correlation = NULL)
{
est <- funpcaEst(DATA,k,correlation)
#est$fi<-est$fi
est$f<-est$f
est$y<-est$y
#est$di<-est$di
#est$fd1<-est$fd1
#est$fd2<-est$fd2
#est$residuals<-est$error
#est$lmm<-est$est
#est$cb.upp<-est$cb.upp
#est$cb.low<-est$cb.low

##editado 8 feb 2021
#est$ef
#est$alpha
    
est$call <- match.call()
class(est) <- "funpca"
est
}
