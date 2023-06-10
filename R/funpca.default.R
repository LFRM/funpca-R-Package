funpca.default <-
function(mat, k)
{

est <- funpcaEst(mat,k)
est$fi <- est$fi
est$f <- est$f
est$fd1 <- est$fd1
est$y <- est$y
est$di <- est$di
est$residuals <-est$error
est$lmm <- est$est
est$ef <- est$ef
est$alpha <- est$alpha
est$f.ucb <- est$f.ucb 
est$f.lcb <- est$f.lcb 
est$fd1.ucb <- est$fd1.ucb 
est$fd1.lcb <- est$fd1.lcb 

est$call <- match.call()
class(est) <- "funpca"
est
}
