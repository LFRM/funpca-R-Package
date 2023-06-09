spectralnorm <-
function(M){result<-sqrt(max(svd(t(M)%*%M)$d));return(result)}
