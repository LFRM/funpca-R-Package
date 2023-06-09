print.funpca <-
function(x,...)
{
q.f<-(x$parameters)[3]
q.di<-(x$parameters)[4]
num.iter<-x$num.iter
  
cat("Call:\n")
print(x$call)
cat("\n")
cat("\n Functional principal component analyis based on the mixed models representation of penalized splines:\n")
cat("\n")
cat("Parameters used:\n")
cat(paste("penalization order in overall trend:",q.f,"\n"))
cat(paste("penalization order in subj spec deviations:",q.di,"\n"))

cat("\n Summary statistics of components:\n")
print(summary(x$fi))

}
