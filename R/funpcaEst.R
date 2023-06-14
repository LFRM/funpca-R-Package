funpcaEst <-
function(mat,k){
    
    # check input
    mat <- data.frame(mat)
    num.times <- nrow(mat)
    num.cases <- ncol(mat)
    nindex.all <- 1:num.times;

    if(num.cases < 2){stop("Error: data matrix must contain at least 2 columns")}
    if(is.null(k)){stop("Error: k must be an integer between 1 and the sample size.")}

    # variables
    threshold <- 5/100
    q <- c(2,2)
    q.trend <- q[1];
    q.ef <- q[2];
    iterations <- num.cases*10
    convergence <- NULL
    MoD.all <- NULL
    ti <- NULL
    di <- NULL
    x <- seq(0,1,length=num.times)
    nd <- floor(num.times/8)
    
    derivv <- c(1,1,1)
    f.d1 <- matrix(NA,num.times,num.cases)
    f.d2 <- matrix(NA,num.times,num.cases)
    alpha.cb <- 0.05

   
#==================================Basis===================================================  
# STEP 0 * Computing Initial Values:
#==================================Basis===================================================
    message('Running initialization * Computing DR Basis...')
    
    if(!exists("Basis")) Basis <- drbasis(nn = num.times, qq = 2, deriv = derivv)
    n = num.times
#==================================Initial===================================================  
# STEP 2 * Computing Initial Values:
#==================================Initial===================================================
    message('Running STEP 1 * Light smoothing...')

    # light smoothing
    q.trend1 <- 2;
    basis <- Basis

    # matrices
    f.X <- basis$eigenvectorsQR[,1:q.trend1]; 
    f.Z <- basis$eigenvectorsQR[,(q.trend1 + 1):n]%*%diag(1/sqrt(basis$eigenvalues[(q.trend1 + 1):n]));
    f.D <- diag(basis$eigenvalues); 
    f.N <- basis$eigenvectorsQR;
    
    # big Matrices
    One.Vec.All <- rep(1,num.times*num.cases)
    X.pop <- kronecker(matrix(1,nrow = num.cases),f.X)
    Z.pop <- kronecker(matrix(1,nrow = num.cases),f.Z)
    y <- c(as.matrix(mat))

    #fit
    est <- lme(y~-1 + X.pop,random = list(One.Vec.All = pdIdent(~Z.pop-1)),correlation = NULL);
    sigma2.e <- est$sigma^2
    sigma2.b.pop <- (est$sigma^2)*exp(2*unlist(est$modelStruct))
    lambda <- sigma2.e/sigma2.b.pop;
    
    f <- matrix(rep(est$fitted[1:num.times,2],num.cases),num.times,num.cases)
    di <- mat-f

#==================================Select k===================================================  
# STEP 2 * Selecting k
#==================================Select k===================================================
    message('Running STEP 2 * Selecting k...');

    svd <- svd(t(di),nu = num.cases,nv = num.times);
    u <- svd$u;
    v <- svd$v;
    d <- svd$d;
    y <- t(di)%*%v;
    ef <- v[,1:k];
    
    if(k == 1&(min(ef) < 0)){ef <- -ef}

#==================================Update===================================================  
# STEP 3.1 * Updating Components
#==================================Update===================================================
    message('Running STEP 3 * Updating Components...')

    # variables
    One.Vec.All <- rep(1,num.times*num.cases)
    time <- rep(x,num.cases)
    casesId <- rep(1:num.cases,each=num.times)
    y <- c(as.matrix(mat))

    # matrices
    basis <- Basis
    f.X <- basis$eigenvectorsQR[,1:q.trend]; 
    f.Z <- basis$eigenvectorsQR[,(q.trend + 1):n]%*%diag(1/sqrt(basis$eigenvalues[(q.trend + 1):n]));
    f.D <- diag(basis$eigenvalues); 
    f.N <- basis$eigenvectorsQR;
    
    # big matrices
    X.pop <- kronecker(matrix(1,nrow = num.cases),f.X)
    Z.pop <- kronecker(matrix(1,nrow = num.cases),f.Z)

    #<-return
    for(iter in 1:iterations){
        message(paste('Iteration ',iter))  
  
        # deviations
        f.Z.ef <- ef
        Z.cases.ef <- kronecker(matrix(1,nrow = num.cases),f.Z.ef)
  
        # fit
        est <- lme(y~ -1 + X.pop,random = list(One.Vec.All = pdIdent(~Z.pop - 1),casesId = pdDiag(~Z.cases.ef - 1)),correlation = NULL)
        sigma2.e <- est$sigma^2
        sigma2.b.pop <- (est$sigma^2)*exp(2*unlist(est$modelStruct)[k+1])
        lambda <- sigma2.e/sigma2.b.pop;
        alpha <- t(est$coefficients$random$casesId)
    
        fi <- matrix(est$fitted[,3],num.times,num.cases)
        di <- ef%*%alpha
        f <- fi-di
        si <- NULL
        ei <- mat - f

#==================================Correct di=================================================== 
# STEP 3.2 * Correction di
#==================================Correct di===================================================
  
        # some variables
        ef.c <- NULL
        ef.d1.c <- NULL
        ef.d2.c <- NULL
        thetaf.c <- NULL
    
        # scale bases by alpha
        if(k == 1){alpha <- matrix(alpha,ncol = 1)}else{alpha <- matrix(alpha,ncol = ncol(alpha)); alpha<-t(alpha)}
        for(i in 1:k){
            if(k == 1){di.m <- 0}else{di.m <- ef[,-i]%*%t(alpha[,-i])}
            r <- mat - f - di.m;
            r <- c(as.matrix(r))
            
            X.ef <- NULL;
            Z.ef <- NULL;
            N.ef <- NULL;
            D.ef <- NULL
    
            for(j in 1:num.cases){
                basisd <- Basis
                alphamat <- diag(rep(alpha[j,i],num.times))
                X.ef.c <- alphamat%*%basisd$eigenvectorsQR[,1:q.ef]; 
                X.ef <- dbind(X.ef,X.ef.c)
                Z.ef.c <- alphamat%*%basisd$eigenvectorsQR[,(q.ef+1):nd]%*%diag(1/sqrt(basisd$eigenvalues[(q.ef+1):nd])); 
                Z.ef <- dbind(Z.ef,Z.ef.c)
                D.ef.c <- diag(basisd$eigenvalues[1:nd]);  
                D.ef <- dbind(D.ef,D.ef.c)
                N.ef.c <- alphamat%*%basisd$eigenvectorsQR[,1:nd]; 
                N.ef <- dbind(N.ef,N.ef.c)
            }
    
            One.Vec.All <- rep(1,num.times*num.cases)
            est.r <- lme(r~ - 1 + X.ef, random=list(One.Vec.All = pdIdent(~ Z.ef-1)), correlation = NULL)
    
            # fit theta
            sigma2.e <- est.r$sigma^2
            sigma2.b <- (est.r$sigma^2)*exp(2*unlist(est.r$modelStruct))
            lambda <- sigma2.e/sigma2.b
    
            add1 <- 0
            add2 <- 0 
    
            X.ef <- basisd$eigenvectorsQR[,1:q.ef]; 
            Z.ef <- basisd$eigenvectorsQR[,(q.ef+1):nd]%*%diag(1/sqrt(basisd$eigenvalues[(q.ef+1):nd])); 
            N.ef <- basisd$eigenvectorsQR[,1:nd]; 
            D.ef <- diag(basisd$eigenvalues[1:nd]);  

            for(j in 1:num.cases){
                r.c <- matrix(r,ncol=num.cases)[,j]
                add1 <- add1+(alpha[j,i]^2)*diag(ncol(N.ef))+(lambda/num.cases)*D.ef
                add2 <- add2+(alpha[j,i])*(t(N.ef)%*%r.c)
            }
            ef.c <- cbind(ef.c,N.ef%*%(solve(add1)%*%add2))
            thetaf.c <- cbind(thetaf.c,solve(add1)%*%add2)
        }
  
        # convergence
        convergence.c <- spectralnorm(ef - ef.c)
        convergence <- c(convergence,convergence.c);
  
        if(convergence[length(convergence)] < (threshold)){
            message('Convergence achieved');
            if(iter == 1){
                ef <- ef.c
                ef.d1 <- ef.d1.c
                ef.d2 <- ef.d2.c
                thetaf <- thetaf.c
                if(k>1){
                    ef <- ef.c
                    ef.d1 <- ef.d1.c
                    ef.d2 <- ef.d2.c
                    thetaf <- thetaf.c
                }
            }
            break;
        }else{
            if(iter == iterations){
                stop("No convergence achieved")
            }
        }
    
        if(k == 1){ef.c <- ef.c/sqrt(sum(ef.c^2))
             }else{
                 thetaf.c <- qr(thetaf.c); thetaf.c<-qr.Q(thetaf.c)
                 ef.c <- N.ef%*%thetaf.c
             }
        
        ef <- ef.c
        thetaf <- thetaf.c
    }# go back ->
    
#==================================Get estimator===================================================
# STEP 4.1 * Get estimations
#==================================Get estimator===================================================
    message('Running STEP 4 * Getting Final Estimators...')

    corcoefs<-"Model fitted with white noise remainder"
    est <- lme(y~ - 1 + X.pop,random = list(One.Vec.All = pdIdent(~Z.pop - 1),casesId = pdDiag(~Z.cases.ef - 1)),correlation = NULL);

    sigma2.e <- est$sigma^2
    sigma2.b.pop <- (est$sigma^2)*exp(2*unlist(est$modelStruct)[k+1])
    lambda <- sigma2.e/sigma2.b.pop;
    cr.Zef <- est$coefficients$random$casesId
    alpha <- cr.Zef

    # report
    fi <- matrix(est$fitted[,3],num.times,num.cases)
    di <- ef%*%t(alpha)
    f <- fi-di
    error <- mat-fi
  
    nf <- ncol(X.pop)
    nr <- ncol(Z.pop)
    beta <- est$coefficients$fixed[1:nf]
    mu <- est$coefficients$random$One.Vec.All[1:nr]

    # confidence bands
    samples.f <- t(mvrnorm(n = 10000, f.X%*%beta+f.Z%*%mu, (sigma2.b)*f.Z%*%t(f.Z)+(sigma2.e)*diag(num.times))) 
    f.ucb <- apply(samples.f,1,quantile,1-alpha.cb/2)
    f.lcb <- apply(samples.f,1,quantile,alpha.cb/2)

    # derivative function
    f.X.der1 <- basis$eigenvectors.d1[,1:q.trend]; 
    f.Z.der1 <- basis$eigenvectors.d1[,(q.trend+1):n]%*%diag(1/sqrt(basis$eigenvalues[(q.trend+1):n]))*sqrt(num.times)
               
    f.X.der2 <- basis$eigenvectors.d2[,1:q.trend]; 
    f.Z.der2 <- basis$eigenvectors.d2[,(q.trend+1):n]%*%diag(1/sqrt(basis$eigenvalues[(q.trend+1):n]))*sqrt(num.times)
               
    fd1 <- matrix(f.X.der1%*%beta+f.Z.der1%*%mu,num.times,num.cases)[,1]
    fd2 <- matrix(f.X.der2%*%beta+f.Z.der2%*%mu,num.times,num.cases)[,1]
    
    samples.fd1 <- t(mvrnorm(n = 10000, f.X.der1%*%beta+f.Z.der1%*%mu, (sigma2.b)*f.Z.der1%*%t(f.Z.der1)+(sigma2.e)*diag(num.times))) 
    fd1.ucb <- apply(samples.fd1,1,quantile,1-alpha.cb/2)
    fd1.lcb <- apply(samples.fd1,1,quantile,alpha.cb/2)

    samples.fd2 <- t(mvrnorm(n = 10000, f.X.der2%*%beta+f.Z.der2%*%mu, (sigma2.b)*f.Z.der2%*%t(f.Z.der2)+(sigma2.e)*diag(num.times))) 
    fd2.ucb <- apply(samples.fd2,1,quantile,1-alpha.cb/2)
    fd2.lcb <- apply(samples.fd2,1,quantile,alpha.cb/2)

    # curvature function
    k <- fd2/(1+fd1^2)^(3/2)
    samples.k <- samples.fd2/(1+samples.fd1^2)^(3/2)
    k.ucb <- apply(samples.k,1,quantile,1-alpha.cb/2)
    k.lcb <- apply(samples.k,1,quantile,alpha.cb/2)

    return(
         list(
            est = est,
            y = mat,
            f = f[,1],
            di = di,
            fi = fi,
            error = error,
            alpha = alpha,
            ef = ef,
            thetaf = thetaf,
            convergence = convergence,
            q = q,
            f.ucb = f.ucb,
            f.lcb = f.lcb,
            fd1 = fd1, 
            fd1.ucb = fd1.ucb,
            fd1.lcb = fd1.lcb,
            fd2 = fd2,
            fd2.ucb = fd2.ucb,
            fd2.lcb = fd2.lcb
         )
   )
}

