drbasis <-
function(nn,qq,deriv){
    
    old <- options() 
    on.exit(options(old)) 

    options(expressions=10000)
    n <- nn
    kset <- seq(1,n);
    tset <- seq(0,1,length.out=n);
    q<-qq
    q.global <- qq
    derivv <- deriv

    E<-exp(1)
    I<-1.i
    Pi<-pi
    Power<-function(x,y){
        if(length(x)<2){
            if(x==E){
                ee<-as.brob(E);
                yRe<-Re(y);
                yIm<-complex(real=0,imaginary=Im(y));
                ans<-(ee^yRe)*(ee^yIm)}else{ans<-x^y}
        }else{
            ans<-x^y
        }
        ans
    }
    Sin<-function(x) sin(x)
    Cos<-function(x) cos(x)
    Sqrt<-function(x) sqrt(x)
    Complex<-function(x,y) complex(real=x,imaginary=y)
    Conjugate<-function(x) Conj(x)
    
    consistent.Mo<-function(M,Mo){
        norm<-function(x){(sum(x^2))^(1/2)}
        sign.change<-function(x) {if(x[1]<=x[2]){-1}else{1}}
        D<-cbind(apply(M+Mo,2,norm),apply(M-Mo,2,norm))
        d<-apply(D,1,sign.change)
        cMo<-Mo%*%diag(d)
        cMo
    }

    if(q==1){
        #polynomials
        nullspace<-rep(1,n)
        nullspace.d1<-rep(0,n)
        nullspace.d2<-rep(0,n)
        
        #ODE solution
        k<-(q.global+1):n
        phi<-function(t){Sqrt(2)*Cos((-1 + k)*Pi*t)}
        phi.d1<-function(t){-(Sqrt(2)*(-1 + k)*Pi*Sin((-1 + k)*Pi*t))}
        phi.d2<-function(t){-(Sqrt(2)*Power(as.complex(-1) + k,2)*Power(Pi,2)*Cos((-1 + k)*Pi*t))}
        
        ev<-c(rep(0,q.global),((((q.global+1):n-(q.global+1)/2)*pi)^(2*q.global)))
        ev.d1<-c(rep(0,q.global),2*Power(Pi,2*q.global)*q.global*Power((-1 - q.global)/2. +
                                                                           (q.global+1):n,-1 + 2*q.global))
        ev.d2<-c(rep(0,q.global),2*Power(Pi,2*q.global)*q.global*(-1 + 2*q.global)*
                     Power((-1 - q.global)/2. + (q.global+1):n,-2 + 2*q.global))

        #output
        if(identical(derivv,c(0,0,0))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-matrix(NA,n,n)
            ev.d1<-rep(NA,n)
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(0,0,1))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-matrix(NA,n,n)
            ev.d1<-rep(NA,n)
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        if(identical(derivv,c(0,1,0))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(0,1,1))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        if(identical(derivv,c(1,0,0))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-matrix(NA,n,n)
            ev.d1<-rep(NA,n)
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(1,0,1))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-matrix(NA,n,n)
            ev.d1<rep(NA,n)
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        if(identical(derivv,c(1,1,0))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(1,1,1))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        
    }
    if(q==2){
        #polynomials
        nullspace<-cbind(rep(1,n),-Sqrt(3) + 2*Sqrt(3)*tset)
        nullspace.d1<-cbind(rep(0,n),2*Sqrt(3))
        nullspace.d2<-cbind(rep(0,n),rep(0,n))
        
        #ODE solution
        k<-(q.global+1):n
        phi<-function(t){
            Re(as.complex(
            Power(as.complex(-1),1 + k)/Power(E,(-1.5 + k)*Pi*(1 - t)) + 
            Power(E,-((-1.5 + k)*Pi*t)) + Sqrt(2)*Cos(Pi/4. + (-1.5 + k)*Pi*t)
            ))
        }
        phi.d1<-function(t){
            Re(as.complex(
            (Power(as.complex(-1),1 + k)*(-1.5 + k)*Pi)/Power(E,(-1.5 + k)*Pi*(1 - t)) - 
            ((-1.5 + k)*Pi)/Power(E,(-1.5 + k)*Pi*t) - 
            Sqrt(2)*(-1.5 + k)*Pi*Sin(Pi/4. + (-1.5 + k)*Pi*t)
            ))
        }
        phi.d2<-function(t){
            Re(as.complex(
            (Power(as.complex(-1),1 + k)*Power(as.complex(-1.5) + k,2)*Power(Pi,2))/
            Power(E,(-1.5 + k)*Pi*(1 - t)) + 
            (Power(as.complex(-1.5) + k,2)*Power(Pi,2))/Power(E,(-1.5 + k)*Pi*t) - 
            Sqrt(2)*Power(as.complex(-1.5) + k,2)*Power(Pi,2)*Cos(Pi/4. + (-1.5 + k)*Pi*t) 
            ))
        }
        
        ev<-c(rep(0,q.global),((((q.global+1):n-(q.global+1)/2)*pi)^(2*q.global)))
        ev.d1<-c(rep(0,q.global),2*Power(Pi,2*q.global)*q.global*Power((-1 - q.global)/2. +
                                                                           (q.global+1):n,-1 + 2*q.global))
        ev.d2<-c(rep(0,q.global),2*Power(Pi,2*q.global)*q.global*(-1 + 2*q.global)*
                     Power((-1 - q.global)/2. + (q.global+1):n,-2 + 2*q.global))
        

        #output
        if(identical(derivv,c(0,0,0))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-matrix(NA,n,n)
            ev.d1<-rep(NA,n)
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(0,0,1))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-matrix(NA,n,n)
            ev.d1<-rep(NA,n)
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        if(identical(derivv,c(0,1,0))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(0,1,1))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        if(identical(derivv,c(1,0,0))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-matrix(NA,n,n)
            ev.d1<-rep(NA,n)
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(1,0,1))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-matrix(NA,n,n)
            ev.d1<rep(NA,n)
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        if(identical(derivv,c(1,1,0))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(1,1,1))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        
    }
    if(q==3){
        #polynomials
        nullspace<-cbind(rep(1,n),-Sqrt(3) + 2*Sqrt(3)*tset,-Sqrt(5) + 6*Sqrt(5)*tset - 6*Sqrt(5)*tset^2)
        nullspace.d1<-cbind(rep(0,n),2*Sqrt(3),6*Sqrt(5) - 12*Sqrt(5)*tset)
        nullspace.d2<-cbind(rep(0,n),rep(0,n),rep(-12*Sqrt(5),n))

        #ODE solution
        k<-(q.global+1):n
        phi<-function(t){
            Re(as.complex(
            (Sqrt(1.5) - Complex(0,1)/Sqrt(2))*
            (Power(as.complex(-1),1 + k)/
            Power(E,Power(as.complex(-1),1/6)*(-2 + k)*Pi*(1 - t)) + 
            Power(E,-(Power(as.complex(-1),1/6)*(-2 + k)*Pi*t))) + 
            (Sqrt(1.5) + Complex(0,1)/Sqrt(2))*
            (Power(as.complex(-1),1 + k)/
            Power(E,(Complex(0,-0.5) + Sqrt(3)/2.)*(-2 + k)*Pi*(1 - t)) + 
            Power(E,-((Complex(0,-0.5) + Sqrt(3)/2.)*(-2 + k)*Pi*t))) - 
            Sqrt(2)*Sin((-2 + k)*Pi*t)
            ))
        }

        phi.d1<-function(t){
            Re(as.complex(
            (Sqrt(1.5) - Complex(0,1)/Sqrt(2))*
            ((Power(as.complex(-1),7/6 + k)*(-2 + k)*Pi)/
            Power(E,Power(as.complex(-1),1/6)*(-2 + k)*Pi*(1 - t)) - 
            (Power(as.complex(-1),1/6)*(-2 + k)*Pi)/
            Power(E,Power(as.complex(-1),1/6)*(-2 + k)*Pi*t)) + 
            (Sqrt(1.5) + Complex(0,1)/Sqrt(2))*
            ((Power(as.complex(-1),1 + k)*(Complex(0,-0.5) + Sqrt(3)/2.)*(-2 + k)*Pi)/
            Power(E,(Complex(0,-0.5) + Sqrt(3)/2.)*(-2 + k)*Pi*(1 - t)) - 
            ((Complex(0,-0.5) + Sqrt(3)/2.)*(-2 + k)*Pi)/
            Power(E,(Complex(0,-0.5) + Sqrt(3)/2.)*(-2 + k)*Pi*t)) - 
            Sqrt(2)*(-2 + k)*Pi*Cos((-2 + k)*Pi*t)
            ))
        }

        phi.d2<-function(t){
            Re(as.complex(
            (Sqrt(1.5) - Complex(0,1)/Sqrt(2))*
            ((Power(as.complex(-1),4/3 + k)*Power(-2 + k,2)*Power(Pi,2))/
            Power(E,Power(as.complex(-1),1/6)*(-2 + k)*Pi*(1 - t)) + 
            (Power(as.complex(-1),1/3)*Power(-2 + k,2)*Power(Pi,2))/
            Power(E,Power(as.complex(-1),1/6)*(-2 + k)*Pi*t)) + 
            (Sqrt(1.5) + Complex(0,1)/Sqrt(2))*
            ((Power(as.complex(-1),1 + k)*Power(Complex(0,-0.5) + Sqrt(3)/2.,2)*Power(-2 + k,2)*Power(Pi,2))/
            Power(E,(Complex(0,-0.5) + Sqrt(3)/2.)*(-2 + k)*Pi*(1 - t)) + 
            (Power(Complex(0,-0.5) + Sqrt(3)/2.,2)*Power(-2 + k,2)*Power(Pi,2))/
            Power(E,(Complex(0,-0.5) + Sqrt(3)/2.)*(-2 + k)*Pi*t)) + 
            Sqrt(2)*Power(-2 + k,2)*Power(Pi,2)*Sin((-2 + k)*Pi*t)   
            ))
        }


        ev<-c(rep(0,q.global),((((q.global+1):n-(q.global+1)/2)*pi)^(2*q.global)))
        ev.d1<-c(rep(0,q.global),2*Power(Pi,2*q.global)*q.global*Power((-1 - q.global)/2. +
                                                                           (q.global+1):n,-1 + 2*q.global))
        ev.d2<-c(rep(0,q.global),2*Power(Pi,2*q.global)*q.global*(-1 + 2*q.global)*
                     Power((-1 - q.global)/2. + (q.global+1):n,-2 + 2*q.global))


        #output
        if(identical(derivv,c(0,0,0))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-matrix(NA,n,n)
            ev.d1<-rep(NA,n)
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(0,0,1))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-matrix(NA,n,n)
            ev.d1<-rep(NA,n)
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        if(identical(derivv,c(0,1,0))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(0,1,1))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        if(identical(derivv,c(1,0,0))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-matrix(NA,n,n)
            ev.d1<-rep(NA,n)
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(1,0,1))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-matrix(NA,n,n)
            ev.d1<rep(NA,n)
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        if(identical(derivv,c(1,1,0))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(1,1,1))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        
    }
    if(q==4){
        #polynomials
        nullspace<-cbind(rep(1,n),
                         -Sqrt(3) + 2*Sqrt(3)*tset,
                         -Sqrt(5) + 6*Sqrt(5)*tset - 6*Sqrt(5)*tset^2,
                         -Sqrt(7) + 12*Sqrt(7)*tset - 30*Sqrt(7)*tset^2 + 20*Sqrt(7)*tset^3)
        nullspace.d1<-cbind(rep(0,n),
                            2*Sqrt(3),
                            6*Sqrt(5) - 12*Sqrt(5)*tset,
                            12*Sqrt(7) - 60*Sqrt(7)*tset + 60*Sqrt(7)*tset^2)
        nullspace.d2<-cbind(rep(0,n),
                            rep(0,n),
                            rep(-12*Sqrt(5),n),
                            -60*Sqrt(7) + 120*Sqrt(7)*tset)

        #ODE solution
        k<-(q.global+1):n
        phi<-function(t){
            Re(as.complex(
            (1 + Sqrt(2))*(Power(as.complex(-1),1 + k)/Power(E,(-2.5 + k)*Pi*(1 - t)) + 
            Power(E,-((-2.5 + k)*Pi*t))) + 
            (1/Sqrt(2) - Complex(0,1)*(1 + 1/Sqrt(2)))*
            (Power(as.complex(-1),1 + k)/Power(E,Power(as.complex(-1),0.25)*(-2.5 + k)*Pi*(1 - t)) + 
            Power(E,-(Power(as.complex(-1),0.25)*(-2.5 + k)*Pi*t))) + 
            (1/Sqrt(2) + Complex(0,1)*(1 + 1/Sqrt(2)))*
            (Power(as.complex(-1),1 + k)/
            Power(E,(Complex(1,-1)*(-2.5 + k)*Pi*(1 - t))/Sqrt(2)) + 
            Power(E,(Complex(-1,1)*(-2.5 + k)*Pi*t)/Sqrt(2))) - 
            Sqrt(2)*Cos(Pi/4. - (-2.5 + k)*Pi*t)
            ))
        }
        phi.d1<-function(t){
            Re(as.complex(
            (1 + Sqrt(2))*(Power(as.complex(-1),1 + k)/Power(E,(-2.5 + k)*Pi*(1 - t)) + 
            Power(E,-((-2.5 + k)*Pi*t))) + 
            (1/Sqrt(2) - Complex(0,1)*(1 + 1/Sqrt(2)))*
            (Power(as.complex(-1),1 + k)/Power(E,Power(as.complex(-1),0.25)*(-2.5 + k)*Pi*(1 - t)) + 
            Power(E,-(Power(as.complex(-1),0.25)*(-2.5 + k)*Pi*t))) + 
            (1/Sqrt(2) + Complex(0,1)*(1 + 1/Sqrt(2)))*
            (Power(as.complex(-1),1 + k)/
            Power(E,(Complex(1,-1)*(-2.5 + k)*Pi*(1 - t))/Sqrt(2)) + 
            Power(E,(Complex(-1,1)*(-2.5 + k)*Pi*t)/Sqrt(2))) - 
            Sqrt(2)*Cos(Pi/4. - (-2.5 + k)*Pi*t)
            ))
        }
        phi.d2<-function(t){
            Re(as.complex(
            (1 + Sqrt(2))*((Power(as.complex(-1),1 + k)*Power(-2.5 + k,2)*Power(Pi,2))/
            Power(E,(-2.5 + k)*Pi*(1 - t)) + 
            (Power(-2.5 + k,2)*Power(Pi,2))/Power(E,(-2.5 + k)*Pi*t)) + 
            (1/Sqrt(2) - Complex(0,1)*(1 + 1/Sqrt(2)))*
            ((Complex(0,1)*Power(as.complex(-1),1 + k)*Power(-2.5 + k,2)*Power(Pi,2))/
            Power(E,Power(as.complex(-1),0.25)*(-2.5 + k)*Pi*(1 - t)) + 
            (Complex(0,1)*Power(-2.5 + k,2)*Power(Pi,2))/
            Power(E,Power(as.complex(-1),0.25)*(-2.5 + k)*Pi*t)) + 
            (1/Sqrt(2) + Complex(0,1)*(1 + 1/Sqrt(2)))*
            ((Complex(0,-1)*Power(as.complex(-1),1 + k)*Power(-2.5 + k,2)*Power(Pi,2))/
            Power(E,(Complex(1,-1)*(-2.5 + k)*Pi*(1 - t))/Sqrt(2)) - 
            (Complex(0,1)*Power(-2.5 + k,2)*Power(Pi,2))/
            Power(E,(Complex(1,-1)*(-2.5 + k)*Pi*t)/Sqrt(2))) + 
            Sqrt(2)*Power(-2.5 + k,2)*Power(Pi,2)*Cos(Pi/4. - (-2.5 + k)*Pi*t)
            ))
        }


        ev<-c(rep(0,q.global),((((q.global+1):n-(q.global+1)/2)*pi)^(2*q.global)))
        ev.d1<-c(rep(0,q.global),2*Power(Pi,2*q.global)*q.global*Power((-1 - q.global)/2. +
                                                                           (q.global+1):n,-1 + 2*q.global))
        ev.d2<-c(rep(0,q.global),2*Power(Pi,2*q.global)*q.global*(-1 + 2*q.global)*
                     Power((-1 - q.global)/2. + (q.global+1):n,-2 + 2*q.global))
        

        #output
        if(identical(derivv,c(0,0,0))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-matrix(NA,n,n)
            ev.d1<-rep(NA,n)
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(0,0,1))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-matrix(NA,n,n)
            ev.d1<-rep(NA,n)
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        if(identical(derivv,c(0,1,0))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(0,1,1))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        if(identical(derivv,c(1,0,0))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-matrix(NA,n,n)
            ev.d1<-rep(NA,n)
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(1,0,1))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-matrix(NA,n,n)
            ev.d1<rep(NA,n)
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        if(identical(derivv,c(1,1,0))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(1,1,1))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        
    }
    if(q==5){
        #polynomials
        nullspace<-cbind(rep(1,n),
                         -Sqrt(3) + 2*Sqrt(3)*tset,
                         -Sqrt(5) + 6*Sqrt(5)*tset - 6*Sqrt(5)*tset^2,
                         -Sqrt(7) + 12*Sqrt(7)*tset - 30*Sqrt(7)*tset^2 + 20*Sqrt(7)*tset^3,
                         3 - 60*tset + 270*tset^2 - 420*tset^3 + 210*tset^4)
        nullspace.d1<-cbind(rep(0,n),
                            2*Sqrt(3),
                            6*Sqrt(5) - 12*Sqrt(5)*tset,
                            12*Sqrt(7) - 60*Sqrt(7)*tset + 60*Sqrt(7)*tset^2,
                            -60 + 540*tset - 1260*tset^2 + 840*tset^3)
        nullspace.d2<-cbind(rep(0,n),
                            rep(0,n),
                            rep(-12*Sqrt(5),n),
                            -60*Sqrt(7) + 120*Sqrt(7)*tset,
                            540 - 2520*tset + 2520*tset^2)

        #ODE solution
        k<-(q.global+1):n
        phi<-function(t){
            Re(as.complex(
            (Sqrt(2)*(1 + Sqrt(5)/2.) - (Complex(0,0.5)*
            (Sqrt(10 - 2*Sqrt(5)) + Sqrt(2*(5 + Sqrt(5)))))/Sqrt(2))*
            (Power(as.complex(-1),1 + k)/Power(E,Power(as.complex(-1),0.1)*(-3 + k)*Pi*(1 - t)) + 
            Power(E,-(Power(as.complex(-1),0.1)*(-3 + k)*Pi*t))) + 
            (-(1/Sqrt(2)) - (Complex(0,0.5)*
            (Sqrt(10 - 2*Sqrt(5)) + Sqrt(2*(5 + Sqrt(5)))))/Sqrt(2))*
            (Power(as.complex(-1),1 + k)/Power(E,Power(as.complex(-1),0.3)*(-3 + k)*Pi*(1 - t)) + 
            Power(E,-(Power(as.complex(-1),0.3)*(-3 + k)*Pi*t))) + 
            (Sqrt(2)*(1 + Sqrt(5)/2.) + 
            (Complex(0,0.5)*(Sqrt(10 - 2*Sqrt(5)) + Sqrt(2*(5 + Sqrt(5)))))/Sqrt(2)
            )*(Power(as.complex(-1),1 + k)/
            Power(E,(Sqrt(0.625 + Sqrt(5)/8.) - Complex(0,0.25)*(-1 + Sqrt(5)))*
            (-3 + k)*Pi*(1 - t)) + 
            Power(E,-((Sqrt(0.625 + Sqrt(5)/8.) - Complex(0,0.25)*(-1 + Sqrt(5)))*
            (-3 + k)*Pi*t))) + (-(1/Sqrt(2)) + 
            (Complex(0,0.5)*(Sqrt(10 - 2*Sqrt(5)) + Sqrt(2*(5 + Sqrt(5)))))/Sqrt(2)
            )*(Power(as.complex(-1),1 + k)/
            Power(E,(Sqrt(0.625 - Sqrt(5)/8.) - Complex(0,0.25)*(1 + Sqrt(5)))*
            (-3 + k)*Pi*(1 - t)) + 
            Power(E,-((Sqrt(0.625 - Sqrt(5)/8.) - Complex(0,0.25)*(1 + Sqrt(5)))*
            (-3 + k)*Pi*t))) - Sqrt(2)*Cos((-3 + k)*Pi*t)
            ))
        }
        
        phi.d1<-function(t){
            Re(as.complex(
            (Sqrt(2)*(1 + Sqrt(5)/2.) - (Complex(0,0.5)*
            (Sqrt(10 - 2*Sqrt(5)) + Sqrt(2*(5 + Sqrt(5)))))/Sqrt(2))*
            ((Power(as.complex(-1),1.1 + k)*(-3 + k)*Pi)/
            Power(E,Power(as.complex(-1),0.1)*(-3 + k)*Pi*(1 - t)) - 
            (Power(as.complex(-1),0.1)*(-3 + k)*Pi)/Power(E,Power(as.complex(-1),0.1)*(-3 + k)*Pi*t)) + 
            (-(1/Sqrt(2)) - (Complex(0,0.5)*
            (Sqrt(10 - 2*Sqrt(5)) + Sqrt(2*(5 + Sqrt(5)))))/Sqrt(2))*
            ((Power(as.complex(-1),1.3 + k)*(-3 + k)*Pi)/
            Power(E,Power(as.complex(-1),0.3)*(-3 + k)*Pi*(1 - t)) - 
            (Power(as.complex(-1),0.3)*(-3 + k)*Pi)/Power(E,Power(as.complex(-1),0.3)*(-3 + k)*Pi*t)) + 
            (Sqrt(2)*(1 + Sqrt(5)/2.) + 
            (Complex(0,0.5)*(Sqrt(10 - 2*Sqrt(5)) + Sqrt(2*(5 + Sqrt(5)))))/Sqrt(2)
            )*((Power(as.complex(-1),1 + k)*(Sqrt(0.625 + Sqrt(5)/8.) - 
            Complex(0,0.25)*(-1 + Sqrt(5)))*(-3 + k)*Pi)/
            Power(E,(Sqrt(0.625 + Sqrt(5)/8.) - Complex(0,0.25)*(-1 + Sqrt(5)))*
            (-3 + k)*Pi*(1 - t)) - 
            ((Sqrt(0.625 + Sqrt(5)/8.) - Complex(0,0.25)*(-1 + Sqrt(5)))*(-3 + k)*
            Pi)/Power(E,(Sqrt(0.625 + Sqrt(5)/8.) - 
            Complex(0,0.25)*(-1 + Sqrt(5)))*(-3 + k)*Pi*t)) + 
            (-(1/Sqrt(2)) + (Complex(0,0.5)*
            (Sqrt(10 - 2*Sqrt(5)) + Sqrt(2*(5 + Sqrt(5)))))/Sqrt(2))*
            ((Power(as.complex(-1),1 + k)*(Sqrt(0.625 - Sqrt(5)/8.) - 
            Complex(0,0.25)*(1 + Sqrt(5)))*(-3 + k)*Pi)/
            Power(E,(Sqrt(0.625 - Sqrt(5)/8.) - Complex(0,0.25)*(1 + Sqrt(5)))*
            (-3 + k)*Pi*(1 - t)) - 
            ((Sqrt(0.625 - Sqrt(5)/8.) - Complex(0,0.25)*(1 + Sqrt(5)))*(-3 + k)*
            Pi)/Power(E,(Sqrt(0.625 - Sqrt(5)/8.) - 
            Complex(0,0.25)*(1 + Sqrt(5)))*(-3 + k)*Pi*t)) + 
            Sqrt(2)*(-3 + k)*Pi*Sin((-3 + k)*Pi*t) 
            ))
        }
        
        phi.d2<-function(t){
            Re(as.complex(
            (Sqrt(2)*(1 + Sqrt(5)/2.) - (Complex(0,0.5)*
            (Sqrt(10 - 2*Sqrt(5)) + Sqrt(2*(5 + Sqrt(5)))))/Sqrt(2))*
            ((Power(as.complex(-1),1.2 + k)*Power(-3 + k,2)*Power(Pi,2))/
            Power(E,Power(as.complex(-1),0.1)*(-3 + k)*Pi*(1 - t)) + 
            (Power(as.complex(-1),0.2)*Power(-3 + k,2)*Power(Pi,2))/
            Power(E,Power(as.complex(-1),0.1)*(-3 + k)*Pi*t)) + 
            (-(1/Sqrt(2)) - (Complex(0,0.5)*
            (Sqrt(10 - 2*Sqrt(5)) + Sqrt(2*(5 + Sqrt(5)))))/Sqrt(2))*
            ((Power(as.complex(-1),1.6 + k)*Power(-3 + k,2)*Power(Pi,2))/
            Power(E,Power(as.complex(-1),0.3)*(-3 + k)*Pi*(1 - t)) + 
            (Power(as.complex(-1),0.6)*Power(-3 + k,2)*Power(Pi,2))/
            Power(E,Power(as.complex(-1),0.3)*(-3 + k)*Pi*t)) + 
            (Sqrt(2)*(1 + Sqrt(5)/2.) + 
            (Complex(0,0.5)*(Sqrt(10 - 2*Sqrt(5)) + Sqrt(2*(5 + Sqrt(5)))))/Sqrt(2)
            )*((Power(as.complex(-1),1 + k)*Power(Sqrt(0.625 + Sqrt(5)/8.) - 
            Complex(0,0.25)*(-1 + Sqrt(5)),2)*Power(-3 + k,2)*Power(Pi,2))/
            Power(E,(Sqrt(0.625 + Sqrt(5)/8.) - Complex(0,0.25)*(-1 + Sqrt(5)))*
            (-3 + k)*Pi*(1 - t)) + 
            (Power(Sqrt(0.625 + Sqrt(5)/8.) - Complex(0,0.25)*(-1 + Sqrt(5)),2)*
            Power(-3 + k,2)*Power(Pi,2))/
            Power(E,(Sqrt(0.625 + Sqrt(5)/8.) - Complex(0,0.25)*(-1 + Sqrt(5)))*
            (-3 + k)*Pi*t)) + (-(1/Sqrt(2)) + 
            (Complex(0,0.5)*(Sqrt(10 - 2*Sqrt(5)) + Sqrt(2*(5 + Sqrt(5)))))/Sqrt(2)
            )*((Power(as.complex(-1),1 + k)*Power(Sqrt(0.625 - Sqrt(5)/8.) - 
            Complex(0,0.25)*(1 + Sqrt(5)),2)*Power(-3 + k,2)*Power(Pi,2))/
            Power(E,(Sqrt(0.625 - Sqrt(5)/8.) - Complex(0,0.25)*(1 + Sqrt(5)))*
            (-3 + k)*Pi*(1 - t)) + 
            (Power(Sqrt(0.625 - Sqrt(5)/8.) - Complex(0,0.25)*(1 + Sqrt(5)),2)*
            Power(-3 + k,2)*Power(Pi,2))/
            Power(E,(Sqrt(0.625 - Sqrt(5)/8.) - Complex(0,0.25)*(1 + Sqrt(5)))*
            (-3 + k)*Pi*t)) + Sqrt(2)*Power(-3 + k,2)*Power(Pi,2)*
            Cos((-3 + k)*Pi*t)
            ))
        }
        
        
        ev<-c(rep(0,q.global),((((q.global+1):n-(q.global+1)/2)*pi)^(2*q.global)))
        ev.d1<-c(rep(0,q.global),2*Power(Pi,2*q.global)*q.global*Power((-1 - q.global)/2. +
                                                                           (q.global+1):n,-1 + 2*q.global))
        ev.d2<-c(rep(0,q.global),2*Power(Pi,2*q.global)*q.global*(-1 + 2*q.global)*
                     Power((-1 - q.global)/2. + (q.global+1):n,-2 + 2*q.global))
        
        
        #output
        if(identical(derivv,c(0,0,0))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-matrix(NA,n,n)
            ev.d1<-rep(NA,n)
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(0,0,1))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-matrix(NA,n,n)
            ev.d1<-rep(NA,n)
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        if(identical(derivv,c(0,1,0))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(0,1,1))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        if(identical(derivv,c(1,0,0))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-matrix(NA,n,n)
            ev.d1<-rep(NA,n)
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(1,0,1))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-matrix(NA,n,n)
            ev.d1<rep(NA,n)
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        if(identical(derivv,c(1,1,0))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(1,1,1))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        
    }
    if(q==6){
        #polynomials
        nullspace<-cbind(rep(1,n),
                         -Sqrt(3) + 2*Sqrt(3)*tset,
                         -Sqrt(5) + 6*Sqrt(5)*tset - 6*Sqrt(5)*tset^2,
                         -Sqrt(7) + 12*Sqrt(7)*tset - 30*Sqrt(7)*tset^2 + 20*Sqrt(7)*tset^3,
                         3 - 60*tset + 270*tset^2 - 420*tset^3 + 210*tset^4,
                         -Sqrt(11) + 30*Sqrt(11)*tset - 210*Sqrt(11)*tset^2 + 560*Sqrt(11)*tset^3 - 
                             630*Sqrt(11)*tset^4 + 252*Sqrt(11)*tset^5
                         )
        nullspace.d1<-cbind(rep(0,n),
                            2*Sqrt(3),
                            6*Sqrt(5) - 12*Sqrt(5)*tset,
                            12*Sqrt(7) - 60*Sqrt(7)*tset + 60*Sqrt(7)*tset^2,
                            -60 + 540*tset - 1260*tset^2 + 840*tset^3,
                            30*Sqrt(11) - 420*Sqrt(11)*tset + 1680*Sqrt(11)*tset^2 - 2520*Sqrt(11)*tset^3 + 
                                1260*Sqrt(11)*tset^4
                            )
        nullspace.d2<-cbind(rep(0,n),
                            rep(0,n),
                            rep(-12*Sqrt(5),n),
                            -60*Sqrt(7) + 120*Sqrt(7)*tset,
                            540 - 2520*tset + 2520*tset^2,
                            -420*Sqrt(11) + 3360*Sqrt(11)*tset - 7560*Sqrt(11)*tset^2 + 5040*Sqrt(11)*tset^3
                            )

        #ODE solution
        k<-(q.global+1):n
        phi<-function(t){
            Re(as.complex(
            (3 + 2*Sqrt(3))*(Power(as.complex(-1),1 + k)/Power(E,(-3.5 + k)*Pi*(1 - t)) + 
            Power(E,-((-3.5 + k)*Pi*t))) + 
            ((1 + Sqrt(3))/2. - Complex(0,0.5)*(5 + 3*Sqrt(3)))*
            (Power(as.complex(-1),1 + k)/
            Power(E,Power(as.complex(-1),1/6)*(-3.5 + k)*Pi*(1 - t)) + 
            Power(E,-(Power(as.complex(-1),1/6)*(-3.5 + k)*Pi*t))) + 
            ((-3 - Sqrt(3))/2. - Complex(0,0.5)*(1 + Sqrt(3)))*
            (Power(as.complex(-1),1 + k)/
            Power(E,Power(as.complex(-1),1/3)*(-3.5 + k)*Pi*(1 - t)) + 
            Power(E,-(Power(as.complex(-1),1/3)*(-3.5 + k)*Pi*t))) + 
            ((-3 - Sqrt(3))/2. + Complex(0,0.5)*(1 + Sqrt(3)))*
            (Power(as.complex(-1),1 + k)/
            Power(E,(0.5 - Complex(0,0.5)*Sqrt(3))*(-3.5 + k)*Pi*(1 - t)) + 
            Power(E,-((0.5 - Complex(0,0.5)*Sqrt(3))*(-3.5 + k)*Pi*t))) + 
            ((1 + Sqrt(3))/2. + Complex(0,0.5)*(5 + 3*Sqrt(3)))*
            (Power(as.complex(-1),1 + k)/
            Power(E,(Complex(0,-0.5) + Sqrt(3)/2.)*(-3.5 + k)*Pi*(1 - t)) + 
            Power(E,-((Complex(0,-0.5) + Sqrt(3)/2.)*(-3.5 + k)*Pi*t))) - 
            Sqrt(2)*Cos(Pi/4. + (-3.5 + k)*Pi*t)
            ))
        }
        
        phi.d1<-function(t){
            Re(as.complex(
            (3 + 2*Sqrt(3))*((Power(as.complex(-1),1 + k)*(-3.5 + k)*Pi)/
            Power(E,(-3.5 + k)*Pi*(1 - t)) - 
            ((-3.5 + k)*Pi)/Power(E,(-3.5 + k)*Pi*t)) + 
            ((1 + Sqrt(3))/2. - Complex(0,0.5)*(5 + 3*Sqrt(3)))*
            ((Power(as.complex(-1),7/6 + k)*(-3.5 + k)*Pi)/
            Power(E,Power(as.complex(-1),1/6)*(-3.5 + k)*Pi*(1 - t)) - 
            (Power(as.complex(-1),1/6)*(-3.5 + k)*Pi)/
            Power(E,Power(as.complex(-1),1/6)*(-3.5 + k)*Pi*t)) + 
            ((-3 - Sqrt(3))/2. - Complex(0,0.5)*(1 + Sqrt(3)))*
            ((Power(as.complex(-1),4/3 + k)*(-3.5 + k)*Pi)/
            Power(E,Power(as.complex(-1),1/3)*(-3.5 + k)*Pi*(1 - t)) - 
            (Power(as.complex(-1),1/3)*(-3.5 + k)*Pi)/
            Power(E,Power(as.complex(-1),1/3)*(-3.5 + k)*Pi*t)) + 
            ((-3 - Sqrt(3))/2. + Complex(0,0.5)*(1 + Sqrt(3)))*
            ((Power(as.complex(-1),1 + k)*(0.5 - Complex(0,0.5)*Sqrt(3))*(-3.5 + k)*Pi)/
            Power(E,(0.5 - Complex(0,0.5)*Sqrt(3))*(-3.5 + k)*Pi*(1 - t)) - 
            ((0.5 - Complex(0,0.5)*Sqrt(3))*(-3.5 + k)*Pi)/
            Power(E,(0.5 - Complex(0,0.5)*Sqrt(3))*(-3.5 + k)*Pi*t)) + 
            ((1 + Sqrt(3))/2. + Complex(0,0.5)*(5 + 3*Sqrt(3)))*
            ((Power(as.complex(-1),1 + k)*(Complex(0,-0.5) + Sqrt(3)/2.)*(-3.5 + k)*Pi)/
            Power(E,(Complex(0,-0.5) + Sqrt(3)/2.)*(-3.5 + k)*Pi*(1 - t)) - 
            ((Complex(0,-0.5) + Sqrt(3)/2.)*(-3.5 + k)*Pi)/
            Power(E,(Complex(0,-0.5) + Sqrt(3)/2.)*(-3.5 + k)*Pi*t)) + 
            Sqrt(2)*(-3.5 + k)*Pi*Sin(Pi/4. + (-3.5 + k)*Pi*t)
            ))
        }
        
        phi.d2<-function(t){
            Re(as.complex(
            (3 + 2*Sqrt(3))*((Power(as.complex(-1),1 + k)*Power(-3.5 + k,2)*Power(Pi,2))/
            Power(E,(-3.5 + k)*Pi*(1 - t)) + 
            (Power(-3.5 + k,2)*Power(Pi,2))/Power(E,(-3.5 + k)*Pi*t)) + 
            ((1 + Sqrt(3))/2. - Complex(0,0.5)*(5 + 3*Sqrt(3)))*
            ((Power(as.complex(-1),4/3 + k)*Power(-3.5 + k,2)*Power(Pi,2))/
            Power(E,Power(as.complex(-1),1/6)*(-3.5 + k)*Pi*(1 - t)) + 
            (Power(as.complex(-1),1/3)*Power(-3.5 + k,2)*Power(Pi,2))/
            Power(E,Power(as.complex(-1),1/6)*(-3.5 + k)*Pi*t)) + 
            ((-3 - Sqrt(3))/2. - Complex(0,0.5)*(1 + Sqrt(3)))*
            ((Power(as.complex(-1),5/3 + k)*Power(-3.5 + k,2)*Power(Pi,2))/
            Power(E,Power(as.complex(-1),1/3)*(-3.5 + k)*Pi*(1 - t)) + 
            (Power(as.complex(-1),2/3)*Power(-3.5 + k,2)*Power(Pi,2))/
            Power(E,Power(as.complex(-1),1/3)*(-3.5 + k)*Pi*t)) + 
            ((-3 - Sqrt(3))/2. + Complex(0,0.5)*(1 + Sqrt(3)))*
            ((Power(as.complex(-1),1 + k)*Power(0.5 - Complex(0,0.5)*Sqrt(3),2)*
            Power(-3.5 + k,2)*Power(Pi,2))/
            Power(E,(0.5 - Complex(0,0.5)*Sqrt(3))*(-3.5 + k)*Pi*(1 - t)) + 
            (Power(0.5 - Complex(0,0.5)*Sqrt(3),2)*Power(-3.5 + k,2)*Power(Pi,2))/
            Power(E,(0.5 - Complex(0,0.5)*Sqrt(3))*(-3.5 + k)*Pi*t)) + 
            ((1 + Sqrt(3))/2. + Complex(0,0.5)*(5 + 3*Sqrt(3)))*
            ((Power(as.complex(-1),1 + k)*Power(Complex(0,-0.5) + Sqrt(3)/2.,2)*
            Power(-3.5 + k,2)*Power(Pi,2))/
            Power(E,(Complex(0,-0.5) + Sqrt(3)/2.)*(-3.5 + k)*Pi*(1 - t)) + 
            (Power(Complex(0,-0.5) + Sqrt(3)/2.,2)*Power(-3.5 + k,2)*Power(Pi,2))/
            Power(E,(Complex(0,-0.5) + Sqrt(3)/2.)*(-3.5 + k)*Pi*t)) + 
            Sqrt(2)*Power(-3.5 + k,2)*Power(Pi,2)*Cos(Pi/4. + (-3.5 + k)*Pi*t)
            ))
        }
        
        ev<-c(rep(0,q.global),((((q.global+1):n-(q.global+1)/2)*pi)^(2*q.global)))
        ev.d1<-c(rep(0,q.global),2*Power(Pi,2*q.global)*q.global*Power((-1 - q.global)/2. +
                                                                           (q.global+1):n,-1 + 2*q.global))
        ev.d2<-c(rep(0,q.global),2*Power(Pi,2*q.global)*q.global*(-1 + 2*q.global)*
                     Power((-1 - q.global)/2. + (q.global+1):n,-2 + 2*q.global))
        
        
        #output
        if(identical(derivv,c(0,0,0))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-matrix(NA,n,n)
            ev.d1<-rep(NA,n)
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(0,0,1))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-matrix(NA,n,n)
            ev.d1<-rep(NA,n)
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        if(identical(derivv,c(0,1,0))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(0,1,1))){
            M<-matrix(NA,n,n)
            Mo<-matrix(NA,n,n)
            ev<-rep(NA,n)
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        if(identical(derivv,c(1,0,0))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-matrix(NA,n,n)
            ev.d1<-rep(NA,n)
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(1,0,1))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-matrix(NA,n,n)
            ev.d1<rep(NA,n)
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        if(identical(derivv,c(1,1,0))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-matrix(NA,n,n)
            ev.d2<-rep(NA,n)
        }
        if(identical(derivv,c(1,1,1))){
            M<-cbind(nullspace,t(apply(as.matrix(tset,nrow=1),1,phi)))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-consistent.Mo(M,Mo)
            ev<-ev
            M.d1<-cbind(nullspace.d1,t(apply(as.matrix(tset,nrow=1),1,phi.d1)))/sqrt(n)/n
            ev.d1<-ev.d1
            M.d2<-cbind(nullspace.d2,t(apply(as.matrix(tset,nrow=1),1,phi.d2)))/sqrt(n)/n^2
            ev.d2<-ev.d2
        }
        
    }
    
    eigenvalues=ev
    eigenvectorsQR=Mo
    eigenvectors=M
    eigenvectors.d1=M.d1;
    eigenvectors.d2=M.d2;
    eigenvalues.d1=ev.d1;
    eigenvalues.d2=ev.d2;
    
    list(
        eigenvectors=eigenvectors,
        eigenvectorsQR=eigenvectorsQR,
        eigenvalues=eigenvalues,
        eigenvalues.d1=eigenvalues.d1,
        eigenvalues.d2=eigenvalues.d2,
        eigenvectors.d1=eigenvectors.d1,
        eigenvectors.d2=eigenvectors.d2
    )
}
