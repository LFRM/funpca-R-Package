dbind <-
function(M1,M2){
    if(is.null(M1)){p<-M2}else{
        if(is.null(M2)){p<-M1}else{
            r1<-dim(M1)[1];	
            r2<-dim(M2)[1];
            c1<-dim(M1)[2];	
            c2<-dim(M2)[2];	
            r<-r1+r2;
            c<-c1+c2;
            p<-matrix(rep(0,r*c),ncol=c);
            p[(1:r1),(1:c1)]<-M1;
            p[((r1+1):r),((c1+1):c)]<-M2
        }}
    return(p)		
}
