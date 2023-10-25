#############################################
########## this is the code #################
########## for binary undirected  networks ##
#############################################

require(rTensor)
require(MASS)
require(logistf)

EPS = 1e-3

logit<-function(x){
	log(x/(1-x))
}


logistic<-function(x){
	1-1/(1+exp(x))
}


norm_vec <- function(x) sqrt(sum(x^2))
scale_vec <- function(x) x/sqrt(sum(x^2))


mytruncate<-function(x,s){
	xtruncate = rep(0, length(x))
	xtruncate[ which(rank(-abs(x))<= s)] = x[which(rank(-abs(x))<= s)]
	xtruncate
}

fun0<-function(A,w,c,X,Z){
    loglikelihood<-sum((A%*%diag(w)%*%t(A)%o%c)*Z@data)-sum(log(1+exp((A%*%diag(w)%*%t(A)%o%c))))
    loglikelihood
}

fun<-function(A,w,c,BB_truncate,X,Z){
    loglikelihood<-sum(((A%*%diag(w)%*%t(A))%o%c+ttm(BB_truncate,X,m = 3))@data*Z@data)-sum(log(1+exp(((A%*%diag(w)%*%t(A))%o%c+ttm(BB_truncate,X,m = 3))@data))) 
    loglikelihood
}

generateTensor<-function(N,n,X,Theta,BB,symmetric=TRUE){
    l<-length(vec(Theta))
    Z<-new("Tensor",3L,c(n,n,N),data=rbinom(n=l,size=1,prob=logistic(vec(Theta+ttm(BB,X,m = 3)))))
    if(symmetric==TRUE){
        for(i in 1:N){
            slice<-matrix(0,n,n)
            slice[upper.tri(slice)]<-Z[,,i]@data[upper.tri(Z[,,i]@data)]
            slice<-slice+t(slice) 
            diag(slice)<-diag(Z[,,i]@data)      
            Z[,,i]@data<-slice
        }
        }else{
            Z<-new("Tensor",3L,c(n,n,N),data=rbinom(n=l,lambda=exp(vec(Theta+ttm(BB,X,m = 3)))))
        }
    Z
}


#######################################################################
## This function estimates the low rank parameter only ################
#######################################################################

mySymGLM0<-function(Z, X, R=null, sparsity=null, niter=50){
	n<-dim(Z)[1]
	N<-dim(Z)[3]
	p<-dim(X)[2]
    ### initialize ####		
    Abar<-apply(Z@data,c(1,2),mean)
    Abar[Abar<=0]<-EPS
    Abar[Abar>=1]<-1-EPS	
	c<-rep(1,N)
    initial<-svd(logit(Abar))
    w<-initial$d[1:R]*(sign(initial$u[1,1:R])*sign(initial$v[1,1:R]))
    A<-initial$u[,1:R] 

    iter<-1; cond<-FALSE; 
    loglikelihood<-fun0(A,w,c,X,Z)

    while(!cond){
        loglikelihood_old<-loglikelihood
    	A_tmp<-A
    	w_tmp<-w
    	vec_a<-ttm(ttm(Z-logistic((A%*%diag(w)%*%t(A))%o%c),t(A),m = 2),t(c),m=3)[,,1]@data
    	direction<-sweep(vec_a,2,w,"*")
    	alpha<-10^-3 
        while(fun0(A_tmp+direction*alpha,w,c,X,Z)<fun0(A,w,c,X,Z)+0.5*alpha*t(as.vector(direction))%*%as.vector(direction)) {
            alpha<-0.5*alpha
        }
        A_tmp<-A_tmp+direction*alpha
    	w<-apply(A_tmp,2,norm_vec)*apply(A_tmp,2,norm_vec)*w_tmp
    	A<-apply(A_tmp,2,scale_vec)
    	Theta<-(A%*%diag(w)%*%t(A))
    	iter<-iter+1
    	loglikelihood<-fun0(A,w,c,X,Z)
    	cond<-(iter > niter)| (loglikelihood-loglikelihood_old <= 0.1)
    }
    BIC<--loglikelihood+log(N*(n+1)*n/2)*(R*n)+2*(R*n)*log(n*n/2)
    return(list(loglikelihood=loglikelihood,BIC=BIC,w=w,A=A))   
}


##################################################################################
## This estimates both the low rank and the sparse parameters ####################
##################################################################################

mySymGLM<-function(Z, X, R=null, sparsity=null, niter=50){
    n<-dim(Z)[1]
    N<-dim(Z)[3]
    p<-dim(X)[2]

    ### initialize ####     
    Abar<-apply(Z@data,c(1,2),mean)
    Abar[Abar<=0]<-EPS
    Abar[Abar>=1]<-1-EPS    
    c<-rep(1,N)
    initial<-svd(logit(Abar))
    w<-initial$d[1:R]*(sign(initial$u[1,1:R])*sign(initial$v[1,1:R]))
    A<-initial$u[,1:R] 
    BB_truncate<-new("Tensor",3L,c(n,n,p),data=0)

    iter<-1; cond<-FALSE; lvec<-c()
    loglikelihood<-fun(A,w,c,BB_truncate,X,Z)

    while(!cond){
        loglikelihood_old<-loglikelihood
        A_tmp<-A
        w_tmp<-w
        vec_a<-ttm(ttm(Z-logistic((A%*%diag(w)%*%t(A))%o%c+ttm(BB_truncate,X,m = 3)@data),t(A),m = 2),t(c),m=3)[,,1]@data*2
        direction<-sweep(vec_a,2,w,"*")
        alpha<-10^-3  
        while(fun(A_tmp+direction*alpha,w,c,BB_truncate,X,Z)<fun(A,w,c,BB_truncate,X,Z)+0.1*alpha*t(as.vector(direction))%*%as.vector(direction)) {
            alpha<-0.5*alpha
        }
        A_tmp<-A_tmp+direction*alpha
        w<-apply(A_tmp,2,norm_vec)*apply(A_tmp,2,norm_vec)*w_tmp
        A<-apply(A_tmp,2,scale_vec)
        Theta<-(A%*%diag(w)%*%t(A))
        tsnr_b<-ttm(Z-logistic((A%*%diag(w)%*%t(A))%o%c+ttm(BB_truncate, X, m = 3)@data),t(X),m=3)@data*2
        beta<-1  
        while (fun(A,w,c,BB_truncate+beta*tsnr_b,X,Z)<fun(A,w,c,BB_truncate,X,Z)+0.1*beta*t(as.vector(tsnr_b))%*%as.vector(tsnr_b)) {
            beta<-0.5*beta
        }
        BB<-BB_truncate+beta*tsnr_b
        BB_truncate<-new("Tensor",3L,c(n,n,p),data=mytruncate(vec(BB),sparsity))
        iter<-iter+1
        loglikelihood<-fun(A,w,c,BB_truncate,X,Z)
        cond<-(iter > niter) | (loglikelihood-loglikelihood_old <= 0.1)
    }
    BIC<--loglikelihood+log(N*(n+1)*n/2)*(R*n+sparsity/2)
    return(list(loglikelihood=loglikelihood,BIC=BIC,w=w,A=A,BB=BB_truncate))  
}


########################################################
###################### tuning  #########################
####### range for R and s can be adjusted   ############
########################################################

best_SymGLM<-function(Z,X){
	print("rank selection")
    n<-dim(Z)[1]
    p<-dim(X)[2]
    result<-list()
    BIC<-c()
    for(R in 2:5){
        output<-mySymGLM0(Z=Z,X=X,R=R,sparsity=0,niter=50)
        BIC[R]<-output$BIC
        plot(BIC,main="BIC for rank selection")
    }
    R<-which.min(BIC)
    print("sparsity level selection")
    BIC2<-c()
    for(i in c(1:10)){
        s<-c(10^seq(-2,0,0.2))[i]
        output<-mySymGLM(Z=Z,X=X,R=R,sparsity=n*n*p*s,niter=50)
        result[[i]]<-output
        BIC2[i]<-output$BIC
        plot(BIC2,main="BIC for sparsity selection")
    }
    result[[which.min(BIC2)]]
}



####################################################
################# Simulation #######################
####################################################

n<-50L
p<-10L
N<-200L
s<-0.1
R<-3


A<-matrix(rnorm(n*R),nrow=n,ncol=R)
m<-A%*%t(A)
Theta<-as.tensor(m%o%rep(1,N))
X<-matrix(rnorm(N*p),nrow=N,ncol=p)
X<-sweep(X,2,apply(X,2,mean),"-")
X<-sweep(X,2,apply(X,2,sd),"/")
slice_vec<-rep(0,n*(n-1)/2)
slice_vec[sample(c(1:(n*(n-1)/2)),n*n*s/2)]<-2
slice<-matrix(0,n,n)
slice[upper.tri(slice)]<-slice_vec
slice<-slice+t(slice)
BB<-as.tensor(slice%o%rep(1,p))
Z<-generateTensor(N,n,X,Theta,BB)

    
result<-best_SymGLM(Z,X)
#result$w
#result$A
#result$BB
