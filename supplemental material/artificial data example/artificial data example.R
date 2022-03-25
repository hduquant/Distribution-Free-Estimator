library(mvtnorm)
library(MASS)
library(matrixcalc)
library(lavaan)
library(expm)
library(doParallel)

method <-c("gammaN.model")

m=3
p=9
ps <- p*(p+1)/2
ms=m*(m-1)/2
q=2*p+ms+1
df <- ps-q
n=45

### read in data
x <- read.csv("artificial data example.csv")
x <- x[,-1]
Scov <-cov(x)
vscov <- vech(Scov)

# first order derivative
sigdsig <- function(p,m,theta0){
  ps <- p*(p+1)/2
  p_h <- p/m
  p2 <- p*2
  
  theta.important<-theta0[4]
  theta0<-  theta0[-4]
  
  lamb <- matrix(0,p,m)
  for (j in 1:m){
    p_hj <- (j-1)*p_h+1
    p_hj1 <- j*p_h
    lamb[p_hj:p_hj1,j] <- theta0[p_hj:p_hj1]
  }
  
  lamb[9,1]<- theta.important
  
  psi_vec <- theta0[(p+1):p2] 
  phi <- matrix(1,m,m)
  k <- 0
  if (m>1){
    for (i in 1:(m-1)){
      for (j in (i+1):m){
        k <- k+1
        phi[j,i] <- theta0[p2+k] 
        phi[i,j] <- theta0[p2+k]
      }
    }
  }
  Sigma0 <- lamb%*%phi%*%t(lamb) + diag(psi_vec) #model-implied covariance
  
  # The following is to evaluate (dSigma0/dtheta) 
  # DL is the derivative with Lambda
  DL <- matrix(0,ps,(p+1)) 
  lambphi <- lamb%*%phi
  lambphi_t <- t(lambphi)
  p_hj1=0
  for (j in 1:m){
    p_hj <-  p_hj1+1
    p_hj1 <-  p_hj+2
    if (j==1){ p_hj1=4}
    for (k in p_hj:p_hj1){
      
      if (j==1){
        tt <- matrix(0,p,m)
        tt[k,j] <- 1
        ttt <- tt%*%lambphi_t+lambphi%*%t(tt)
        DL[,k] <- vech(ttt)} else {
          tt <- matrix(0,p,m)
          tt[k-1,j] <- 1
          ttt <- tt%*%lambphi_t+lambphi%*%t(tt)
          DL[,k] <- vech(ttt)
        }
      
      if (k==4){
        tt <- matrix(0,p,m)
        tt[9,1] <- 1
        ttt <- tt%*%lambphi_t+lambphi%*%t(tt)
        DL[,k] <- vech(ttt)
      }
      # print(list(j,k,tt))
    }
  }
  #Dps is the derivative with Psi
  Dps <- matrix(0,ps,p)
  for (j in 1:p){
    tt <- matrix(0,p,p)
    tt[j,j] <- 1
    Dps[,j] <- vech(tt)
  }
  
  #Dphi is the derivative with phi
  if (m>1){
    ms <- m*(m-1)/2 
    Dphi <- matrix(0,ps,ms)  
    k <- 0 
    for (i in 1:(m-1)){
      for (j in (i+1):m){
        k <- k+1 
        tt <- matrix(0,m,m) 
        tt[j,i] <- 1 ; tt[i,j] <- 1 
        ttt <- lamb%*%tt%*%t(lamb)
        Dphi[,k] <- vech(ttt)
      }
    } 
    vdsig <- cbind(DL,Dps,Dphi)
  } else {  vdsig <- cbind(DL,Dps)}
  
  out <- list('Sigma0'=Sigma0, 'vdsig'=vdsig)
  return(out)
}

  ###### gamma matrices
## gammaadfu 
mean_xt <- apply(x,2,mean)
x_c <- x-matrix(1,n,1)%*%mean_xt
k=1
sigmaele<-matrix(NA,nrow=ps,ncol=n)
for(i in 1:p){
  for(j in i:p){
    # print(c(i,j))
    sigmaele[k,]<- x_c[,i]*x_c[,j]
    k=k+1
  }
}

kk=1
index<-matrix(NA,nrow=p,ncol=p)
for(i in 1:p){
  for(j in i:p){
    index[i,j]=index[j,i]=kk
    kk=kk+1
  }}

gammaadfu<-matrix(NA,nrow=ps,ncol=ps)
for(i in 1:p){
  for(j in i:p){
    for (k in  1:p){
      for (l in 1:p){
        if ( index[k,l]>=index[i,j] & l>=k){
          #    print(c(i,j,k,l))
          gammaadfu[index[k,l],index[i,j]]<- n*(n-1)/(n-2)/(n-3)*(sum(sigmaele[index[k,l],]*sigmaele[index[i,j],])/n-
                                                                   sum(sigmaele[index[k,l],])*sum(sigmaele[index[i,j],])/n^2)-
            n/(n-2)/(n-3)*(sum(sigmaele[index[k,i],])*sum(sigmaele[index[l,j],])/n^2+
                             sum(sigmaele[index[k,j],])*sum(sigmaele[index[i,l],])/n^2-
                             2/(n-1)*sum(sigmaele[index[k,l],])*sum(sigmaele[index[i,j],])/n^2 )
          gammaadfu[index[i,j],index[k,l]]<-gammaadfu[index[k,l],index[i,j]]
        }
      }
    }
  }
}

## gammaadf 
sigmaijkl=c()
gammaadf=matrix(NA,nrow=ps,ncol=ps)
k=1
for(i in 1:ps){
  for(j in i:ps){
    sigmaijkl[k]<-sum(sigmaele[i,]*sigmaele[j,])/n #n
    gammaadf[i,j]<-sigmaijkl[k]-sum(sigmaele[i,])*sum(sigmaele[j,])/(n)^2
    gammaadf[j,i]<- gammaadf[i,j]
    k=k+1
  }
}

    a<-0
    #  starting values
    #################
    theta0=c(0.82, 0.54, 0.69, 0.46, 0.97, 0.96, 0.93, 0.71, 0.90, 0.45, 0.65, 0.93, 0.60,
             0.48, 0.31, 0.42, 0.41, 0.56, 0.29,0.55, 0.39, 0.24)
    
    sig <- sigdsig(p,m, theta0)  
    vdsig<-sig$vdsig  
    Sigma0<-sig$Sigma0  
    vsig0 <- vech(Sigma0)
    
    diconverg=0
    
    # Stop when the update is smaller than ep 
    ep <- 0.0001 
    
    for(t in 1:300){
      ######## gammaN.model; updated with theta0 and sigma0
      if (method=="gammaN.model"){
        
        gammaNm<-matrix(NA,nrow=ps,ncol=ps)
        for(i in 1:p){
          for(j in i:p){
            for (k in  1:p){
              for (l in 1:p){
                if ( index[k,l]>=index[i,j] & l>=k){
                  #  print(c(i,j,k,l))
                  gammaNm[index[k,l],index[i,j]]<-  Sigma0[i,k]*Sigma0[j,l]+Sigma0[i,l]*Sigma0[j,k]
                  gammaNm[index[i,j],index[k,l]]= gammaNm[index[k,l],index[i,j]]
                }
              }
            }
          }
        }
        
        weight<-try(solve(a* gammaadf+(1-a)*gammaNm))
        if( class( weight)=="try-error"){
          diconverg <- 1   # not converge
          break
        }
      }
      
      stdi <- try(solve(t(vdsig) %*% weight%*%vdsig))
      if( class(  stdi)=="try-error"){
        diconverg <- 1   # not converge
        break
      }
      eresdu <- vscov-vsig0 #s-sigma
      dtheta <- t( eresdu) %*% weight %*%  vdsig%*%  stdi
      theta0 <- theta0 + dtheta
      delta <- max(abs(dtheta))
      
      sig <- sigdsig(p,m, theta0)  # gives a list of 2: sig$Sigma0; sig$vdsig
      vdsig<-sig$vdsig 
      Sigma0<-sig$Sigma0
      vsig0 <- vech(Sigma0)
      
      if(delta<=ep) {
        #      sig <- sigdsig(p,m, theta0)  # gives a list of 2: sig$Sigma0; sig$vdsig
        #### ## gammaN.model; updated with theta0 and sigma0
        if (method=="gammaN.model"){
          gammaNm<-matrix(NA,nrow=ps,ncol=ps)
          for(i in 1:p){
            for(j in i:p){
              for (k in  1:p){
                for (l in 1:p){
                  if ( index[k,l]>=index[i,j] & l>=k){
                    #  print(c(i,j,k,l))
                    gammaNm[index[k,l],index[i,j]]<-  Sigma0[i,k]*Sigma0[j,l]+Sigma0[i,l]*Sigma0[j,k]
                    gammaNm[index[i,j],index[k,l]]= gammaNm[index[k,l],index[i,j]]
                  }
                }
              }
            }
          }
          
          weight<-try(solve(a* gammaadf+(1-a)*gammaNm))
          if( class( weight)=="try-error"){
            diconverg <- 1   # not converge
            break
          }
        }#get final weight (gamma.N.modelbased)
        #  vdsig<-sig$vdsig 
        # Sigma0<-sig$Sigma0
        # vsig0 <- vech(Sigma0)#first order
        break};
    }
    
if (t<300 & diconverg==0) {
        # Test start here;
        dtd<-try(solve(t(vdsig) %*% weight %*%vdsig))
        r.Var<-dtd%*%t(vdsig)%*%weight%*%gammaadf%*%weight%*%vdsig%*%dtd
        r.SE <- sqrt(diag(r.Var)/(n-1))
        
        r.Varu<-dtd%*%t(vdsig)%*%weight%*%gammaadfu%*%weight%*%vdsig%*%dtd
        r.SEu <- sqrt(diag(r.Varu)/(n-1))
        
        result<-cbind( c(theta0),r.SE, c(theta0)/r.SE,r.SEu, c(theta0)/r.SEu)
       # write.csv(result,"/Users/handu/Desktop/sem projects/test project/simulation/real data/real4.csv")
       
        dwe<-t(vdsig)%*%weight
        U<- weight-t(dwe)%*%dtd%*%dwe
        ugamma <-  U %*% gammaadf
        ugammau <-  U %*% gammaadfu
      
        
        Fstats<-(t(vscov-vsig0)%*% weight %*%(vscov-vsig0))
        
        Tstats<-(t(vscov-vsig0)%*% weight %*%(vscov-vsig0))*(n-1)
        
        c2u<-sum(diag(ugammau ))/qr(ugammau)$rank
        rTstats2u <- Tstats/c2u
        
        c2<-sum(diag(ugamma ))/qr(ugamma)$rank
        rTstats2 <- Tstats/c2
        
      ### T
       glspvalue<-1-pchisq(Tstats, df)
       glsrTstats<-Tstats
       print(c(glsrTstats,glspvalue)) #6.222278e+01 1.811610e-05
      ### T_SB=T_cor1
       glspvalue2<-1-pchisq(rTstats2, qr(ugamma)$rank)
       glsrTstats2<-rTstats2
       print(c(glsrTstats2,glspvalue2)) # 36.26355625  0.03876749
      ### T_SB^U=T_cor1^U
       glspvalue2u<-1-pchisq(rTstats2u, qr(ugammau)$rank)
       glsrTstats2u<-rTstats2u
       print(c(glsrTstats2u,glspvalue2u)) #33.52629717  0.07232489

      ###### T_MVA2^U
       si=apply(x, 1, function(z) vech((z-colMeans(x))%*%t(z-colMeans(x))))
       U5=sqrtm(U)
       wi=apply(si, 2, function(si)  U5%*%si)
       wi=Re(wi)
       wbar=rowMeans(wi)
       yi=wi-wbar
       Dmatrix<-function(i){
         t(yi[,i])%*%yi[,i]
       } 
       D<-lapply(1:n, Dmatrix)
       D<-unlist(D)
       a2c.u<-1/(n*(n-1)*(n-2)*(n-3))*((n-2)*(n-1)^3*sum(diag(U%*%gammaadfu%*%U%*%gammaadfu))-
                                         n*(n-1)*sum(D^2)+sum(diag((n-1)*U%*%gammaadfu))^2 )
       c13 <-a2c.u/sum(diag(ugammau )) 
       df2.2.u<-(sum(diag(ugammau )))^2/ a2c.u
       rTstats13 <- Tstats/c13
       pvalue13<-pchisq(rTstats13,df2.2.u,lower.tail =FALSE)
       print(c(rTstats13,pvalue13)) # 30.32705196  0.08121776
}
       