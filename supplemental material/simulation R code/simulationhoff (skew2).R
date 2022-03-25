library(mvtnorm)
library(MASS)
library(matrixcalc)
library(doParallel)
library(CompQuadForm)
library(expm)

RUNVARS <- commandArgs(T)
RUNMETHOD <- as.integer(RUNVARS[1])
id <- as.integer(RUNVARS[2]) 
ncores <- as.integer(RUNVARS[3])
FILEPATH <- RUNVARS[4]

#----------------------------------------------------------------
# Computing transformed skewness and kurtosis
#----------------------------------------------------------------
Mardia.tran1 <- function(x){
  n <- nrow(x)
  p <- ncol(x)
  x.c <- x-matrix(apply(x,2,mean),nrow=n,ncol=p,byrow=T)
  
  S <- cov(x)*(n-1)/n
  S_inv <- solve(S)
  D <- x.c%*%S_inv%*%t(x.c)
  
  b1 <- sum(D^3)
  b2 <- sum(diag(D^2))/n   #Mardia's kurtosis
  
  df <- p*(p+1)*(p+2)
  sk <- b1/(n*df)
  kt <- b2/(p*(p+2))
  
  return(list('skew'=sk,'kurt'=kt))
}

################
################
# ----------------------------------- 
# Hao Wu's method
HW.sim.skew <- function(seed,n,p,m,lambphi_12.2, psi12){
  set.seed(seed)
  # one_p_mat <- matrix(rep(1,p*n), ncol=n)
  sqrt2 <- sqrt(2)
  f <- matrix(rnorm(m*n),ncol=n)
  f_s <- (f*f-1)/sqrt2;   #z_{\xi}, standardized chisquare;f~N(0,1),E(f*f)=1
  e <- matrix(rnorm(p*n),ncol=n)
  e_s <- (e*e-1)/sqrt2; 
  # generating random z1 and z2
  chis1 <- rchisq(n=n,df=1)
  z1 <- sqrt(chis1/1)
  eta<-rchisq(n=n,df=2-1)
  z2<-sqrt((1*chis1+eta)/2)
  x_i <- lambphi_12.2%*%f_s + matrix(rep(psi12,n),ncol=n)*e_s
  z1_mat <- matrix(rep(z1,p/m),ncol=p/m)
  z2_mat <- matrix(rep(z2,p/m),ncol=p/m)
  x_i<- t(x_i)
  x_i[,1:(p/m)] <-  x_i[,1:(p/m)]*  z1_mat  # First 5 variables (factor 1)
  x_i[,(p/m+1):(2*p/m)] <-   x_i[,(p/m+1):(2*p/m)]*z2_mat  # Second 5 variables (factor 2)
  return(x_i)
}
###########


################################################################################
cond <- NULL
cond[[1]] <- expand.grid(n_vec=c(40,60,100,200,300,500,1000),
                         m_vec=c(3),
                         p=15)
cond[[2]] <- expand.grid(n_vec=c(100,200,300,500,1000),
                         m_vec=c(3),
                         p=30)

condition<-do.call(rbind, cond)
index.l<-c(1:nrow(condition))
################################################################################

#alist<-seq(0,1,by=0.02)
methodlist<-c("gammaN.model")
dislist<-c("nskew")
it<-c(1:20)

#grid[(grid$method=="ridgeD"&grid$simdist=="enorm"),]

grid <- expand.grid(it,methodlist,dislist,index.l)
grid<-as.data.frame(grid)
colnames(grid)<-c("it","method","simdist","index")

args <- grid[id, ]

do_one <- function(it,method,simdist,index){
  
  m=condition[index,2]
  p=condition[index,3]
  n=condition[index,1]
  ms=m*(m-1)/2
  q=2*p+ms
  ps=p*(p+1)/2
  p2=p*2
  df=ps-q
  p_m=p/m

################
################

set.seed(which(simdist==dislist)*13^2+n*17)
lamb_uni<-sort(sample(seq(70,95,5),p_m,replace=TRUE))/100
lamb <- matrix(0,p,m)
for (j in 1:m){
  lamb[(p_m*(j-1)+1):(p_m*j),j] <- lamb_uni
}

lamb_t <- t(lamb)
phi <- 0.5*diag(m)+matrix(0.5,m,m) #desired phi
phi2 <- phi #rescaled phi for HW's method
phi2[1,2]=phi2[2,1]=0.5/0.9 #based on simulation, E(z1z2)=0.9004993
phi2[1,3]=phi2[3,1]=0.5/0.798013  #based on simulation, E(z1)=0.798013
phi2[2,3]=phi2[3,2]=0.5/0.8863193 #based on simulation, E(z2)=0.8863193

eigen_phi <- eigen(phi)
phi_val <- eigen_phi$values
phi_vec <- eigen_phi$vectors
phi_h <- phi_vec%*%diag(sqrt(phi_val))%*%t(phi_vec)
lambphi_12 <- lamb%*%phi_h   # to simulate data from skewed distribution; phi^0.5

eigen_phi2 <- eigen(phi2)
phi_val2 <- eigen_phi2$values
phi_vec2 <- eigen_phi2$vectors
phi_h2 <- phi_vec2%*%diag(sqrt(phi_val2))%*%t(phi_vec2)
lambphi_12.2 <- lamb%*%phi_h2   # to simulate data from skewed distribution; phi^0.5

lamblamb <- lamb%*%phi%*%lamb_t
psi <- diag(p)-diag(diag(lamblamb)) # psi is specified to ensure sigma's diagonals are 1
psi_vec <- diag(psi)
psi12 <- sqrt(psi_vec) # to simulate data from skewed distribution

sig <- lamblamb+psi    #BIG SIGMA
eigen_sig <- eigen(sig)
sig_vec <- eigen_sig$vectors
sig_val <- eigen_sig$values
sig12 <-  sig_vec%*%diag(sqrt(sig_val))%*%t(sig_vec)

theta_p <- c(rep(lamb_uni,m),psi_vec,phi[lower.tri(phi, diag = FALSE)])

result.sum=matrix(NA,ncol=83+q*5, nrow=50)

#simulate<-function(rep){
for (rep in 1:50){
############################################# 
###############simulate x
#### Below is the SAME for each method 
  seed=(it-1)*50+rep 

  x<-HW.sim.skew(seed,n,p,m,lambphi_12.2, psi12)

#get sample covariance
  ps <- p*(p+1)/2
  Scov <-cov(x) 
  vscov <- vech(Scov)
  
  si=apply(x, 1, function(z) vech((z-colMeans(x))%*%t(z-colMeans(x))))

#get model-implied covariance and first order derivative
sigdsig <- function(p,m,theta0){
  ps <- p*(p+1)/2
  p_h <- p/m
  p2 <- p*2
  lamb <- matrix(0,p,m)
  for (j in 1:m){
    p_hj <- (j-1)*p_h+1
    p_hj1 <- j*p_h
    lamb[p_hj:p_hj1,j] <- theta0[p_hj:p_hj1]
  }
  
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
  DL <- matrix(0,ps,p) 
  lambphi <- lamb%*%phi
  lambphi_t <- t(lambphi)
  for (j in 1:m){
    p_hj <- (j-1)*p_h+1
    p_hj1 <- j*p_h
    for (k in p_hj:p_hj1){
      tt <- matrix(0,p,m)
      tt[k,j] <- 1
      ttt <- tt%*%lambphi_t+lambphi%*%t(tt)
      DL[,k] <- vech(ttt)
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

#######GLM
  ep <- 0.0001
  ps <- p*(p+1)/2
  q <- length(theta_p)
  df <- ps-q

  ###### gamma matrices
  ## gammaadf 
  mean_xt <- apply(x,2,mean)
  x_c <- x-matrix(1,n,1)%*%mean_xt
  k=1
  sigmaele<-matrix(NA,nrow=ps,ncol=n)
  for(i in 1:p){
    for(j in i:p){
      sigmaele[k,]<- x_c[,i]*x_c[,j]
      k=k+1
    }
  }
  sigmaij=sigmakl=rowSums( sigmaele)/n
  
  sigmaijkl=c()
  gammaadf=matrix(NA,nrow=ps,ncol=ps)
  k=1
  for(i in 1:ps){
    for(j in i:ps){
      sigmaijkl[k]<-sum(sigmaele[i,]*sigmaele[j,])/n
      gammaadf[i,j]<-sigmaijkl[k]-sum(sigmaele[i,])*sum(sigmaele[j,])/n^2
      gammaadf[j,i]<- gammaadf[i,j]
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
  
  gammaadf.u<-matrix(NA,nrow=ps,ncol=ps)
  for(i in 1:p){
    for(j in i:p){
      for (k in  1:p){
        for (l in 1:p){
          if ( index[k,l]>=index[i,j] & l>=k){
            #    print(c(i,j,k,l))
            gammaadf.u[index[k,l],index[i,j]]<- n*(n-1)/(n-2)/(n-3)*(sum(sigmaele[index[k,l],]*sigmaele[index[i,j],])/n-
                                                                       sum(sigmaele[index[k,l],])*sum(sigmaele[index[i,j],])/n^2)-
              n/(n-2)/(n-3)*(sum(sigmaele[index[k,i],])*sum(sigmaele[index[l,j],])/n^2+
                               sum(sigmaele[index[k,j],])*sum(sigmaele[index[i,l],])/n^2-
                               2/(n-1)*sum(sigmaele[index[k,l],])*sum(sigmaele[index[i,j],])/n^2 )
            gammaadf.u[index[i,j],index[k,l]]<-gammaadf.u[index[k,l],index[i,j]]
          }
        }
      }
    }
  }
  
#### Above is the SAME for each method  
############################################# 
  
##### Try different a
   a<-0
   
  ######
   
   #starting value
   
   set.seed(which(simdist==dislist)*13^2+n*17)
   
   theta0=theta_p-runif(q,0.1,0.2)
   if (sum(theta0<=0 | theta0>=1)!=0){
     theta0[which(theta0<=0 | theta0>=1)]=theta_p[which(theta0<=0 | theta0>=1)]}
   
   sig <- sigdsig(p,m, theta0)  # gives a list of 2: sig$Sigma0; sig$vdsig
   vdsig<-sig$vdsig  #first order, sigma-dot
   Sigma0<-sig$Sigma0 # model_implied covariance matrix
   vsig0 <- vech(Sigma0)
   
   diconverg=0
   ###### 
   for(t in 1:200){
     ###### ## gammaN.model; updated with theta0 and sigma0
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
       break};
   }
  
  if(t<200){
    if(diconverg==0){  #converge within 120 iteration 
    # Test start here;
    dtd<-try(solve(t(vdsig) %*% weight %*%vdsig))
    if( class( dtd)=="try-error"){
      break
    }
    
    dwe<-t(vdsig)%*%weight
    U<- weight-t(dwe)%*%dtd%*%dwe
    
    r.Var<-dtd%*%t(vdsig)%*%weight%*%gammaadf%*%weight%*%vdsig%*%dtd
    r.SE <- sqrt(diag(r.Var)/(n-1))
    r.Var.u<-dtd%*%t(vdsig)%*%weight%*%gammaadf.u%*%weight%*%vdsig%*%dtd
    r.SE.u <- sqrt(diag(r.Var.u)/(n-1))
    Var<-dtd
    SE <- sqrt(diag(Var)/(n-1))
   
    ###############################
    
    ugamma <- U %*% gammaadf
    rank<-qr(ugamma)$rank
    
    ###### calculate a2    
    U5=sqrtm(U)
    wi=apply(si, 2, function(si)  U5%*%si)
    wi=Re(wi)
    wbar=rowMeans(wi)
    yi=wi-wbar
    H= yi%*%t(yi)
    
    Dmatrix<-function(i){
      t(yi[,i])%*%yi[,i]
    } 
    D<-lapply(1:n, Dmatrix)
    D<-unlist(D)
    a2c<-1/(n*(n-1)*(n-2)*(n-3))*((n-2)*(n-1)*sum(diag(H%*%H))-
                                    n*(n-1)*sum(D^2)+sum(diag(H))^2 )
    
    ######
    # 7 T ststistics
    Fstats<-(t(vscov-vsig0)%*% weight %*%(vscov-vsig0))
    Tstats<-(t(vscov-vsig0)%*% weight %*%(vscov-vsig0))*(n-1)
    
    Tstats.b<-(n-(2*p+11)/6-2*m/3)* Fstats
    hq=((1+8*q)^0.5-1)/2
    Tstats.s<-(n-1-(p*(2*p^2+3*p-1)-hq*(2*hq^2+3*hq-1))/12/df)* Fstats
    Tstats.s2<-(n-1-(p*(2*p^2+3*p-1)-hq*(2*hq^2+3*hq-1))/12/rank)* Fstats
    q1=q-2*p
    Tstats.e1<-(n-(2.381+0.367*p+0.003*q1))*Fstats
    Tstats.e2<-(n-(2.299+0.365*p+0.038*m))*Fstats
    Tstats.e3<-(n-(2.262+0.365*p+0.052*m-0.002*q1))*Fstats
    
    #noncorrect
    pvalue1<-c()
    pvalue1[1]<-pchisq(Tstats,df,lower.tail =FALSE)
    pvalue1[2]<-pchisq(Tstats.b,df,lower.tail =FALSE)
    pvalue1[3]<-pchisq(Tstats.s,df,lower.tail =FALSE)
    pvalue1[4]<-pchisq(Tstats.s2,df,lower.tail =FALSE)
    pvalue1[5]<-pchisq(Tstats.e1,df,lower.tail =FALSE)
    pvalue1[6]<-pchisq(Tstats.e2,df,lower.tail =FALSE)
    pvalue1[7]<-pchisq(Tstats.e3,df,lower.tail =FALSE)
    
    #recaled correction
    c1 <- sum(diag(ugamma ))/df
    rTstats1 <- Tstats/c1
    rTstats1.b <- Tstats.b/c1
    rTstats1.s <- Tstats.s/c1
    rTstats1.s2 <- Tstats.s2/c1
    rTstats1.e1 <- Tstats.e1/c1
    rTstats1.e2 <- Tstats.e2/c1
    rTstats1.e3 <- Tstats.e3/c1
    pvalue2<-c()
    pvalue2[1]<-pchisq(rTstats1,df,lower.tail =FALSE)
    pvalue2[2]<-pchisq(rTstats1.b,df,lower.tail =FALSE)
    pvalue2[3]<-pchisq(rTstats1.s,df,lower.tail =FALSE)
    pvalue2[4]<-pchisq(rTstats1.s2,df,lower.tail =FALSE)
    pvalue2[5]<-pchisq(rTstats1.e1,df,lower.tail =FALSE)
    pvalue2[6]<-pchisq(rTstats1.e2,df,lower.tail =FALSE)
    pvalue2[7]<-pchisq(rTstats1.e3,df,lower.tail =FALSE)
    
    #rank deficient correction; rank
    c2<-sum(diag(ugamma ))/rank
    rTstats2 <- Tstats/c2
    rTstats2.b <- Tstats.b/c2
    rTstats2.s <- Tstats.s/c2
    rTstats2.s2 <- Tstats.s2/c2
    rTstats2.e1 <- Tstats.e1/c2
    rTstats2.e2 <- Tstats.e2/c2
    rTstats2.e3 <- Tstats.e3/c2
    pvalue3<-c()
    pvalue3[1]<-pchisq(rTstats2,rank,lower.tail =FALSE)
    pvalue3[2]<-pchisq(rTstats2.b,rank,lower.tail =FALSE)
    pvalue3[3]<-pchisq(rTstats2.s,rank,lower.tail =FALSE)
    pvalue3[4]<-pchisq(rTstats2.s2,rank,lower.tail =FALSE)
    pvalue3[5]<-pchisq(rTstats2.e1,rank,lower.tail =FALSE)
    pvalue3[6]<-pchisq(rTstats2.e2,rank,lower.tail =FALSE)
    pvalue3[7]<-pchisq(rTstats2.e3,rank,lower.tail =FALSE)
    
    #mean-variance correction
    ugamma2 <- ugamma %*% ugamma
    c3 <-sum(diag( ugamma2))/sum(diag(ugamma )) 
    df2<-(sum(diag(ugamma )))^2/ sum(diag( ugamma2))
    rTstats3 <- Tstats/c3
    rTstats3.b <- Tstats.b/c3
    rTstats3.s <- Tstats.s/c3
    rTstats3.s2 <- Tstats.s2/c3
    rTstats3.e1 <- Tstats.e1/c3
    rTstats3.e2 <- Tstats.e2/c3
    rTstats3.e3 <- Tstats.e3/c3
    pvalue4<-c()
    pvalue4[1]<-pchisq(rTstats3,df2,lower.tail =FALSE)
    pvalue4[2]<-pchisq(rTstats3.b,df2,lower.tail =FALSE)
    pvalue4[3]<-pchisq(rTstats3.s,df2,lower.tail =FALSE)
    pvalue4[4]<-pchisq(rTstats3.s2,df2,lower.tail =FALSE)
    pvalue4[5]<-pchisq(rTstats3.e1,df2,lower.tail =FALSE)
    pvalue4[6]<-pchisq(rTstats3.e2,df2,lower.tail =FALSE)
    pvalue4[7]<-pchisq(rTstats3.e3,df2,lower.tail =FALSE)
    
    #mean-variance correction + correct a2
    c4 <-a2c/sum(diag(ugamma )) 
    df2.2<-(sum(diag(ugamma )))^2/ a2c
    rTstats4 <- Tstats/c4
    rTstats4.b <- Tstats.b/c4
    rTstats4.s <- Tstats.s/c4
    rTstats4.s2 <- Tstats.s2/c4
    rTstats4.e1 <- Tstats.e1/c4
    rTstats4.e2 <- Tstats.e2/c4
    rTstats4.e3 <- Tstats.e3/c4
    pvalue5<-c()
    pvalue5[1]<-pchisq(rTstats4,df2.2,lower.tail =FALSE)
    pvalue5[2]<-pchisq(rTstats4.b,df2.2,lower.tail =FALSE)
    pvalue5[3]<-pchisq(rTstats4.s,df2.2,lower.tail =FALSE)
    pvalue5[4]<-pchisq(rTstats4.s2,df2.2,lower.tail =FALSE)
    pvalue5[5]<-pchisq(rTstats4.e1,df2.2,lower.tail =FALSE)
    pvalue5[6]<-pchisq(rTstats4.e2,df2.2,lower.tail =FALSE)
    pvalue5[7]<-pchisq(rTstats4.e3,df2.2,lower.tail =FALSE)
    
    #cor2; df
    c5<-0.5*(c1+c2)
    rTstats5 <- Tstats/c5
    rTstats5.b <- Tstats.b/c5
    rTstats5.s <- Tstats.s/c5
    rTstats5.s2 <- Tstats.s2/c5
    rTstats5.e1 <- Tstats.e1/c5
    rTstats5.e2 <- Tstats.e2/c5
    rTstats5.e3 <- Tstats.e3/c5
    pvalue6<-c()
    pvalue6[1]<-pchisq(rTstats5,df,lower.tail =FALSE)
    pvalue6[2]<-pchisq(rTstats5.b,df,lower.tail =FALSE)
    pvalue6[3]<-pchisq(rTstats5.s,df,lower.tail =FALSE)
    pvalue6[4]<-pchisq(rTstats5.s2,df,lower.tail =FALSE)
    pvalue6[5]<-pchisq(rTstats5.e1,df,lower.tail =FALSE)
    pvalue6[6]<-pchisq(rTstats5.e2,df,lower.tail =FALSE)
    pvalue6[7]<-pchisq(rTstats5.e3,df,lower.tail =FALSE)
    
    #cor3; df
    c6<-0.5*(1/c1+1/c2)
    rTstats6 <- Tstats*c6
    rTstats6.b <- Tstats.b*c6
    rTstats6.s <- Tstats.s*c6
    rTstats6.s2 <- Tstats.s2*c6
    rTstats6.e1 <- Tstats.e1*c6
    rTstats6.e2 <- Tstats.e2*c6
    rTstats6.e3 <- Tstats.e3*c6
    pvalue7<-c()
    pvalue7[1]<-pchisq(rTstats6,df,lower.tail =FALSE)
    pvalue7[2]<-pchisq(rTstats6.b,df,lower.tail =FALSE)
    pvalue7[3]<-pchisq(rTstats6.s,df,lower.tail =FALSE)
    pvalue7[4]<-pchisq(rTstats6.s2,df,lower.tail =FALSE)
    pvalue7[5]<-pchisq(rTstats6.e1,df,lower.tail =FALSE)
    pvalue7[6]<-pchisq(rTstats6.e2,df,lower.tail =FALSE)
    pvalue7[7]<-pchisq(rTstats6.e3,df,lower.tail =FALSE)
    
    #cor4
    pvalue8<-(pvalue2+pvalue4)/2
    
    #cor5
    pvalue9<-(pvalue2+pvalue5)/2
    
    ms<-Mardia.tran1(x)$skew
    mk<-Mardia.tran1(x)$kurt
    
    ######## gamma2; unbiased gamma
    ugamma.u <- U %*% gammaadf.u
    rank.u <-qr(ugamma.u)$rank
    #recaled correction
    c10 <- sum(diag(ugamma.u ))/df
    rTstats10 <- Tstats/c10
    pvalue10<-pchisq(rTstats10,df,lower.tail =FALSE)
    
    #rank deficient correction
    c11<-sum(diag(ugamma.u ))/rank.u
    rTstats11 <- Tstats/c11
    pvalue11<-pchisq(rTstats11,rank.u,lower.tail =FALSE)
    
    #mean-variance correction
    ugamma2.u <- ugamma.u %*% ugamma.u
    c12 <-sum(diag( ugamma2.u))/sum(diag(ugamma.u )) 
    df2.u<-(sum(diag(ugamma.u )))^2/ sum(diag( ugamma2.u))
    rTstats12 <- Tstats/c12
    pvalue12<-pchisq(rTstats12,df2.u,lower.tail =FALSE)
    
    ###### mean-variance correction a2    
    a2c.u<-1/(n*(n-1)*(n-2)*(n-3))*((n-2)*(n-1)^3*sum(diag(U%*%gammaadf.u%*%U%*%gammaadf.u))-
                                      n*(n-1)*sum(D^2)+sum(diag((n-1)*U%*%gammaadf.u))^2 )
    c13 <-a2c.u/sum(diag(ugamma.u )) 
    df2.2.u<-(sum(diag(ugamma.u )))^2/ a2c.u
    rTstats13 <- Tstats/c13
    pvalue13<-pchisq(rTstats13,df2.2.u,lower.tail =FALSE)
    
    result.sum[rep,]<-c(t,p,m,n,q,df,df2,df2.2,df2.u,df2.2.u,rank,rank.u,ms,mk,simdist,
                  Fstats,
                  pvalue1,pvalue2,pvalue3, pvalue4,pvalue5,pvalue6,
                  pvalue7,pvalue8,pvalue9,pvalue10,pvalue11,pvalue12,pvalue13,
                  theta_p,c(theta0),
                  c(r.SE),c(r.SE.u),c(SE))
      }
    } 
  }#replication

################################################
################################################

fn <- paste0(FILEPATH,"/result/",simdist,"/p",p,"m",m,"n",n,"result",it,".csv")
write.table(result.sum, file = fn,  sep = ",", col.names = FALSE,append = TRUE)

}
########################
do.call(do_one, args)
