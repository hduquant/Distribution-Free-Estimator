library(lavaan)
library(mvtnorm)
vech <- lav_matrix_vech
method <- "gammaN.model"

sim.norm <- function(seed,n,p,sig12){
  set.seed(seed)
  z <- rmvnorm(n=n, mean=rep(0,p))
  x <- z%*%sig12
  return(x)
}

  m=1
  p=5
  n=300
  ms=m*(m-1)/2
  q=2*p+ms
  ps=p*(p+1)/2
  p2=p*2
  df=ps-q
  p_m=p/m

set.seed(1)
lamb_uni<-sort(sample(seq(70,95,5),p_m,replace=TRUE))/100
lamb <- matrix(0,p,m)
for (j in 1:m){
  lamb[(p_m*(j-1)+1):(p_m*j),j] <- lamb_uni
}

lamb_t <- t(lamb)
phi <- 0.5*diag(m)+matrix(0.5,m,m)
eigen_phi <- eigen(phi)
phi_val <- eigen_phi$values
phi_vec <- eigen_phi$vectors
phi_h <- phi_vec%*%diag(sqrt(phi_val))%*%t(phi_vec)

lambphi_12 <- lamb%*%phi_h   # to simulate data from skewed distribution; phi^0.5
lamblamb <- lamb%*%phi%*%lamb_t
psi <- diag(p)-diag(diag(lamblamb))
psi_vec <- diag(psi)
psi12 <- sqrt(psi_vec) # to simulate data from skewed distribution

sig <- lamblamb+psi    #BIG SIGMA
eigen_sig <- eigen(sig)
sig_vec <- eigen_sig$vectors
sig_val <- eigen_sig$values
sig12 <-  sig_vec%*%diag(sqrt(sig_val))%*%t(sig_vec)
theta_p <- c(rep(lamb_uni,m),psi_vec,phi[lower.tri(phi, diag = FALSE)])

# simulate data
x <- sim.norm(seed=1,n,p,sig12)

#get sample covariance
  ps <- p*(p+1)/2
  Scov <-cov(x) #*(n-1)/n # consistent with  sample.cov.rescale = FALSE
  vscov <- vech(Scov)

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


# GLS based on DLS code
  ep <- 0.0001
  ps <- p*(p+1)/2
  q <- length(theta_p)
  df <- ps-q

  ###### W matrix
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
  
##### 
   a<-1
   set.seed(1)
    theta0=rep(1, length(theta_p)) # to mimic lavaan with start = "simple"
   
   sig <- sigdsig(p,m, theta0)  # gives a list of 2: sig$Sigma0; sig$vdsig
   vdsig<-sig$vdsig  #first order, sigma-dot
   Sigma0<-sig$Sigma0 # model_implied covariance matrix
   vsig0 <- vech(Sigma0)
   
   diconverg=0
######  iterations
for(t in 1:200) {

     ###### ## gammaN.model; updated with theta0 and sigma0
     gammaNm <- 2 * lavaan:::lav_matrix_duplication_ginv_pre_post(Sigma0 %x% Sigma0)
       
     weight<-try(solve((1-a)* gammaadf+a*gammaNm))
       if( class(weight)[1] =="try-error"){
         diconverg <- 1   # not converge
         break
       }
     
     stdi <- try(solve(t(vdsig) %*% weight%*%vdsig))
     if( class(stdi)[1]=="try-error"){
       diconverg <- 1   # not converge
       break
     }
     eresdu <- vscov-vsig0 #s-sigma
     dtheta <- t( eresdu) %*% weight %*%  vdsig%*%  stdi
     theta0.old <- theta0
     theta0 <- theta0 + dtheta
     delta <- max(abs(dtheta))
     rms.x <- sqrt(mean((theta0.old - theta0)^2))

     cat("rms.x = ", rms.x, "\n")
     
     sig <- sigdsig(p,m, theta0)  # gives a list of 2: sig$Sigma0; sig$vdsig
     vdsig<-sig$vdsig 
     Sigma0<-sig$Sigma0
     vsig0 <- vech(Sigma0)
     
     if(delta<=ep) {
         # we can stop
         gammaNm <- 2 * lavaan:::lav_matrix_duplication_ginv_pre_post(Sigma0 %x% Sigma0)
         weight<-try(solve((1-a)* gammaadf+a*gammaNm))
         if( class(weight)[1] =="try-error"){
           diconverg <- 1   # not converge
           break
         }
    
       break
     } # delta <= ep
}
   
   if(t<200){
     if(diconverg==0){  
       # Test start here;
       dtd<-try(solve(t(vdsig) %*% weight %*%vdsig))
       if( class(dtd)[1]=="try-error"){
         break
       }
       r.Var<-dtd%*%t(vdsig)%*%weight%*%gammaadf%*%weight%*%vdsig%*%dtd
       r.SE <- sqrt(diag(r.Var)/(n-1))
       Var<-dtd
       SE <- sqrt(diag(Var)/(n-1))
       dwe<-t(vdsig)%*%weight
       U<- weight-t(dwe)%*%dtd%*%dwe
       ugamma <-  U %*% gammaadf
       
       rank<-qr(ugamma)$rank
     }
   }


##########################
#lavaan
###########################
   i.com=f.com=ff.com=NULL
   for (mi in 1:m){
     p.m<-p/m
     begin<-paste0("NA*x",((mi-1)*p/m+1))
     i.com[[mi]]      <- paste(begin,"+",paste("x",((mi-1)*p/m+2):(mi*p/m),
                                               sep = "", 
                                               collapse = " + "))
     f.com[[mi]]      <-paste("f", mi,"=~ ", i.com[[mi]]  , sep = "")
     ff.com[[mi]]      <-paste("f", mi,"~~ ", "1*f" , mi, sep = "")
   }
   
   model <- paste(f.com,ff.com, sep = "\n")
   
   datalavaan<-x
   colnames(datalavaan)<-paste0("x",c(1:p))
   
   # Please compare to this one with model-implied covariance matrix and gauss-newton (gn) iteration
   fit <- cfa(model, data =   datalavaan, estimator = "DLS", verbose = TRUE,
              start = "simple", 
              optim.gn.tol.x = 1e-04,
              optim.gn.stephalf.max = 0L,
              estimator.args = list(dls.a = 1.0,
                                    dls.GammaNT = "model",
                                    dls.FtimesNminus1 = TRUE),
              sample.cov.rescale = FALSE)
   lavInspect(fit, "options")$optim.method
   
   # The one that you perhaps use is with sample covariance matrix and quasi-newton (nlminb) iteration
   fit2 <- cfa(model, data =   datalavaan,estimator="GLS")
   lavInspect(fit2, "options")$optim.method

   
   Delta <- lavTech(fit, "Delta")[[1]]
   Gamma <- lavTech(fit, "gamma")[[1]]
   WLS.V <- lavTech(fit, "WLS.V")[[1]]
   Sigma <- lavTech(fit, "implied")[[1]]$cov
   GammaNT <- 2 * lavaan:::lav_matrix_duplication_ginv_pre_post(Sigma %x% Sigma)
   lavaan.U = lavTech(fit,"UfromUGamma")[[1]]
   UGamma=lavTech(fit,"UGamma") 
   
   
##########################
#compare my results with lavaan
###########################  
# compare coefficients: no difference
cat("theta:\n")
print(range(parTable(fit)$est[-6]-theta0))

# compare sigma   
cat("Sigma:\n")
print(range(Sigma-Sigma0))

# compare gamma_ADF : tiny difference
cat("Gamma:\n")
print(range(Gamma-gammaadf))

# compare gamma_N : tiny difference
cat("GammaNT:\n")
print(range(GammaNT-gammaNm))

# compare delta
cat("Delta:\n")
print(range(Delta-vdsig))

# compare U
cat("U:\n")
print(range(lavaan.U-U))

# compare Ugamma
cat("UGamma:\n")
print(range(UGamma-ugamma))
   

   
   
