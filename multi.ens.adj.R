##' @title multi.ens.adj
##' @name  multi.ens.adj
##' @author Michael Dietze \email{dietze@@bu.edu}
##' 
##' 
##' @param Xf Dataframe or matrix of forecast state variables for different ensembles.
##' @param cf vector assigning forecast disturbance class to each ensemble member
##' @param mu.f A vector with forecast mean estimates of state variables.
##' @param Pf  A cov matrix of forecast state variables.  
##' @param Xa A matrix of posterior samples of state variables.
##' @param ca A vector of posterior samples of disturbance class
##' 
##' @return Returns a matrix of adjusted analysis mean estimates of state variables and class assignments
##' @export

multi.ens.adj<-function(Xf, cf, mu.f, Pf, Xa, ca){
  
  if(FALSE){  ## set params
    ## single
    Xf = Nf
    mu.f = rbind(priors$mufN,priors$mufD)
    Pf   = abind::abind(solve(priors$pfN),
                 solve(priors$pfD),
                 along=3)
    cf = d.prior
    Xa = dat[,nsel]
    ca = dat[,"D"]
    ## multi class
    Xf   = Nf
    cf   = d.class
    mu.f = t(priors$muf)
    Pf   = abind::abind(    ## find a way to do this with apply(priors$pf,3,solve)
      solve(priors$pf[,,1]),
      solve(priors$pf[,,2]),
      solve(priors$pf[,,3]),
      along=3)
    Xa = dat[,Xsel]
    ca = d.class.post
  }
  if(FALSE){ ## figures
    plot(density(Xf[,1]),main="")
    lines(density(Xa[,1]),col=2)
    points(Xf[,1],rep(0,nrow(Xf)),pch="|")
    plot(ecdf(Xf[,1]))
    lines(ecdf(Xa[,1]),col=2)
  }
  
  ## reassign classes
  ff = table(cf)/length(cf) ## class frequency in the forecast
  uc = sort(unique(cf))
  fa = table_by(ca,uc)/length(ca) ## class frequency in the posterior
  df = pmax(fa-ff,0) ## difference in frequency
  df = df/sum(df)    ## posterior reassignment frequency
  cA = cf            ## class assigned in the Analysis
  for(i in seq_along(uc)){
    if(fa[i] < ff[i]){ ## if a class decreases in frequency in the analysis
      sel.c = which(cf == uc[i])
      cA[sel.c] = sample(uc,length(sel.c),prob=df,replace = TRUE)
    }
  }
  cases = table(cf,cA) ## summary of reassignments
  
  Z <- Xf*0
  
  for(k in seq_along(uc)){
    sel.c = which(cf == uc[k])
    
    ## SVD of forecast covariances
    S_f  <- svd(Pf[,,k])
    L_f  <- S_f$d
    V_f  <- S_f$v
    
    ## normalize
    for(i in sel.c){
      
      Z[i,] <- 1/sqrt(L_f) * t(V_f)%*%(Xf[i,]-mu.f[k,])
      
    }
    Z[is.na(Z)]<-0
    Z[is.infinite(Z)] <- 0
    if(FALSE){  ## graphical check on rescaling
      for(j in 1:3){
        plot(Xf[sel.c,j],Z[sel.c,j])
      }
    }    
    
  }
  
  ### ANALYSIS
  X_a <- Xf*0
  for(k in seq_along(uc)){
    sel.c = which(cA == uc[k])
    if(length(sel.c) == 0) next
    
    ## analysis (by class)
    mu.a <- colMeans(Xa[ca == uc[k],])
    Pa   <- cov(Xa[ca == uc[k],]) ## calculate analysis cov for each class
    S_a  <- svd(Pa)
    L_a  <- S_a$d
    V_a  <- S_a$v
  
    ## analysis ensemble 
    for(i in sel.c){
      # she decomposed Pa - then it's putting it back together but with a different Z which comes from the likelihood of that ens    
      X_a[i,] <- V_a %*%diag(sqrt(L_a))%*%Z[i,] + mu.a
    }
    
    if(sum(abs(mu.a - colMeans(X_a))) > 1)
      warning('Problem with ensemble adjustment (1)')
    
  }
  
  if(FALSE){ ## figures
    for(i in 1:3){
    plot(density(Xf[,i]),main="")
    lines(density(Xa[,i]),col=2)
    points(X_a[,i],rep(0,nrow(X_a)),pch="|",col=cf+1)
    }
    pairs(X_a)
  }
  
#  if(sum(diag(Pa) - diag(cov(X_a))) > 5 | sum(diag(Pa) - diag(cov(X_a))) < -5) logger.warn('Problem with ensemble adjustment (2)')
  
#  analysis <- as.data.frame(X_a)
  
  return(cbind(cA,X_a))
}


table_by <- function(x,cat){
  cnt = rep(0,length(cat))
  names(cnt) = cat
  for(i in seq_along(cat)){
    cnt[i] = sum(x==cat[i])
  }
  cnt
}
