tobit <- function(x){
  x[x < 0] <- 0
  x
}

tobit2 <- function(x,upper){
  x[x < 0] <- 0
  x[x>upper] <- upper
  x
}


disturbance <- function(x,                         ## x = Bleaf, Bsoil, Bstem
                        mu0 = c(1,1),              ## mean leaf & stem biomass after disturbance
                        V0 = diag(c(0.15,0.25))^2, ## var in post disturbance leaf and stem biomass, 
                        alloc.soil = 0.25){        ## fraction of removed C that goes in the soil
  old.biomass <- x[c(1,3)]
  new.biomass <- tobit(mvtnorm::rmvnorm(1,mu0,V0)) ## draw disturbed leaf and stem
  check.max <- which(new.biomass > old.biomass)
  if(length(check.max)>0) new.biomass[check.max] = old.biomass[check.max]
  residual = sum(old.biomass-new.biomass)
  x[c(1,3)] <- new.biomass
  x[2] <- x[2] + residual*alloc.soil
  removal <- residual*(1-alloc.soil)
  return(x)
}

disturbance.p <- function(x,                         ## x = Bleaf, Bsoil, Bstem
                        mu0 = c(1,1),              ## mean leaf & stem % retention after disturbance
                        V0 = diag(c(0.15,0.25))^2, ## var in post disturbance leaf and stem biomass, 
                        alloc.soil = 0.25){        ## fraction of removed C that goes in the soil
  old.biomass <- x[c(1,3)]
  new.biomass <- tobit(mvtnorm::rmvnorm(1,mu0*old.biomass,tcrossprod(old.biomass,old.biomass)*V0)) ## draw disturbed leaf and stem
  check.max <- which(new.biomass > old.biomass)
  if(length(check.max)>0) new.biomass[check.max] = old.biomass[check.max]
  residual = sum(old.biomass-new.biomass)
  x[c(1,3)] <- new.biomass
  x[2] <- x[2] + residual*alloc.soil
  removal <- residual*(1-alloc.soil)
  return(x)
}

disturbance2 <- function(x,type){ ## x = Bleaf, Bsoil, Bstem
  old.biomass <- x[c(1,3)]
  new.biomass <- tobit(mvtnorm::rmvnorm(1,mu0[,type]*old.biomass,V0[,,type])) ## draw disturbed leaf and stem
  check.max <- which(new.biomass > old.biomass)
  if(length(check.max)>0) new.biomass[check.max] = old.biomass[check.max]
  residual = sum(old.biomass-new.biomass)
  x[c(1,3)] <- new.biomass
  x[2] <- x[2] + residual*alloc.soil
  removal <- residual*(1-alloc.soil)
  return(x)
}

# random correlated bernouli
rcbern <- function(p,rho,mu){
  
  ## calculate cov
  D = as.matrix(dist(seq_along(mu),diag = TRUE,upper = TRUE))
  SIGMA <- 1/(1-rho^2)*rho^D
  
  ## draw rMVN
  x = mvtnorm::rmvnorm(1,rep(0,length(p)),SIGMA)
  
  ## inverse distribution transform
  y = pnorm(x)
  qbinom(y,1,p)
  
}