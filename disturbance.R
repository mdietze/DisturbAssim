tobit <- function(x){
  x[x < 0] <- 0
  x
}

disturbance <- function(x){ ## x = Bleaf, Bsoil, Bstem
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