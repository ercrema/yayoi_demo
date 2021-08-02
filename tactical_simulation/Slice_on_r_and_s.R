library(here)
library(nimbleCarbon)
load(here('tactical_simulation','simdata.RData'))

## Nimble Functions and Models ----
cov_GPL2 <- nimbleFunction(
    run = function(dists = double(2), rhosq = double(0), etasq = double(0), sigmasq = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    deltaij <- matrix(nrow = n, ncol = n,init = TRUE)
    diag(deltaij) <- 1
    for(i in 1:n)
      for(j in 1:n)
        result[i, j] <- etasq*exp(-rhosq*dists[i,j]^2)+sigmasq*deltaij[i,j]
    return(result)
  })
Ccov_GPL2 <- compileNimble(cov_GPL2)

model <- nimbleCode({
  for (i in 1:N){
    # Model
    local_r[i] <- r + s[i]
    theta[i] ~ dExponentialGrowth(a=a,b=b,r=local_r[i])
    
    c14age[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
    sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
    sigmaDate[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
    cra[i] ~ dnorm(mean=c14age[i],sd=sigmaDate[i]);
  }
  #priors
  r ~ dnorm(0,sd=0.001)
  etasq ~ dexp(1000);
  rhosq ~ dexp(1000);
  mu_s[1:N] <- 0;
  cov_s[1:N, 1:N] <- cov_GPL2(dists[1:N, 1:N], rhosq, etasq, 0.000001)
  s[1:N] ~ dmnorm(mu_s[1:N], cov = cov_s[1:N, 1:N])
})

## Compile Model and Setup MCMC
model <- nimbleModel(model,constants = constants,data=data,inits=inits)
cModel <- compileNimble(model)
conf <- configureMCMC(model)
conf$addMonitors('s')
conf$removeSamplers('s[1:500]')
conf$removeSampler('r')
conf$addSampler(c('r','s[1:500]'), type='AF_slice') 

MCMC <- buildMCMC(conf)
cMCMC <- compileNimble(MCMC)
results <- runMCMC(cMCMC,inits=inits, nchain=1,niter = 10000, thin=5,nburnin = 5000, progressBar = TRUE) 
save(results,file=here('tactical_simulation','slice_sr_res.RData')) 