#Full tutorial
#https://darrenjw.wordpress.com/2010/08/15/metropolis-hastings-mcmc-algorithms/


#what does dnorm() do
#http://ww2.coastal.edu/kingw/statistics/R-tutorials/prob.html
metrop1=function(n=1000,eps=0.5)
{
  vec=vector("numeric", n)
  x=0
  vec[1]=x
  for (i in 2:n) {
    innov=runif(1,-eps,eps)
    can=x+innov
    aprob=min(1,dnorm(can)/dnorm(x))
    u=runif(1)
    if (u < aprob)
      x=can
    vec[i]=x
  }
  vec
}

#First a candidate value (can) is constructed by perturbing the current state of the 
#chain with the uniform innovation. Then the acceptance probability is computed, and 
#the chain is updated appropriately depending on whether the proposed new value is 
#accepted. Note the standard trick of picking an event with probability p by checking 
#to see if u<p, where u is a draw from a U(0,1). 

###improved version of metrop1####
#The next version of the code (below) makes use of the fact that min(1,A) is redundant 
#if you are just going to compare to a uniform (since values of the ratio larger than 
#1 will be accepted with probability 1, as they should), and also that it is unnecessary 
#to recalculate the likelihood of the current value of the chain at every iteration – better 
#to store and re-use the value. That obviously doesn’t make much difference here, but 
#for real problems likelihood computations are often the main computational bottleneck. 

metrop2=function(n=1000,eps=0.5)
{
  vec=vector("numeric", n)
  x=0
  oldlik=dnorm(x)
  vec[1]=x
  for (i in 2:n) {
    innov=runif(1,-eps,eps)
    can=x+innov
    lik=dnorm(can)
    a=lik/oldlik
    u=runif(1)
    if (u < a) {
      x=can
      oldlik=lik
    }
    vec[i]=x
  }
  vec
}

### more improvement ###
#However, this code still has a very serious flaw. It computes likelihoods! 
#For this problem it isn’t a major issue, but in general likelihoods are the 
#product of a very large number of small numbers, and numerical underflow 
#is a constant hazard. For this reason (and others), no-one ever computes 
#likelihoods in code if they can possibly avoid it, but instead log-likelihoods 
#(which are the sum of a large number of reasonably-sized numbers, and therefore 
#numerically very stable). We can use these log-likelihoods to calculate a 
#log-acceptance ratio, which can then be compared to the log of a uniform for 
#the accept-reject decision. We end up with the code below, which now doesn’t 
#look too different to the kind of code one might actually write… 

metrop3=function(n=1000,eps=0.5)
{
  vec=vector("numeric", n)
  x=0
  oldll=dnorm(x,log=TRUE)
  vec[1]=x
  for (i in 2:n) {
    can=x+runif(1,-eps,eps)
    loglik=dnorm(can,log=TRUE)
    loga=loglik-oldll
    if (log(runif(1)) < loga) {
      x=can
      oldll=loglik
    }
    vec[i]=x
  }
  vec
}

#output drawing functionalities
plot.mcmc<-function(mcmc.out)
{
  op=par(mfrow=c(2,2))
  plot(ts(mcmc.out),col=2)
  hist(mcmc.out,30,col=3)
  qqnorm(mcmc.out,col=4)
  abline(0,1,col=2)
  acf(mcmc.out,col=2,lag.max=100)
  par(op)
}

metrop.out<-metrop1(10000,1)
plot.mcmc(metrop.out)