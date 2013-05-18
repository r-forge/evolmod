library(robustDist)

## set the seed
set.seed(942387484)

## define a stationary distribution of the Markov model
stat.dist = c(0.2,0.2,0.3,0.3)

## define a GTR Markov chain 
my.model = gtr.mc(c(0.5,0.3,0.6,0.2,0.3,0.2), stat.dist, T)

## compute eigen-decomposition of the generator
my.model.eigen = as.eigen(my.model)

## define evolutionary times for simulation
my.time = 0.1

## get transition probabilities  
my.prob = matexp(my.model.eigen,my.time)


## matrix that defines transition labeling
## count only 1<->2,1<->3,1<->4
reg.matrix = matrix(c(0,1,1,1,1,0,0,0,1,0,0,0,1,0,0,0), nrow = 4, ncol = 4)


## True distance according to my.model
true.dist = pair.conv.dist(rescale.mc(my.model,my.time), reg.matrix)


## set alignment length
sim.sites = 2000

## set number of simulations
sim.num = 200


## allocate memory of two types of estimators and true distance
conv.gtr.dist = numeric(sim.num)
robust.gtr.dist = numeric(sim.num)
conv.f84.dist = numeric(sim.num)
robust.f84.dist = numeric(sim.num)


for (i in 1:sim.num){
  
  ## simulate a pairwise alignment
  my.align = sim.pair.align(my.prob, my.model$stationary, sim.sites)
    
  ## estimate F84 and GTR models
  my.stat = emp.freq(my.align,4)
  my.f84 = fit.f84(my.stat,my.align)
  my.gtr = fit.mc(gtr.mc, c(0,0,0,0,0,0),my.stat,my.align)
    
  ## calculate eigen-decomposition of the estimated generators
  my.f84.eigen = as.eigen(my.f84)
  my.gtr.eigen = as.eigen(my.gtr)
  
  ## compute conventional distances
  conv.f84.dist[i] = pair.conv.dist(my.f84, reg.matrix)
  conv.gtr.dist[i] = pair.conv.dist(my.gtr, reg.matrix)
  
  ## compute robust distances
  robust.f84.dist[i] = pair.robust.dist(my.f84.eigen, reg.matrix, my.align)
  robust.gtr.dist[i] = pair.robust.dist(my.gtr.eigen, reg.matrix, my.align)
  
}


## plot simulation results

par(mfrow = c(2,2), cex.axis=1.5,cex.lab=1.5, mar=c(5,5,4,2)+0.1, cex.main=1.5)

min.x = min(conv.f84.dist,conv.gtr.dist,robust.f84.dist,robust.gtr.dist, true.dist) - 0.1*true.dist
max.x = max(conv.f84.dist,conv.gtr.dist,robust.f84.dist,robust.gtr.dist, true.dist) + 0.1*true.dist

hist(conv.f84.dist, main = "Conventional F84 Distances", col="royalblue", xlim=c(min.x,max.x),
     xlab="Mean Number of A<->X mutations")
abline(v=true.dist, lty=2, lwd=3, col="red")
legend(0.063, 150, legend=c("true\nvalue"), lty=2,lwd=3,col="red", cex=1.25)
box()

hist(conv.gtr.dist, main = "Conventional GTR Distances",col="royalblue", xlim=c(min.x,max.x),
          xlab="Mean Number of A<->X mutations")
abline(v=true.dist, lty=2,lwd=3, col="red")
box()

hist(robust.f84.dist,main = "Robust F84 Distances",col="royalblue", xlim=c(min.x,max.x),
     xlab="Mean Number of A<->X mutations")
abline(v=true.dist, lty=2, lwd=3, col="red")
box()

hist(robust.gtr.dist,main = "Robust GTR Distances",col="royalblue", xlim=c(min.x,max.x),
     xlab="Mean Number of A<->X mutations")
abline(v=true.dist, lty=2, lwd=3, col="red")
box()

