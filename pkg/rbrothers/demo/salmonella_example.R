library(rbrothers)

### Salmonella fimH Analysis ########

my.fimH.align = read.dna(file=system.file("extdata/Salmonella/Salmonella_fimH8.phy",package="rbrothers"))
write.dna(my.fimH.align,"Salmonella_fimH8.phy")

fimH.db<-dualbrothers(round(runif(1,1,100000),0),"Salmonella_fimH8", top_lambda=0.001, par_lambda=1, window.size=c(300,400), length=51000000, burnin=1000000, subsample=5000)


plot(fimH.db,seetrees="TRUE",threshold=0.5)
plotmcmc(fimH.db)
summary(fimH.db)
breakpointCI(fimH.db,350,500)
breakpointCI(fimH.db,570,700)


### Salmonella fimA Analysis ########

my.fimA.align = read.dna(file=system.file("extdata/Salmonella/Salmonella_fimA8.phy",package="rbrothers"))
write.dna(my.fimA.align,"Salmonella_fimA8.phy")

fimA.db<-dualbrothers(round(runif(1,1,10000),0),"Salmonella_fimA8", top_lambda=0.001, par_lambda=1, window.size=250)


plot(fimA.db,seetrees="TRUE",threshold=0.5)
plotmcmc(fimA.db)
summary(fimA.db)






