####|Testing for phylo signal in residuals of 

# packages and data
    library(ape)
    library(arm)
    library(coda)
    library(dplyr)
    library(geiger)
    library(MCMCglmm)
    library(parallel)
    library(phangorn)
    library(phylobase)

    trees =  read.tree("Data/trees.tre")
    load('Data/DAT_res.Rdata')
    d = d_
    s = s_

# MCC tree and covariance matrix
    namesd<-d %>% distinct(Species)
    row.names(namesd) <- namesd$Species
    namecheck_fidd <- name.check(trees[[1]], namesd)
    trees_fidd <- lapply(trees, drop.tip, namecheck_fidd$tree_not_data)
    class(trees_fidd) <- "multiPhylo"
    length(trees_fidd[[1]]$tip.label)
    tree_fidd<-maxCladeCred(trees_fidd)
    tree_fidd$node.label<-NULL
    inv.phylo_d <- inverseA(tree_fidd, nodes="ALL", scale=TRUE)
    phyloresd <-solve(inv.phylo_d$Ainv)
    rownames(phyloresd) <-rownames(inv.phylo_d$Ainv)

    namess<-s %>% distinct(Species)
    row.names(namess) <- namess$Species
    namecheck_fids <- name.check(trees[[1]], namess)
    trees_fids <- lapply(trees, drop.tip, namecheck_fids$tree_not_data)
    class(trees_fids) <- "multiPhylo"
    length(trees_fids[[1]]$tip.label)
    tree_fids<-maxCladeCred(trees_fids)
    tree_fids$node.label<-NULL
    inv.phylo_s <- inverseA(tree_fids, nodes="ALL", scale=TRUE)
    phyloress <-solve(inv.phylo_s$Ainv)
    rownames(phyloress) <-rownames(inv.phylo_s$Ainv)

#MCMC phylogenetic models
  nitt = 1000000
  burnin=20000
  thin=1000
  n=(nitt-burnin)/thin # number of simulations per tree model
  p=0 # number of fixed parameters

  l = list()
  # stringency index
    mcmc1<-MCMCglmm(res ~ 1, random=~Species, family="gaussian", ginverse=list(Species = inv.phylo_s$Ainv), 
                    verbose=F, nitt=nitt, burnin=burnin, thin=thin, data=s)
    
    plot(mcmc1$Sol)
    plot(mcmc1$VCV)

    autocorr(mcmc1$Sol)
    autocorr.diag(mcmc1$VCV)

    #pr=length(unique(s$animal))+p+1 
    #sol_r=data.frame(mcmc1$Sol[1:n,(p+2):pr ])
                  
    #summary(mcmc1)
    lambda <- mcmc1$VCV[,"Species"] / (mcmc1$VCV[,"Species"] + mcmc1$VCV[,"units"])
    #summary(lambda)

    l[['s']] = data.table(model = 'Table S2 - 1c', lambda = posterior.mode(lambda), lwr = HPDinterval(lambda, 0.95)[1], lwr = HPDinterval(lambda, 0.95)[2])


  # period
    mcmc2<-MCMCglmm(res ~ 1, random=~Species, family="gaussian", ginverse=list(Species = inv.phylo_d$Ainv), 
                    verbose=F, nitt=nitt, burnin=burnin, thin=thin, data=d)
    
    plot(mcmc2$Sol)
    plot(mcmc2$VCV)

    autocorr(mcmc2$Sol)
    autocorr.diag(mcmc2$VCV)

    #summary(mcmc2)
    lambda <- mcmc2$VCV[,"Species"] / (mcmc2$VCV[,"Species"] + mcmc2$VCV[,"units"])
    #summary(lambda2)
    l[['d']] = data.table(model = 'Table S1 - 1d', lambda = posterior.mode(lambda), lwr = HPDinterval(lambda, 0.95)[1], lwr = HPDinterval(lambda, 0.95)[2])

  do.call(rbind, l)  