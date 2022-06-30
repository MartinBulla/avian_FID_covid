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
mcmc1<-MCMCglmm(res ~ 1, random=~Species, family="gaussian", ginverse=list(Species = inv.phylo_s$Ainv), 
                verbose=F, nitt=1000000, burnin=20000, thin=100, data=s)
summary(mcmc1)
lambda <- mcmc1$VCV[,"Species"] / (mcmc1$VCV[,"Species"] + mcmc1$VCV[,"units"])
lambda <- posterior.mode(lambda)
summary(lambda)

mcmc2<-MCMCglmm(res ~ 1, random=~Species, family="gaussian", ginverse=list(Species = inv.phylo_d$Ainv), 
                verbose=F, nitt=1000000, burnin=20000, thin=100, data=d)
summary(mcmc2)
lambda2 <- mcmc2$VCV[,"Species"] / (mcmc2$VCV[,"Species"] + mcmc2$VCV[,"units"])
lambda2 <- posterior.mode(lambda2)
summary(lambda2)