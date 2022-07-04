####|Testing for phylo signal in residuals of period (Table S1 - 1d) and stringency index models (Table S2 - 1c)

# packages and data
    library(ape)
    library(arm)
    library(brms)
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
 
# MCMC phylogenetic models
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
    # using Jarrod's prior - improves posterior distribution of Species
        # an inverse-gamma prior for the residual variance with shape and scale equal to 0.001 was used for the residual variance (that is, variance set to 1 and the degree of belief parameter to 0.002).
        # the parameter-expanded priors for Species to give scaled F-distributions with numerator and denominator degrees of freedom set to 1 and a scale parameter of 1,000.
        prior_Jarrod<-list(R=list(V=1,nu=0.002),
                                    G=list(
                                        G0=list(V=1,nu=1, alpha.mu=0, alpha.V=1000)#,
                                        #G3=list(V=1,nu=0.002)                                      
                                        ))  
        mcmc2_2<-MCMCglmm(res ~ 1, random=~Species, family="gaussian", ginverse=list(Species = inv.phylo_d$Ainv), 
                    verbose=F, nitt=nitt, burnin=burnin, thin=thin, data=d, prior = prior_Jarrod)
        plot(mcmc2_2$Sol)
        plot(mcmc2_2$VCV)
        plot(log(mcmc2_2$VCV))

        autocorr(mcmc2_2$Sol)
        autocorr.diag(mcmc2_2$VCV)

        lambda <- mcmc2_2$VCV[,"Species"] / (mcmc2_2$VCV[,"Species"] + mcmc2_2$VCV[,"units"])
        #summary(lambda2)
        l[['d']] = data.table(model = 'Table S1 - 1d', lambda = posterior.mode(lambda), lwr = HPDinterval(lambda, 0.95)[1], upr = HPDinterval(lambda, 0.95)[2])
    
    # using in-build prior leads to poor posterior distribution of Species
        mcmc2<-MCMCglmm(res ~ 1, random=~Species, family="gaussian", ginverse=list(Species = inv.phylo_d$Ainv), 
                        verbose=F, nitt=nitt, burnin=burnin, thin=thin, data=d)

        plot(mcmc2$Sol)
        plot(mcmc2$VCV)
        plot(log(mcmc2$VCV))

        autocorr(mcmc2$Sol)
        autocorr.diag(mcmc2$VCV)

        #summary(mcmc2)
        #lambda <- mcmc2$VCV[,"Species"] / (mcmc2$VCV[,"Species"] + mcmc2$VCV[,"units"])
        #summary(lambda2)
        #l[['d']] = data.table(model = 'Table S1 - 1d', lambda = posterior.mode(lambda), lwr = HPDinterval(lambda, 0.95)[1], lwr = HPDinterval(lambda, 0.95)[2])
  do.call(rbind, l)  

# compare model with and without phylogeny
    mcmc2_no = MCMCglmm(res ~ 1, random=~Species, family="gaussian",
                    verbose=F, nitt=nitt, burnin=burnin, thin=thin, data=d)
    mcmc2_no$DIC-mcmc2$DIC # 0.02710599

    mcmc1_no = MCMCglmm(res ~ 1, random=~Species, family="gaussian",  
                    verbose=F, nitt=nitt, burnin=burnin, thin=thin, data=s)
    mcmc1_no$DIC-mcmc1$DIC # -10.94093
# export
  save(file = 'Data/DAT_mcmc.Rdata', mcmc1_no, mcmc1, mcmc2_no, mcmc2)

# End   
# END