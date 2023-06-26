####|Testing for phylo signal in residuals of period (Table S1 - 1d) and stringency index models (Table S2 - 1c)

# packages and data
    library(ape)
    library(arm)
    library(brms)
    library(coda)
    library(dplyr)
    library(geiger)
    require(ggplot2)
    library(MCMCglmm)
    library(parallel)
    library(phangorn)
    library(phylobase)

    sapply(c('magrittr','data.table','brms','arm','lme4','car',
           'ape','phytools','geiger','stringr','plyr'),
         require, character.only = TRUE)

    trees =  read.tree("Data/trees.tre")
  
    d[, scinam := Species]

    s[, scinam := Species]
    

    # get residuals
      # period
        ms <- lmer(scale(log(FID)) ~
            scale(log(SD)) +
            scale(log(FlockSize)) +
            scale(log(BodyMass)) +
            scale(sin(rad)) + scale(cos(rad)) +
            # scale(Day)+
            scale(Temp) +
            scale(Covid) +
            (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (scale(Covid) | Country) + (scale(Covid) | IDLocality) + (1 | sp_loc),
        data = d, REML = FALSE,
        control <- lmerControl(
            optimizer = "optimx", optCtrl = list(method = "nlminb")
        )
        )
        d[, res_period := resid(ms)]
      # stringency
       m01a <- lmer(scale(log(FID)) ~
        scale(Year) +
        scale(log(SD)) +
        scale(log(FlockSize)) +
        scale(log(BodyMass)) +
        scale(sin(rad)) + scale(cos(rad)) +
        # scale(Day)+
        scale(Temp) +
        scale(StringencyIndex) +
        (1|weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (scale(StringencyIndex) | Country) + (1 | IDLocality) + (1 | sp_loc),
        data = s, REML = FALSE)
        s[, res := resid(m01a)]
      # google    
        ss = s[!is.na(parks_percent_change_from_baseline)]
        ss[, country_year := paste(Country, Year)] # table(paste(s$Country, s$Year))
        ss[parks_percent_change_from_baseline < 0, google := "before_zero"]
        ss[parks_percent_change_from_baseline > 0, google := "after_zero"]
        ss[, sp_country_google := paste(sp_country, google)]
        mg01a = lmer(scale(log(FID)) ~
            scale(Year) +
            scale(log(SD)) +
            scale(log(FlockSize)) +
            scale(log(BodyMass)) +
            scale(sin(rad)) + scale(cos(rad)) +
            # scale(Day)+
            scale(Temp) +
            scale(parks_percent_change_from_baseline) +
            (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (scale(parks_percent_change_from_baseline) | Country) + (1 | IDLocality) + (1 | sp_loc),
        data = ss, REML = FALSE
        )
        ss[, res := resid(mg01a)]
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
    
    namesss <- ss %>% distinct(Species)
    row.names(namesss) <- namesss$Species
    namecheck_fidss <- name.check(trees[[1]], namesss)
    trees_fidss <- lapply(trees, drop.tip, namecheck_fidss$tree_not_data)
    class(trees_fidss) <- "multiPhylo"
    length(trees_fidss[[1]]$tip.label)
    tree_fidss <- maxCladeCred(trees_fidss)
    tree_fidss$node.label <- NULL
    inv.phylo_ss <- inverseA(tree_fidss, nodes = "ALL", scale = TRUE)
    phyloresss <- solve(inv.phylo_ss$Ainv)
    rownames(phyloresss) <- rownames(inv.phylo_ss$Ainv)

# brmc
  load('Data/DAT_brms.Rdata')
  # period - no - need to run (takes time), if you load the above DAT_brms.Rdata  
    # without phylo
      priors <- get_prior(res_period ~ 0 + Intercept + (1|Species) , data=d)    #
      mP_no = brm(form    = res_period ~ 0 + Intercept + (1|Species), 
                    data    = d,   
                    cores   = 2,
                    chains  = 5,
                    control = list(adapt_delta = 0.99),
                    iter    = 5000,
                    thin = 5,
                    sample_prior="yes",
                    save_pars = save_pars(all = TRUE),
                    prior   = priors,
                    seed = 5 
                    )
       plot(mP_no, ask = FALSE)
       pp_check(mP_no, ndraws = 100)
       mcmc_plot(mP_no, type = "acf")
       summary(mP_no)   
    # with phylo
       priors <- get_prior(res_period ~ 0 + Intercept + (1|Species) + (1|gr(scinam, cov = A)), data=d, data2   = list(A = phyloresd))   
       mP_yes = brm(form    = res_period ~ 0 + Intercept + (1|Species) + (1|gr(scinam, cov = A)), 
                    data    = d,  
                    data2   = list(A = phyloresd), 
                    cores   = 2,
                    chains  = 5,
                    control = list(adapt_delta = 0.99),
                    iter    = 5000,
                    thin = 5,
                    sample_prior="yes",
                    save_pars = save_pars(all = TRUE),
                    prior   = priors,
                    seed = 5  
                    ) 
       plot(mP_yes, ask = FALSE)
       pp_check(mP_yes, ndraws = 100)
       #pp_check(m_yes, type = "scatter_avg_grouped", group = "Species") +  geom_abline(intercept = 0, slope = 1 , color = "red", lty = 2)
       mcmc_plot(mP_yes, type = "acf")
       summary(mP_yes)

       # lambda
          hypothesis(mP_yes, "sd_scinam__Intercept^2 / (sd_scinam__Intercept^2 + sd_Species__Intercept^2) = 0", class = NULL)
          hypothesis(mP_yes, "sd_scinam__Intercept^2 / (sd_scinam__Intercept^2 + sd_Species__Intercept^2 + sigma^2) = 0", class = NULL)


          v_sc <- (VarCorr(mP_yes, summary = FALSE)$scinam$sd)^2
          v_sp <- (VarCorr(mP_yes, summary = FALSE)$Species$sd)^2
          v_r <- (VarCorr(mP_yes, summary = FALSE)$residual$sd)^2
          summary(as.mcmc(v_sc / (v_sc + v_sp + v_r)))

  # compare period models
        # 1. LOOic
          loo(mP_yes)    
          loo(mP_no)   
  
  
        # 2. Bayes factor
          bayes_factor(mP_no, mP_yes)
  
  
        # 3. Poterior probability  
          post_prob(mP_yes, mP_no)
  
  # stringency index  - no - need to run (takes time), if you load the above DAT_brms.Rdata  
    # without phylo
      priors <- get_prior(res ~ 0 + Intercept + (1|Species) , data=s)    #
      m_no = brm(form    = res ~ 0 + Intercept + (1|Species), 
                    data    = s,   
                    cores   = 2,
                    chains  = 5,
                    control = list(adapt_delta = 0.99),
                    iter    = 5000,
                    thin = 5, 
                    sample_prior="yes",
                    save_pars = save_pars(all = TRUE),
                    prior   = priors,
                    seed = 5 
                    )
       plot(m_no,  ask = FALSE)
       pp_check(m_no, ndraws = 100)
       mcmc_plot(m_no, type = "acf")
       summary(m_no)    
    # with phylo
       priors <- get_prior(res ~ 0 + Intercept + (1|Species) + (1|gr(scinam, cov = A)), data=s, data2   = list(A = phyloress))   
       m_yes = brm(form    = res ~ 0 + Intercept + (1|Species) + (1|gr(scinam, cov = A)), 
                    data    = s,  
                    data2   = list(A = phyloress), 
                    cores   = 2,
                    chains  = 5,
                    control = list(adapt_delta = 0.99),
                    iter    = 5000,
                    thin = 5,
                    sample_prior="yes",
                    save_pars = save_pars(all = TRUE),
                    prior   = priors,
                    seed = 5 
                    ) 
       plot(m_yes, ask = FALSE)
       pp_check(m_yes, ndraws = 100)
       #pp_check(m_yes, type = "scatter_avg_grouped", group = "Species") +  geom_abline(intercept = 0, slope = 1 , color = "red", lty = 2)
       mcmc_plot(m_yes, type = "acf")
       summary(m_yes)

       # lambda
          hypothesis(m_yes, "sd_scinam__Intercept^2 / (sd_scinam__Intercept^2 + sd_Species__Intercept^2) = 0", class = NULL) # lambda
          hypothesis(m_yes, "sd_scinam__Intercept^2 / (sd_scinam__Intercept^2 + sd_Species__Intercept^2 + sigma^2) = 0", class = NULL)    

          v_sc <- (VarCorr(m_yes, summary = FALSE)$scinam$sd)^2
          v_sp <- (VarCorr(m_yes, summary = FALSE)$Species$sd)^2
          v_r <- (VarCorr(m_yes, summary = FALSE)$residual$sd)^2
          summary(as.mcmc(v_sc / (v_sc + v_sp + v_r)))
         

          Mode <- function(x) {
                  ux <- unique(x)
                  ux[which.max(tabulate(match(x, ux)))]
                }

          Mode(as.mcmc(v_sc / (v_sc + v_sp + v_r)))       

  # compare stringency models
        # 1. LOOic
          loo(m_yes)    
          loo(m_no)  
  
        # 2. Bayes factor
          bayes_factor(m_no, m_yes)
  
  
        # 3. Poterior probability  
          post_prob(m_yes, m_no)       
# Google Mobility  - no - need to run (takes time), if you load the above DAT_brms.Rdata
# without phylo
priors <- get_prior(res ~ 0 + Intercept + (1 | Species), data = ss) #
mG_no = brm(
    form = res ~ 0 + Intercept + (1 | Species),
    data = ss,
    cores = 2,
    chains = 5,
    control = list(adapt_delta = 0.99),
    iter = 5000,
    thin = 5,
    sample_prior = "yes",
    save_pars = save_pars(all = TRUE),
    prior = priors,
    seed = 5
)
plot(mG_no, ask = FALSE)
pp_check(mG_no, ndraws = 100)
mcmc_plot(mG_no, type = "acf")
summary(mG_no)
# with phylo
priors <- get_prior(res ~ 0 + Intercept + (1 | Species) + (1 | gr(scinam, cov = A)), data = ss, data2 = list(A = phyloress))
mG_yes = brm(
    form = res ~ 0 + Intercept + (1 | Species) + (1 | gr(scinam, cov = A)),
    data = ss,
    data2 = list(A = phyloress),
    cores = 2,
    chains = 5,
    control = list(adapt_delta = 0.99),
    iter = 5000,
    thin = 5,
    sample_prior = "yes",
    save_pars = save_pars(all = TRUE),
    prior = priors,
    seed = 5
)
plot(mG_yes, ask = FALSE)
pp_check(mG_yes, ndraws = 100)
# pp_check(m_yes, type = "scatter_avg_grouped", group = "Species") +  geom_abline(intercept = 0, slope = 1 , color = "red", lty = 2)
mcmc_plot(mG_yes, type = "acf")
summary(mG_yes)

# lambda
hypothesis(mG_yes, "sd_scinam__Intercept^2 / (sd_scinam__Intercept^2 + sd_Species__Intercept^2) = 0", class = NULL) # lambda
hypothesis(mG_yes, "sd_scinam__Intercept^2 / (sd_scinam__Intercept^2 + sd_Species__Intercept^2 + sigma^2) = 0", class = NULL)

v_sc <- (VarCorr(m_yes, summary = FALSE)$scinam$sd)^2
v_sp <- (VarCorr(m_yes, summary = FALSE)$Species$sd)^2
v_r <- (VarCorr(m_yes, summary = FALSE)$residual$sd)^2
summary(as.mcmc(v_sc / (v_sc + v_sp + v_r)))


Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

Mode(as.mcmc(v_sc / (v_sc + v_sp + v_r)))

# compare stringency models
# 1. LOOic
loo(mG_yes)
loo(mG_no)

# 2. Bayes factor
bayes_factor(mG_no, mG_yes)


# 3. Poterior probability
post_prob(mG_yes, mG_no)
  # export models
    # save(file = 'Data/rev_DAT_brms.Rdata', mP_no, mP_yes, m_no, m_yes, mG_no, mG_yes)

# END