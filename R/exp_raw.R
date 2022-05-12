#' ---
#' title: "Some plots for avian FIDs before & during COVID using raw data"
#' author: "Martin Bulla"
#' date: "`r Sys.time()`"
#' output: 
#'     html_document:
#'         toc: true
#'         toc_float: true
#'         toc_depth: 5
#'         code_folding: hide
#'         bibliography: shorebird_bipInc_sex.bibtex
#'         link-citations: yes
#' ---

#+ r setup, include=FALSE 
knitr::opts_chunk$set(message = FALSE, warning = FALSE, cache = TRUE)

#' ##### Code to load tools & data
  # tools
    require(arm)
    require(data.table)
    require(effects)
    require(here)
    require(ggimage)
    require(ggplot2)
    require(ggpubr)
    require(grid)
    require(gtable)
    require(MASS)
    require(multcomp)
    require(optimx)
    require(PerformanceAnalytics)
    require(rphylopic)
    require(viridis)
  
    # constants
    round_ = 3 # number of decimal places to round model coefficients
    nsim = 5000 # number of simulations to extract estimates and 95%CrI
    ax_lines = "grey60" # defines color of the axis lines
    colors <- c("#999999", "#E69F00", "#56B4E9") #viridis(3)

    # function to remove ggplot components
      gtable_filter_remove <- function (x, name, trim = TRUE){
        matches <- !(x$layout$name %in% name)
        x$layout <- x$layout[matches, , drop = FALSE]
        x$grobs <- x$grobs[matches]
        if (trim) 
          x <- gtable_trim(x)
        x
      }
    # customized ggplot theme
      theme_MB = theme(  
                title = element_text(size=8, colour="grey30"),
                axis.line = element_blank(),
                #axis.line = element_line(colour="grey70", size=0.25),
                axis.title = element_text(size=7, colour="grey30"),
                axis.title.y = element_text(vjust=3.5),
                axis.title.x = element_text(vjust=1),
                axis.text = element_text(size=6),#, vjust = 0.5, hjust=1),# margin=units(0.5,"mm")),
                axis.ticks.length=unit(0.5,"mm"),
                axis.ticks = element_line(colour = "grey70", size = 0.1),
                #axis.ticks.margin,
                
                strip.text.x = element_text(size = 6, color="grey30",  margin=margin(1,1,1,1,"mm")), #grey50
                strip.text.y = element_text(size = 6, color="grey30",  margin=margin(1,1,1,1,"mm")), #grey50
                strip.background = element_rect(fill="grey99",colour="grey70", size=0.25),
                  #strip.background = element_blank(), 
                  #strip.text = element_blank(),
                panel.spacing = unit(0, "mm"),
                panel.background=element_blank(),
                panel.border = element_rect(colour="grey70", size=0.1, fill = NA), #panel.border=element_blank(),
                panel.grid = element_blank(),

                legend.text=element_text(size=6),
                legend.title=element_text(size=6),
                legend.key = element_rect(colour = NA, fill = NA),
                legend.key.height= unit(0.5,"line"),
                legend.key.width = unit(0.25, "cm"),
                legend.margin = margin(0,0,0,0, unit="cm"),
                legend.box.margin = margin(l = -6), #legend.justification = c(-1,0),
                legend.background = element_blank()
                )  
    # function for estimates
      est_out =function(model = m, label = "", nsim = 5000){
          bsim = sim(model, n.sim=5000) 
          v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
          ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
          data.table(predictor=rownames(coef(summary(model))),estimate=v, lwr=ci[1,], upr=ci[2,], model = paste(label, "N =", nobs(model)))
        }
  # data
    o  =  fread('Data/phylopic.txt')
    setnames(o, old = c('Name', 'Code'), new = c('genus2', 'uid'))
    t = fread(here::here('Data/taxonomy.txt'))
    d = fread(here::here('Data/data.txt'))
    # adjust correct assignment of season (Year)
    d[Country == 'Australia' & Year == 2020 & Covid == 0, Year:=2019]
    d[Country == 'Australia' & Year == 2021 & Day>139, Year:=2020]
    d[Country == 'Australia' & Year == 2022 & Day>139, Year:=2021]

    d = d[order(Year, IDLocality, Day, Hour)]
    d[Country%in%'Czech_Republic', Country:='Czech Republic']
    d[, genus := sub("_.*", "", Species)]
    d[, sp_day_year := paste(Year, Species, Day, sep="_")]
    d[, sp_loc := paste(Species, IDLocality, sep="_")]
    d[, sp_country := paste(Species, Country, sep="\n")]
    d[, rad:=(2*pi*Hour) / 24]
    d[, Day_:= Day]
    #d[, FID_z := scale(FID), by = Species]
    #d[, SD_z := scale(SD), by = Species]
    d[, FID_ln := log(FID)]
    d[, SD_ln := log(SD)]
    d[, body_ln := log(BodyMass)]
    d[, flock_ln := log(FlockSize)]
    #d[Country == 'Australia', Day_:= abs(Day - 189)]
    
    d1 = d[Covid == 1, .N, by = Species]
    d2 = d[Covid == 0, .N, by = Species]
    setnames(d1, old = 'N', new ='N_during')
    setnames(d2, old = 'N', new ='N_before')
    dd = merge(d1,d2) # species with data before and during
    da = merge(d1,d2, all = TRUE)

    d1p = d[Country!='Poland' & Covid == 1, .N, by = Species]
    d2p = d[Country!='Poland' & Covid == 0, .N, by = Species]
    setnames(d1p, old = 'N', new ='N_during')
    setnames(d2p, old = 'N', new ='N_before')
    ddp = merge(d1p,d2p) # species with data before and during, but without Poland data

    p1 = d[Covid == 1, .N, by = .(IDLocality, Species)]
    p2 = d[Covid == 0, .N, by = .(IDLocality, Species)]
    setnames(p1, old = 'N', new ='N_during')
    setnames(p2, old = 'N', new ='N_before')
    pp = merge(p1,p2)  # species-localities with data before and during
    pa = merge(p1,p2, all = TRUE)

    p1p = d[Country!='Poland' & Covid == 1, .N, by = .(IDLocality, Species)]
    p2p = d[Country!='Poland' & Covid == 0, .N, by = .(IDLocality, Species)]
    setnames(p1p, old = 'N', new ='N_during')
    setnames(p2p, old = 'N', new ='N_before')
    ppp = merge(p1p,p2p)  # species-localities with data before and during,but without Poland data 
   
    p1p4 = d[Year!=2014 & Country!='Poland' & Covid == 1, .N, by = .(IDLocality, Species)]
    p2p4 = d[Year!=2014 &Country!='Poland' & Covid == 0, .N, by = .(IDLocality, Species)]
    setnames(p1p4, old = 'N', new ='N_during')
    setnames(p2p4, old = 'N', new ='N_before')
    ppp4 = merge(p1p4,p2p4) # species-localities with data before and during, but without Poland & 2014 data

#' ### Models
#' To control for pseudo-replication we fitted year, genus, species, country, species at specific location and species at a given day of a given year as random intercepts and period (before lockdown and during lockdown) as random slope to all random intercepts. Given the singular fit, we removed random slope from intercepts where correlations were -1 and kept the model with random slope for country and species at sampling location.
  
  # prepare estimates 
    # 01a all data, all random slopes - singularity 
      mf_max=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) + (Covid|genus)+(Covid|Species)+(1|sp_day_year) + (Covid|Country) + (Covid|sp_loc),
                  data = d)# # (Covid|IDLocality) +
      est_mf_01a = est_out(mf_max, '01a) (Year) + (Covid|genus) + (Covid|Species) + (1|sp_day_year) + (Covid|Country) + (Covid|sp_loc)')  
    # 01b all data, all random slopes, but some without cor to avoid singularity
      mf_max_=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) +(0+Covid|genus)+(0+Covid|Species)+(1|sp_day_year) + (Covid|Country) + (0+Covid|IDLocality) +(Covid|sp_loc),
                  data = d, REML = FALSE) # (Covid|IDLocality) +
      est_mf_01b = est_out(mf_max, '01b) (1|Year) + (0+Covid|genus) + (0+Covid|Species) + (0+Covid|sp_day_year) + (Covid|Country) + (0+Covid|IDLocality) + (Covid|sp_loc)')   
    # 01c all data, random slopes that allow for non-singular fit 
      mf=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) +(1|genus)+(1|Species)+(1|sp_day_year) + (Covid|Country) + (Covid|sp_loc),
                  data = d, REML = FALSE, control = lmerControl(
                           optimizer ='optimx', optCtrl=list(method='nlminb')))# (Covid|IDLocality) +
      #d[,res := resid(mf)]
      est_mf_01c = est_out(mf, '01c) (1|Year) +(1|genus) + (1|Species) + (1|sp_day_year) + (Covid|Country) + (Covid|sp_loc)')
    # 01d no 2014, random slopes that allow for non-singular fit 
      mf1d=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) +(1|genus)+(1|Species)+(1|sp_day_year) + (Covid|Country) + (Covid|sp_loc),
                  data = d[Year!=2014], REML = FALSE, control = lmerControl(
                           optimizer ='optimx', optCtrl=list(method='nlminb')))# (Covid|IDLocality) +
      est_mf_01d = est_out(mf1d, '01d) (1|Year) +(1|genus) + (1|Species) + (1|sp_day_year) + (Covid|Country) + (Covid|sp_loc); no 2014')
    # 01e no PL, random slopes that allow for non-singular fit 
      mf1e=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) +(1|genus)+(1|Species)+(1|sp_day_year) + (Covid|Country) + (Covid|sp_loc),
                  data = d[Country!='Poland'], REML = FALSE, control = lmerControl(
                           optimizer ='optimx', optCtrl=list(method='nlminb')))# (Covid|IDLocality) +
      est_mf_01e = est_out(mf1e, '01e) (1|Year) +(1|genus) + (1|Species) + (1|sp_day_year) + (Covid|Country) + (Covid|sp_loc); no Poland')
    # 01e no 2014 '& PL, random slopes that allow for non-singular fit 
      mf1f=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) +(1|genus)+(1|Species)+(1|sp_day_year) + (Covid|Country) + (Covid|IDLocality) + (1|sp_loc),
                  data = d[Year!=2014 & Country!='Poland'], REML = FALSE, control = lmerControl(
                           optimizer ='optimx', optCtrl=list(method='nlminb')))# (Covid|IDLocality) +
      est_mf_01f = est_out(mf1f, '01f) (1|Year) +(1|genus) + (1|Species) + (1|sp_day_year) + (Covid|Country) + (Covid|IDLocality)+ (1|sp_loc); no Poland & no 2014')  
    
    # 02a) >9 per species
      d[, Nsp := .N, Species]
      mf_10=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) + (1|genus)+(1|Species)+(1|sp_day_year) + (Covid|Country) + (Covid|IDLocality) + (Covid|sp_loc),
                  data = d[Nsp>9], REML = FALSE,
                  control = lmerControl(
                           optimizer ='optimx', optCtrl=list(method='nlminb'))) # (Covid|IDLocality) +
      est_mf_02a = est_out(mf_10, '02a) (1|Year) + (1|genus) + (1|Species) + (1|sp_day_year) + (Covid|Country) + (Covid|IDLocality) + (Covid|sp_loc); >9/species')  
    # 02b) >9 per species, no P
      dp = d[Country!='Poland']
      dp[, Nsp := .N, Species]
      mf_10p=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) + (1|genus)+(1|Species)+(1|sp_day_year) + (Covid|Country) + (Covid|IDLocality) + (Covid|sp_loc),
                  data = dp[Nsp>9], REML = FALSE) # (Covid|IDLocality) +
      est_mf_02b = est_out(mf_10p, '02b) (1|Year) + (1|genus) + (1|Species) + (1|sp_day_year) + (Covid|Country) + (Covid|IDLocality) + (Covid|sp_loc); >9/species, no Poland')    
    
    # 03a) before & during > 4/species, random slopes that allow for non-singular fit 
      dx = dd[N_during>4 & N_before >4]
      #dxx = d[Species %in% dx$Species]
      #dx2 = dd[N_during>9 & N_before >9]
      mf_5ba=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) + (1|genus)+(1|Species)+(1|sp_day_year) + (Covid|Country) + (Covid|sp_loc),
                  data = d[Species %in% dx$Species],
                  REML = FALSE) # (Covid|IDLocality) +
      est_mf_03a = est_out(mf_5ba, '03a) (1|Year) +(1|genus) + (1|Species) + (1|sp_day_year) + (Covid|Country) + (Covid|sp_loc); >4/species/period')  
    # 03b) before & during > 4/species, random structure that allow for non-singular fit no Poland 
      dx = ddp[N_during>4 & N_before >4]
      #dxx = d[Species %in% dx$Species]
      #dx2 = dd[N_during>9 & N_before >9]
      mf_5baP=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) + (1|genus)+(1|Species)+(1|sp_day_year) + (Covid|Country) + (Covid|sp_loc),
                  data = d[Species %in% dx$Species],
                  REML = FALSE,
                  control = lmerControl(
                      optimizer ='optimx', optCtrl=list(method='nlminb'))) # (Covid|IDLocality) +
      est_mf_03b = est_out(mf_5baP, '03b) (1|Year) +(1|genus) + (1|Species) + (1|sp_day_year) + (Covid|Country) + (Covid|sp_loc); >4/species/period & no Poland')  

    # 04a) before & during > 9/species, random slopes that allow for non-singular fit 
      dx = dd[N_during>9 & N_before >9]
      #dxx = d[Species %in% dx$Species]
      #dx2 = dd[N_during>9 & N_before >9]
      mf_9ba=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) + (1|Species)+(1|sp_day_year) + (Covid|Country) + (Covid|sp_loc),
                  data = d[Species %in% dx$Species],
                  REML = FALSE) # (1|genus) explained 0 so taken out, but for reporting could stay(Covid|IDLocality) +
      est_mf_04a = est_out(mf_9ba, '04a) (1|Year) +  (1|Species) + (1|sp_day_year) + (Covid|Country) + (Covid|sp_loc); >9/species/period')  
    # 04b) before & during > 9/species, random slopes that allow for non-singular fit no Poland
      dx = ddp[N_during>9 & N_before >9]
      #dxx = d[Species %in% dx$Species]
      #dx2 = dd[N_during>9 & N_before >9]
      mf_9baP=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) + (1|genus) + (1|Species)+(1|sp_day_year) + (Covid|Country) + (Covid|sp_loc),
                  data = d[Species %in% dx$Species],
                  REML = FALSE) # (1|genus) explained 0 so taken out, but for reporting could stay(Covid|IDLocality) +
      est_mf_04b = est_out(mf_9baP, '04b) (1|Year) + (1|genus) + (1|Species) + (1|sp_day_year) + (Covid|Country) + (Covid|sp_loc); >9/species/period & no Poland')  

    # 05a) before & during > 4/species/locality/covid, random structure that allow for non-singular fit
      dx = pp[N_during>4 & N_before >4]
      dxx = d[paste(IDLocality, Species) %in% paste(dx$IDLocality, dx$Species)]
      #dxx = d[Species %in% dx$Species]
      #dx2 = dd[N_during>9 & N_before >9]
      mf_05a=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) + (Covid|Species)+ (1|sp_day_year) +(0+Covid|Country) + (Covid|IDLocality) + (1|sp_loc),
                  data = dxx,
                  REML = FALSE, control = lmerControl(
                  optimizer ='optimx', optCtrl=list(method='nlminb'))) 
                  # (1|genus) explained 0 so taken out, but for reporting could stay
                  # could be also with (1|Species) + (0+Covid|sp_loc)
      est_mf_05a = est_out(mf_05a, '05a) (1|Year) + (1|genus)+(Covid|Species) + (1|sp_day_year)+ (0+Covid|Country) + (Covid|IDLocality) + (1|sp_loc); >4/species/locality/period')   
    # 05b) before & during > 4/species/locality/covid, random structure that allow for non-singular fit no poland
      dx = ppp[N_during>4 & N_before >4]
      dxx = d[paste(IDLocality, Species) %in% paste(dx$IDLocality, dx$Species)]
      #dxx = d[Species %in% dx$Species]
      #dx2 = dd[N_during>9 & N_before >9]
      mf_05b=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) + (Covid|Species)+ (1|sp_day_year) +(0+Covid|Country) + (Covid|IDLocality) + (1|sp_loc),
                  data = dxx,
                  REML = FALSE, control = lmerControl(
                  optimizer ='optimx', optCtrl=list(method='nlminb'))) 
                  # (1|genus) explained 0 so taken out, but for reporting could stay
                  # could be also with (1|Species) + (0+Covid|sp_loc)
      est_mf_05b = est_out(mf_05b, '05b) (1|Year) + (1|genus)+(Covid|Species) + (1|sp_day_year)+ (0+Covid|Country) + (Covid|IDLocality) + (1|sp_loc); >4/species/locality/period & no Poland')  
    # 05c) diff model structures for previous
      dx = ppp[N_during>4 & N_before >4]
      dxx = d[paste(IDLocality, Species) %in% paste(dx$IDLocality, dx$Species)]
      #dxx = d[Species %in% dx$Species]
      #dx2 = dd[N_during>9 & N_before >9]
      mf_05c=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) + (1|Species)+ (1|sp_day_year) +(0+Covid|Country) + (Covid|IDLocality) + (0+Covid|sp_loc),
                  data = dxx,
                  REML = FALSE,
                  control = lmerControl(
                  optimizer ='optimx', optCtrl=list(method='nlminb'))) # (Covid|IDLocality) +
      est_mf_05c = est_out(mf_05c, '05c) (1|Year) + (1|Species) + (1|sp_day_year) + (0+Covid|Country) + (Covid|IDLocality) +(0+Covid|sp_loc); >4/species/locality/period & no Poland')  
    # 05d) & without 2014
      dx = ppp4[N_during>4 & N_before >4]
      dxx = d[paste(IDLocality, Species) %in% paste(dx$IDLocality, dx$Species)]
      #dxx = d[Species %in% dx$Species]
      #dx2 = dd[N_during>9 & N_before >9]
      mf_05d=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) + (Covid|Species)+ (1|sp_day_year) +(0+Covid|Country) + (Covid|IDLocality) + (0+Covid|sp_loc),
                  data = dxx,
                  REML = FALSE, control = lmerControl(
                  optimizer ='optimx', optCtrl=list(method='nlminb'))) 

      est_mf_05d = est_out(mf_05d, '05d) (1|Year) + (Covid|Species) + (0+Covid|Country) + (Covid|IDLocality) +(0+Covid|sp_loc); >4/species/locality/period no Poland & no 2014')  
    
    # 06a) 1 obs/species/locality/covid
      dxx = d[paste(IDLocality, Species) %in% paste(ppp$IDLocality, ppp$Species)]
      dxx[, pk:=1:nrow(dxx)]
      #length(dxx[, unique(sp_loc)])
      #length(dxx[, unique(paste(IDLocality, Species, Covid))])
      set.seed(42)
      dxx = dxx[,.SD[sample(.N, min(1,.N))],by = .(IDLocality, Species, Covid)]

      m06a=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (Covid|Species)+ (1|sp_day_year) +(1|IDLocality) + (1|sp_loc), # runs also with (Covid|IDLocality), but seed 1
                  data = dxx,
                  REML = FALSE) # (1|Year)  Country & Genus can stay (but explains 0 var)
        est_m06a = est_out(m06a, '06a) (Covid|Species)+ (1|sp_day_year) +(1|IDLocality) + (1|sp_loc); 1/species/locality/period no Poland')  
    # 06b) 1 obs/species/locality/covid from localities with >4 before/after
      dx = ppp[N_during>4 & N_before >4]
      dxx = d[paste(IDLocality, Species) %in% paste(dx$IDLocality, dx$Species)]
      dxx[, pk:=1:nrow(dxx)]
      set.seed(42)
      dxx = dxx[,.SD[sample(.N, min(1,.N))],by = .(IDLocality, Species, Covid)]

      m06b=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Species), #similar if sp_loc used
                  data = dxx,
                  REML = FALSE
                  ) 
                  # year/genus/Country/IDLocality can stay (but explains 0 var)
                  # also possible (1|Country) or (0+Covid|IDLocality) or (0+Covid|Country) + (0+Covid|IDLocality)
      # singular fit but just because some random effects estimated as null
      est_m06b = est_out(m06b, '06b) (1|Species); 1/species/locality/period no Poland & before/during >4')  

    # 07a) avg obs/species/locality/covid from localities
      dxx = d[paste(IDLocality, Species) %in% paste(ppp$IDLocality, ppp$Species)]
      #length(dxx[, unique(paste(sp_loc, Covid))])
      m = lm(log(FID) ~ log(SD),dxx)
      dxx[, resid_FID := resid(m)]
      a = dxx[, mean(resid_FID) ,by = .(Country, IDLocality, genus, Species, sp_loc, Covid)]
      setnames(a, old = 'V1', new = 'resid_FID_avg')
      #ggplot(a, aes(x = resid_FID_avg)) +geom_histogram()
     
      m07a=lmer(scale(resid_FID_avg)~
                  scale(Covid)+
                   (Covid|Species)+ (Covid|Country) + (Covid|IDLocality) + (1|sp_loc),
                  data = a,
                  REML = FALSE) # (1|genus) explains 0, so could stay
      est_m07a = est_out(m07a, '07a) (1|genus) + (Covid|Species)+ (Covid|Country) + (Covid|IDLocality) + (1|sp_loc); 1 average resid_FID/species/locality/period no Poland')  
    # 07b) avg obs/species/locality/covid from localities with >4 before/afteer
      dx = ppp[N_during>4 & N_before >4]
      dxx = d[paste(IDLocality, Species) %in% paste(dx$IDLocality, dx$Species)]
      dxx[, pk:=1:nrow(dxx)]
      m = lm(log(FID) ~ log(SD),dxx)
      dxx[, resid_FID := resid(m)]
      a = dxx[, mean(resid_FID) ,by = .(Country, IDLocality, genus, Species, sp_loc, Covid)]
      setnames(a, old = 'V1', new = 'resid_FID_avg')
      #ggplot(a, aes(x = resid_FID_avg)) +geom_histogram()
     
      m07b=lmer(scale(resid_FID_avg)~
                  scale(Covid)+
                  (1|genus) + (Covid|Species)+ (0+Covid|Country) + (Covid|IDLocality) + (1|sp_loc),
                  data = a,
                  REML = FALSE,
                  control = lmerControl(
                  optimizer ='optimx', optCtrl=list(method='nlminb')))  # (Covid|IDLocality) +
      est_m07b = est_out(m07b, '07b) (1|genus) + (Covid|Species)+ (0+Covid|Country) + (Covid|IDLocality) + (1|sp_loc); 1 average resid_FID/species/locality/period no Poland & before/during >4')  
    # 08) avg per species and country from above
      dxx = d[paste(IDLocality, Species) %in% paste(ppp$IDLocality, ppp$Species)]
      #length(dxx[, unique(paste(sp_loc, Covid))])
      m = lm(log(FID) ~ log(SD),dxx)
      dxx[, resid_FID := resid(m)]
      a = dxx[, mean(resid_FID) ,by = .(Country, IDLocality, genus, Species, sp_loc, Covid)]
      setnames(a, old = 'V1', new = 'resid_FID_avg')
      aa=a[,mean(resid_FID_avg),.(genus, Country, Species, Covid)]
      setnames(aa, old = 'V1', new = 'resid_FID_avg')
      #ggplot(a, aes(x = resid_FID_avg)) +geom_histogram()
     
      m08=lmer(scale(resid_FID_avg)~
                  scale(Covid)+
                  (1|genus) + (1|Species)+ (1|Country), #(Covid|Country) gives -1 cor
                  data = aa,
                  REML = FALSE)
      
      est_m08 = est_out(m08, '08) (1|genus) + (1|Species)+ (1|Country) ; 1 average resid_FID/species/country/period no Poland')  
    # 09) avg per species from above
      dxx = d[paste(IDLocality, Species) %in% paste(pp$IDLocality, pp$Species) & Country!='Poland']
      m = lm(log(FID) ~ log(SD),dxx)
      dxx[, resid_FID := resid(m)]
      a = dxx[, mean(resid_FID) ,by = .(Country, IDLocality, genus, Species, sp_loc, Covid)]
      setnames(a, old = 'V1', new = 'resid_FID_avg')
      aa=a[,mean(resid_FID_avg),.(genus, Country, Species, Covid)]
      setnames(aa, old = 'V1', new = 'resid_FID_avg')
      aaa=aa[,mean(resid_FID_avg),.(genus, Species, Covid)]
      setnames(aaa, old = 'V1', new = 'resid_FID_avg')
      #ggplot(a, aes(x = resid_FID_avg)) +geom_histogram()
     
      m09=lmer(scale(resid_FID_avg)~
                  scale(Covid)+
                  (1|genus) + (1|Species), # similar as if only (Covid|genus) used
                  data = aaa,
                  REML = FALSE) 
      est_m09 = est_out(m09, '09) (1|genus) + (1|Species); 1 average resid_FID/species/period no Poland')  
    
      t.test(aaa$resid_FID_avg[aaa$Covid == 0], aaa$resid_FID_avg[aaa$Covid == 1],paired = TRUE, alternative = "two.sided")
      
  # plot
      x = rbind(est_mf_01a,est_mf_01b, est_mf_01c, est_mf_01d, est_mf_01e, est_mf_01f, est_mf_02a,est_mf_02b, est_mf_03a,est_mf_03b,est_mf_04b,est_mf_04a, est_mf_05a,est_mf_05b,est_mf_05c,est_mf_05d,est_m06a,est_m06b,est_m07a,est_m07b,est_m08,est_m09)

      g = 
      ggplot(x[predictor == 'scale(Covid)'], aes(y = model, x = estimate, col = model)) +
          geom_vline(xintercept = 0, col = "grey30", lty =3)+
          geom_errorbar(aes(xmin = lwr, xmax = upr, col = model), width = 0, position = position_dodge(width = 0.01) ) +
          #ggtitle ("Sim based")+
          geom_point(position = position_dodge(width = 0.01)) +
          #scale_colour_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          #scale_fill_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
          scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE)) + 
          scale_y_discrete(limits=rev)+
          coord_fixed(ratio = 0.05)+
          #scale_shape(guide = guide_legend(reverse = TRUE)) + 
          #scale_x_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2,2, by = 1), labels = seq(-2,2, by = 1)) +
          labs(y = NULL ,x = "Standardized effect size") +
          #ylim(c(0,100))+
          #coord_flip()+
          theme_bw() +
          theme( legend.position ="none",
                plot.title = element_text(size=7),
                legend.title=element_text(size=7), 
                legend.text=element_text(size=6),
                ##legend.spacing.y = unit(0.1, 'cm'), 
                legend.key.height= unit(0.5,"line"),
                #plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r =0.5, unit =  "pt"),
                panel.grid = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = ax_lines, size = 0.25),
                axis.line.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x= element_line( colour = ax_lines, size = 0.25),
                axis.ticks.length = unit(1, "pt"),
                axis.text.x = element_text(colour="black", size = 6),
                axis.text.y=element_text(colour="black", size = 7),
                axis.title=element_text(size=9)
                )
      g
      ggsave(here::here('Outputs/effect_sizes_2022-05-12.png'),g, width = 30, height =8, units = 'cm')

  # model ass
    acf(resid(m16), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))

#' plot before and during
  # prepare data
     dxx = d[paste(IDLocality, Species) %in% paste(pp$IDLocality, pp$Species)]
        #length(dxx[, unique(paste(sp_loc, Covid))])
     m = lm(log(FID) ~ log(SD),dxx)
     dxx[, resid_FID := resid(m)]
     a = dxx[, .(mean(resid_FID),sd(resid_FID),mean(FID),.N) ,by = .(Country, IDLocality, genus, Species, sp_loc, Covid)]
     setnames(a, old = c('V1','V2','V3'), new = c('resid_FID_avg', 'SD','FID_avg'))
     a[is.na(SD), SD := 0]   

     aw = reshape(a, idvar = c('Country', 'IDLocality', 'genus', 'Species', 'sp_loc'), timevar = 'Covid', direction = "wide")  
     aw[, Species := gsub('[_]', ' ', Species)]
     aw = merge(aw,t, all.x =TRUE)
     table(aw$Family)
   
     x = aw[, .N, by = Species]
     x[order(Species)]
     aw[, genus2 := genus]
     aw[Species %in% x[N%in%c(1,2), Species], genus2:='other']

  # check
     ggplot(a, aes(x = resid_FID_avg)) + geom_histogram()
     ggplot(a, aes(x=FID_avg, y = resid_FID_avg)) + 
        stat_smooth() + 
        geom_point() +
        stat_cor(aes(label = ..r.label..),  label.x = 0.3, size = 2) +
        scale_x_continuous(trans = 'log10')
    ggplot(aw, aes(x = N.0-N.1)) +geom_histogram()
    nrow(aw[abs(N.0-N.1)>2])
    nrow(aw[!abs(N.0-N.1)>2])
  # USE removed x-axis
    o[genus2=='Motacilla' | uid%in%c('67a9ecfd-58ba-44a4-9986-243b6e610419'), uid:='cf522e02-35cc-44f5-841c-0e642987c2e4']
    o[,size:=0.2]
    o[ genus2%in%c('Anas','Columba','Dendrocopos','Sturnus'), size := c(0.25, 0.25, 0.15, 0.1)]
    o[, FID_avg.0 := 1.5]
    o[, FID_avg.1 := 20]
    o[ genus2%in%c('Anas','Columba'), FID_avg.0 := c(1.7, 1.7)]

    g = 
    ggplot(aw, aes(x = FID_avg.0, y = FID_avg.1)) + 
      #geom_errorbar(aes(ymin = FID_avg.1-SD.1, ymax = FID_avg.1+SD.1, col = Country), width = 0) +
      #geom_errorbar(aes(xmin = FID_avg.0-SD.0, xmax = FID_avg.0+SD.0, col = Country), width = 0) +
      geom_point(pch = 21, alpha = 0.7, aes(col = Country)) + 
        #ggtitle ("Sim based")+
      geom_abline(intercept = 0, slope = 1, lty =3, col = "grey80")+
      facet_wrap(~genus2) +
      geom_phylopic(data = o, aes(image = uid),  color = "grey80", size = o$size) + # ,
      scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
      scale_x_continuous("Before COVID flight initiation distance", expand = c(0, 0), trans = 'log10') +
      scale_y_continuous("During COVID flight initiation distance", expand = c(0, 0), trans = 'log10') +
      labs(title = "Species means per sampling location")+
      theme_MB  +
      theme(
              plot.title = element_text(size=7),
              strip.background = element_blank(),
              #panel.spacing = unit(1, "mm"),
              legend.position = c(1, 0),
              legend.justification = c(1, 0)
              )  
    gg <- ggplotGrob(g) #gg$layout$name
    ggx <- gtable_filter_remove(gg, name = paste0("axis-b-", c(2, 4), "-4"),
                                     trim = FALSE)
    grid.draw(ggx)
    ggsave('Outputs/before_after_Genus.png',ggx, width=4.5,height=4,dpi=600)
   
  # resid per species
     ggplot(aw, aes(x = resid_FID_avg.0, y = resid_FID_avg.1)) + 
      geom_point(pch = 21, alpha = 0.7, aes(col = Country)) + 
      geom_abline(intercept = 0, slope = 1, lty =3, col = "grey80")+
      facet_wrap(~Species) +
      scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
      scale_x_continuous("residual FID before COVID", expand = c(0, 0)) +
      scale_y_continuous("residual FID during COVID", expand = c(0, 0)) +
      labs(title = "Species means per sampling location")+
      theme_MB +
      theme(
            plot.title = element_text(size=7),
            strip.background = element_blank()
            )
  # resid per genus2
     ggplot(aw, aes(x = resid_FID_avg.0, y = resid_FID_avg.1)) + 
      geom_point(pch = 21, alpha = 0.7, aes(col = Country)) + 
      geom_abline(intercept = 0, slope = 1, lty =3, col = "grey80")+
      facet_wrap(~genus2) +
      scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
      scale_x_continuous("residual FID before COVID", expand = c(0, 0)) +
      scale_y_continuous("residual FID during COVID", expand = c(0, 0)) +
      labs(title = "Species means per sampling location")+
      theme_MB  +
      theme(
            plot.title = element_text(size=7),
            strip.background = element_blank(),
            legend.position = c(1, 0),
            legend.justification = c(1, 0)
            )
  # FID per genus2
     ggplot(aw, aes(x = FID_avg.0, y = FID_avg.1)) + 
      geom_point(pch = 21, alpha = 0.7, aes(col = Country)) + 
      geom_abline(intercept = 0, slope = 1, lty =3, col = "grey80")+
      facet_wrap(~genus2) +
      scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
      scale_x_continuous("FID before COVID", expand = c(0, 0)) +
      scale_y_continuous("FID during COVID]", expand = c(0, 0)) +
      labs(title = "Species means per sampling location")+
      theme_MB +
      theme(
            plot.title = element_text(size=7),
            strip.background = element_blank(),
            legend.position = c(1, 0),
            legend.justification = c(1, 0)
            )
  # FID per genus2 - flexible x
     ggplot(aw, aes(x = FID_avg.0, y = FID_avg.1)) + 
      geom_point(pch = 21, alpha = 0.7, aes(col = Country)) + 
      geom_abline(intercept = 0, slope = 1, lty =3, col = "grey80")+
      facet_wrap(~genus2, scales = 'free') +
      scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
      scale_x_continuous("FID before COVID", expand = c(0, 0)) +
      scale_y_continuous("FID during COVID]", expand = c(0, 0)) +
      labs(title = "Species means per sampling location")+
      theme_MB +
      theme(
            plot.title = element_text(size=7),
            strip.background = element_blank(),
            legend.position = c(1, 0),
            legend.justification = c(1, 0)
            )
  # log FID per genus2 - moved x-axis
    g = 
    ggplot(aw, aes(x = FID_avg.0, y = FID_avg.1)) + 
      geom_point(pch = 21, alpha = 0.7, aes(col = Country)) + 
      geom_abline(intercept = 0, slope = 1, lty =3, col = "grey80")+
      facet_wrap(~genus2) +
      scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
      scale_x_continuous("FID before COVID", expand = expansion(mult =c(0,0.075)), trans = 'log10') +
      scale_y_continuous("FID during COVID", expand = c(0, 0), trans = 'log10') +
      labs(title = "Species means per sampling location")+
      theme_MB  +
      theme(
              plot.title = element_text(size=7),
              strip.background = element_blank(),
              #panel.spacing = unit(1, "mm"),
              legend.position = c(1, 0),
              legend.justification = c(1, 0)
              )  



#' show no difference graphically
#' 
#' ### exploration
#' #### sample sizes
    nrow(dd) # only 55 species scored during and after covid - is that an issue?)
    nrow(da[!Species%in%dd$Species & is.na(N_before)]) # the rest 80 scored  during covid, and only 13 before
    nrow(da[!Species%in%dd$Species & !is.na(N_before)]) # only 13 before


    data.frame(table(d$Species,d$Covid))
    table(d$Year)
    summary(as.factor(data.frame(table(d$sp_day_year))$Freq))
    length(unique(d$sp_day_year))
  

    nrow(dd) # only 55 species scored during and after covid - is that an issue?)
    nrow(da[!Species%in%dd$Species & is.na(N_before)]) # the rest 80 scored  during covid, and only 13 before
    nrow(da[!Species%in%dd$Species & !is.na(N_before)]) # only 13 before

    length(unique(pp$Species))
    nrow(pp[N_during>4 & N_before>4])

     table(dx$IDLocality)
     table(dx$IDLocality, dx$Year)   
#' #### distributions & correlations  
   
   chart.Correlation(d[, c('FID_ln', 'SD_ln', 'flock_ln', 'body_ln', 'Hour','Temp', 'Day')], histogram=TRUE, pch=19)
   mtext("Single observations", side=3, line=3)
    #ggplot(d, aes(x = FID))+geom_histogram()
    #ggplot(d, aes(x = log(FID)))+geom_histogram()
    #ggplot(d, aes(x = FID_z))+geom_histogram()
    #ggplot(d, aes(x = log(FID)))+geom_histogram()
    #ggplot(d, aes(x = SD))+geom_histogram()
    #ggplot(d, aes(x = log(SD)))+geom_histogram()
    #ggplot(d, aes(x = FlockSize))+geom_histogram()
    #ggplot(d, aes(x = BodyMass))+geom_histogram()
    #ggplot(d, aes(x = Hour))+geom_histogram()
    #ggplot(d, aes(x = Temp))+geom_histogram()
    #ggplot(d, aes(x = Day))+geom_histogram()
    ggplot(d, aes(x = Day, y = Temp, col = Country))+geom_point(pch = 21, alpha =0.5)
    #ggplot(d, aes(x = Day_, y = Temp, col = Country))+geom_point(pch = 21, alpha =0.5)



   ggplot(d, aes(x = Year, y = FID)) + geom_smooth() + scale_y_continuous(trans = 'log10')
   ggplot(d, aes(x = Day, y = FID)) + geom_smooth() + scale_y_continuous(trans = 'log10')
   ggplot(d, aes(x = Day, y = FID)) + geom_point() + geom_smooth() + scale_y_continuous(trans = 'log10')
   
   ggplot(d, aes(x = Day, y = FID)) +  geom_smooth(aes(col = as.factor(Year))) + scale_y_continuous(trans = 'log10') +  scale_color_viridis(discrete=TRUE)
   ggplot(d, aes(x = Day, y = FID)) + geom_point(aes(col = Country), alpha = 0.8) + geom_smooth(col ='red', fill = 'red') + scale_y_continuous(trans = 'log10') + facet_wrap(~Year, nrow = 1) + scale_color_viridis(discrete=TRUE) 
   
   ggplot(d[!Country%in%'Australia'], aes(x = Day, y = FID)) + geom_point(aes(col = Country), alpha = 0.8) + geom_smooth(col ='red', fill = 'red') + scale_y_continuous(trans = 'log10') + facet_wrap(~Year, nrow = 1) + scale_color_viridis(discrete=TRUE) 
   ggplot(d[!Country%in%c('Australia','Poland')], aes(x = Day, y = FID)) + geom_point(aes(col = Country), alpha = 0.8) + geom_smooth(col ='red', fill = 'red') + scale_y_continuous(trans = 'log10') + facet_wrap(~Year, nrow = 1) + scale_color_viridis(discrete=TRUE) 
   
   ggplot(d, aes(x = Day, y = FID)) + geom_point(aes(col = Country), alpha = 0.2) + geom_smooth(aes(col =Country)) + scale_y_continuous(trans = 'log10') + facet_wrap(~Year, nrow = 1) + scale_color_viridis(discrete=TRUE) 


   ggplot(d, aes(x = Day, y = FID)) + geom_point(aes(col = Country, shape = as.factor(Covid)), alpha = 0.8) + geom_smooth(col ='red', fill = 'red') + scale_y_continuous(trans = 'log10') + facet_wrap(~Year, nrow = 1) + scale_color_viridis(discrete=TRUE) 
   
   ggplot(d, aes(x = paste(Year,Day), y = FID)) + geom_point(aes(col = Country), alpha = 0.8) + geom_smooth(col ='red', fill = 'red') + scale_y_continuous(trans = 'log10') + scale_color_viridis(discrete=TRUE) 

   px = pp[N_during>4 & N_before>4]
   dxx = d[paste(IDLocality, Species) %in% paste(px$IDLocality, px$Species)]
    table(dxx$IDLocality, dxx$Year)
   g = 
   ggplot(dxx, aes(x = as.factor(Year), y = FID, col = Year)) + geom_boxplot() + facet_wrap(~sp_loc) + scale_y_continuous(trans = 'log10') + scale_x_discrete(guide = guide_axis(angle = 45))

   ggsave(here::here('Outputs/year_trend.png'),g, width = 28, height =25, units = 'cm')
#' #### stringency 
   s = d[Covid == 1]
   s[, Nsp := .N, by ='Species']
   s[, sp := gsub('[_]', ' ', Species)]
  ##### Explore
     ggplot(s, aes(x=StringencyIndex)) + geom_histogram()
     ggplot(s, aes(x=StringencyIndex)) + geom_histogram()+ facet_wrap(~Country)
     ggplot(s, aes(x = Day, y = FID)) + geom_point(aes(col = Country), alpha = 0.8) + geom_smooth(col ='red', fill = 'red') + scale_y_continuous(trans = 'log10') + facet_wrap(~Year, nrow = 1) + scale_color_viridis(discrete=TRUE) 

     ggplot(s, aes(x = Day, y = StringencyIndex)) + geom_point(aes(col = Country), alpha = 0.8) + geom_smooth(col ='red', fill = 'red') +  facet_wrap(~Year, nrow = 1) + scale_color_viridis(discrete=TRUE) 
    
     table(s$Species)
  ##### MODEL outputs
    # estimates 
     m01a=lmer(scale(log(FID))~
          scale(log(SD))+
          scale(log(FlockSize))+
          scale(log(BodyMass))+
          scale(sin(rad)) + scale(cos(rad)) + 
          #scale(Day)+
          scale(Temp)+
          scale(StringencyIndex)+
          (scale(StringencyIndex)|Species) + (1|Country) + (1|IDLocality),  
          data = s, REML = FALSE, 
          control = lmerControl(
              optimizer ='optimx', optCtrl=list(method='L-BFGS-B'))
          ) 
          # (1|Year), (1|genus), (1|sp_day_year), (1|sp_loc),  explain nothing - could stay
     est_m01a = est_out(m01a, '01a) (scale(StringencyIndex)|Species) + (1|Country) + (1|IDLocality)')  


     m01b=lmer(scale(log(FID))~ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(StringencyIndex)+ 
        (scale(StringencyIndex)|genus)+(1|Species)+(1|sp_day_year) + 
        (1|Country) + (scale(StringencyIndex)|IDLocality) +(1|sp_loc),   
        data = s, REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
     s[, res := resid(m01b)]
        # (1|Year) explains nothing - could stay 
     est_m01b = est_out(m01b, '01b) (scale(StringencyIndex)|genus)+(1|Species)+(1|sp_day_year) + (1|Country) + (scale(StringencyIndex)|IDLocality) +(1|sp_loc)') 

     m02a=lmer(scale(log(FID))~
          scale(log(SD))+
          scale(log(FlockSize))+
          scale(log(BodyMass))+
          scale(sin(rad)) + scale(cos(rad)) + 
          #scale(Day)+
          scale(Temp)+
          scale(StringencyIndex)+
          (scale(StringencyIndex)|Species) + (1|Country) + (1|IDLocality),  
          data = s[Nsp>4], REML = FALSE, 
          control = lmerControl(
              optimizer ='optimx', optCtrl=list(method='L-BFGS-B'))
          ) 
     est_m02a = est_out(m02a, '02a)  (scale(StringencyIndex)|Species) + (1|Country) + (1|IDLocality); >4/species') 

     m02b=lmer(scale(log(FID))~ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(StringencyIndex)+ 
        (scale(StringencyIndex)|genus)+(1|Species)+(1|sp_day_year) + 
        (1|Country) + (scale(StringencyIndex)|IDLocality) +(1|sp_loc),   
        data = s[Nsp>4], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
        # (1|Year) explains nothing - could stay 
     est_m02b = est_out(m02b, '02b) (scale(StringencyIndex)|genus)+(1|Species)+(1|sp_day_year) + (1|Country) + (scale(StringencyIndex)|IDLocality) +(1|sp_loc);>4/species') 
    # plot
      xs = rbind(est_m01a,est_m01b, est_m02a, est_m02b)

      g = 
      ggplot(xs[predictor == 'scale(StringencyIndex)'], aes(y = model, x = estimate, col = model)) +
          geom_vline(xintercept = 0, col = "grey30", lty =3)+
          geom_errorbar(aes(xmin = lwr, xmax = upr, col = model), width = 0, position = position_dodge(width = 0.01) ) +
          #ggtitle ("Sim based")+
          geom_point(position = position_dodge(width = 0.01)) +
          #scale_colour_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          #scale_fill_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
          scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE)) + 
          scale_y_discrete(limits=rev)+
          coord_fixed(ratio = 0.05)+
          #scale_shape(guide = guide_legend(reverse = TRUE)) + 
          #scale_x_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2,2, by = 1), labels = seq(-2,2, by = 1)) +
          labs(y = NULL ,x = "Standardized effect size") +
          #ylim(c(0,100))+
          #coord_flip()+
          theme_bw() +
          theme( legend.position ="none",
                plot.title = element_text(size=7),
                legend.title=element_text(size=7), 
                legend.text=element_text(size=6),
                ##legend.spacing.y = unit(0.1, 'cm'), 
                legend.key.height= unit(0.5,"line"),
                #plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r =0.5, unit =  "pt"),
                panel.grid = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = ax_lines, size = 0.25),
                axis.line.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x= element_line( colour = ax_lines, size = 0.25),
                axis.ticks.length = unit(1, "pt"),
                axis.text.x = element_text(colour="black", size = 6),
                axis.text.y=element_text(colour="black", size = 7),
                axis.title=element_text(size=9)
                )
      g
      ggsave(here::here('Outputs/effect_sizes_stringency.png'),g, width = 30, height =3, units = 'cm')
  ##### visualise
    g = 
    ggplot(s[Nsp>9], aes(x = StringencyIndex, y = FID)) +
      stat_smooth(se = FALSE, aes(colour = 'LOESS smooth'), lwd = 0.5)+ # show_guide=TRUE
      #stat_smooth(method = 'rlm', se = FALSE, col = 'black', lwd = 0.5)+
      geom_point(pch = 21, alpha = 0.7, aes(fill = Country), col = 'white') + 
      facet_wrap(~sp, ncol = 6) +
      scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
      scale_x_continuous("COVID measure stringency index", expand = c(0, 0)) +
      scale_y_continuous("Flight initiation distance [m]", expand = c(0, 0), trans = 'log10') +
      #annotate("text", x = 1, y = 1, label = c(rep("", 52),"Observation"), hjust = -0.08, size = 1) +
      #labs(title = "Species means per sampling location")+
      scale_colour_manual(values=c('grey60'))+
      #scale_color_manual(name = 'try', values = c('LOESS smoothed = "grey60"'))+
      theme_MB  +
      theme(
          plot.title = element_text(size=7),
          strip.background = element_blank(),
          #panel.spacing = unit(1, "mm"),
          legend.position = c(1, 0),
          legend.justification = c(1, 0),
          legend.title = element_blank(),
          legend.spacing.y = unit(-0.78, "cm")
          ) 

    gg <- ggplotGrob(g) #gg$layout$name
    ggx <- gtable_filter_remove(gg, name = paste0("axis-b-", c(2, 4), "-9"), trim = FALSE)
    grid.draw(ggx)
    ggsave('Outputs/raw_stringency_loess.png',ggx, width=6,height=7,dpi=600)
   

#ggdraw(g)+draw_label("Observation", x = 0.95, y = 0.01, size =5)
   g = 
    ggplot(s[Nsp>9], aes(x = StringencyIndex, y = FID, group = Species)) +
      #stat_smooth(se = FALSE, col = 'black', lwd = 0.5)+
      stat_smooth(method = 'rlm', se = FALSE, col = 'black', lwd = 0.5)+
      #scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
      scale_x_continuous("COVID measure stringency index", expand = c(0, 0)) +
      scale_y_continuous("Flight initiation distance [m]", expand = c(0, 0), trans = 'log10') +
      #labs(title = "Species means per sampling location")+
      theme_MB  +
      theme(
          plot.title = element_text(size=7),
          strip.background = element_blank(),
          #panel.spacing = unit(1, "mm"),
          legend.position = c(1, 0),
          legend.justification = c(1, 0),
          legend.title = element_blank()
          )  

#' not used
  save(d,s, file = 'Data/for_Peto.Rdata')
  # get pictures
    t = ggimage::phylopic_uid('Turdus merula')
    t = ggimage::phylopic_uid(x$Species)

    cat= image_data (t$uid, size = 128)[[1]]


     
      


  test1 <- data.frame("x" = rnorm(4),
                     "y" = rnorm(4),
                     "genus2" = rep(c(1, 2), 2))
  test2 <- data.frame("x" = rnorm(4),
                     "y" = rnorm(4),
                     "genus2" = rep(c(1, 2), 2),
                     "uia" = rep(c(o$uid[1],
                                        o$uid[2]), 2))
  ggplot(test1, aes(x, y)) +
    geom_point() +
    facet_wrap(~ genus2) +
    geom_phylopic(data = test2, aes(image = uia))

     

  pics <- data.frame("FID_avg.0" = 1.5,
                     "FID_avg.1" = 20,
                     "genus2" = o$genus2,
                     "uia" = o$uid)


  ggplot(test, aes(x, y)) +
    geom_point() +
    facet_wrap(~ genus2) +
    geom_phylopic(pics,aes(x = FID_avg.0, y = FID_avg.1, image = uia))+


# end      