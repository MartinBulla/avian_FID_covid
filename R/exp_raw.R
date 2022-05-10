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
    require(ggplot2)
    require(MASS)
    require(multcomp)
    require(optimx)
    require(PerformanceAnalytics)
    require(viridis)
  
    # constants
    round_ = 3 # number of decimal places to round model coefficients
    nsim = 5000 # number of simulations to extract estimates and 95%CrI
    ax_lines = "grey60" # defines color of the axis lines
    colors <- c("#999999", "#E69F00", "#56B4E9") #viridis(3)

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
    d = fread(here::here('Data/data.txt'))
    d = d[order(Year, IDLocality, Day, Hour)]
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
    dd = merge(d1,d2)  
    da = merge(d1,d2, all = TRUE)

    p1 = d[Covid == 1, .N, by = .(IDLocality, Species)]
    p2 = d[Covid == 0, .N, by = .(IDLocality, Species)]
    setnames(p1, old = 'N', new ='N_during')
    setnames(p2, old = 'N', new ='N_before')
    pp = merge(p1,p2)  
    pa = merge(p1,p2, all = TRUE)
    
#' ### Models
#' To control for pseudo-replication we fitted year, genus, species, country, species at specific location and species at a given day of a given year as random intercepts and period (before lockdown and during lockdown) as random slope to all random intercepts. Given the singular fit, we removed random slope from intercepts where correlations were -1 and kept the model with random slope for country and species at sampling location.
  
# CHECK SINGULARITIES and ADD average model a model with 1 sample per species and locality
  # prepare estimates 
    # all data, all random slopes - singularity 
      mf_max=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) + (Covid|genus)+(Covid|Species)+(Covid|sp_day_year) + (Covid|Country) + (Covid|sp_loc),
                  data = d)# # (Covid|IDLocality) +
      est_mf_max = est_out(mf_max, '01) (Year) + (Covid|genus) + (Covid|Species) + (CovidCovid|sp_day_year) + (Covid|Country) + (Covid|sp_loc)')  
    # all data, all random slopes, but some without cor to avoid singularity
      mf_max_=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) +(0+Covid|genus)+(0+Covid|Species)+(0+Covid|sp_day_year) + (Covid|Country) + (0+Covid|IDLocality) +(Covid|sp_loc),
                  data = d, REML = FALSE) # (Covid|IDLocality) +
      est_mf_max_ = est_out(mf_max, '02) (1|Year) + (0+Covid|genus) + (0+Covid|Species) + (0+Covid|sp_day_year) + (Covid|Country) + (0+Covid|IDLocality) + (Covid|sp_loc)')   
    # all data, random slopes that allow for non-singular fit 
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
      est_mf = est_out(mf, '03) (1|Year) +(1|genus) + (1|Species) + (1|sp_day_year) + (Covid|Country) + (Covid|sp_loc)')
    
    # >9 per species
      d[, Nsp := .N, Species]
      mf_10=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) + (1|genus)+(1|Species)+(1|sp_day_year) + (Covid|Country) + (Covid|sp_loc),
                  data = d[Nsp>9], REML = FALSE) # (Covid|IDLocality) +
      est_mf_10 = est_out(mf_10, '04) (1|Year) + (1|genus) + (1|Species) + (1|sp_day_year) + (Covid|Country) + (Covid|sp_loc); >9/species')  
    
    # before & during > 4/species, random slopes that allow for non-singular fit 
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
      est_mf_5ba = est_out(mf_5ba, '05) (1|Year) +(1|genus) + (1|Species) + (1|sp_day_year) + (Covid|Country) + (Covid|sp_loc); >4/species/period')  
    # before & during > 9/species, random slopes that allow for non-singular fit 
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
      est_mf_9ba = est_out(mf_9ba, '06) (1|Year) +  (1|Species) + (1|sp_day_year) + (Covid|Country) + (Covid|sp_loc); >9/species/period')  
    # no Poland before & during > 4/species, random structure that allow for non-singular fit 
      dx = dd[N_during>4 & N_before >4]
      dxx = d[Species %in% dx$Species & Country!='Poland']
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
            (1|Year) + (1|genus) + (1|Species)+ (Covid|Country) + (Covid|sp_loc),
            data = d[Species %in% dx$Species & Country!='Poland'],
            REML = FALSE, control = lmerControl(
                  optimizer ='optimx', optCtrl=list(method='nlminb'))) 
                  # (Covid|IDLocality) +
      est_mf_5baP = est_out(mf_5baP, '07)  (1|Year) + (1|genus) + (1|Species) + (Covid|Country) + (Covid|sp_loc); >4/species/period without PL')  
    # no Poland before & during > 4/species/locality/covid, random structure that allow for non-singular fit
      dx = pp[N_during>4 & N_before >4]
      dxx = d[paste(IDLocality, Species) %in% paste(dx$IDLocality, dx$Species) & Country!='Poland']
      #dxx = d[Species %in% dx$Species]
      #dx2 = dd[N_during>9 & N_before >9]
      mf_5ba_loc_P=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) + (1|Species)+ (1|sp_day_year) +(0+Covid|Country) + (Covid|IDLocality) + (0+Covid|sp_loc),
                  data = dxx,
                  REML = FALSE, control = lmerControl(
                  optimizer ='optimx', optCtrl=list(method='nlminb'))) 
                  # (1|genus) explained 0 so taken out, but for reporting could stay
                  # (Covid|IDLocality) +
      est_mf_5ba_loc_P = est_out(mf_5ba_loc_P, '08) (1|Year) + (1|genus)+(1|Species) + (1|sp_day_year)+ (0+Covid|Country) + (Covid|IDLocality) + (0+Covid|sp_loc); >4/species/locality/period without PL')  
    # diff model structures for previous
      dx = pp[N_during>4 & N_before >4]
      dxx = d[paste(IDLocality, Species) %in% paste(dx$IDLocality, dx$Species) & Country!='Poland']
      #dxx = d[Species %in% dx$Species]
      #dx2 = dd[N_during>9 & N_before >9]
      mf_5ba_loc_P_2=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) + (1|Species)+ (1|sp_day_year) +(0+Covid|Country) + (Covid|IDLocality) + (1|sp_loc),
                  data = dxx,
                  REML = FALSE,
                  control = lmerControl(
                  optimizer ='optimx', optCtrl=list(method='nlminb'))) # (Covid|IDLocality) +
      est_mf_5ba_loc_P_2 = est_out(mf_5ba_loc_P_2, '09) (1|Year) + (1|Species) + (1|sp_day_year) + (0+Covid|Country) + (Covid|IDLocality) +(1|sp_loc); >4/species/locality/period without PL')  
    # year diff structures Poland before & during > 4/species/locality, random structure that allow for non-singular fit
      dx = pp[N_during>4 & N_before >4]
      dxx = d[paste(IDLocality, Species) %in% paste(dx$IDLocality, dx$Species) & Country!='Poland']
      #dxx = d[Species %in% dx$Species]
      #dx2 = dd[N_during>9 & N_before >9]
      mf_5ba_loc_P_2_year=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year)+  (1|genus)+(1|Species)+ (1|sp_day_year) +(0+Covid|Country) + (Covid|IDLocality) + (0+Covid|sp_loc),
                  data = dxx,
                  REML = FALSE,
                  control = lmerControl(
                  optimizer ='optimx', optCtrl=list(method='nlminb'))) # (Covid|IDLocality) +
      est_mf_5ba_loc_P_2_year = est_out(mf_5ba_loc_P_2_year, '10) (1|Year) + (1|genus) + (1|Species) + (1|sp_day_year)+(0+Covid|Country) + (Covid|IDLocality) +(0+Covid|sp_loc); >4/species/locality/period without PL')  
    # year diff structures Poland before & during > 4/species/locality, random structure that allow for non-singular fit; without 2014
      dx = pp[N_during>4 & N_before >4]
      dxx = d[paste(IDLocality, Species) %in% paste(dx$IDLocality, dx$Species) & Country!='Poland' & !IDLocality%in%d[Year == 2014, unique(IDLocality)]]
      #dxx = d[Species %in% dx$Species]
      #dx2 = dd[N_during>9 & N_before >9]
      mf_5ba_loc_P_2_year14=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year)+  (1|Species)+ (1|sp_day_year) +(0+Covid|Country) + (Covid|IDLocality) + (0+Covid|sp_loc),
                  data = dxx,
                  REML = FALSE) # (Covid|IDLocality) +
      est_mf_5ba_loc_P_2_year14 = est_out(mf_5ba_loc_P_2_year14, '11) (1|Year) + (1|Species) + (0+Covid|Country) + (Covid|IDLocality) +(0+Covid|sp_loc); >4/species/locality/period without PL & 2014')  
    
    # 1 obs/species/locality/covid
      dxx = d[paste(IDLocality, Species) %in% paste(pp$IDLocality, pp$Species) & Country!='Poland']
      dxx[, pk:=1:nrow(dxx)]
      #length(dxx[, unique(sp_loc)])
      #length(dxx[, unique(paste(IDLocality, Species, Covid))])
      set.seed(42)
      pkk = dxx[, sample(pk, size = 1), by = .(IDLocality, Species, Covid)]
      d11 = dxx[pk %in% pkk$V1]

      m11=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year)+(0+Covid|Species)+ (1|sp_day_year) +(0+Covid|Country) + (0+Covid|IDLocality) + (1|sp_loc),
                  data = d11,
                  REML = FALSE) # (Covid|IDLocality) +
      est_m11 = est_out(m11, '12) (1|Year)+(0+Covid|Species)+ (1|sp_day_year) +(0+Covid|Country) + (0+Covid|IDLocality) + (1|sp_loc); 1/species/locality/period without PL')  
    # 1 obs/species/locality/covid from localities with >4 before/after
      dx = pp[N_during>4 & N_before >4]
      dxx = d[paste(IDLocality, Species) %in% paste(pp$IDLocality, pp$Species) & Country!='Poland']
      dxx[, pk:=1:nrow(dxx)]
     
      dxx=dxx[paste(IDLocality, Species) %in% paste(dx$IDLocality, dx$Species) ]
      #length(dxx[, unique(paste(IDLocality, Species, Covid))])
      #table(dxx$sp_loc, dxx$Covid)
      set.seed(42)
      pkk = dxx[, sample(pk, size = 1), by = .(IDLocality, Species, Covid)]
      d11 = dxx[pk %in% pkk$V1]

     m13=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                 (1|Species)+ (0+Covid|Country) + (0+Covid|IDLocality) + (1|sp_loc),
                  data = d11,
                  REML = FALSE)
      # singular fit but just because some random effects estimated as null
      est_m13 = est_out(m13, '13) 1|Species) +  (0+Covid|Country) + (0+Covid|IDLocality) + (1|sp_loc); 1/species/locality/period without PL & before/during >4')  
    # 1 obs/species/locality/covid from localities with >4 before/after dif random structure
      dx = pp[N_during>4 & N_before >4]
      dxx = d[paste(IDLocality, Species) %in% paste(pp$IDLocality, pp$Species) & Country!='Poland']
      dxx[, pk:=1:nrow(dxx)]
      
      dxx=dxx[paste(IDLocality, Species) %in% paste(dx$IDLocality, dx$Species) ]
      #length(dxx[, unique(paste(IDLocality, Species, Covid))])
      #table(dxx$sp_loc, dxx$Covid)
      set.seed(42)
      pkk = dxx[, sample(pk, size = 1), by = .(IDLocality, Species, Covid)]
      d11 = dxx[pk %in% pkk$V1]

      m14=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Species), # similar if sp_loc useed
                  data = d11,
                  REML = FALSE) # (Covid|IDLocality) +
      est_m14 = est_out(m14, '14) (1|Species); 1/species/locality/period without PL & before/during >4')  
    
    # avg obs/species/locality/covid from localities
      # 1 obs/species/locality/covid
      dxx = d[paste(IDLocality, Species) %in% paste(pp$IDLocality, pp$Species) & Country!='Poland']
      length(dxx[, unique(paste(sp_loc, Covid))])
      m = lm(log(FID) ~ log(SD),dxx)
      dxx[, resid_FID := resid(m)]
      a = dxx[, mean(resid_FID) ,by = .(Country, IDLocality, genus, Species, sp_loc, Covid)]
      setnames(a, old = 'V1', new = 'resid_FID_avg')
      #ggplot(a, aes(x = resid_FID_avg)) +geom_histogram()
     
      m15=lmer(scale(resid_FID_avg)~
                  scale(Covid)+
                  (1|genus) + (Covid|Species)+ (Covid|Country) + (Covid|IDLocality) + (1|sp_loc),
                  data = a,
                  REML = FALSE) # (Covid|IDLocality) +
      est_m15 = est_out(m15, '15) (1|genus) + (Covid|Species)+ (Covid|Country) + (Covid|IDLocality) + (1|sp_loc); 1 average/species/locality/period without PL; on means from resid_FID')  
    
      
      
  # plot
      x = rbind(est_mf_max, est_mf_max_, est_mf, est_mf_10, est_mf_5ba,est_mf_9ba,est_mf_5baP,est_mf_5ba_loc_P, est_mf_5ba_loc_P_2,est_mf_5ba_loc_P_2_year,est_mf_5ba_loc_P_2_year14,est_m11,est_m13,est_m14,est_m15)

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
      ggsave(here::here('Outputs/effect_sizes.png'),g, width = 20, height =10, units = 'cm')



    acf(resid(m11), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))


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
     dxx = d[paste(IDLocality, Species) %in% paste(dx$IDLocality, dx$Species) & Country!='Poland']
    table(dxx$IDLocality, dxx$Year)
    ggplot(dxx, aes(x = as.factor(Year), y = FID, col = Year)) + geom_boxplot() + facet_wrap(~sp_loc) + scale_y_continuous(trans = 'log10')
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

#' ### MODELS   
  # FID  
    m0=lmer(FID~
                scale(log(SD))+
                scale(log(FlockSize))+
                scale(log(BodyMass))+
                scale(sin(rad)) + scale(cos(rad)) + 
                scale(Temp)+
                scale(Covid)+
                (Covid|genus)+(Covid|Species)+(Covid|sp_day_year) + (Covid|Country) + (Covid|sp_loc),
                data = d) # (Covid|IDLocality) +

    d[, sp_day_year_covid := paste(sp_day_year, Covid)]
    m=lmer(FID~
                scale(log(SD))+
                scale(log(FlockSize))+
                scale(log(BodyMass))+
                scale(sin(rad)) + scale(cos(rad)) + 
                scale(Temp)+
                scale(Covid)+
                (Covid|genus)+(Covid|Species)+(1|sp_day_year_covid) + (Covid|Country) + (Covid|sp_loc),
                data = d)

    m=lmer(FID~
                scale(log(SD))+
                scale(log(FlockSize))+
                scale(log(BodyMass))+
                scale(sin(rad)) + scale(cos(rad)) + 
                scale(Temp)+
                scale(Covid)+
                (Covid|genus)+(Covid|Species)+(1|sp_day_year_covid) + (Covid|Country) + (Covid|sp_loc),
                data = d[Species%in%dd$Species])
    summary(glht(m))
    
    plot(allEffects(m))
  
    print(summary(m))

     acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))

#'***
  # log FID
    ml=lmer(log(FID)~scale(Covid)+
                scale(log(SD))+
                scale(log(FlockSize))+
                scale(log(BodyMass))+
                scale(sin(rad)) + scale(cos(rad)) + 
                scale(Temp)+
                (Covid|genus)+(Covid|Species)+(Covid|Country) + (Covid|IDLocality) +(Covid|sp_loc),
                data = d)

    summary(glht(ml))
    
    plot(allEffects(ml))

    print(summary(ml))
#'***
  # scaled within species
    mz=lmer(FID_z~scale(Covid)+
                scale((SD_z))+
                scale(log(FlockSize))+
                scale(log(BodyMass))+
                scale(sin(rad)) + scale(cos(rad)) + 
                scale(Temp)+
                (Covid|genus)+(Covid|Species)+(Covid|Country) + (Covid|IDLocality) +(Covid|sp_loc),
                data = d)

    summary(glht(mz))
    
    plot(allEffects(mz))

    print(summary(mz))  

#' ### test randome slope structures
    mt=lmer(FID~
      scale(Covid)+
      scale(log(SD))+
      scale(log(FlockSize))+
      scale(log(BodyMass))+
      scale(sin(rad)) + scale(cos(rad)) + 
      scale(Temp)+
      (Covid|genus)+(1|Species)+(1|Country) + (1|IDLocality) +(1|sp_loc),
      data = d)

    
    summary(glht(mt))
    print(summary(mt))

#' species specific for 10 most common species
  l = list()
  lp = list()
  for(i in unique(ds$Species)){
      #i = "Passer_domesticus"  
      dsi = ds[Species%in%i]
      dsi[, Country := as.factor(Country)]
      m = lmer(scale(FID)~
        scale(Covid)+
        scale(log(SD))+
        scale(log(FlockSize))+
        scale(sin(rad)) + scale(cos(rad)) + 
        scale(Temp)+
        (1|Country) + (1|IDLocality),
        data = dsi)

      bsim = sim(m, n.sim=2000) 
      v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
      ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))  
      l[[i]]=data.frame(species=i,effect=rownames(coef(summary(m)))[2],estimate=v[2], lwr=ci[1,2], upr=ci[2,2])
      print(i) 
  }
  ll = data.table(do.call(rbind,l) ) 
  
  g = 
  ggplot(ll, aes(y = effect, x = estimate, col = species)) +
    geom_vline(xintercept = 0, col = "grey30", lty =3)+
    geom_errorbar(aes(xmin = lwr, xmax = upr, col = species), width = 0.1, position = position_dodge(width = 0.4) ) +
    ggtitle ("Including Poland")+
    geom_point(position = position_dodge(width = 0.4)) +
    #scale_colour_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
    #scale_fill_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
    scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
    scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE)) + 
    #scale_shape(guide = guide_legend(reverse = TRUE)) + 
    #scale_x_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2,2, by = 1), labels = seq(-2,2, by = 1)) +
    labs(y = NULL ,x = "Standardized effect size",title = 'Species-specific models') +
    #ylim(c(0,100))+
    #coord_flip()+
    theme_bw() +
    theme( #legend.position ="right",
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
          axis.text.x = element_text(colour="black", size = 7),
          axis.text.y=element_text(colour="black", size = 7),
          axis.title=element_text(size=9)
          )
  g
#' species specific for 10 most common species, without POLAND
  l = list()
  lp = list()
  x = unique(ds$Species)
  x = x[x!="Larus_ridibundus"]
  for(i in x){
      #i = "Larus_ridibundus"  
      dsi = ds[!Country %in% 'Poland' & Species%in%i]
      m = lmer(scale(FID)~
        scale(Covid)+
        scale(log(SD))+
        scale(log(FlockSize))+
        scale(sin(rad)) + scale(cos(rad)) + 
        scale(Temp)+
        (1|Country) + (1|IDLocality),
        data = dsi)

      bsim = sim(m, n.sim=2000) 
      v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
      ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))  
      l[[i]]=data.frame(species=i,effect=rownames(coef(summary(m)))[2],estimate=v[2], lwr=ci[1,2], upr=ci[2,2])
      print(i) 
  }
  ll = data.table(do.call(rbind,l) ) 
  
  g = 
  ggplot(ll, aes(y = effect, x = estimate, col = species)) +
    geom_vline(xintercept = 0, col = "grey30", lty =3)+
    geom_errorbar(aes(xmin = lwr, xmax = upr, col = species), width = 0.1, position = position_dodge(width = 0.4) ) +
     ggtitle ("Without Poland")+
    geom_point(position = position_dodge(width = 0.4)) +
    #scale_colour_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
    #scale_fill_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
    scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
    scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE)) + 
    #scale_shape(guide = guide_legend(reverse = TRUE)) + 
    #scale_x_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2,2, by = 1), labels = seq(-2,2, by = 1)) +
    labs(y = NULL ,x = "Standardized effect size",title = 'Species-specific models') +
    #ylim(c(0,100))+
    #coord_flip()+
    theme_bw() +
    theme( #legend.position ="right",
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
          axis.text.x = element_text(colour="black", size = 7),
          axis.text.y=element_text(colour="black", size = 7),
          axis.title=element_text(size=9)
          )
  g
#' species specific for 10 most common species
  summary(factor(d$Species))
  db = d[Covid %in% 0, .N, by=Species]
  setnames(db, old = 'N', new = 'N_before')

  da=d[Covid %in% 1, .N, by=Species]
  setnames(da, old = 'N', new = 'N_after')

  dab = merge(da,db, all.x = TRUE)
  dab = dab[!is.na(N_before)]

  dab_ = dab[N_before>10 & N_after>10]

  d_= d[Species%in%unique(dab_$Species)]

  d_c = d_[, length(unique(Country)), by = .(Species)]

  dab__ = dab_[Species%in%d_c[!V1%in%c(1,2), Species]]
  
  l = list()
  lp = list()
  for(i in unique(dab__$Species)){
      #i = "Passer_domesticus"  
      dsi = d[Species%in%i]
      dsi[, Country := as.factor(Country)]
      m = lmer(scale(FID)~
        scale(Covid)+
        scale(log(SD))+
        scale(log(FlockSize))+
        scale(sin(rad)) + scale(cos(rad)) + 
        scale(Temp)+
        (1|Country) + (1|IDLocality),
        data = dsi)

      bsim = sim(m, n.sim=2000) 
      v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
      ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))  
      l[[i]]=data.frame(species=i,effect=rownames(coef(summary(m)))[2],estimate=v[2], lwr=ci[1,2], upr=ci[2,2])
      print(i) 
  }
  ll = data.table(do.call(rbind,l) ) 
  
  g = 
  ggplot(ll, aes(y = effect, x = estimate, col = species)) +
    geom_vline(xintercept = 0, col = "grey30", lty =3)+
    geom_errorbar(aes(xmin = lwr, xmax = upr, col = species), width = 0.1, position = position_dodge(width = 0.4) ) +
    #ggtitle ("Including Poland & species with >10 for before & during from more than one country")+
    geom_point(position = position_dodge(width = 0.4)) +
    #scale_colour_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
    #scale_fill_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
    scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
    scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE)) + 
    #scale_shape(guide = guide_legend(reverse = TRUE)) + 
    #scale_x_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2,2, by = 1), labels = seq(-2,2, by = 1)) +
    labs(y = NULL ,x = "Standardized effect size",title = 'Species-specific models\nData including Poland & species with >10 for before & during from more than one country') +
    #ylim(c(0,100))+
    #coord_flip()+
    theme_bw() +
    theme( #legend.position ="right",
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
          axis.text.x = element_text(colour="black", size = 7),
          axis.text.y=element_text(colour="black", size = 7),
          axis.title=element_text(size=9)
          )
  g

#' 
# end      