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
      require(PerformanceAnalytics)
      require(viridis)
    
    # constants
      round_ = 3 # number of decimal places to round model coefficients
      nsim = 5000 # number of simulations to extract estimates and 95%CrI
      ax_lines = "grey60" # defines color of the axis lines
      colors <- c("#999999", "#E69F00", "#56B4E9") #viridis(3)

    # Customized ggplot theme
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

    # data
      d = fread(here::here('Data/data.txt'))
      d[, genus := sub("_.*", "", Species)]
      d[, sp_loc := paste(Species, IDLocality, sep="_")]
      d[, sp_country := paste(Species, Country, sep="\n")]
      d[, rad:=(2*pi*Hour) / 24]
      d[, Day_:= Day]
      d[, FID_z := scale(FID), by = Species]
      d[, SD_z := scale(SD), by = Species]
      d[, FID_ln := log(FID)]
      d[, SD_ln := log(SD)]
      d[, body_ln := log(BodyMass)]
      d[, flock_ln := log(FlockSize)]
      #d[Country == 'Australia', Day_:= abs(Day - 189)]
    
      s = fread(here::here('Data/datatop10.txt'))
      s[, genus := sub("_.*", "", Species)]
      s[, sp_loc := paste(Species, IDLocality, sep="_")]
      s[, sp_country := paste(Species, Country, sep="\n")]
      
      ds = d[Species%in%unique(s$Species)]
      #d[Country == 'Australia', Day_:= abs(Day - 189)]
#' ### exploration
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
    m=lmer(FID~scale(Covid)+
                scale(log(SD))+
                scale(log(FlockSize))+
                scale(log(BodyMass))+
                scale(sin(rad)) + scale(cos(rad)) + 
                scale(Temp)+
                (Covid|genus)+(Covid|Species)+(Covid|Country) + (Covid|IDLocality) +(Covid|sp_loc),
                data = d)

    
    summary(glht(m))
    
    plot(allEffects(m))
  
    print(summary(m))

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