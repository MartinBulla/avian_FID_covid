#' ---
#' title: "Some plots for avian FIDs before & during COVID"
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
    require(data.table)
    require(ggplot2)
    require(MASS)
    require(arm)
    require(multcomp)
    require(effects)

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
    d = fread('Data/datafull.txt')
    setnames(d,old = c('mean(FID)'), new = c('mean_FID'))
    d[Covid == 0, when:='before']
    d[Covid == 1, when:='during']
    d$Covid = NULL
    dw = reshape(d, idvar = c('Species','Country','IDLocality'), timevar = 'when', direction = "wide")  
    setnames(dw,old = c('mean_FID.before', 'n.before','mean_FID.during','n.during'), new = c('FID_before', 'n_before','FID_during', 'n_during'))

    dw[, sp_loc := paste(Species, IDLocality, sep="_")]
    dw[, sp_country := paste(Species, Country, sep="\n")]

#' ### PLOTS for paired dataset
#+ fig.width=6,fig.height=6
ggplot(dw, aes(x = FID_before, y = FID_during)) + 
    geom_point(pch = 21, alpha = 0.7) + 
    geom_abline(intercept = 0, slope = 1, lty =3, col = "grey80")+
    facet_wrap(~sp_country) +
    #scale_color_manual(values=c(male, female))+ 
    scale_x_continuous("FID before COVID [m]", expand = c(0, 0)) +
    scale_y_continuous("FID during COVID [m]]", expand = c(0, 0)) +
    theme_MB

#+ fig.width=5,fig.height=5
ggplot(dw, aes(x = FID_before, y = FID_during, col = Country)) + 
    geom_point(pch = 21, alpha = 0.7) + 
    geom_abline(intercept = 0, slope = 1, lty =3, col = "grey80")+
    facet_wrap(~Species) +
    #scale_color_manual(values=c(male, female))+ 
    scale_x_continuous("FID before COVID [m]", expand = c(0, 0)) +
    scale_y_continuous("FID during COVID [m]]", expand = c(0, 0)) +
    theme_MB    

#+ fig.width=5,fig.height=5
ggplot(dw, aes(x = FID_before, y = FID_during)) + 
    geom_point(pch = 21, alpha = 0.7, aes(col = Country)) + 
    geom_smooth(method = 'lm', se=FALSE, aes(x = FID_before, y = FID_during), col = 'black')+
    geom_abline(intercept = 0, slope = 1, lty =3, col = "grey80")+
    facet_wrap(~Species) +
    #scale_color_manual(values=c(male, female))+ 
    scale_x_continuous("FID before COVID [m]", expand = c(0, 0)) +
    scale_y_continuous("FID during COVID [m]]", expand = c(0, 0)) +
    labs(title = "Fit for species")+
    theme_MB      

#+ fig.width=5,fig.height=5
ggplot(dw, aes(x = FID_before, y = FID_during)) + 
    geom_point(pch = 21, alpha = 0.7, aes(col = Country)) + 
    geom_smooth(method = 'rlm', se=FALSE, aes(x = FID_before, y = FID_during), col = 'black')+
    geom_abline(intercept = 0, slope = 1, lty =3, col = "grey80")+
    facet_wrap(~Species) +
    #scale_color_manual(values=c(male, female))+ 
    scale_x_continuous("FID before COVID [log10(m)]", expand = c(0, 0), trans = 'log10') +
    scale_y_continuous("FID during COVID [log10(m)]]", expand = c(0, 0), trans = 'log10') +
    labs(title = "Fit for species")+
    theme_MB    


#+ fig.width=6,fig.height=3
 ggplot(dw, aes(x = FID_before, y = FID_during)) + 
    geom_point(pch = 21, alpha = 0.7, aes(col = Species)) + 
    geom_smooth(method = 'lm', se=FALSE, aes(x = FID_before, y = FID_during), col = 'black')+
    geom_abline(intercept = 0, slope = 1, lty =3, col = "grey80")+
    facet_wrap(~Country) +
    #scale_color_manual(values=c(male, female))+ 
    scale_x_continuous("FID before COVID [m]", expand = c(0, 0)) +
    scale_y_continuous("FID during COVID [m]]", expand = c(0, 0)) +
    labs(title = "Fit for country")+
    theme_MB      
#+ fig.width=6,fig.height=3
 ggplot(dw, aes(x = FID_before, y = FID_during)) + 
    geom_point(pch = 21, alpha = 0.7, aes(col = Species)) + 
    geom_smooth(method = 'rlm', se=FALSE, aes(x = FID_before, y = FID_during), col = 'black')+
    geom_abline(intercept = 0, slope = 1, lty =3, col = "grey80")+
    facet_wrap(~Country) +
    #scale_color_manual(values=c(male, female))+ 
    scale_x_continuous("FID before COVID [log10(m)]", expand = c(0, 0), trans = 'log10') +
    scale_y_continuous("FID during COVID [log10(m)]]", expand = c(0, 0), trans = 'log10') +
    labs(title = "Fit for country")+
    theme_MB  

#' Efekt covidu uplne zmizi ked kontrolujes na ranom slope len u country a je slaby ked len u species a zadny ked len u genus... cize daco tam v tych random slopoch je.
#' ### MODELS
    dd = d
    dd[when == 'before', when01:= 0]
    dd[when == 'during', when01:= 1]
    dd[, sp_loc := paste(Species, IDLocality, sep="_")]
    dd[, sp_country := paste(Species, Country, sep="\n")]
    
    m<-lmer(mean_FID~scale(when01)+
               scale(log10(n)) + 
               (when01|Species)+(when01|Country) + (when01|IDLocality) +(1|sp_loc),
             data = dd)

    
    summary(glht(m))
    plot(allEffects(m))
  
    print(summary(m))

#'***

    m<-lmer(log(mean_FID)~scale(when01)+
               scale(log10(n)) + 
               (when01|Species)+(when01|Country) + (when01|IDLocality) +(1|sp_loc),
             data = dd)

    summary(glht(m))
    plot(allEffects(m))

    print(summary(m))
# end      