START WITH FIGURE SY and finish with C
Make tables and check model ass

#' ##### Code to load tools & data
  # packages
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
    require(png)
    require(rphylopic)
    require(viridis)
  # constants
    round_ = 3 # number of decimal places to round model coefficients
    nsim = 5000 # number of simulations to extract estimates and 95%CrI
    ax_lines = "grey60" # defines color of the axis lines
    colors <- c("#999999", "#E69F00", "#56B4E9") #viridis(3)
  # functions
    # to add images to panels
      annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) {
        layer(data = data, stat = StatIdentity, position = PositionIdentity, 
            geom = ggplot2:::GeomCustomAnn,
            inherit.aes = TRUE, params = list(grob = grob, 
                                              xmin = xmin, xmax = xmax, 
                                              ymin = ymin, ymax = ymax))
      }
    # to remove ggplot components
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
    # for estimates
        est_out =function(model = m, label = "", nsim = 5000){
            bsim = sim(model, n.sim=5000) 
            v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
            data.table(predictor=rownames(coef(summary(model))),estimate=v, lwr=ci[1,], upr=ci[2,], model = paste(label, "N =", nobs(model)))
          }
    # change color
      change_col = function(replace_black, theimg) {
        r_b = col2rgb(replace_black) / 255
        #theimg[theimg == 1] <- 2
        for (i in 1:3) {
            theimg[,,i][theimg[,,i] == 0] <- r_b[i]
        }
        return(theimg)
        }
    # for Supplementary Table output based on sim
        R_out = function(name = "define", model = m, nsim = 5000){
         bsim <- sim(model, n.sim=nsim)  
         l=data.frame(summary(model)$varcor)
         l = l[is.na(l$var2),]
         l$var1 = ifelse(is.na(l$var1),"",l$var1)
         l$pred = paste(l$grp,l$var1)

         q50={}
         q025={}
         q975={}
         pred={}
         
         # variance of random effects
         for (ran in names(bsim@ranef)) {
           #ran =names(bsim@ranef)[1]
           ran_type = l$var1[l$grp == ran]
           for(i in ran_type){
              # i = ran_type[2]
            q50=c(q50,quantile(apply(bsim@ranef[[ran]][,,i], 1, var), prob=c(0.5)))
            q025=c(q025,quantile(apply(bsim@ranef[[ran]][,,i], 1, var), prob=c(0.025)))
            q975=c(q975,quantile(apply(bsim@ranef[[ran]][,,i], 1, var), prob=c(0.975)))
            pred= c(pred,paste(ran, i))
            }
           }
         # residual variance
         q50=c(q50,quantile(bsim@sigma^2, prob=c(0.5)))
         q025=c(q025,quantile(bsim@sigma^2, prob=c(0.025)))
         q975=c(q975,quantile(bsim@sigma^2, prob=c(0.975)))
         pred= c(pred,'Residual')

         ci = c(round(100*q025/sum(q025))[1], round(100*q975/sum(q975))[1])
         ci = ci[order(ci)]
         
         ri=data.table(model = name, repeatability=paste0(round(100*q50/sum(q50)),'%')[1], CI = paste0(paste(ci[1], ci[2], sep ="-"), '%'))
         
         
         return(ri)
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

    s = d[Covid == 1]
    s[, Nsp := .N, by ='Species']
    s[, sp := gsub('[_]', ' ', Species)]
#' ## Figures
# Figure A
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

          aw[genus2=='Phoenicurus', unique(Species)]

        o[genus2=='Motacilla' | uid%in%c('67a9ecfd-58ba-44a4-9986-243b6e610419'), uid:='cf522e02-35cc-44f5-841c-0e642987c2e4']
        o[genus2=='Sylvia', uid:='67a9ecfd-58ba-44a4-9986-243b6e610419']

        o[,size:=0.2]
        o[ genus2%in%c('Anas','Columba','Dendrocopos','Sturnus'), size := c(0.25, 0.25, 0.15, 0.1)]
        o[, FID_avg.0 := 1.5]
        o[, FID_avg.1 := 20]
        o[ genus2%in%c('Anas','Columba'), FID_avg.0 := c(1.7, 1.7)]
        
        o[, resid_FID_avg.0 := -1.7]
        o[, resid_FID_avg.1 := 0.7]
        o[ genus2%in%c('Anas','Columba'), resid_FID_avg.0 := c(-1.6, -1.6)]

        o[, genus2 := factor(genus2, levels = c('Anas', 'Larus', 'Columba', 'Dendrocopos','Picus', 'Motacilla','Erithacus','Phoenicurus','Turdus', 'Sylvia','Parus','Sitta','Pica','Garrulus','Corvus','Sturnus','Passer','Fringilla','other'))]
        
        aw[, genus2 := factor(genus2, levels = c('Anas', 'Larus', 'Columba', 'Dendrocopos', 'Picus', 'Motacilla','Erithacus','Phoenicurus','Turdus', 'Sylvia','Parus','Sitta','Pica','Garrulus','Corvus','Sturnus','Passer','Fringilla','other'))]

        aw[, genus2 := factor(genus2, levels = c('Anas', 'Larus', 'Columba', 'Dendrocopos', 'Picus', 'Motacilla','Erithacus','Phoenicurus','Turdus', 'Sylvia','Parus','Sitta','Pica','Garrulus','Corvus','Sturnus','Passer','Fringilla','other'))]
    # plot from web
        g = 
        ggplot(aw, aes(x = FID_avg.0, y = FID_avg.1)) + 
          #geom_errorbar(aes(ymin = FID_avg.1-SD.1, ymax = FID_avg.1+SD.1, col = Country), width = 0) +
          #geom_errorbar(aes(xmin = FID_avg.0-SD.0, xmax = FID_avg.0+SD.0, col = Country), width = 0) +
          #geom_point(pch = 21, alpha = 0.7, aes(col = Country)) + 
          geom_point(pch = 21, alpha = 0.7, aes(fill = Country), col = 'white') + 
            #ggtitle ("Sim based")+
          geom_abline(intercept = 0, slope = 1, lty =3, col = "grey80")+
          facet_wrap(~genus2) +
          geom_phylopic(data = o, aes(image = uid),  color = "grey80", size = o$size) + # ,
          scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = FALSE))  +
          scale_x_continuous("Before shutdown - flight initiation distance [m]", expand = c(0, 0), trans = 'log10') +
          scale_y_continuous("During shutdown - flight initiation distance [m]", expand = c(0, 0), trans = 'log10') +
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
        ggsave('Outputs/before_after_Genus_fill.png',ggx, width=4.5,height=4,dpi=600)
    # USE plot from files
        columba = rasterGrob(change_col('#CCCCCC',readPNG('Data/PhyloPic/Columba.png')))
        Dendrocopos = rasterGrob(change_col('#CCCCCC',readPNG('Data/PhyloPic/Dendrocopos.png')))
        Larus = rasterGrob(change_col('#CCCCCC',readPNG('Data/PhyloPic/Larus_flip.png')))
        Picus = rasterGrob(change_col('#CCCCCC',readPNG('Data/PhyloPic/Picus.png')))
        Motacilla = rasterGrob(change_col('#CCCCCC',readPNG('Data/PhyloPic/Motacilla.png')))
        Erithacus = rasterGrob(change_col('#CCCCCC',readPNG('Data/PhyloPic/Erithacus.png')))
        Phoenicurus = rasterGrob(change_col('#CCCCCC',readPNG('Data/PhyloPic/Phoenicurus.png')))
        Turdus = rasterGrob(change_col('#CCCCCC',readPNG('Data/PhyloPic/Turdus.png')))
        Sylvia = rasterGrob(change_col('#CCCCCC',readPNG('Data/PhyloPic/Sylvia.png')))
        Parus = rasterGrob(change_col('#CCCCCC',readPNG('Data/PhyloPic/Parus_flip.png')))
        Sitta = rasterGrob(change_col('#CCCCCC',readPNG('Data/PhyloPic/Sitta.png')))
        Pica = rasterGrob(change_col('#CCCCCC',readPNG('Data/PhyloPic/Pica.png')))
        Garrulus = rasterGrob(change_col('#CCCCCC',readPNG('Data/PhyloPic/Garrulus.png')))
        Corvus = rasterGrob(change_col('#CCCCCC',readPNG('Data/PhyloPic/Corvus.png')))
        Sturnus = rasterGrob(change_col('#CCCCCC',readPNG('Data/PhyloPic/Sturnus.png')))
        Passer = rasterGrob(change_col('#CCCCCC',readPNG('Data/PhyloPic/Passer_flip.png')))
        Fringilla = rasterGrob(change_col('#CCCCCC',readPNG('Data/PhyloPic/Fringilla.png')))
        other = rasterGrob(change_col('#CCCCCC',readPNG('Data/PhyloPic/other.png')))

        g = 
        ggplot(aw, aes(x = FID_avg.0, y = FID_avg.1)) + 
          #geom_errorbar(aes(ymin = FID_avg.1-SD.1, ymax = FID_avg.1+SD.1, col = Country), width = 0) +
          #geom_errorbar(aes(xmin = FID_avg.0-SD.0, xmax = FID_avg.0+SD.0, col = Country), width = 0) +
          #geom_point(pch = 21, alpha = 0.7, aes(col = Country)) + 
          annotation_custom2(anas, data=o[genus2 == 'Anas'], xmin = 0.05, xmax =0.5, ymax = 2.6)+
          annotation_custom2(Larus, data=o[genus2 == 'Larus'], xmin = 0.05, xmax =0.5, ymax = 2.6)+
          annotation_custom2(columba, data=o[genus2 == 'Columba'], xmin = 0.05, xmax =0.4, ymax = 2.7)+
          annotation_custom2(Dendrocopos, data=o[genus2 == 'Dendrocopos'], xmin = 0.05, xmax =0.25, ymax = 2.6)+
          annotation_custom2(Picus, data=o[genus2 == 'Picus'], xmin = 0.05, xmax =0.4, ymax = 2.7)+
          annotation_custom2(Motacilla, data=o[genus2 == 'Motacilla'], xmin = 0.05, xmax =0.5, ymax = 2.7)+
          annotation_custom2(Erithacus, data=o[genus2 == 'Erithacus'], xmin = 0.05, xmax =0.35, ymax = 2.7)+
          annotation_custom2(Phoenicurus, data=o[genus2 == 'Phoenicurus'], xmin = 0.05, xmax =0.46, ymax = 2.7)+
          annotation_custom2(Turdus, data=o[genus2 == 'Turdus'], xmin = 0.05, xmax =0.5, ymax = 2.7)+
          annotation_custom2(Sylvia, data=o[genus2 == 'Sylvia'], xmin = 0.05, xmax =0.5, ymax = 2.7)+
          annotation_custom2(Parus, data=o[genus2 == 'Parus'], xmin = 0.05, xmax =0.42, ymax = 2.7)+
          annotation_custom2(Sitta, data=o[genus2 == 'Sitta'], xmin = 0.05, xmax =0.5, ymax = 2.7)+
          annotation_custom2(Pica, data=o[genus2 == 'Pica'], xmin = 0.05, xmax =0.5, ymax = 2.5)+
          annotation_custom2(Garrulus, data=o[genus2 == 'Garrulus'], xmin = 0.05, xmax =0.6, ymax = 2.7)+
          annotation_custom2(Corvus, data=o[genus2 == 'Corvus'], xmin = 0.05, xmax =0.4, ymax = 2.55)+
          annotation_custom2(Sturnus, data=o[genus2 == 'Sturnus'], xmin = 0.05, xmax =0.24, ymax = 2.65)+
          annotation_custom2(Passer, data=o[genus2 == 'Passer'], xmin = 0.05, xmax =0.36, ymax = 2.65)+
          annotation_custom2(Fringilla, data=o[genus2 == 'Fringilla'], xmin = 0.05, xmax =0.5, ymax = 2.7)+
          annotation_custom2(other, data=o[genus2 == 'other'], xmin = 0.05, xmax =0.45, ymax = 2.4)+
          geom_point(pch = 21, alpha = 0.7, aes(fill = Country), col = 'white') + 
            #ggtitle ("Sim based")+
          geom_abline(intercept = 0, slope = 1, lty =3, col = "grey80")+
          facet_wrap(~genus2) +
          #geom_phylopic(data = o, aes(image = uid),  color = "grey80", size = o$size) + # ,
          scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = FALSE))  +
          scale_x_continuous("Before shutdown - flight initiation distance [m]", expand = c(0, 0), trans = 'log10') +
          scale_y_continuous("During shutdown - flight initiation distance [m]", expand = c(0, 0), trans = 'log10') +
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
        ggsave('Outputs/before_after_Genus_fill_files_.png',ggx, width=4.5,height=4,dpi=600)
    # legend
        # correlation between mean fid and residual fid
            ggplot(a, aes(x = resid_FID_avg)) + geom_histogram()
            ggplot(a, aes(x=FID_avg, y = resid_FID_avg)) + 
                stat_smooth() + 
                geom_point() +
                stat_cor(aes(label = ..r.label..),  label.x = 0.3, size = 2) +
                scale_x_continuous(trans = 'log10')
            ggplot(aw, aes(x = N.0-N.1)) +geom_histogram()
            nrow(aw[abs(N.0-N.1)>2])
            nrow(aw[!abs(N.0-N.1)>2])
    # supplement
        g2 =     
        ggplot(aw, aes(x = resid_FID_avg.0, y = resid_FID_avg.1)) + 
          geom_point(pch = 21, alpha = 0.7, aes(fill = Country), col = 'white') + 
          geom_abline(intercept = 0, slope = 1, lty =3, col = "grey80")+
          facet_wrap(~genus2) +
          #geom_phylopic(data = o, aes(image = uid),  color = "grey80", size = o$size) + 
          scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
          scale_x_continuous("Before COVID residual escape distance", expand = c(0, 0)) +
          scale_y_continuous("During COVID residual escape distance", expand = c(0, 0)) +
          labs(title = "Species means per sampling location")+
          theme_MB  +
          theme(
                plot.title = element_text(size=7),
                strip.background = element_blank(),
                legend.position = "none",
                #legend.position = c(1, 0),
                legend.justification = c(1, 0)
                )     

        gg2 <- ggplotGrob(g2) #gg$layout$name
        ggx2 <- gtable_filter_remove(gg2, name = paste0("axis-b-", c(2, 4), "-4"),
                                         trim = FALSE)

        grid.draw(cbind(cbind(ggx,ggx2 ,  size = "last")))
        ggsave('Outputs/Fig_Sw1_2.png',cbind(ggx,ggx2 , size = "last"), width=4.5*2,height=4,dpi=600)

# Figure B
  
  g = 
    ggplot(s[Nsp>9], aes(x = StringencyIndex, y = FID)) +
      stat_smooth(se = FALSE, aes(colour = 'LOESS smooth'), lwd = 0.5)+ # show_guide=TRUE
      #stat_smooth(method = 'rlm', se = FALSE, col = 'black', lwd = 0.5)+
      geom_point(pch = 21, alpha = 0.7, aes(fill = Country), col = 'white') + 
      facet_wrap(~sp, ncol = 6) +
      scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = FALSE))  +
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
          strip.text.x = element_text(size = 4.5, color="grey30",  margin=margin(1,1,1,1,"mm")),
          #panel.spacing = unit(1, "mm"),
          legend.position = c(1, -0.01),
          legend.justification = c(1, 0),
          legend.title = element_blank(),
          #legend.spacing.y = unit(-0.78, "cm")
          legend.spacing.y = unit(0.02, "cm")
          ) 

    gg <- ggplotGrob(g) #gg$layout$name
    ggx <- gtable_filter_remove(gg, name = paste0("axis-b-", c(2, 4), "-9"), trim = FALSE)
    grid.draw(ggx)
    ggsave('Outputs/raw_stringency_loess_2.png',ggx, width=6,height=7,dpi=600)
   
# Figure C


# Figure Sx
   d[, sin_rad:=sin(rad)]
   d[, cos_rad:=cos(rad)]

   dp = d[,c('StringencyIndex','FID_ln','SD_ln', 'flock_ln', 'body_ln', 'sin_rad', 'cos_rad','Temp', 'Day')]
   setnames(dp,old = c('StringencyIndex','FID_ln','SD_ln', 'flock_ln', 'body_ln', 'sin_rad', 'cos_rad','Temp', 'Day'), new = c('Stringency\nindex','Escape distance\nln(m)','Starting distance\nln(m)', 'Flock size\nln (m)', 'Body mass\nln(m)', 'Sinus\n of radians', 'Cosinus\nof radians','Temperature\n°C', 'Day'))
   
   png("Outputs/Fig_Sx.png", width =17, height = 17, units = "cm", bg = "transparent", res = 600)
   chart.Correlation(dp, histogram=TRUE, pch=19, alpha = 0.5)
   mtext("Single observations", side=3, line=3)
   dev.off()

# Figure Sy
  # estimates 
    # 01a all data, all random slopes - singularity 
      m1a=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) + (scale(Covid)|genus)+(scale(Covid)|Species)+(1|sp_day_year) + (scale(Covid)|Country) +(scale(Covid)|IDLocality) + (scale(Covid)|sp_loc),
                  data = d)# # (Covid|IDLocality) +
      est_m1a = est_out(m1a, '01a) (1|Year) + (scale(Covid)|genus) + (scale(Covid)|Species) + (1|sp_day_year) + (scale(Covid)|Country) + (scale(Covid)|IDLocality) + (scale(Covid)|sp_loc)')  
    # 01b all data, all random slopes, but some without cor to avoid singularity
      m1b=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  #(1|Year) +(0+Covid|genus)+(0+Covid|Species)+(1|sp_day_year) + (Covid|Country) + (0+Covid|IDLocality) +(Covid|sp_loc)
                  (1|Year) + (0+scale(Covid)|genus)+(0+scale(Covid)|Species)+(1|sp_day_year) + (scale(Covid)|Country) +(scale(Covid)|IDLocality) + (scale(Covid)|sp_loc),
                  data = d, REML = FALSE) # (0+scale(Covid)|genus) is zero 
      est_m1b = est_out(m1b, '01b) (1|Year) + (0+scale(Covid)|genus) + (0+scale(Covid)|Species) + (1|sp_day_year) + (scale(Covid)|Country) + (scale(Covid)|IDLocality) + (scale(Covid)|sp_loc)')   
    # 01c all data, random slopes that allow for non-singular fit 
      m1c=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) +(1|genus)+(1|Species)+(1|sp_day_year) + (scale(Covid)|Country) + (scale(Covid)|IDLocality) + (scale(Covid)|sp_loc),
                  data = d, REML = FALSE, control = lmerControl(
                           optimizer ='optimx', optCtrl=list(method='nlminb')))# (Covid|IDLocality) +
      #d[,res := resid(mf)]
      est_m1c = est_out(m1c, '01c) (1|Year) + (1|genus) + (1|Species) + (1|sp_day_year) + (scale(Covid)|Country) + (scale(Covid)|IDLocality) + (scale(Covid)|sp_loc)')
    # 01d all data, random slopes that allow for non-singular fit with simple struccture
      m1d=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) +(1|genus)+(1|Species)+(1|sp_day_year) + (scale(Covid)|Country) + (1|IDLocality) + (1|sp_loc),
                  data = d, REML = FALSE, control = lmerControl(
                           optimizer ='optimx', optCtrl=list(method='nlminb')))# (Covid|IDLocality) +
      #d[,res := resid(mf)]
      est_m1d = est_out(m1d, '01d) (1|Year) + (1|genus) + (1|Species) + (1|sp_day_year) + (scale(Covid)|Country) + (1|IDLocality) + (1|sp_loc)')

    # 02a) before & during > 4/species - singularity 
      dx = dd[N_during>4 & N_before >4]
      #dxx = d[Species %in% dx$Species]
      #dx2 = dd[N_during>9 & N_before >9]
      m2a=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) + (scale(Covid)|genus)+(scale(Covid)|Species)+(1|sp_day_year) + (scale(Covid)|Country) +(scale(Covid)|IDLocality) + (scale(Covid)|sp_loc),
                  #(1|Year) + (1|genus)+(1|Species)+(1|sp_day_year) + (Covid|Country) + (Covid|sp_loc),
                  data = d[Species %in% dx$Species],
                  REML = FALSE) # (Covid|IDLocality) +
      est_m2a = est_out(m2a, '02a) (1|Year) + (scale(Covid)|genus)+(scale(Covid)|Species) + (1|sp_day_year) + (scale(Covid)|Country) + (scale(Covid)|IDLocality) + (scale(Covid)|sp_loc); >4/species/period')  
    # 02b) before & during > 4/species, random slopes that allow for non-singular fit and simple structure
      dx = dd[N_during>4 & N_before >4]
      #dxx = d[Species %in% dx$Species]
      #dx2 = dd[N_during>9 & N_before >9]
      m2b=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) +(1|genus)+(1|Species)+(1|sp_day_year) + (scale(Covid)|Country) + (1|IDLocality) + (1|sp_loc),
                  #(1|Year) + (1|genus)+(1|Species)+(1|sp_day_year) + (Covid|Country) + (Covid|sp_loc),
                  data = d[Species %in% dx$Species],
                  REML = FALSE) # (Covid|IDLocality) +
      est_m2b = est_out(m2b, '02b) (1|Year) + (1|genus) + (1|Species)+(1|sp_day_year) + (scale(Covid)|Country) + (1|IDLocality) + (1|sp_loc); >4/species/period')  
       
    # 03a) before & during > 9/species - singularity 
      dx = dd[N_during>9 & N_before >9]
      #dxx = d[Species %in% dx$Species]
      #dx2 = dd[N_during>9 & N_before >9]
      m3a=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) + (scale(Covid)|genus)+(scale(Covid)|Species)+(1|sp_day_year) + (scale(Covid)|Country) +(scale(Covid)|IDLocality) + (scale(Covid)|sp_loc),
                  #(1|Year) + (1|genus)+(1|Species)+(1|sp_day_year) + (Covid|Country) + (Covid|sp_loc),
                  data = d[Species %in% dx$Species],
                  REML = FALSE) # (Covid|IDLocality) +
      est_m3a = est_out(m3a, '03a) (1|Year) + (scale(Covid)|genus)+(scale(Covid)|Species) + (1|sp_day_year) + (scale(Covid)|Country) + (scale(Covid)|IDLocality) + (scale(Covid)|sp_loc); >9/species/period')  
    # 03b) before & during > 9/species, random slopes that allow for non-singular fit and simple structure
      dx = dd[N_during>9 & N_before >9]
      #dxx = d[Species %in% dx$Species]
      #dx2 = dd[N_during>9 & N_before >9]
      m3b=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) +(1|genus)+(1|Species)+(1|sp_day_year) + (scale(Covid)|Country) + (1|IDLocality) + (1|sp_loc),
                  #(1|Year) + (1|genus)+(1|Species)+(1|sp_day_year) + (Covid|Country) + (Covid|sp_loc),
                  data = d[Species %in% dx$Species],
                  REML = FALSE) # (Covid|IDLocality) +
      est_m3b = est_out(m3b, '03b) (1|Year) + (1|genus) + (1|Species)+(1|sp_day_year) + (scale(Covid)|Country) + (1|IDLocality) + (1|sp_loc); >9/species/period')  
  # plot
      xs = rbind(est_m1a,est_m1b,est_m1c, est_m2a, est_m2b,est_m3a, est_m3b)

      g = 
      ggplot(xs[predictor == 'scale(Covid)'], aes(y = model, x = estimate, col = model)) +
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
          labs(y = NULL ,x = "Standardized effect size", title = "a) Effect of Period (before/during shutdown)") +
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
                axis.title=element_text(size=7)
                )
      g
      ggsave(here::here('Outputs/Figure_Sy.png'),g, width = 30, height =5, units = 'cm')
        
   
# Figure Sz   
    # estimates 
     m01a=lmer(scale(log(FID))~ 
        Year+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(StringencyIndex)+ 
        (scale(StringencyIndex)|genus)+(1|Species)+(1|sp_day_year) + (1|Country) + (scale(StringencyIndex)|IDLocality) +(1|sp_loc),
        data = s, REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
     s[, res := resid(m01a)]
        # (1|Year) explains nothing - could stay 
     est_m01a = est_out(m01a, '01a) (scale(StringencyIndex)|genus)+(1|Species)+(1|sp_day_year) + (1|Country) + (scale(StringencyIndex)|IDLocality) +(1|sp_loc)')
     m01b=lmer(scale(log(FID))~
          Year+       
          scale(log(SD))+
          scale(log(FlockSize))+
          scale(log(BodyMass))+
          scale(sin(rad)) + scale(cos(rad)) + 
          #scale(Day)+
          scale(Temp)+
          scale(StringencyIndex)+
          (1|genus) +(1|Species) + (1|sp_day_year) + (1|Country) + (scale(StringencyIndex)|IDLocality)+(1|sp_loc),  
          data = s, REML = FALSE, 
          control = lmerControl(
              optimizer ='optimx', optCtrl=list(method='nlminb'))
          ) 
          # (1|Year), (1|genus), (1|sp_day_year), (1|sp_loc),  explain nothing - could stay
     est_m01b = est_out(m01b, '01b) (1|genus) +(1|Species) + (1|sp_day_year) + (1|Country) + (scale(StringencyIndex)|IDLocality)+(1|sp_loc)')  
     m01c=lmer(scale(log(FID))~
          Year+       
          scale(log(SD))+
          scale(log(FlockSize))+
          scale(log(BodyMass))+
          scale(sin(rad)) + scale(cos(rad)) + 
          #scale(Day)+
          scale(Temp)+
          scale(StringencyIndex)+
          (1|Country) + (scale(StringencyIndex)|IDLocality),  
          data = s, REML = FALSE, 
          control = lmerControl(
              optimizer ='optimx', optCtrl=list(method='nlminb'))
          ) 
          # (1|Year), (1|genus), (1|sp_day_year), (1|sp_loc),  explain nothing - could stay
     est_m01c = est_out(m01c, '01c) (1|Country) + (scale(StringencyIndex)|IDLocality)')  

     m02a=lmer(scale(log(FID))~ 
        Year+ 
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
     est_m02a = est_out(m02a, '02a) (scale(StringencyIndex)|genus)+(1|Species)+(1|sp_day_year) + (1|Country) + (scale(StringencyIndex)|IDLocality) +(1|sp_loc);>4/species')   
     m02b=lmer(scale(log(FID))~
          Year+ 
          scale(log(SD))+
          scale(log(FlockSize))+
          scale(log(BodyMass))+
          scale(sin(rad)) + scale(cos(rad)) + 
          #scale(Day)+
          scale(Temp)+
          scale(StringencyIndex)+
          (1|genus) +(1|Species) + (1|sp_day_year) + (1|Country) + (scale(StringencyIndex)|IDLocality)+(1|sp_loc),  
          data = s[Nsp>4], REML = FALSE, 
          control = lmerControl(
              optimizer ='optimx', optCtrl=list(method='nlminb'))
          ) 
     est_m02b = est_out(m02b, '02b)  (1|genus) +(1|Species) + (1|sp_day_year) + (1|Country) + (scale(StringencyIndex)|IDLocality)+(1|sp_loc); >4/species') 
     m02c=lmer(scale(log(FID))~
          Year+ 
          scale(log(SD))+
          scale(log(FlockSize))+
          scale(log(BodyMass))+
          scale(sin(rad)) + scale(cos(rad)) + 
          #scale(Day)+
          scale(Temp)+
          scale(StringencyIndex)+
          (1|Country) + (scale(StringencyIndex)|IDLocality),  
          data = s[Nsp>4], REML = FALSE, 
          control = lmerControl(
              optimizer ='optimx', optCtrl=list(method='nlminb'))
          ) 
     est_m02c = est_out(m02c, '02c)  (1|Country) + (scale(StringencyIndex)|IDLocality); >4/species') 
     
     m03a=lmer(scale(log(FID))~ 
        Year+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(StringencyIndex)+ 
        (scale(StringencyIndex)|genus)+(1|Species)+(1|sp_day_year) + 
        (1|Country) + (scale(StringencyIndex)|IDLocality) +(1|sp_loc),   
        data = s[Nsp>9], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
        # (1|Year) explains nothing - could stay 
     est_m03a = est_out(m03a, '03a) (scale(StringencyIndex)|genus)+(1|Species)+(1|sp_day_year) + (1|Country) + (scale(StringencyIndex)|IDLocality) +(1|sp_loc);>9/species') 
     m03b=lmer(scale(log(FID))~
          Year+ 
          scale(log(SD))+
          scale(log(FlockSize))+
          scale(log(BodyMass))+
          scale(sin(rad)) + scale(cos(rad)) + 
          #scale(Day)+
          scale(Temp)+
          scale(StringencyIndex)+
          (1|genus) +(1|Species) + (1|sp_day_year) + (1|Country) + (scale(StringencyIndex)|IDLocality)+(1|sp_loc),  
          data = s[Nsp>9], REML = FALSE, 
          control = lmerControl(
              optimizer ='optimx', optCtrl=list(method='nlminb'))
          ) 
     est_m03b = est_out(m03b, '03b)  (1|genus) +(1|Species) + (1|sp_day_year) + (1|Country) + (scale(StringencyIndex)|IDLocality)+(1|sp_loc); >9/species') 
     m03c=lmer(scale(log(FID))~
          Year+ 
          scale(log(SD))+
          scale(log(FlockSize))+
          scale(log(BodyMass))+
          scale(sin(rad)) + scale(cos(rad)) + 
          #scale(Day)+
          scale(Temp)+
          scale(StringencyIndex)+
          (1|Country) + (scale(StringencyIndex)|IDLocality),  
          data = s[Nsp>9], REML = FALSE, 
          control = lmerControl(
              optimizer ='optimx', optCtrl=list(method='nlminb'))
          ) 
     est_m03c = est_out(m03c, '03c)  (1|Country) + (scale(StringencyIndex)|IDLocality); >9/species') 

    # plot
      xs0 = rbind(est_m01a,est_m01b,est_m01c, est_m02a, est_m02b,est_m02c, est_m03a, est_m03b, est_m03c)

      r g0
      ggsave(here::here('Outputs/Figure_Sz.png'),g0, width = 30, height =5, units = 'cm')
  

# Figure Szy   
      ggsave(here::here('Outputs/Figure_Szb.png'),rbind(ggplotGrob(g),ggplotGrob(g0)), width = 30, height =10, units = 'cm')


# END