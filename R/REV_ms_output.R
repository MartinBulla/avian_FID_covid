#' ---
#' title: <font size="2">Supporting code, figures and tables for</font><br><font size="5">Urban birds’ tolerance towards humans was largely unaffected by COVID-19 shutdown-induced variation in human presence</font>
#' author: <font size="2">Peter Mikula, Martin Bulla, Daniel T. Blumstein, Yanina Benedetti, Kristina Floigl, Jukka Jokimäki, Marja-Liisa Kaisanlahti-Jokimäki, Gábor Markó, Federico Morelli, Anders Pape Møller, Anastasiia Siretckaia, Sára Szakony, Michael A. Weston, Farah Abou Zeid, Piotr Tryjanowski & Tomáš Albrecht</font><br><br><font size="2">created by Martin Bulla</font><br>
#' date: <font size="1.5">`r Sys.time()`</font>
#' bibliography: _bib.bib
#' output:
#'     html_document:
#'         toc: true
#'         toc_float: true
#'         toc_depth: 4
#'         code_folding: hide
#'         link-citations: yes
#' base:  href="/[avian_FID_covid]/"
#' ---

#' <style> body {text-align: justify}</style>

#+ r setup, include=FALSE
knitr::opts_chunk$set(message = FALSE, warning = FALSE, cache = TRUE)

#' ***

#' ###### Code to load tools and data
#+ start, echo = T, results = 'hide', warning=FALSE
 # packages
    require(ape)
    require(arm)
    require(brms)
    require(car)
    require(coda)
    require(data.table)
    require(dplyr)
    require(effects)
    require(emo)
    require(foreach)
    require(geiger)
    require(ggimage)
    require(ggplot2)
    require(ggpubr)
    require(ggsci)
    require(ggtext)
    require(grid)
    require(gtable)
    require(here)
    require(kableExtra)
    require(magrittr)
    require(MASS)
    require(MCMCglmm)
    require(multcomp)
    require(optimx)
    require(pander)
    require(parallel)
    require(patchwork)
    require(performance)  
    require(PerformanceAnalytics)
    require(phangorn)
    require(phylobase)
    require(phytools)
    require(plyr)
    require(png)
    require(RColorBrewer)
    require(rmeta)
    require(rphylopic)
    require(scales)
    require(stringr)
    require(viridis)

 # constants
    save_plot = TRUE
    round_ = 3 # number of decimal places to round model coefficients
    nsim = 5000 # number of simulations to extract estimates and 95%CrI
    ax_lines = "grey60" # defines color of the axis lines
    #colors <- c("#999999", "#E69F00", "#56B4E9") #viridis(3)
    col_fit = '#f2e030'# line color of loes smoothed fit '#fdff00'#'#8cff32'#'#abff32'#'orange'
    set.seed(42)
    #width_ = .7 # spacing between error bars
    #col_ = c(brewer.pal(n =12, name = "Paired"), 'grey30','grey80')
 # functions
    # mode
      Mode <- function(x) {
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
      }
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
                  axis.ticks = element_line(colour = "grey70", linewidth = 0.1),
                  #axis.ticks.margin,
                  
                  strip.text.x = element_text(size = 6, color="grey30",  margin=margin(1,1,1,1,"mm")), #grey50
                  strip.text.y = element_text(size = 6, color="grey30",  margin=margin(1,1,1,1,"mm")), #grey50
                  strip.background = element_rect(fill="grey99",colour="grey70", linewidth=0.25),
                    #strip.background = element_blank(), 
                    #strip.text = element_blank(),
                  panel.spacing = unit(0, "mm"),
                  panel.background=element_blank(),
                  panel.border = element_rect(colour="grey70", linewidth=0.1, fill = NA), #panel.border=element_blank(),
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
        theme_MB2 =  theme(
            title = element_text(size=8, colour="grey30"),
            #axis.title.y = element_text(vjust=3.5),
            #axis.title.x = element_text(vjust=1),
            axis.title = element_text(size = 7, color = 'black'),
            axis.text = element_text(size=6),
            axis.line = element_blank(),#axis.line = element_line(colour = ax_lines, linewidth = 0.25),
            axis.ticks = element_line(colour = "grey70", linewidth = 0.1),
            axis.ticks.length = unit(1, "pt"),
        
            legend.position = c(1, 0.06),
            legend.justification = c(1, 0),
            legend.title = element_blank(),
            legend.spacing.y = unit(-0.9, "cm"),
            legend.text=element_text(size=6),
            legend.key = element_rect(colour = NA, fill = NA),
            legend.key.height= unit(0.5,"line"),
            legend.key.width = unit(0.25, "cm"),
            legend.margin = margin(0,0,0,0, unit="cm"),
            legend.box.margin = margin(l = -6), #legend.justification = c(-1,0),
            legend.background = element_blank(),

            strip.background = element_blank(),
            strip.text.x = element_text(size = 5, color="grey30",  margin=margin(1,1,1,1,"mm")), 
            strip.text.y = element_text(size = 5, color="grey30",  margin=margin(1,1,1,1,"mm")), 

            panel.spacing = unit(0, "mm"),
            panel.background=element_blank(),
            panel.border = element_rect(colour="grey70", linewidth=0.1, fill = NA),
            panel.grid = element_blank()
        ) 
    # for estimates
        est_out =function(model = m, label = "", nsim = 5000){
            bsim = sim(model, n.sim=5000) 
            v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
            ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 
            sd = apply(bsim@fixef, 2, sd)
            data.table(predictor=rownames(coef(summary(model))),estimate=v, lwr=ci[1,], upr=ci[2,], sd = sd, model = paste(label, "N =", nobs(model)))
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
      m_out = function(model = m, type = "mixed", 
        name = "define", dep = "define", fam = 'Gaussian',
        round_ = 3, nsim = 5000, aic = FALSE, save_sim = here::here('Data/model_sim/'), back_tran = FALSE, perc_ = 1, R2 = FALSE){
          # perc_ 1 = proportion or 100%
        bsim = sim(model, n.sim=nsim)  
        
        if(save_sim!=FALSE){save(bsim, file = paste0(save_sim, name,'.RData'))}
       
        if(type != "mixed"){
          v = apply(bsim@coef, 2, quantile, prob=c(0.5))
          ci = apply(bsim@coef, 2, quantile, prob=c(0.025,0.975)) 

          if(back_tran == TRUE & fam == "binomial"){
           v = perc_*plogis(v)
           ci = perc_*plogis(ci)
           }
          if(back_tran == TRUE & fam == "binomial_logExp"){
                v = perc_*(1-plogis(v))
                ci = perc_*(1-plogis(ci))
                ci = rbind(ci[2,],ci[1,])
               }

          if(back_tran == TRUE & fam == "Poisson"){
           v = exp(v)
           ci = exp(ci)
          }

         oi=data.frame(type='fixed',effect=rownames(coef(summary(model))),estimate=v, lwr=ci[1,], upr=ci[2,])
          rownames(oi) = NULL
          oi$estimate_r=round(oi$estimate,round_)
          oi$lwr_r=round(oi$lwr,round_)
          oi$upr_r=round(oi$upr,round_)
          if(perc_ == 100){
           oi$estimate_r = paste0(oi$estimate_r,"%")
           oi$lwr_r = paste0(oi$lwr_r,"%")
           oi$upr_r = paste0(oi$upr_r,"%")
          }
         x=data.table(oi[c('type',"effect", "estimate_r","lwr_r",'upr_r')]) 
       
        }else{
         v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
         ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975)) 

         if(back_tran == TRUE & fam == "binomial"){
          v = perc_*plogis(v)
          ci = perc_*plogis(ci)
         }
          if(back_tran == TRUE & fam == "binomial_logExp"){
                v = perc_*(1-plogis(v))
                ci = perc_*(1-plogis(ci))
                ci = rbind(ci[2,],ci[1,])
               }

          if(back_tran == TRUE & fam == "Poisson"){
            v = exp(v)
            ci = exp(ci)
         }

        oi=data.table(type='fixed',effect=rownames(coef(summary(model))),estimate=v, lwr=ci[1,], upr=ci[2,])
            rownames(oi) = NULL
            oi[,estimate_r := round(estimate,round_)]
            oi[,lwr_r := round(lwr,round_)]
            oi[,upr_r :=round(upr,round_)]
            if(perc_ == 100){
             oi[,estimate_r := paste0(estimate_r,"%")]
             oi[,lwr_r := paste0(lwr_r,"%")]
             oi[,upr_r := paste0(upr_r,"%")]
            }
         oii=oi[,c('type',"effect", "estimate_r","lwr_r",'upr_r')] 
        
         l=data.frame(summary(model)$varcor)
         l = l[is.na(l$var2),]
         l$var1 = ifelse(is.na(l$var1),"",l$var1)
         l$pred = paste(l$grp,l$var1)

         q050={}
         q025={}
         q975={}
         pred={}
         
         # variance of random effects
         for (ran in names(bsim@ranef)) {
           ran_type = l$var1[l$grp == ran]
           for(i in ran_type){
            q050=c(q050,quantile(apply(bsim@ranef[[ran]][,,ran_type], 1, var), prob=c(0.5)))
            q025=c(q025,quantile(apply(bsim@ranef[[ran]][,,ran_type], 1, var), prob=c(0.025)))
            q975=c(q975,quantile(apply(bsim@ranef[[ran]][,,ran_type], 1, var), prob=c(0.975)))
            pred= c(pred,paste(ran, i))
            }
           }
         # residual variance
         q050=c(q050,quantile(bsim@sigma^2, prob=c(0.5)))
         q025=c(q025,quantile(bsim@sigma^2, prob=c(0.025)))
         q975=c(q975,quantile(bsim@sigma^2, prob=c(0.975)))
         pred= c(pred,'Residual')

         ri=data.table(type='random',effect=pred, estimate_r=round(100*q050/sum(q050)), lwr_r=round(100*q025/sum(q025)), upr_r=round(100*q975/sum(q975)))
           
         ri[lwr_r>upr_r, lwr_rt := upr_r]
         ri[lwr_r>upr_r, upr_rt := lwr_r]
         ri[!is.na(lwr_rt), lwr_r := lwr_rt]
         ri[!is.na(upr_rt), upr_r := upr_rt]
         ri$lwr_rt = ri$upr_rt = NULL

         ri[,estimate_r := paste0(estimate_r,'%')]
         ri[,lwr_r := paste0(lwr_r,'%')]
         ri[,upr_r := paste0(upr_r,'%')]
        
        x = data.table(rbind(oii,ri))
        }
        
        x[1, model := name]                                                                
        x[1, response := dep]                                                                
        x[1, error_structure := fam]      
        N = length(resid(model))                                                          
        x[1, N := N ]                                                                

        x=x[ , c('model', 'response', 'error_structure', 'N', 'type',"effect", "estimate_r","lwr_r",'upr_r')] 

        if (aic == TRUE){   
            x[1, AIC := AIC(update(model,REML = FALSE))] 
            }
        if (aic == "AICc"){
            aicc = AICc(model)
            x[1, AICc := aicc] 
        }
        if(type == "mixed" & nrow(x[type=='random' & estimate_r =='0%'])==0 & R2 == TRUE){
          R2_mar = as.numeric(r2_nakagawa(model)$R2_marginal)
          R2_con = as.numeric(r2_nakagawa(model)$R2_conditional)
          x[1, R2_mar := R2_mar]
          x[1, R2_con := R2_con]
         }
        x[is.na(x)] = ""
        return(x)
      } 
    # model assumption function
      m_ass = function(name = 'define', mo = m0, dat = d, fixed = NULL, categ = NULL, trans = "none", spatial = TRUE, temporal = TRUE, PNG = TRUE, outdir = 'outdir', n_col=8, width_ = 10, height_ = 5.5){
       l=data.frame(summary(mo)$varcor)
       l = l[is.na(l$var2),]
       nt = if(temporal==TRUE){1}else{0}
       ns = if(spatial==TRUE){7}else{0}
       n = 3+nrow(l)-1+length(fixed)+length(categ) +  nt +  ns
     
       if(PNG == TRUE){
        png(paste(outdir,name, ".png", sep=""), width=width_,height=height_,units="in",res=300) # width = 6; res = 150 ok for html
        par(mfrow=c(4, n_col),tcl = -0.08, cex = 0.5, cex.main = 0.9,#ceiling(n/n_col),n_col)
            oma = c(1,1,4,1),
            mar = c(2, 2, 2, 1), mgp=c(1,0,0)
            )
         }else{
          dev.new(width=width_,height=height_)
          par(mfrow=c(4,n_col), tcl = -0.08, cex = 0.5, cex.main = 0.9,#ceiling(n/n_col),n_col)
            oma = c(1,1,2,1),
            mar = c(2, 2, 2, 1), mgp=c(1,0,0)
            )
        }
       
       scatter.smooth(fitted(mo),resid(mo),col='grey');abline(h=0, lty=2, col ='red')
       scatter.smooth(fitted(mo),sqrt(abs(resid(mo))), col='grey')
       qqnorm(resid(mo), main=list("Normal Q-Q Plot: residuals"),col='grey');qqline(resid(mo), col = 'red')
       #unique(l$grp[l$grp!="Residual"])
       for(i in unique(l$grp[l$grp!="Residual"])){
        #i = "mean_year"
        ll=ranef(mo)[names(ranef(mo))==i][[1]]
        if(ncol(ll)==1){
         qqnorm(ll[,1], main = paste(i,names(ll)[1]),col='grey',);qqline(ll[,1], col ='red')
         }else{
          qqnorm(ll[,1], main = paste(i,names(ll)[1]),col='grey');qqline(ll[,1], col ='red')
          qqnorm(ll[,2], main = paste(i,names(ll)[2]),col='grey');qqline(ll[,2], col ='red')
         }
        }
        
       # variables
         scatter={} 
         for (i in rownames(summary(mo)$coef)) {
              # i = "lat_abs" #i = rownames(summary(mo)$coef)[9]
            j=sub("\\).*", "", sub(".*\\(", "",i)) 
            scatter[length(scatter)+1]=j
          }
          x = data.frame(scatter=unique(scatter)[2:length(unique(scatter))],
                          log_ = grepl("log",rownames(summary(mo)$coef)[2:length(unique(scatter))]), stringsAsFactors = FALSE)
          for (i in 1:length(fixed)){
              jj =fixed[i] # jj = fixed[1]
              #print(jj)
              variable=dat[, ..jj][[1]]
              if(trans[i]=='log'){
              scatter.smooth(resid(mo)~log(variable),xlab=paste('log(',jj,')',sep=''), col = 'grey');abline(h=0, lwd=1, lty = 2, col ='red')
              }else if(trans[i]=='abs'){
              scatter.smooth(resid(mo)~abs(variable),xlab=paste('abs(',jj,')',sep=''), col = 'grey');abline(h=0, lwd=1, lty = 2, col ='red')
              }else if(trans[i]=='sin'){scatter.smooth(resid(mo)~sin(variable),xlab=paste('sin(',jj,')',sep=''), col = 'grey');abline(h=0, lwd=1, lty = 2, col ='red')
              }else if(trans[i]=='cos'){scatter.smooth(resid(mo)~cos(variable),xlab=paste('cos(',jj,')',sep=''), col = 'grey');abline(h=0, lwd=1, lty = 2, col ='red')
              }else{
              scatter.smooth(resid(mo)~variable,xlab=jj,col = 'grey');abline(h=0, lwd=1, lty = 2, col ='red')
            }
           }
          
          if(length(categ)>0){
            for(i in categ){
               variable=dat[, ..i][[1]]
                boxplot(resid(mo)~variable, medcol='grey', whiskcol='grey', staplecol='grey', boxcol='grey', outcol='grey', xlab = i);abline(h=0, lty=3, lwd=1, col = 'red')
               }
          }     
              
        if(temporal == TRUE){
            acf(resid(mo), type="p", main=list("Temporal autocorrelation:\npartial series residual"))
            }
        if(spatial == TRUE){    
          spdata=data.frame(resid=resid(mo), x=dat$Lon, y=dat$Lat)
            spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
            #cex_=c(1,2,3,3.5,4)
            cex_=c(1,1.5,2,2.5,3)
            spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
          plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8), xlab = 'Longitude', ylab = 'Latitude')
            legend("topright", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)))
          plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('residuals <0'), xlab = 'Longitude', ylab = 'Latitude')
          plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('residual >=0'), xlab = 'Longitude', ylab = 'Latitude')
          
          if('Australia'%in%unique(dat$Country) & length(unique(dat$Country))>1){
          # EU
          dat$res = resid(mo)
          spdata=data.frame(resid = dat$res[dat$Country!='Australia'], x=dat$Lon[dat$Country!='Australia'], y=dat$Lat[dat$Country!='Australia'])
          spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
            #cex_=c(1,2,3,3.5,4)
            cex_=c(1,1.5,2,2.5,3)
            spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
          plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('EU -  residuals <0'), xlab = 'Longitude', ylab = 'Latitude')
          plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('EU residuals >=0)'), xlab = 'Longitude', ylab = 'Latitude')

          # Australia
          spdata=data.frame(resid = dat$res[dat$Country=='Australia'], x=dat$Lon[dat$Country=='Australia'], y=dat$Lat[dat$Country=='Australia'])
          spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
            #cex_=c(1,2,3,3.5,4)
            cex_=c(1,1.5,2,2.5,3)
            spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
          plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Australia residuals <0'), xlab = 'Longitude', ylab = 'Latitude')
          plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Australia residuals >=0'), xlab = 'Longitude', ylab = 'Latitude')
          }
        }
  
       mtext(stringr::str_wrap(paste(paste0(name," model: "), slot(mo,"call")[1],'(',slot(mo,"call")[2],sep=''), width = ceiling(nchar(paste(slot(mo,"call")[1],'(',slot(mo,"call")[2],sep=''))/2)+10), side = 3, line = 1, cex=0.5,outer = TRUE, col = 'darkblue') #ceiling(nchar(paste(slot(mo,"call")[1],'(',slot(mo,"call")[2],sep=''))/2)
       if(PNG==TRUE){dev.off()}
      }
    
 # data
    t = fread(here::here("Data/taxonomy.txt"))
    
    ph  =  fread(here::here('Data/phylopic.txt'))
    setnames(ph, old = c('Name', 'Code'), new = c('genus2', 'uid'))

    g = fread(here::here('Data/google_mobility.txt')) #fwrite(d, here::here('Data/data.txt'), sep ='\t')
    g[, Year := as.integer(substring(date, nchar(date)-3, nchar(date)))]
    g[nchar(date)==9, date:=paste0('0',date)]
    g[, date_ :=as.Date(date, format = '%d.%m.%Y')]
    g[, Day :=yday(date_)]
    setnames(g, old = "country_region", new = "Country")

    # adjust correct assignment of breeding season (Year) for Australiato the the year of when a given breeding season started
    g[Country != "Australia", Day := Day - 92 + 1] # 1 April = start of breeding season (1st day) = 92 day of the year
    g[Country == "Australia", Day := Day - 228 + 1] # 15 Augusst = start of breeding season (1st day) = 228 day of the year
    
    d = fread(here::here('Data/data.txt')) #fwrite(d, here::here('Data/data.txt'), sep ='\t')
    d[Human == 0, humans := 0]
    d[Human > 0, humans := 1] # dh = dh[!(Country%in%'Czechia' & Year == 2018)]
  
    # adjust correct assignment of breeding season (Year) for Australiato the the year of when a given breeding season started
    d[Country == 'Australia' & Year == 2020 & Covid == 0, Year:=2019]
    d[Country == 'Australia' & Year == 2021 & Day>139, Year:=2020]
    d[Country == 'Australia' & Year == 2022 & Day>139, Year:=2021]

    d = d[order(Year, IDLocality, Day, Hour)]
    d[, year_ := as.character(Year)]
    d[Country %in% c("Czech_Republic", "Czech Republic"), Country := "Czechia"]
    d[, genus := sub("_.*", "", Species)]
    d[, sp := gsub("[_]", " ", Species)]
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
    d[, weekday := weekdays(date_)]
    #d[Country == 'Australia', Day_:= abs(Day - 189)]
    
    # add species abbreviations
    
    si = data.table(Species = unique(d$Species), nch = nchar(unique(d$Species)))
    si[, scinam2:=Species]
    si[nch>15, scinam2:=paste0(substr(Species, 1,3),'. ',sub('.*_', '', Species))]
    si[is.na(scinam2), scinam2 := Species]
    si[, nch2 := nchar(scinam2)]
    si[, scinam_abb:=scinam2]
    si[nch2>15, scinam_abb:=paste0(substr(scinam2, 1, nch2 - (nch2-14)),'.')]
    si[, scinam_abb := gsub("[_]", " ", scinam_abb)]
    si[scinam_abb%in%'Ana. platyrhyn.', scinam_abb:= 'Anas platyrhy.']
    d = merge(d, si[,.(Species, scinam_abb)], by = 'Species')

    # species with data before and during
    d1 = d[Covid == 1, .N, by = Species]
    d2 = d[Covid == 0, .N, by = Species]
    setnames(d1, old = 'N', new ='N_during')
    setnames(d2, old = 'N', new ='N_before')
    dd = merge(d1, d2) # species with data before and during
    da = merge(d1,d2, all = TRUE)

    # species with data before and during
    d1p = d[Country!='Poland' & Covid == 1, .N, by = Species]
    d2p = d[Country!='Poland' & Covid == 0, .N, by = Species]
    setnames(d1p, old = 'N', new ='N_during')
    setnames(d2p, old = 'N', new ='N_before')
    ddp = merge(d1p,d2p) # species with data before and during, but without Poland data

    # species-localities with data before and during
    p1 = d[Covid == 1, .N, by = .(IDLocality, Species)]
    p2 = d[Covid == 0, .N, by = .(IDLocality, Species)]
    setnames(p1, old = 'N', new ='N_during')
    setnames(p2, old = 'N', new ='N_before')
    pp = merge(p1,p2)  # species-localities with data before and during
    pa = merge(p1,p2, all = TRUE)

    # species-localities with data before and during,but without Poland data
    p1p = d[Country!='Poland' & Covid == 1, .N, by = .(IDLocality, Species)]
    p2p = d[Country!='Poland' & Covid == 0, .N, by = .(IDLocality, Species)]
    setnames(p1p, old = 'N', new ='N_during')
    setnames(p2p, old = 'N', new ='N_before')
    ppp = merge(p1p,p2p)  # species-localities with data before and during,but without Poland data 

    # species-localities with data before and during, but without Poland & 2014 data
    p1p4 = d[Year!=2014 & Country!='Poland' & Covid == 1, .N, by = .(IDLocality, Species)]
    p2p4 = d[Year!=2014 & Country!='Poland' & Covid == 0, .N, by = .(IDLocality, Species)]
    setnames(p1p4, old = 'N', new ='N_during')
    setnames(p2p4, old = 'N', new ='N_before')
    ppp4 = merge(p1p4,p2p4) # species-localities with data before and during, but without Poland & 2014 data

    # add google mobility
    d = merge(d, g[,.(Country,  date_, parks_percent_change_from_baseline)], all.x = TRUE, by = c('Country', 'date_'))
     
    # limit to data with # of humans
    dh = d[!is.na(Human)] # summary(factor(dh$Year))
    dh[, Nsp := .N, by = "Species"]
    dh[, Country := factor(Country, levels = (c("Finland", "Poland", "Czechia", "Hungary")))]
    dh[Covid == 0, Period := "Before COVID-19 shutdown"]
    dh[Covid == 1, Period := "During COVID-19 shutdown"]

    # limit to data with # of humans > 0
    dhh <- dh[Human > 0]
    dhh[, Country := factor(Country, levels = (c("Finland", "Poland", "Czechia", "Hungary")))]
    dhh[Covid == 0, Period := "Before COVID-19 shutdown"]
    dhh[Covid == 1, Period := "During COVID-19 shutdown"]
    
    # limit to COVID-19 period (Stringeny data)
    s = d[Covid == 1]
    s[, Nsp := .N, by = "Species"]
    s[, year_weekday := paste(Year, weekday)]

    # limit to data with Google Mobility
    ss = s[!is.na(parks_percent_change_from_baseline)]
    ss[, country_year := paste(Country, Year)] #table(paste(s$Country, s$Year))   
    ss[, GoogleMobility := parks_percent_change_from_baseline]
    ss[parks_percent_change_from_baseline<0, google := 'before_zero']
    ss[parks_percent_change_from_baseline>0, google := 'after_zero']
    ss[, sp_country_google:= paste(sp_country, google)]

    g[, weekday := weekdays(date_)]

    # limit stringency data to those with # of humans
    sh <- s[!is.na(Human)]
    sh[, year_day := paste(Year, Day)]
    sh[, year_weekday := paste(Year, weekday)]

    # limit Google data # limit to data with # of humans
    ssh <- ss[!is.na(Human)]
    ssh[, year_day := paste(Year, Day)]
    ssh[, year_weekday := paste(Year, weekday)]

    # add scinam
    d[, scinam := Species]
    s[, scinam := Species]
    ss[, scinam := Species]
    dh[, scinam := Species]
#'
#' ***
#' 
#' ### Introduction
#' To facilitate transparency, the following document contains example code used to generate the results of the manuscript. Thus, apart from the supplementary figures and tables, below we display also main text figures and plot of model assumptions. The figures and tables are ordered according to their appearance in the main text. 
#' 
#' The code is displayed upon clicking the "code" icon at the top right, above each display item. Continuous variables were  standardised by subtracting the mean and dividing by the standard deviation. For descriptions of variables see Methods of the [paper](https://doi.org/10.1101/2022.07.15.500232). 'Residual' in tables indicates residual variance.
#' 
#' When using this content **PLEASE CITE** the [paper](https://doi.org/10.1101/2022.07.15.500232) and this repository[<sup>1</sup>](#1)
#' 
#' Questions can be directed to bulla.mar@gmail.com & petomikula158@gmail.com
#' 
#' <a name="1"></a>(1) Martin Bulla, Peter Mikula, Daniel T. Blumstein, Yanina Benedetti, Kristina Floigl, Jukka Jokimäki, Marja-Liisa Kaisanlahti-Jokimäki, Gábor Markó, Federico Morelli, Anders Pape Møller, Anastasiia Siretckaia, Sára Szakony, Michael A. Weston, Farah Abou Zeid, Piotr Tryjanowski & Tomáš Albrecht (2024). *Supporting information for "Urban birds’ tolerance towards humans was largely unaffected by COVID-19 shutdown-induced variation in human presence"*, GitHub, [https://martinbulla.github.io/avian_FID_covid/](https://martinbulla.github.io/avian_FID_covid/); DOI: [https://doi.org/10.5281/zenodo.11633804](https://doi.org/10.5281/zenodo.11633804).
#' 
#' 
#' ***
#' 
#' ### Repository: files & folders
#' [**Supplementary information**, including code](https://martinbulla.github.io/avian_FID_covid/): the current html document with supplementary informatiion, figures and tables and plot of [model assumptions](#MS_).
#'  
#' [**R**](https://github.com/MartinBulla/avian_FID_covid/tree/main/R/)-scripts used in the analysis:  
#' - [_runRmarkdown.R](https://github.com/MartinBulla/avian_FID_covid/tree/main/R/_runRmarkdown.R) generates htmls from the following R-script:  
#' - [REV_ms_output.R](https://github.com/MartinBulla/avian_FID_covid/tree/main/R/) used to generate the [Supplement](https://martinbulla.github.io/avian_FID_covid/), contains all scripts used to generate the paper outputs, including the display items  
#' <br />
#' [**Data**](https://github.com/MartinBulla/avian_FID_covid/tree/main/Data):  
#' - *raw data* (for their description see [READ_ME](https://github.com/MartinBulla/avian_FID_covid/tree/main/Data/READ_ME.md))   
#' - *manipulated data* (starting with 'dat_') generated with R-scripts and used in the further analyses  
#' - [*posterior model simulations*](https://osf.io/at85z/) - located at OSF
#' - [*phylopic pictures*](https://github.com/MartinBulla/avian_FID_covid/tree/main/Data/Pics/) used in the graphs
#'  
#' [**Outputs**](https://github.com/MartinBulla/avian_FID_covid/tree/main/Outputs/): separate files of all outputs used in the manuscript and this Supplement
#'  
#' [LICENSE](https://github.com/MartinBulla/avian_FID_covid/tree/main/LICENSE.txt): terms of reuse - applicable only after this work is published as a preprint or in a scientific journal, until then the data are not available for reuse.
#'  
#' ***
#' 
#' ### Effect sizes
# HUMAN #s
## predictions
# full
mhs <- lmer(scale(log(FID)) ~
  scale(log(SD)) +
  scale(log(FlockSize)) +
  scale(log(BodyMass)) +
  scale(sin(rad)) + scale(cos(rad)) +
  # scale(Day)+
  scale(Temp) +
  scale(Human) +
  (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) +
  (scale(Human) | Country) + (1 | IDLocality) + (1 | sp_loc),
data = dh, REML = FALSE
)

est_mhs <- est_out(mhs, "ALL: (1|Year) + (1|weekday) + (1|genus) + (1|Species)  + (1|sp_day_year) + (scale(Human) |Country) + (1|IDLocality) +(1|sp_loc)")
est_mhs[, control_for_starting_distance := "yes"]

mhx <- lmer(scale(log(FID)) ~
  # scale(log(SD)) +
  scale(log(FlockSize)) +
  scale(log(BodyMass)) +
  scale(sin(rad)) + scale(cos(rad)) +
  # scale(Day)+
  scale(Temp) +
  scale(Human) +
  (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) +
  (scale(Human) | Country) + (1 | IDLocality) + (1 | sp_loc),
data = dh, REML = FALSE
)
est_mhx <- est_out(mhx, "ALL: (1|Year) + (1|weekday) + (1|genus) + (1|Species) + (1|sp_day_year) + (scale(Human) |Country) + (1|IDLocality) +(1|sp_loc)")
est_mhx[, control_for_starting_distance := "no"]


# CZ 
chs <- lmer(scale(log(FID)) ~
  scale(Year) +
  scale(log(SD)) +
  scale(log(FlockSize)) +
  scale(log(BodyMass)) +
  scale(sin(rad)) + scale(cos(rad)) +
  # scale(Day)+
  scale(Temp) +
  scale(Human) +
  (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = dh[Country == "Czechia"], REML = FALSE
)
est_chs <- est_out(chs, "Czechia:(1|weekday)+(1|genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
est_chs[, control_for_starting_distance := "yes"]

chx <- lmer(scale(log(FID)) ~
  scale(Year) +
  scale(log(SD)) +
  scale(log(FlockSize)) +
  scale(log(BodyMass)) +
  scale(sin(rad)) + scale(cos(rad)) +
  # scale(Day)+
  scale(Temp) +
  scale(Human) +
  (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = dh[Country == "Czechia"], REML = FALSE
)
est_chx <- est_out(chx, "Czechia: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
est_chx[, control_for_starting_distance := "no"]

# FI
fhs <- lmer(scale(log(FID)) ~
  scale(Year) +
  scale(log(SD)) +
  scale(log(FlockSize)) +
  scale(log(BodyMass)) +
  scale(sin(rad)) + scale(cos(rad)) +
  # scale(Day)+
  scale(Temp) +
  scale(Human) +
  (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
data = dh[Country == "Finland"], REML = FALSE
)
est_fhs <- est_out(fhs, "Finland: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
est_fhs[, control_for_starting_distance := "yes"]

fhx <- lmer(scale(log(FID)) ~
  scale(Year) +
  # scale(log(SD)) +
  scale(log(FlockSize)) +
  scale(log(BodyMass)) +
  scale(sin(rad)) + scale(cos(rad)) +
  # scale(Day)+
  scale(Temp) +
  scale(Human) +
  (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
data = dh[Country == "Finland"], REML = FALSE
)
est_fhx <- est_out(fhx, "Finland: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
est_fhx[, control_for_starting_distance := "no"]

# HU - singular fits only due to sp_loc estimated as zero (removing it changes no results)
hhs <- lmer(scale(log(FID)) ~
  scale(Year) +
  scale(log(SD)) +
  scale(log(FlockSize)) +
  scale(log(BodyMass)) +
  scale(sin(rad)) + scale(cos(rad)) +
  # scale(Day)+
  scale(Temp) +
  scale(Human) +
  (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = dh[Country == "Hungary"], REML = FALSE
)
est_hhs <- est_out(hhs, "Hungary: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
est_hhs[, control_for_starting_distance := "yes"]

hhx <- lmer(scale(log(FID)) ~
  scale(Year) +
  # scale(log(SD)) +
  scale(log(FlockSize)) +
  scale(log(BodyMass)) +
  scale(sin(rad)) + scale(cos(rad)) +
  # scale(Day)+
  scale(Temp) +
  scale(Human) +
  (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
# (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = dh[Country == "Hungary"], REML = FALSE
)
est_hhx <- est_out(hhx, "Hungary: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
est_hhx[, control_for_starting_distance := "no"]

# PL
phs <- lmer(scale(log(FID)) ~
  scale(log(SD)) +
  scale(log(FlockSize)) +
  scale(log(BodyMass)) +
  scale(sin(rad)) + scale(cos(rad)) +
  # scale(Day)+
  scale(Temp) +
  scale(Human) +
  (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = dh[Country == "Poland"], REML = FALSE
)
est_phs <- est_out(phs, "Poland: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)")
est_phs[, control_for_starting_distance := "yes"]

phx <- lmer(scale(log(FID)) ~
  # scale(log(SD)) +
  scale(log(FlockSize)) +
  scale(log(BodyMass)) +
  scale(sin(rad)) + scale(cos(rad)) +
  # scale(Day)+
  scale(Temp) +
  scale(Human) +
  (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year),
# (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = dh[Country == "Poland"], REML = FALSE
)
est_phx <- est_out(phx, "Poland: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)")
est_phx[, control_for_starting_distance := "no"]

# combine
est_mhs[, Country := "All\n(mixed model)"]
est_mhx[, Country := "All\n(mixed model)"]
est_chs[, Country := "Czechia"]
est_chx[, Country := "Czechia"]
est_hhs[, Country := "Hungary"]
est_hhx[, Country := "Hungary"]
est_phs[, Country := "Poland"]
est_phx[, Country := "Poland"]
est_fhs[, Country := "Finland"]
est_fhx[, Country := "Finland"]

oh <- rbind(
  est_mhs, est_mhx,
  est_chs, est_chx,
  est_hhs, est_hhx,
  est_phs, est_phx,
  est_fhs, est_fhx
)
save(oh, file = here::here("Data/dat_est_numberOFhumans_rev.Rdata")) # save(oh, file = here::here("Data/dat_est_humans_rev_no_2018-CZ.Rdata"))

# estimatees for table
mhs_out <- m_out(name = "Table S1a - full a", dep = "Escape distance", model = mhs, nsim = 5000, R2 = FALSE)
mhx_out <- m_out(name = "Table S1a - full b", dep = "Escape distance", model = mhx, nsim = 5000, R2 = FALSE)
chs_out <- m_out(name = "Table S1a - CZ a", dep = "Escape distance", model = chs, nsim = 5000, R2 = FALSE)
chx_out <- m_out(name = "Table S1a - CZ b", dep = "Escape distance", model = chx, nsim = 5000, R2 = FALSE)
fhs_out <- m_out(name = "Table S1a - FI a", dep = "Escape distance", model = fhs, nsim = 5000, R2 = FALSE)
fhx_out <- m_out(name = "Table S1a - FI b", dep = "Escape distance", model = fhx, nsim = 5000, R2 = FALSE)
hhs_out <- m_out(name = "Table S1a - HU a", dep = "Escape distance", model = hhs, nsim = 5000, R2 = FALSE)
hhx_out <- m_out(name = "Table S1a - HU b", dep = "Escape distance", model = hhx, nsim = 5000, R2 = FALSE)
phs_out <- m_out(name = "Table S1a - PL a", dep = "Escape distance", model = phs, nsim = 5000, R2 = FALSE)
phx_out <- m_out(name = "Table S1a - PL b", dep = "Escape distancey", model = phx, nsim = 5000, R2 = FALSE)

out_FID_h <- rbind(mhs_out, mhx_out, fhs_out, fhx_out, phs_out, phx_out, chs_out, chx_out, hhs_out, hhx_out, fill = TRUE)
out_FID_h[is.na(out_FID_h)] <- ""
out_FID_h[, effect := gsub("scale\\(Human\\)", "# of humans present", effect)]
out_FID_h[, effect := gsub("scale\\(Year\\)", "year", effect)]
out_FID_h[, effect := gsub("scale\\(log\\(SD\\)", "starting distance (ln)", effect)]
out_FID_h[, effect := gsub("scale\\(Temp\\)", "temperature", effect)]
out_FID_h[, effect := gsub("scale\\(log\\(FlockSize\\)\\)", "flock size (ln)", effect)]
out_FID_h[, effect := gsub("scale\\(log\\(BodyMass\\)\\)", "body mass (ln)", effect)]
out_FID_h[, effect := gsub("scale\\(sin\\(rad\\)\\)", "time (sine of radians)", effect)]
out_FID_h[, effect := gsub("scale\\(cos\\(rad\\)\\)", "time (cosine of radians)", effect)]
out_FID_h[, effect := gsub("Species", "species", effect)]
out_FID_h[, effect := gsub("Country", "country", effect)]
out_FID_h[, effect := gsub("sp_day_year", "species within day & year", effect)]
out_FID_h[, effect := gsub("IDLocality", "site", effect)]
out_FID_h[, effect := gsub("sp_loc", "species within site", effect)]
out_FID_h[, effect := gsub("country # of humans present", "# of humans (slope) | country", effect)]
fwrite(file = here::here("Outputs/Table_S1a.csv"), out_FID_h)

# GOOGLE
## predictions
 # full
 mgs <- lmer(scale(log(FID)) ~
     scale(Year) +
     scale(log(SD)) +
     scale(log(FlockSize)) +
     scale(log(BodyMass)) +
     scale(sin(rad)) + scale(cos(rad)) +
     # scale(Day)+
     scale(Temp) +
     scale(parks_percent_change_from_baseline) +
     (1|weekday) +  (1| genus) + (1 | Species) + (1 | sp_day_year) + 
     (scale(parks_percent_change_from_baseline)| Country) + (1 | IDLocality) + (1 | sp_loc),
 data = ss, REML = FALSE,
 control = lmerControl(
     optimizer = "optimx", optCtrl = list(method = "nlminb")
 )
 )
 est_mgs <- est_out(mgs, "ALL: (1|weekday) + (1|genus) + (1|Species)  + (1|sp_day_year) + (scale(Google)|Country) + (1|IDLocality) +(1|sp_loc)")
 est_mgs[, control_for_starting_distance := "yes"]

 mgx <- lmer(scale(log(FID)) ~
    scale(Year) +
    #scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(parks_percent_change_from_baseline) +
     (1|weekday) + (1| genus) + (1 | Species) + (1 | sp_day_year) + 
     (scale(parks_percent_change_from_baseline) | Country) + (1 | IDLocality) + (1 | sp_loc),
 data = ss, REML = FALSE,
 control = lmerControl(
    optimizer = "optimx", optCtrl = list(method = "nlminb")
)
)
est_mgx <- est_out(mgx, "ALL: (1|weekday) + (1|genus) + (1|Species) + (1|sp_day_year) + (scale(Google)|Country) + (1|IDLocality) +(1|sp_loc)")
est_mgx[, control_for_starting_distance := "no"]


# CZ - singular fits only due to genera estimated as zero (removing it changes no results)
  cgs <- lmer(scale(log(FID)) ~
    #scale(Year) +
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(parks_percent_change_from_baseline) +
    (1|weekday) + (1 | genus) + (1| Species) + (1 | sp_day_year) + (1| IDLocality) + (1|sp_loc),
    #(1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
    data = ss[Country == "Czechia"], REML = FALSE
    )
  est_cgs <- est_out(cgs, "Czechia:(1|weekday)+(1|genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
  est_cgs[, control_for_starting_distance := "yes"]
  
  cgx <- lmer(scale(log(FID)) ~
    #scale(Year) +
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(parks_percent_change_from_baseline) +
    (1|weekday) + (1 | genus) + (1| Species) + (1 | sp_day_year) + (1| IDLocality) + (1|sp_loc),
    #(1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
    data = ss[Country == "Czechia"], REML = FALSE
    )
  est_cgx <- est_out(cgx, "Czechia: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
  est_cgx[, control_for_starting_distance := "no"]

# FI
fgs <- lmer(scale(log(FID)) ~
    scale(Year) +
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(parks_percent_change_from_baseline) +
    (1|weekday) + (1 | genus) + (1| Species) + (1 | sp_day_year) + (1| IDLocality) + (1|sp_loc),
data = ss[Country == "Finland"], REML = FALSE
)
est_fgs <- est_out(fgs, "Finland: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
est_fgs[, control_for_starting_distance := "yes"]

fgx <- lmer(scale(log(FID)) ~
    scale(Year) +
    # scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(parks_percent_change_from_baseline) +
    (1|weekday) + (1 | genus) + (1| Species) + (1 | sp_day_year) + (1| IDLocality) + (1|sp_loc),
data = ss[Country == "Finland"], REML = FALSE
)
est_fgx <- est_out(fgx, "Finland: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
est_fgx[, control_for_starting_distance := "no"]

# HU - singular fits only due to sp_loc estimated as zero (removing it changes no results)
hgs <- lmer(scale(log(FID)) ~
    scale(Year) +
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(parks_percent_change_from_baseline) +
     (1|weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = ss[Country == "Hungary"], REML = FALSE
)
est_hgs <- est_out(hgs, "Hungary: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
est_hgs[, control_for_starting_distance := "yes"]

hgx <- lmer(scale(log(FID)) ~
    scale(Year) +
    # scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(parks_percent_change_from_baseline) +
    (1|weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
# (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = ss[Country == "Hungary"], REML = FALSE
)
est_hgx <- est_out(hgx, "Hungary: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
est_hgx[, control_for_starting_distance := "no"]

# AU - singular fits only due to Year and random slope estimated as zero (removing those changes no results)
ags <- lmer(scale(log(FID)) ~
    scale(Year) +
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(parks_percent_change_from_baseline) +
      (1|weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = ss[Country == "Australia"], REML = FALSE
)
est_ags <- est_out(ags, "Australia: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
est_ags[, control_for_starting_distance := "yes"]

agx <- lmer(scale(log(FID)) ~
    scale(Year) +
    # scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(parks_percent_change_from_baseline) +
    (1|weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
# (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = ss[Country == "Australia"], , REML = FALSE,
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
)
est_agx <- est_out(agx, "Australia: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
est_agx[, control_for_starting_distance := "no"]

# PL 
pgs <- lmer(scale(log(FID)) ~
    scale(Year) +
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(parks_percent_change_from_baseline) +
    (1|weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = ss[Country == "Poland"], REML = FALSE
)
est_pgs <- est_out(pgs, "Poland: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)")
est_pgs[, control_for_starting_distance := "yes"]

pgx <- lmer(scale(log(FID)) ~
    scale(Year) +
    # scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(parks_percent_change_from_baseline) +
    (1|weekday) + (1 | genus) + (1 | Species)+ (1 | sp_day_year),
# (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = ss[Country == "Poland"], REML = FALSE
)
est_pgx <- est_out(pgx, "Poland: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)")
est_pgx[, control_for_starting_distance := "no"]

  # combine
    est_mgs[, Country := 'All\n(mixed model)']
    est_mgx[, Country := "All\n(mixed model)"]
    est_ags[, Country := "Australia"]
    est_agx[, Country := "Australia"]
    est_cgs[, Country := "Czechia"]
    est_cgx[, Country := "Czechia"]
    est_hgs[, Country := "Hungary"]
    est_hgx[, Country := "Hungary"]
    est_pgs[, Country := "Poland"]
    est_pgx[, Country := "Poland"]
    est_fgs[, Country := "Finland"]
    est_fgx[, Country := "Finland"]

    og = rbind(est_mgs, est_mgx, 
            est_ags, est_agx, 
            est_cgs, est_cgx, 
            est_hgs, est_hgx,
            est_pgs, est_pgx, 
            est_fgs, est_fgx)
    save(og, file = here::here('Data/dat_est_Google_rev.Rdata'))
  # estimatees for table
  mgs_out <- m_out(name = "Table S1b - full a", dep = "Escape distance", model = mgs, nsim = 5000, R2 = FALSE)
  mgx_out <- m_out(name = "Table S1b - full b", dep = "Escape distance", model = mgx, nsim = 5000, R2 = FALSE)
  cgs_out <- m_out(name = "Table S1b - CZ a", dep = "Escape distance", model = cgs, nsim = 5000, R2 = FALSE)
  cgx_out <- m_out(name = "Table S1b - CZ b", dep = "Escape distance", model = cgx, nsim = 5000, R2 = FALSE)
  fgs_out <- m_out(name = "Table S1b - FI a", dep = "Escape distance", model = fgs, nsim = 5000, R2 = FALSE)
  fgx_out <- m_out(name = "Table S1b - FI b", dep = "Escape distance", model = fgx, nsim = 5000, R2 = FALSE)
  hgs_out <- m_out(name = "Table S1b - HU a", dep = "Escape distance", model = hgs, nsim = 5000, R2 = FALSE)
  hgx_out <- m_out(name = "Table S1b - HU b", dep = "Escape distance", model = hgx, nsim = 5000, R2 = FALSE)
  ags_out <- m_out(name = "Table S1b - AU a", dep = "Escape distance", model = ags, nsim = 5000, R2 = FALSE)
  agx_out <- m_out(name = "Table S1b - AU b", dep = "Escape distancey", model = agx, nsim = 5000, R2 = FALSE)
  pgs_out <- m_out(name = "Table S1b - PL a", dep = "Escape distance", model = pgs, nsim = 5000, R2 = FALSE)
  pgx_out <- m_out(name = "Table S1b - PL b", dep = "Escape distancey", model = pgx, nsim = 5000, R2 = FALSE)

  out_FID_g <- rbind(mgs_out, mgx_out, fgs_out, fgx_out, pgs_out, pgx_out, cgs_out, cgx_out, hgs_out, hgx_out, ags_out, agx_out, fill = TRUE)
  out_FID_g[is.na(out_FID_g)] <- ""
  out_FID_g[, effect := gsub("scale\\(parks_percent_change_from_baseline\\)", "Google Mobility", effect)]
  out_FID_g[, effect := gsub("scale\\(Year\\)", "year", effect)]
  out_FID_g[, effect := gsub("scale\\(log\\(SD\\)", "starting distance (ln)", effect)]
  out_FID_g[, effect := gsub("scale\\(Temp\\)", "temperature", effect)]
  out_FID_g[, effect := gsub("scale\\(log\\(FlockSize\\)\\)", "flock size (ln)", effect)]
  out_FID_g[, effect := gsub("scale\\(log\\(BodyMass\\)\\)", "body mass (ln)", effect)]
  out_FID_g[, effect := gsub("scale\\(sin\\(rad\\)\\)", "time (sine of radians)", effect)]
  out_FID_g[, effect := gsub("scale\\(cos\\(rad\\)\\)", "time (cosine of radians)", effect)]
  out_FID_g[, effect := gsub("Species", "species", effect)]
  out_FID_g[, effect := gsub("Country", "country", effect)]
  out_FID_g[, effect := gsub("sp_day_year", "species within day & year", effect)]
  out_FID_g[, effect := gsub("IDLocality", "site", effect)]
  out_FID_g[, effect := gsub("sp_loc", "species within site", effect)]
  out_FID_g[, effect := gsub("Google Mobility (slope) | country", "country Google Mobility", effect)]
  fwrite(file = here::here("Outputs/Table_S1b.csv"), out_FID_g)

# STRINGENCY
# predictions for fig and table for stringency
 # full
 mss <- lmer(scale(log(FID)) ~
     scale(Year) +
     scale(log(SD)) +
     scale(log(FlockSize)) +
     scale(log(BodyMass)) +
     scale(sin(rad)) + scale(cos(rad)) +
     # scale(Day)+
     scale(Temp) +
     scale(StringencyIndex) +
     (1|weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (scale(StringencyIndex) | Country) + (1 | IDLocality) + (1 | sp_loc),
 data = s, REML = FALSE,
 control = lmerControl(
     optimizer = "optimx", optCtrl = list(method = "nlminb")
 )
 )
 est_mss <- est_out(mss, "ALL: (1|weekday) + (1|genus) + (1|Species) + (1|sp_day_year) + (scale(StringencyIndex)|Country) + (1|IDLocality) +(1|sp_loc)")
 est_mss[, control_for_starting_distance := "yes"]

 msx <- lmer(scale(log(FID)) ~
    scale(Year) +
    #scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(StringencyIndex) +
    (1|weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (scale(StringencyIndex) | Country) + (1| IDLocality) + (1 | sp_loc),
 data = s, REML = FALSE,
 control = lmerControl(
    optimizer = "optimx", optCtrl = list(method = "nlminb")
)
)
est_msx <- est_out(msx, "ALL: (1|weekday) + (1|genus) + (1|Species) + (1|sp_day_year) + (scale(StringencyIndex)|Country) + (1|IDLocality) +(1|sp_loc)")
est_msx[, control_for_starting_distance := "no"]


# CZ - singular fits only due to genera estimated as zero (removing it changes no results)
  css <- lmer(scale(log(FID)) ~
    #scale(Year) +
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(StringencyIndex) +
    (1|weekday) + (1 | genus) + (1| Species) + (1 | sp_day_year) + (1| IDLocality) + (1|sp_loc),
    #(1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
    data = s[Country == "Czechia"], REML = FALSE
    )
  est_css <- est_out(css, "Czechia: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
  est_css[, control_for_starting_distance := "yes"]
  
  csx <- lmer(scale(log(FID)) ~
    #scale(Year) +
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(StringencyIndex) +
    (1|weekday) + (1 | genus) + (1| Species) + (1 | sp_day_year) + (1| IDLocality) + (1|sp_loc),
    #(1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
    data = s[Country == "Czechia"], REML = FALSE
    )
  est_csx <- est_out(csx, "Czechia: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
  est_csx[, control_for_starting_distance := "no"]

# FI
fss <- lmer(scale(log(FID)) ~
    scale(Year) +
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(StringencyIndex) +
    (1|weekday) + (1 | genus) + (1| Species) + (1 | sp_day_year) + (1| IDLocality) + (1|sp_loc),
data = s[Country == "Finland"], REML = FALSE
)
est_fss <- est_out(fss, "Finland: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)+(scale(StringencyIndex)|IDLocality)+(1|sp_loc)")
est_fss[, control_for_starting_distance := "yes"]

fsx <- lmer(scale(log(FID)) ~
    scale(Year) +
    # scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(StringencyIndex) +
    (1|weekday) + (1 | genus) + (1| Species) + (1 | sp_day_year) + (1| IDLocality) + (1|sp_loc),
data = s[Country == "Finland"], REML = FALSE
)
est_fsx <- est_out(fsx, "Finland: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)+(scale(StringencyIndex)|IDLocality)+(1|sp_loc)")
est_fsx[, control_for_starting_distance := "no"]

# HU - singular fits only due to sp_loc estimated as zero (removing it changes no results)
hss <- lmer(scale(log(FID)) ~
    scale(Year) +
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(StringencyIndex) +
    (1|weekday) +  (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = s[Country == "Hungary"], REML = FALSE
)
est_hss <- est_out(hss, "Hungary: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
est_hss[, control_for_starting_distance := "yes"]

hsx <- lmer(scale(log(FID)) ~
    scale(Year) +
    # scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(StringencyIndex) +
    (1|weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
# (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = s[Country == "Hungary"], REML = FALSE
)
est_hsx <- est_out(hsx, "Hungary: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
est_hsx[, control_for_starting_distance := "no"]

# AU - singular fits only due to Year and random slope estimated as zero (removing those changes no results)
ass <- lmer(scale(log(FID)) ~
    scale(Year) +
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(StringencyIndex) +
     (1|weekday) +  (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = s[Country == "Australia"], REML = FALSE
)
est_ass <- est_out(ass, "Australia: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
est_ass[, control_for_starting_distance := "yes"]

asx <- lmer(scale(log(FID)) ~
    scale(Year) +
    # scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(StringencyIndex) +
    (1|weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality) + (1 | sp_loc),
# (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = s[Country == "Australia"], , REML = FALSE,
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
)
est_asx <- est_out(asx, "Australia: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)+(1|IDLocality)+(1|sp_loc)")
est_asx[, control_for_starting_distance := "no"]

# PL 
pss <- lmer(scale(log(FID)) ~
    scale(Year) +
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(StringencyIndex) +
    (1|weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = s[Country == "Poland"], REML = FALSE
)
est_pss <- est_out(pss, "Poland: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)")
est_pss[, control_for_starting_distance := "yes"]

psx <- lmer(scale(log(FID)) ~
    scale(Year) +
    # scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(StringencyIndex) +
    (1|weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year),
# (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = s[Country == "Poland"], REML = FALSE
)
est_psx <- est_out(psx, "Poland: (1|weekday) + (1|genus)+(1|Species)+(1|sp_day_year)")
est_psx[, control_for_starting_distance := "no"]

  # combine
    est_mss[, Country := 'All\n(mixed model)']
    est_msx[, Country := "All\n(mixed model)"]
    est_ass[, Country := "Australia"]
    est_asx[, Country := "Australia"]
    est_css[, Country := "Czechia"]
    est_csx[, Country := "Czechia"]
    est_hss[, Country := "Hungary"]
    est_hsx[, Country := "Hungary"]
    est_pss[, Country := "Poland"]
    est_psx[, Country := "Poland"]
    est_fss[, Country := "Finland"]
    est_fsx[, Country := "Finland"]

    os = rbind(est_mss, est_msx, 
            est_ass, est_asx, 
            est_css, est_csx, 
            est_hss, est_hsx,
            est_pss, est_psx, 
            est_fss, est_fsx)
    save(os, file = here::here('Data/dat_est_Stringency_rev.Rdata'))

# estimates for table
  mss_out <- m_out(name = "Table S1c - full a", dep = "Escape distance", model = mss, nsim = 5000, R2 = FALSE)
  msx_out <- m_out(name = "Table S1c  - full b", dep = "Escape distance", model = msx, nsim = 5000, R2 = FALSE)
  css_out <- m_out(name = "Table S1c  - CZ a", dep = "Escape distance", model = css, nsim = 5000, R2 = FALSE)
  csx_out <- m_out(name = "Table S1c  - CZ b", dep = "Escape distance", model = csx, nsim = 5000, R2 = FALSE)
  fss_out <- m_out(name = "Table S1c  - FI a", dep = "Escape distance", model = fss, nsim = 5000, R2 = FALSE)
  fsx_out <- m_out(name = "Table S1c  - FI b", dep = "Escape distance", model = fsx, nsim = 5000, R2 = FALSE)
  hss_out <- m_out(name = "Table S1c  - HU a", dep = "Escape distance", model = hss, nsim = 5000, R2 = FALSE)
  hsx_out <- m_out(name = "Table S1c  - HU b", dep = "Escape distance", model = hsx, nsim = 5000, R2 = FALSE)
  ass_out <- m_out(name = "Table S1c  - AU a", dep = "Escape distance", model = ass, nsim = 5000, R2 = FALSE)
  asx_out <- m_out(name = "Table S1c  - AU b", dep = "Escape distancey", model = asx, nsim = 5000, R2 = FALSE)
  pss_out <- m_out(name = "Table S1c  - PL a", dep = "Escape distance", model = pss, nsim = 5000, R2 = FALSE)
  psx_out <- m_out(name = "Table S1c  - PL b", dep = "Escape distancey", model = psx, nsim = 5000, R2 = FALSE)

  out_FID_s <- rbind(mss_out, msx_out, fss_out, fsx_out, pss_out, psx_out, css_out, csx_out, hss_out, hsx_out, ass_out, asx_out, fill = TRUE)
  out_FID_s[is.na(out_FID_s)] <- ""
  out_FID_s$R2_mar = out_FID_s$R2_con = NULL
  out_FID_s[, effect := gsub("scale\\(StringencyIndex\\)", "Stringency index", effect)]
  out_FID_s[, effect := gsub("scale\\(Year\\)", "year", effect)]
  out_FID_s[, effect := gsub("scale\\(log\\(SD\\)", "starting distance (ln)", effect)]
  out_FID_s[, effect := gsub("scale\\(Temp\\)", "temperature", effect)]
  out_FID_s[, effect := gsub("scale\\(log\\(FlockSize\\)\\)", "flock size (ln)", effect)]
  out_FID_s[, effect := gsub("scale\\(log\\(BodyMass\\)\\)", "body mass (ln)", effect)]
  out_FID_s[, effect := gsub("scale\\(sin\\(rad\\)\\)", "time (sine of radians)", effect)]
  out_FID_s[, effect := gsub("scale\\(cos\\(rad\\)\\)", "time (cosine of radians)", effect)]
  out_FID_s[, effect := gsub("Species", "species", effect)]
  out_FID_s[, effect := gsub("Country", "country", effect)]
  out_FID_s[, effect := gsub("sp_day_year", "species within day & year", effect)]
  out_FID_s[, effect := gsub("IDLocality", "site", effect)]
  out_FID_s[, effect := gsub("sp_loc", "species within site", effect)]
  out_FID_s[type == "random" & grepl("Stringency index", effect, fixed = TRUE), effect := paste("Stringency index (slope) |", gsub(" Stringency index", "", effect))]
  fwrite(file = here::here("Outputs/Table_S1c.csv"), out_FID_s)

# PERIOD
## predictions
### full
ms <- lmer(scale(log(FID)) ~
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(Covid) +
    (1 | Year) + (1 | weekday) + (1| genus) + (1| Species) + (1 | sp_day_year) + (scale(Covid) | Country) + (scale(Covid) | IDLocality) + (1 | sp_loc),
data = d, REML = FALSE,
control <- lmerControl(
    optimizer = "optimx", optCtrl = list(method = "nlminb")
)
)
est_ms <- est_out(ms, "All: (1|Year) + (1|genus) + (1|Species) + (1|sp_day_year) + (scale(Covid)|Country) + (scale(Covid)|IDLocality) + (1|sp_loc)")
est_ms[, control_for_starting_distance := "yes"]

mx <- lmer(scale(log(FID)) ~
    #scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(Covid) +
    (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (scale(Covid) | Country) + (scale(Covid) | IDLocality) + (1 | sp_loc),
    data = d, REML = FALSE,
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
# (Covid|IDLocality) +
est_mx <- est_out(mx, "All: (1|Year) + (1|genus) + (1|Species) + (1|sp_day_year) + (scale(Covid)|Country) + (scale(Covid)|IDLocality) + (1|sp_loc)")
est_mx[, control_for_starting_distance := "no"]

### CZ - singular fits only due to genera estimated as zero (removing it changes no results)
  cs <- lmer(scale(log(FID)) ~
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(Covid) +
    (1 | Year) + (1 | weekday) + (1 | genus) + (1| Species) + (1 | sp_day_year) + (scale(Covid)| IDLocality) + (1|sp_loc),
    #(1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
    data = d[Country == "Czechia"], REML = FALSE
    )
  est_cs <- est_out(cs, "Czechia: (1|weekday)+(1|genus)+(1|Species)+(1|sp_day_year)+(scale(Covid)|IDLocality)+(1|sp_loc)")
  est_cs[, control_for_starting_distance := "yes"]
  
  cx <- lmer(scale(log(FID)) ~
    #scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(Covid) +
    (1 | Year) + (1 | weekday) + (1 | genus) + (1| Species) + (1 | sp_day_year) + (scale(Covid)| IDLocality) + (1|sp_loc),
    #(1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
    data = d[Country == "Czechia"], REML = FALSE
    )
  est_cx <- est_out(cx, "Czechia: (1|weekday)+(1|genus)+(1|Species)+(1|sp_day_year)+(scale(Covid)|IDLocality)+(1|sp_loc)")
  est_cx[, control_for_starting_distance := "no"]

### FI
fs <- lmer(scale(log(FID)) ~
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(Covid) +
    (1 | Year) + (1 | weekday) + (1 | genus) + (1| Species) + (1 | sp_day_year) + (scale(Covid)| IDLocality) + (1|sp_loc),
data = d[Country == "Finland"], REML = FALSE
)
est_fs <- est_out(fs, "Finland: (1|weekday)+(1|genus)+(1|Species)+(1|sp_day_year)+(scale(Covid)|IDLocality)+(1|sp_loc)")
est_fs[, control_for_starting_distance := "yes"]

fx <- lmer(scale(log(FID)) ~
    # scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(Covid) +
     (1 | Year) + (1 | weekday) + (1 | genus) + (1| Species) + (1 | sp_day_year) + (scale(Covid)| IDLocality) + (1|sp_loc),
data = d[Country == "Finland"], REML = FALSE
)
est_fx <- est_out(fx, "Finland: (1|weekday)+(1|genus)+(1|Species)+(1|sp_day_year)+(scale(Covid)|IDLocality)+(1|sp_loc)")
est_fx[, control_for_starting_distance := "no"]

### HU - singular fits only due to sp_loc estimated as zero (removing it changes no results)
hs <- lmer(scale(log(FID)) ~
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(Covid) +
    (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (scale(Covid) | IDLocality) + (1 | sp_loc),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = d[Country == "Hungary"], REML = FALSE
)
est_hs <- est_out(hs, "Hungary: (1|weekday)+(1|genus)+(1|Species)+(1|sp_day_year)+(scale(Covid)|IDLocality)+(1|sp_loc)")
est_hs[, control_for_starting_distance := "yes"]

hx <- lmer(scale(log(FID)) ~
    # scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(Covid) +
    (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (scale(Covid) | IDLocality) + (1 | sp_loc),
# (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = d[Country == "Hungary"], REML = FALSE
)
est_hx <- est_out(hx, "Hungary: (1|weekday)+(1|genus)+(1|Species)+(1|sp_day_year)+(scale(Covid)|IDLocality)+(1|sp_loc)")
est_hx[, control_for_starting_distance := "no"]

### AU - singular fits only due to Year and random slope estimated as zero (removing those changes no results)
as <- lmer(scale(log(FID)) ~
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(Covid) +
     (1 | Year) +(1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (scale(Covid) | IDLocality) + (1 | sp_loc),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = d[Country == "Australia"], REML = FALSE
)
est_as <- est_out(as, "Australia: (1|weekday)+(1|genus)+(1|Species)+(1|sp_day_year)+(scale(Covid)|IDLocality)+(1|sp_loc)")
est_as[, control_for_starting_distance := "yes"]

ax <- lmer(scale(log(FID)) ~
    # scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(Covid) +
    (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (scale(Covid) | IDLocality) + (1 | sp_loc),
# (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = d[Country == "Australia"], , REML = FALSE,
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
)
est_ax <- est_out(ax, "Australia: (1|weekday)+(1|genus)+(1|Species)+(1|sp_day_year)+(scale(Covid)|IDLocality)+(1|sp_loc)")
est_ax[, control_for_starting_distance := "no"]

### PL 
ps <- lmer(scale(log(FID)) ~
    scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(Covid) +
    (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = d[Country == "Poland"], REML = FALSE
)
est_ps <- est_out(ps, "Poland: (1|weekday)+(1|genus)+(1|Species)+(1|sp_day_year)")
est_ps[, control_for_starting_distance := "yes"]

px <- lmer(scale(log(FID)) ~
    # scale(log(SD)) +
    scale(log(FlockSize)) +
    scale(log(BodyMass)) +
    scale(sin(rad)) + scale(cos(rad)) +
    # scale(Day)+
    scale(Temp) +
    scale(Covid) +
    (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year),
# (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = d[Country == "Poland"], REML = FALSE
)
est_px <- est_out(px, "Poland: (1|weekday)+(1|genus)+(1|Species)+(1|sp_day_year)")
est_px[, control_for_starting_distance := "no"]

  # combine
    est_ms[, Country := 'All\n(mixed model)']
    est_mx[, Country := "All\n(mixed model)"]
    est_as[, Country := "Australia"]
    est_ax[, Country := "Australia"]
    est_cs[, Country := "Czechia"]
    est_cx[, Country := "Czechia"]
    est_hs[, Country := "Hungary"]
    est_hx[, Country := "Hungary"]
    est_ps[, Country := "Poland"]
    est_px[, Country := "Poland"]
    est_fs[, Country := "Finland"]
    est_fx[, Country := "Finland"]

    o = rbind(est_ms, est_mx, 
            est_as, est_ax, 
            est_cs, est_cx, 
            est_hs, est_hx,
            est_ps, est_px, 
            est_fs, est_fx)
    save(o, file = here::here('Data/dat_est_rev.Rdata'))
 
  #AIC(hs, hx)
  #AIC(cs, cx)
  #0AIC(ps, px)
  #AIC(as, ax)
  #AIC(fs, fx)

# supplementary Table S1d - fid_period  
  ms_out <- m_out(name = "Table S1d - full a", dep = "Escape distance", model = ms, nsim = 5000)
  mx_out <- m_out(name = "Table S1d - full b", dep = "Escape distance", model = mx, nsim = 5000)
  cs_out <- m_out(name = "Table S1d - CZ a", dep = "Escape distancey", model = cs, nsim = 5000)
  cx_out <- m_out(name = "Table S1d - CZ b", dep = "Escape distance", model = cx, nsim = 5000)
  fs_out <- m_out(name = "Table S1d - FI a", dep = "Escape distance", model = fs, nsim = 5000)
  fx_out <- m_out(name = "Table S1d - FI b", dep = "Escape distance", model = fx, nsim = 5000)
  hs_out <- m_out(name = "Table S1d - HU a", dep = "Escape distance", model = hs, nsim = 5000)
  hx_out <- m_out(name = "Table S1d - HU b", dep = "Escape distance", model = hx, nsim = 5000)
  as_out <- m_out(name = "Table S1d - AU a", dep = "Escape distance", model = as, nsim = 5000)
  ax_out <- m_out(name = "Table S1d - AU b", dep = "Escape distance", model = ax, nsim = 5000)
  ps_out <- m_out(name = "Table S1d - PL a", dep = "Escape distance", model = ps, nsim = 5000)
  px_out <- m_out(name = "Table S1d - PL b", dep = "Escape distance", model = px, nsim = 5000)
  
  out_FID_c <- rbind(fs_out, fx_out, ps_out, px_out, cs_out, cx_out, hs_out, hx_out, as_out, ax_out, fill = TRUE)
  out_FID_c[is.na(out_FID_c)] <- ""
  out_FID_c$R2_mar = out_FID_c$R2_con = NULL
  out_FID_c[, effect := gsub("scale\\(Covid\\)", "Period", effect)]
  out_FID_c[, effect := gsub("scale\\(Year\\)", "year", effect)]
  out_FID_c[, effect := gsub("scale\\(log\\(SD\\)\\)", "starting distance (ln)", effect)]
  out_FID_c[, effect := gsub("scale\\(Temp\\)", "temperature", effect)]
  out_FID_c[, effect := gsub("scale\\(log\\(FlockSize\\)\\)", "flock size (ln)", effect)]
  out_FID_c[, effect := gsub("scale\\(log\\(BodyMass\\)\\)", "body mass (ln)", effect)]
  out_FID_c[, effect := gsub("scale\\(sin\\(rad\\)\\)", "time (sine of radians)", effect)]
  out_FID_c[, effect := gsub("scale\\(cos\\(rad\\)\\)", "time (cosine of radians)", effect)]
  out_FID_c[, effect := gsub("Species", "species", effect)]
  out_FID_c[, effect := gsub("Year", "year", effect)]
  out_FID_c[, effect := gsub("sp_day_year", "species within day & year", effect)]
  out_FID_c[, effect := gsub("IDLocality", "site", effect)]
  out_FID_c[, effect := gsub("sp_loc", "species within site", effect)]
  out_FID_c[, effect := gsub("site Period", "Period (slope) | site", effect)]
  fwrite(file = here::here("Outputs/Table_S1d.csv"), out_FID_c)

#+ F_1, fig.width=16/2.5, fig.height = 6.5/2.5
width_ <- .5 # spacing between error bars
col_ <- c("#357EBDFF", "#D43F3AFF", "#46B8DAFF", "#5CB85CFF", "#EEA236FF", "#9632B8FF", "#9632B8FF")[7:1] # col_ <- c("#357EBDFF", "#D43F3AFF", "#46B8DAFF", "#5CB85CFF", "#9632B8FF", "#9632B8FF")[6:1] # "#D43F3AFF", "#B8B8B8FF" # show_col(col_)
tx_legend_tit = 6
## a)
load(here::here("Data/dat_est_numberOFhumans_rev.Rdata")) # load(here::here("Data/dat_est_humans_rev_no_2018-CZ.Rdata"))
oh[predictor %in% c("scale(Human)"), predictor := "# of humans"]
oho <- oh[predictor %in% c("# of humans")]
oho[, N := as.numeric(sub(".*N = ", "", model))]
# add meta-analytical mean
oho_s <- oho[control_for_starting_distance == "yes"]
met <- summary(meta.summaries(d = oho_s$estimate, se = oho_s$sd, method = "fixed", weights = oho_s$N))$summci
oho_met <- data.table(predictor = "# of humans", estimate = met[2], lwr = met[1], upr = met[3], sd = NA, model = NA, control_for_starting_distance = "yes", Country = "Combined\n(meta-analytical)", N = NA)

oho_sx <- oho[control_for_starting_distance == "no"]
metx <- summary(meta.summaries(d = oho_sx$estimate, se = oho_sx$sd, method = "fixed", weights = oho_sx$N))$summci
oho_metx <- data.table(predictor = "# of humans", estimate = metx[2], lwr = metx[1], upr = metx[3], sd = NA, model = NA, control_for_starting_distance = "no", Country = "Combined\n(meta-analytical)", N = NA)

oho_au <- data.table(predictor = "# of humans", estimate = NA, lwr = NA, upr = NA, sd = NA, model = NA, control_for_starting_distance = c("no", "yes"), Country = "Australia", N = 0) # dummy for australia
oho <- rbind(oho, oho_met, oho_metx, oho_au)

oho[, Country := factor(Country, levels = rev(c("Finland", "Poland", "Czechia", "Hungary", "Australia", "Combined\n(meta-analytical)", "All\n(mixed model)")))]

# prepare for adding N
oho[, N := as.character(N)]
oho[control_for_starting_distance == "no" | is.na(N), N := ""]
oho[, n_pos := .1]

g1a <-
  ggplot(oho, aes(x = estimate, y = Country, col = Country, shape = control_for_starting_distance)) +
  geom_text(aes(x = n_pos, label = N), vjust = 1, size = 1.75, position = ggstance::position_dodgev(width_)) +
  geom_vline(xintercept = 0, color = "grey", linetype = "dotted") +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0, position = ggstance::position_dodgev(width_)) +
  # geom_point(position = ggstance::position_dodgev(.6)) +
  geom_point(position = position_dodge(width = width_), bg = "white", size = 1.1) +
  # scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
  # scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) +
  # geom_text( aes(x = n_pos,label = N), vjust = 0, size = 1.75, position = ggstance::position_dodgev(width_))+ # 3 positions for 3 bars
  # annotate("text", x=log10(3), y=85, label= "Used", col = "grey30", size = 2.5)+

  scale_shape_manual(name = "Controlled for\nstarting distance", guide = guide_legend(reverse = TRUE), values = c(21, 19)) +
  # scale_color_jama(guide = "none")+ #, palette = 'light'
  scale_color_manual(guide = "none", values = col_) + # guide_legend(reverse = TRUE)
  scale_x_continuous(breaks = round(seq(-0.3, 0.1, by = 0.1), 1), limits = c(-0.3, 0.132)) +
  ylab("") +
  xlab("Hourly human levels\n[# of humans]") +
  labs(tag = "a") +
  # coord_cartesian(xlim = c(-.15, .15)) +
  # scale_x_continuous(breaks = round(seq(-.15, .15, by = 0.05),2)) +
  theme_bw() +
  theme(
    plot.tag = element_text(size = 7,face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = tx_legend_tit),
    legend.text = element_text(size = 6),
    # legend.spacing.y = unit(0.1, 'cm'),
    legend.key.height = unit(0.5, "line"),
    legend.margin = margin(0, 0, 0, 0),
    # legend.position=c(0.5,1.6),
    plot.title = element_text(color = "grey", size = 7),
    plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r = 1, unit = "pt"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = ax_lines, size = 0.25),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(colour = ax_lines, size = 0.25),
    # axis.text.x = element_text()
    axis.ticks.length = unit(1, "pt"),
    axis.text.x = element_text(, size = 6),
    axis.text.y = element_text(colour = "black", size = 7),
    axis.title = element_text(size = 7)
  )

## b
load(here::here("Data/dat_est_Google_rev.Rdata"))
og[predictor %in% c("scale(parks_percent_change_from_baseline)"), predictor := "Google Mobility"]
ogo <- og[predictor %in% c("Google Mobility")]
ogo[, N := as.numeric(sub(".*N = ", "", model))]
# add meta-analytical mean
ogo_s <- ogo[control_for_starting_distance == "yes"]
met <- summary(meta.summaries(d = ogo_s$estimate, se = ogo_s$sd, method = "fixed", weights = ogo_s$N))$summci
ogo_met <- data.table(predictor = "Period", estimate = met[2], lwr = met[1], upr = met[3], sd = NA, model = NA, control_for_starting_distance = "yes", Country = "Combined\n(meta-analytical)", N = NA)

ogo_sx <- ogo[control_for_starting_distance == "no"]
metx <- summary(meta.summaries(d = ogo_sx$estimate, se = ogo_sx$sd, method = "fixed", weights = ogo_sx$N))$summci
ogo_metx <- data.table(predictor = "Period", estimate = metx[2], lwr = metx[1], upr = metx[3], sd = NA, model = NA, control_for_starting_distance = "no", Country = "Combined\n(meta-analytical)", N = NA)

ogo <- rbind(ogo, ogo_met, ogo_metx)

ogo[, Country := factor(Country, levels = rev(c("Finland", "Poland", "Czechia", "Hungary", "Australia", "Combined\n(meta-analytical)", "All\n(mixed model)")))]

# prepare for adding N
ogo[, N := as.character(N)]
ogo[control_for_starting_distance == "no" | is.na(N), N := ""]
ogo[, n_pos := .15]

g1b <-
  ggplot(ogo, aes(x = estimate, y = Country, col = Country, shape = control_for_starting_distance)) +
  geom_vline(xintercept = 0, color = "grey", linetype = "dotted") +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0, position = ggstance::position_dodgev(width_)) +
  # geom_point(position = ggstance::position_dodgev(.6)) +
  geom_point(position = position_dodge(width = width_), bg = "white", size = 1.1) +
  # scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
  # scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) +
  # geom_text( aes(x = n_pos,label = N), vjust = 0, size = 1.75, position = ggstance::position_dodgev(width_))+ # 3 positions for 3 bars
  # annotate("text", x=log10(3), y=85, label= "Used", col = "grey30", size = 2.5)+
  geom_text(aes(x = n_pos, label = N), vjust = 1, size = 1.75, position = ggstance::position_dodgev(width_)) +
  scale_shape_manual(name = "Controlled for\nstarting distance", guide = guide_legend(reverse = TRUE), values = c(21, 19)) +
  # scale_color_jama(guide = "none")+ #, palette = 'light'
  scale_color_manual(guide = "none", values = col_) + # guide_legend(reverse = TRUE)
  scale_x_continuous(breaks = round(seq(-0.3, 0.2, by = 0.1), 1)) +
  ylab("") +
  xlab("Daily human levels\n[Google Mobility]") +
  labs(tag = "b") +
  # coord_cartesian(xlim = c(-.15, .15)) +
  # scale_x_continuous(breaks = round(seq(-.15, .15, by = 0.05),2)) +
  theme_bw() +
  theme(
    plot.tag = element_text(size = 7,face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = tx_legend_tit),
    legend.text = element_text(size = 6),
    # legend.spacing.y = unit(0.1, 'cm'),
    legend.key.height = unit(0.5, "line"),
    legend.margin = margin(0, 0, 0, 0),
    # legend.position=c(0.5,1.6),
    plot.title = element_text(color = "grey", size = 7),
    plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r = 1, unit = "pt"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = ax_lines, size = 0.25),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(colour = ax_lines, size = 0.25),
    # axis.text.x = element_text()
    axis.ticks.length = unit(1, "pt"),
    axis.text.x = element_text(, size = 6),
    axis.text.y = element_text(colour = "black", size = 7),
    axis.title = element_text(size = 7)
  )

## c)
load(here::here("Data/dat_est_Stringency_rev.Rdata"))
os[predictor %in% c("scale(StringencyIndex)"), predictor := "Stringency Index"]
oso <- os[predictor %in% c("Stringency Index")]
oso[, N := as.numeric(sub(".*N = ", "", model))]
# add meta-analytical mean
oso_s <- oso[control_for_starting_distance == "yes"]
met <- summary(meta.summaries(d = oso_s$estimate, se = oso_s$sd, method = "fixed", weights = oso_s$N))$summci
oso_met <- data.table(predictor = "Period", estimate = met[2], lwr = met[1], upr = met[3], sd = NA, model = NA, control_for_starting_distance = "yes", Country = "Combined\n(meta-analytical)", N = NA)

oso_sx <- oso[control_for_starting_distance == "no"]
metx <- summary(meta.summaries(d = oso_sx$estimate, se = oso_sx$sd, method = "fixed", weights = oso_sx$N))$summci
oso_metx <- data.table(predictor = "Period", estimate = metx[2], lwr = metx[1], upr = metx[3], sd = NA, model = NA, control_for_starting_distance = "no", Country = "Combined\n(meta-analytical)", N = NA)

oso <- rbind(oso, oso_met, oso_metx)

oso[, Country := factor(Country, levels = rev(c("Finland", "Poland", "Czechia", "Hungary", "Australia", "Combined\n(meta-analytical)", "All\n(mixed model)")))]

# prepare for adding N
oso[, N := as.character(N)]
oso[control_for_starting_distance == "no" | is.na(N), N := ""]
oso[, n_pos := 0.35]

g1c <-
  ggplot(oso, aes(x = estimate, y = Country, col = Country, shape = control_for_starting_distance)) +
  geom_vline(xintercept = 0, color = "grey", linetype = "dotted") +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0, position = ggstance::position_dodgev(width_)) +
  # geom_point(position = ggstance::position_dodgev(.6)) +
  geom_point(position = position_dodge(width = width_), bg = "white", size = 1.1) +
  # scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
  # scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) +
  # geom_text( aes(x = n_pos,label = N), vjust = 0, size = 1.75, position = ggstance::position_dodgev(width_))+ # 3 positions for 3 bars
  # annotate("text", x=log10(3), y=85, label= "Used", col = "grey30", size = 2.5)+
  geom_text(aes(x = n_pos, label = N), vjust = 1, size = 1.75, position = ggstance::position_dodgev(width_)) +
  scale_shape_manual(name = "Controlled for\nstarting distance", guide = guide_legend(reverse = TRUE), values = c(21, 19)) +
  # scale_color_jama(guide = "none")+ #, palette = 'light'
  scale_color_manual(guide = "none", values = col_) + # guide_legend(reverse = TRUE)
  scale_x_continuous(breaks = round(seq(-0.2, 0.4, by = 0.2), 1), , expand = c(0.04, 0.04)) +
  ylab("") +
  xlab("Weekly human levels\n[Stringency index]") +
  labs(tag = "c") +
  # coord_cartesian(xlim = c(-.15, .15)) +
  # scale_x_continuous(breaks = round(seq(-.15, .15, by = 0.05),2)) +
  theme_bw() +
  theme(
    legend.position = "right",
    plot.tag = element_text(size = 7,face = "bold"),
    legend.title = element_text(size = tx_legend_tit),
    legend.text = element_text(size = 6),
    # legend.spacing.y = unit(0.1, 'cm'),
    legend.key.height = unit(0.5, "line"),
    legend.margin = margin(0, 0, 0, 0),
    # legend.position=c(0.5,1.6),
    plot.title = element_text(color = "grey", size = 7),
    plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r = 1, unit = "pt"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = ax_lines, size = 0.25),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(colour = ax_lines, size = 0.25),
    # axis.text.x = element_text()
    axis.ticks.length = unit(1, "pt"),
    axis.text.x = element_text(, size = 6),
    axis.text.y = element_text(colour = "black", size = 7),
    axis.title = element_text(size = 7)
  )

## d)
load(here::here("Data/dat_est_rev.Rdata"))
o[predictor %in% c("scale(Covid)"), predictor := "Period"]
oo <- o[predictor %in% c("Period")]
oo[, N:=as.numeric(sub('.*N = ', '', model))]
# add meta-analytical mean
  oo_s = oo[control_for_starting_distance == 'yes']
  met = summary(meta.summaries(d = oo_s$estimate, se = oo_s$sd, method = "fixed", weights = oo_s$N))$summci
  oo_met = data.table(predictor = "Period", estimate = met[2], lwr = met[1], upr = met[3], sd = NA, model = NA, control_for_starting_distance = "yes", Country = "Combined\n(meta-analytical)", N = NA)

  oo_sx = oo[control_for_starting_distance == "no"]
  metx = summary(meta.summaries(d = oo_sx$estimate, se = oo_sx$sd, method = "fixed", weights = oo_sx$N))$summci
  oo_metx = data.table(predictor = "Period", estimate = metx[2], lwr = metx[1], upr = metx[3], sd = NA, model = NA, control_for_starting_distance = "no", Country = "Combined\n(meta-analytical)", N = NA)
  
  oo = rbind(oo, oo_met, oo_metx)
    
oo[, Country := factor(Country, levels = rev(c("Finland", "Poland", "Czechia", "Hungary", "Australia", "Combined\n(meta-analytical)", "All\n(mixed model)")))]

# prepare for adding N
oo[control_for_starting_distance == "no" | is.na(N), N := ""]
oo[, n_pos := 0.9]

g1d = 
ggplot(oo, aes(x = estimate, y = Country, col = Country, shape = control_for_starting_distance)) +
    geom_vline(xintercept = 0, color = "grey", linetype = "dotted") +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0, position = ggstance::position_dodgev(width_)) +
    # geom_point(position = ggstance::position_dodgev(.6)) +
    geom_point(position = position_dodge(width = width_), bg = "white", size = 1.1) +
    # scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
    # scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) +
    # geom_text( aes(x = n_pos,label = N), vjust = 0, size = 1.75, position = ggstance::position_dodgev(width_))+ # 3 positions for 3 bars
    # annotate("text", x=log10(3), y=85, label= "Used", col = "grey30", size = 2.5)+
    geom_text( aes(x = n_pos,label = N), vjust = 1, size = 1.75, position = ggstance::position_dodgev(width_))+
    scale_shape_manual(name = "Controlled for\nstarting distance", guide = guide_legend(reverse = TRUE), values = c(21, 19)) +
    #scale_color_jama(guide = "none")+ #, palette = 'light'
    scale_color_manual(guide = "none", values = col_) + #guide_legend(reverse = TRUE)
    scale_x_continuous(breaks = round(seq(-0.6, 0.9, by = 0.3), 1)) +
    ylab("") +
    xlab("Yearly human levels\n[Period]") +
    labs(tag = "d") +
    # coord_cartesian(xlim = c(-.15, .15)) +
    # scale_x_continuous(breaks = round(seq(-.15, .15, by = 0.05),2)) +
    theme_bw() +
    theme(
        legend.position = "right",
        plot.tag = element_text(size = 7,face = "bold"),
        legend.title = element_text(size = tx_legend_tit),
        legend.text = element_text(size = 6),
        # legend.spacing.y = unit(0.1, 'cm'),
        legend.key.height = unit(0.5, "line"),
        legend.margin = margin(0, 0, 0, 0),
        # legend.position=c(0.5,1.6),
        plot.title = element_text(color = "grey", size = 7),
        plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r = 1, unit = "pt"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = ax_lines, size = 0.25),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(colour = ax_lines, size = 0.25),
        # axis.text.x = element_text()
        axis.ticks.length = unit(1, "pt"),
        axis.text.x = element_text(, size = 6),
        axis.text.y = element_text(colour = "black", size = 7),
        axis.title = element_text(size = 7)
    )

## combine

# gg_S6 <- ggarrange(
#  gs6a, NULL, gs6b, widths = c(1,0.05,1),
#  nrow = 1, common.legend = TRUE, legend = 'right'
# )

gg_f1 <- g1a + g1b + g1c + g1d + 
  plot_layout(guides = "collect", nrow = 1) +
  plot_annotation(caption = 
    "Standardized effect sizes on flight initiation distance", 
    theme = theme(plot.caption = element_text(size = 7, hjust=0.5))
  )
    
### adjust tag
gg_f1[[1]] <- gg_f1[[1]] + theme(
  # plot.margin = margin(r = 2, l = 2),
  plot.tag.position = c(.44, .98)
)

### remove y-axis from second subplot & adjust tag
gg_f1[[2]] <- gg_f1[[2]] + theme(
  plot.margin = margin(r = 2, l = 2),
  plot.tag.position = c(.035, .98),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y = element_blank()
)

### remove y-axis from third subplot & adjust tag
gg_f1[[3]] <- gg_f1[[3]] + theme(
  plot.tag.position = c(.035, .98),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y = element_blank()
)

### remove y-axis from fourth subplot & adjust tag
gg_f1[[4]] <- gg_f1[[4]] + theme(
  plot.tag.position = c(.035, .98),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title.y = element_blank()
)

if (save_plot == TRUE) {
  ggsave(here::here("Outputs/Fig_1_width-160mm_r2.png"), gg_f1, width = 16, height = 6.5, unit = "cm", dpi = 600)
}

gg_f1

#' <a name="F_1">
#' **Figure 1</a> | Avian tolerance towards humans across four temporal scales.** Avian tolerance according to human levels during (a) the escape distance trial (hourly scale) as number of humans within a 50 meter radius, (b) the day of the escape tiral as proxied by Google Mobility, (c) week of the escape trial as proxied by the stringency of governmental measures, and (d) the period of the escape trial, i.e. before vs during the COVID-19 shutdowns. (a-d) The dots with horizontal lines represent estimated standardised effect size and their 95% confidence intervals, the numbers sample sizes. For the **countries** and **'All'**, the effect sizes and 95% confidence intervals come from the joint posterior distribution of 5000 simulated values generated by the *sim* function from the *arm* package ([Gelman et al 2016](https://www.cambridge.org/highereducation/books/data-analysis-using-regression-and-multilevel-hierarchical-models/32A29531C7FD730C3A68951A17C9D983#overview)) using the mixed model outputs controlled for starting distance of the observer (`r knitr::asis_output("\U25CF")`) or not (`r knitr::asis_output("\U25CB")`; Table [S1a](T_S1a), [S1b](#T_S4), [S1c](#T_S5) and [S1d](#T_S6)). The models were further controlled for year (in case of Period fitted as a random intercept), flock size (ln-transfomred), body size (ln-transformed), temperature (also a proxy for a day within the breeding season: r~Pearson~ = `r round(cor(d$Temp,d$Day_),2)`; Fig. [S6](#F_S6)), and time of a day, as well as for the non-independence of data points by fitting random intercepts of weekday, genus, species, species at a given day and year, country (in "All" - a global mixed model), site, and species within a site, while fitting, in case of "All" (global model on all countries) Number of humans, Google Mobility, Stringency index, or Period as a random slope within country (i.e. allowing for different effect in each country) and in case of Period models (country-specific or 'All') also within site (i.e. allowing for different Period effect at each site). Fitting Number of humans, Google Mobility, Stringency index or Period as random slope at other random intercepts, simplifying the random effects structures, or using conservative datasets produces similar results (Fig. [S1](#F_S1), Table [S2a](#T_S2a), [S2b](#T_S2b), [S2c](#T_S2c) & [S2d](#T_S2d)). The multicollinearity was small as correlations between predictors were weak (Fig. [S6](#F_S6)). For the **'Combined'**, the estimates and 95% confidence intervals represent the meta-analytical means based on the country-specific estimates and their standard deviation (from the country-specific models), and sample size per country. Number of humans data are from before and during COVID-19 shutdowns. Google Mobility data were sourced from [Google Mobility Reports](https://www.google.com/covid19/mobility) and Stringency index data from [Hale et al. 2021](https://ourworldindata.org/covid-stringency-index), and their analyses were restricted to data collected in the Period during the COVID-19 shutdowns. Note that the effect sizes are small and estimates tend to centre around zero.
#' 
#' <a name="T_S1a">
#' **Table S1a | Escape distance in relation to # of humans**</a>
out_FID_h = fread(here::here("Outputs/Table_S1a.csv"))
out_FID_h$R2_mar <- out_FID_h$R2_con <- out_FID_h$error_structure <- out_FID_h$response <- NULL
out_FID_h[model != "", model := c(
  "All countries", "All countries, without starting distance",
  "Finland", "Finland, without starting distance",
  "Poland", "Poland, without starting distance",
  "Czechia", "Czechia, without starting distance",
  "Hungary", "Hungary, without starting distance"
)]
setnames(out_FID_h, old = c("estimate_r", "lwr_r", "upr_r"), new = c("estimate", "lower", "upper"))
out_FID_h[is.na(out_FID_h)] <- ""
out_FID_h %>%
  kbl() %>%
  kable_paper("hover", full_width = F) %>%
  scroll_box(width = "100%", height = "650px")
#' For model descriptions see legend of Fig. [1](#F_1) above, for descriptions of variables Methods of the [paper](https://doi.org/10.1101/2022.07.15.500232).
#'
#' <a name="T_S1b">
#' **Table S1b | Escape distance in relations to Google Mobility**</a>
out_FID_g <- fread(here::here("Outputs/Table_S1b.csv"))
out_FID_g$R2_mar = out_FID_g$R2_con = out_FID_g$error_structure = out_FID_g$response = NULL
out_FID_g[model != "", model := c(
  "All countries", "All countries, without starting distance",
  "Finland", "Finland, without starting distance",
  "Poland", "Poland, without starting distance",
  "Czechia", "Czechia, without starting distance",
  "Hungary", "Hungary, without starting distance",
  "Australia", "Australia, without starting distance"
)]
setnames(out_FID_g, old = c("estimate_r", "lwr_r", "upr_r"), new = c("estimate", "lower", "upper"))
out_FID_g[is.na(out_FID_g)] <- ""
out_FID_g %>%
  kbl() %>%
  kable_paper("hover", full_width = F) %>%
  scroll_box(width = "100%", height = "650px")
#' For model descriptions see legend of Fig. [1](#F_1) above, for descriptions of variables Methods of the [paper](https://doi.org/10.1101/2022.07.15.500232).
#'
#' <a name="T_S1c">
#' **Table S1c | Escape distance in relation to Stringency index**</a>
out_FID_s <- fread(here::here("Outputs/Table_S1c.csv"))
out_FID_s$error_structure = out_FID_s$response = NULL
out_FID_s[model != "", model := c(
  "All countries", "All countries, without starting distance",
  "Finland", "Finland, without starting distance",
  "Poland", "Poland, without starting distance",
  "Czechia", "Czechia, without starting distance",
  "Hungary", "Hungary, without starting distance",
  "Australia", "Australia, without starting distance"
)]
setnames(out_FID_s, old = c("estimate_r", "lwr_r", "upr_r"), new = c("estimate", "lower", "upper"))
out_FID_s[is.na(out_FID_s)] <- ""
out_FID_s %>%
  kbl() %>%
  kable_paper("hover", full_width = F) %>%
  scroll_box(width = "100%", height = "650px")
#' For model descriptions see legend of Fig. [1](#F_1) above, for descriptions of variables Methods of the [paper](https://doi.org/10.1101/2022.07.15.500232).
#'
#' <a name="T_S1d"> 
#' **Table S1d | Escape distance in relations to Period**</a>
out_FID_c <- fread(here::here("Outputs/Table_S1d.csv"))
out_FID_c$error_structure = out_FID_c$response = NULL
out_FID_c[model!="", model:=c('Finland', 'Finland, without starting distance', 
                              'Poland', 'Poland, without starting distance',
                              'Czechia', 'Czechia, without starting distance',
                              'Hungary', 'Hungary, without starting distance',
                              'Australia', 'Australia, without starting distance')]
setnames(out_FID_c, old = c("estimate_r", "lwr_r", "upr_r"), new = c("estimate", "lower", "upper"))
out_FID_c[is.na(out_FID_c)] <- ""
out_FID_c %>%
  kbl() %>%
  kable_paper("hover", full_width = F)  %>%
  scroll_box(width = "100%", height = "650px")
#' For model descriptions see legend of Fig. [1](#F_1) above, for descriptions of variables Methods of the [paper](https://doi.org/10.1101/2022.07.15.500232).
#' 
#' ***
#' 
#' ### Alternative models give similar results
#+ Fig_S1_alt_mod, fig.width=11, fig.height = 8
  # width = 30, height = 20, units = "cm"
  # prepare estimates Period
    # 01a all data, main text
      m1a <- lmer(scale(log(FID)) ~
        scale(log(SD)) +
        scale(log(FlockSize)) +
        scale(log(BodyMass)) +
        scale(sin(rad)) + scale(cos(rad)) +
        # scale(Day)+
        scale(Temp) +
        scale(Covid) +
        (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (scale(Covid) | Country) + (scale(Covid) | IDLocality) + (1 | sp_loc),
      data = d, REML = FALSE,
       control = lmerControl(
         optimizer = "optimx", optCtrl = list(method = "nlminb")
       )
      )
     est_m1a = est_out(m1a, "01a) (1|Year) + (1|weekday) + (1|genus) + (1|Species) + (1|sp_day_year) + (scale(Covid)|Country) + (scale(Covid)|IDLocality) + (1|sp_loc)")
    # 01b all data, all random slopes - singularity 
      m1b=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) + (1 | weekday) + (scale(Covid)|genus)+(scale(Covid)|Species)+(1|sp_day_year) + (scale(Covid)|Country) +(scale(Covid)|IDLocality) + (scale(Covid)|sp_loc),
                  data = d, REML = FALSE
                )
                  # # (Covid|IDLocality) +
      est_m1b = est_out(m1b, "01b) (1|Year) + (1|weekday) + (scale(Covid)|genus) + (scale(Covid)|Species) + (1|sp_day_year) + (scale(Covid)|Country) + (scale(Covid)|IDLocality) + (scale(Covid)|sp_loc)")
    # 01c all data, all random slopes, but some without cor to avoid singularity
      m1c=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  #(1|Year) +(0+Covid|genus)+(0+Covid|Species)+(1|sp_day_year) + (Covid|Country) + (0+Covid|IDLocality) +(Covid|sp_loc)
                  (1|Year)+ (1|weekday) + (0+scale(Covid)|genus)+(0+scale(Covid)|Species)+(1|sp_day_year) + (scale(Covid)|Country) +(scale(Covid)|IDLocality) + (scale(Covid)|sp_loc),
                  data = d, REML = FALSE) # (0+scale(Covid)|genus) is zero 
      est_m1c = est_out(m1c, "01c) (1|Year) + (1|weekday) + (0+scale(Covid)|genus) + (0+scale(Covid)|Species) + (1|sp_day_year) + (scale(Covid)|Country) + (scale(Covid)|IDLocality) + (scale(Covid)|sp_loc)")
    # 01d all data, random slopes that allow for non-singular fit 
      m1d=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) + (1|weekday) +(1|genus)+(1|Species)+(1|sp_day_year) + (scale(Covid)|Country) + (scale(Covid)|IDLocality) + (scale(Covid)|sp_loc),
                  data = d, REML = FALSE)# (Covid|IDLocality) +
      #d[,res := resid(mf)]
      est_m1d = est_out(m1d, "01d) (1|Year) + (1|weekday)+ (1|genus) + (1|Species) + (1|sp_day_year) + (scale(Covid)|Country) + (scale(Covid)|IDLocality) + (scale(Covid)|sp_loc)")
    # 01e all data, random slopes that allow for non-singular fit with simple structure
      m1e=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) + (1|weekday)+(1|genus)+(1|Species)+(1|sp_day_year) + (scale(Covid)|Country) + (1|IDLocality) + (1|sp_loc),
                  data = d, REML = FALSE)# (Covid|IDLocality) +
      est_m1e = est_out(m1e, "01e) (1|Year) + (1|weekday) + (1|genus) + (1|Species) + (1|sp_day_year) + (scale(Covid)|Country) + (1|IDLocality) + (1|sp_loc)")

    # 02a) before & during > 4/species - main text 
      dx = dd[N_during>4 & N_before >4]
      #dxx = d[Species %in% dx$Species]
      #dx2 = dd[N_during>9 & N_before >9]
      m2a <- lmer(scale(log(FID)) ~
        scale(log(SD)) +
        scale(log(FlockSize)) +
        scale(log(BodyMass)) +
        scale(sin(rad)) + scale(cos(rad)) +
        # scale(Day)+
        scale(Temp) +
        scale(Covid) +
        (1 | Year) + (1| weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (scale(Covid) | Country) + (scale(Covid) | IDLocality) + (1 | sp_loc),
      data = d[Species %in% dx$Species], REML = FALSE,
       control = lmerControl(
         optimizer = "optimx", optCtrl = list(method = "nlminb")
       )
      )
      est_m2a = est_out(m2a, "02a) (1|Year) + (1|weekday) + (1|genus) + (1|Species) + (1|sp_day_year) + (scale(Covid)|Country) + (scale(Covid)|IDLocality) + (1|sp_loc); >4/species/period")
    
    # 02b) before & during > 4/species - singularity
      dx = dd[N_during > 4 & N_before > 4]
      m2b=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) + (1| weekday)+ (scale(Covid)|genus)+(scale(Covid)|Species)+(1|sp_day_year) + (scale(Covid)|Country) +(scale(Covid)|IDLocality) + (scale(Covid)|sp_loc),
                  #(1|Year) + (1|genus)+(1|Species)+(1|sp_day_year) + (Covid|Country) + (Covid|sp_loc),
                  data = d[Species %in% dx$Species],
                  REML = FALSE) # (Covid|IDLocality) +
      est_m2b = est_out(m2b, "02b) (1|Year) + (1|weekday)+ (scale(Covid)|genus)+(scale(Covid)|Species) + (1|sp_day_year) + (scale(Covid)|Country) + (scale(Covid)|IDLocality) + (scale(Covid)|sp_loc); >4/species/period")
    # 02c) before & during > 4/species, random slopes that allow for non-singular fit and simple structure
      dx = dd[N_during>4 & N_before >4]
      #dxx = d[Species %in% dx$Species]
      #dx2 = dd[N_during>9 & N_before >9]
      m2c=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year)+ (1| weekday) +(1|genus)+(1|Species)+(1|sp_day_year) + (scale(Covid)|Country) + (1|IDLocality) + (1|sp_loc),
                  #(1|Year) + (1|genus)+(1|Species)+(1|sp_day_year) + (Covid|Country) + (Covid|sp_loc),
                  data = d[Species %in% dx$Species],
                  REML = FALSE) # (Covid|IDLocality) +
      est_m2c = est_out(m2c, "02c) (1|Year) + (1|weekday)+ (1|genus) + (1|Species)+(1|sp_day_year) + (scale(Covid)|Country) + (1|IDLocality) + (1|sp_loc); >4/species/period")
    # 03a) before & during > 9/species - main text
      dx = dd[N_during > 9 & N_before > 9]
      m3a <- lmer(scale(log(FID)) ~
       scale(log(SD)) +
       scale(log(FlockSize)) +
       scale(log(BodyMass)) +
       scale(sin(rad)) + scale(cos(rad)) +
       # scale(Day)+
       scale(Temp) +
       scale(Covid) +
       (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (scale(Covid) | Country) + (scale(Covid) | IDLocality) + (1 | sp_loc),
      data = d[Species %in% dx$Species], REML = FALSE,
      )
     est_m3a = est_out(m3a, "03a) (1|Year) + (1|weekday) + (1|genus) + (1|Species) + (1|sp_day_year) + (scale(Covid)|Country) + (scale(Covid)|IDLocality) + (1|sp_loc); >9/species/period")
    # 03b) before & during > 9/species - singularity 
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
                  (1|Year) + (1| weekday) + (scale(Covid)|genus)+(scale(Covid)|Species)+(1|sp_day_year) + (scale(Covid)|Country) +(scale(Covid)|IDLocality) + (scale(Covid)|sp_loc),
                  #(1|Year) + (1|genus)+(1|Species)+(1|sp_day_year) + (Covid|Country) + (Covid|sp_loc),
                  data = d[Species %in% dx$Species],
                  REML = FALSE) # (Covid|IDLocality) +
      est_m3b = est_out(m3b, "03b) (1|Year) + (1|weekday) + (scale(Covid)|genus)+(scale(Covid)|Species) + (1|sp_day_year) + (scale(Covid)|Country) + (scale(Covid)|IDLocality) + (scale(Covid)|sp_loc); >9/species/period")
    # 03c) before & during > 9/species, random slopes that allow for non-singular fit and simple structure
      dx = dd[N_during>9 & N_before >9]
      #dxx = d[Species %in% dx$Species]
      #dx2 = dd[N_during>9 & N_before >9]
      m3c=lmer(scale(log(FID))~
                  scale(log(SD))+
                  scale(log(FlockSize))+
                  scale(log(BodyMass))+
                  scale(sin(rad)) + scale(cos(rad)) + 
                  #scale(Day)+
                  scale(Temp)+
                  scale(Covid)+
                  (1|Year) + (1| weekday) +(1|genus)+(1|Species)+(1|sp_day_year) + (scale(Covid)|Country) + (1|IDLocality) + (1|sp_loc),
                  #(1|Year) + (1|genus)+(1|Species)+(1|sp_day_year) + (Covid|Country) + (Covid|sp_loc),
                  data = d[Species %in% dx$Species],
                  REML = FALSE) # (Covid|IDLocality) +
      est_m3c = est_out(m3c, '03c) (1|Year) + (1|genus) + (1|Species)+(1|sp_day_year) + (scale(Covid)|Country) + (1|IDLocality) + (1|sp_loc); >9/species/period')  
  
  # prepare estimates Stringency 
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
    est_m01a = est_out(m01a, "01a) (1|weekday) + (1|genus) + (1|Species) + (1|sp_day_year) + (scale(StringencyIndex)|Country) + (1|IDLocality) +(1|sp_loc)")

    m01b=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(StringencyIndex)+ 
        (1|weekday) + (scale(StringencyIndex)|genus)+(1|Species)+(1|sp_day_year) + (scale(StringencyIndex)|Country) + (1|IDLocality) +(1|sp_loc),
        data = s, REML = FALSE
        )  
        # (1|Year) explains nothing - could stay 
     est_m01b = est_out(m01b, "01b) (1|weekday) + (scale(StringencyIndex)|genus) + (1|Species) + (1|sp_day_year) + (scale(StringencyIndex)|Country) + (1|IDLocality) +(1|sp_loc)")
     
     m01c=lmer(scale(log(FID))~
          scale(Year)+       
          scale(log(SD))+
          scale(log(FlockSize))+
          scale(log(BodyMass))+
          scale(sin(rad)) + scale(cos(rad)) + 
          #scale(Day)+
          scale(Temp)+
          scale(StringencyIndex)+
          (scale(StringencyIndex)|Country) + (1|IDLocality),  
          data = s, REML = FALSE) 
          # (1|Year), (1|genus), (1|sp_day_year), (1|sp_loc),  explain nothing - could stay
     est_m01c = est_out(m01c, "01c) (scale(StringencyIndex)|Country) + (1|IDLocality)")

     m02a <- lmer(scale(log(FID)) ~
       scale(Year) +
       scale(log(SD)) +
       scale(log(FlockSize)) +
       scale(log(BodyMass)) +
       scale(sin(rad)) + scale(cos(rad)) +
       # scale(Day)+
       scale(Temp) +
       scale(StringencyIndex) +
       (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (scale(StringencyIndex) | Country) + (1 | IDLocality) + (1 | sp_loc),
     data = s[Nsp>4], REML = FALSE
     )

     est_m02a = est_out(m02a, "02a) (1|weekday) + (1|genus) + (1|Species) + (1|sp_day_year) + (scale(StringencyIndex)|Country) + (1|IDLocality) +(1|sp_loc); >4/species")
    m02b=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(StringencyIndex)+ 
        (1|weekday) + (scale(StringencyIndex)|genus)+(1|Species)+(1|sp_day_year) + (scale(StringencyIndex)|Country) + (1|IDLocality) +(1|sp_loc),
        data = s[Nsp>4], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
        # (1|Year) explains nothing - could stay 
     est_m02b = est_out(m02b, "02b) (1|weekday) + (scale(StringencyIndex)|genus) + (1|Species) + (1|sp_day_year) + (scale(StringencyIndex)|Country) + (1|IDLocality) +(1|sp_loc); >4/species")
     
     m02c=lmer(scale(log(FID))~
          scale(Year)+       
          scale(log(SD))+
          scale(log(FlockSize))+
          scale(log(BodyMass))+
          scale(sin(rad)) + scale(cos(rad)) + 
          #scale(Day)+
          scale(Temp)+
          scale(StringencyIndex)+
          (scale(StringencyIndex)|Country) + (1|IDLocality),  
          data = s[Nsp>4], REML = FALSE) 
          # (1|Year), (1|genus), (1|sp_day_year), (1|sp_loc),  explain nothing - could stay
      est_m02c = est_out(m02c, "02c) (scale(StringencyIndex)|Country) + (1|IDLocality); >4/species")

      m03a <- lmer(scale(log(FID)) ~
           scale(Year) +
           scale(log(SD)) +
           scale(log(FlockSize)) +
           scale(log(BodyMass)) +
           scale(sin(rad)) + scale(cos(rad)) +
           # scale(Day)+
           scale(Temp) +
           scale(StringencyIndex) +
           (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (scale(StringencyIndex) | Country) + (1 | IDLocality) + (1 | sp_loc),
         data = s[Nsp > 9], REML = FALSE
         )

         est_m03a = est_out(m03a, "03a) (1|weekday) + (1|genus) + (1|Species) + (1|sp_day_year) + (scale(StringencyIndex)|Country) + (1|IDLocality) +(1|sp_loc); >9/species")
         m03b = lmer(scale(log(FID)) ~
           scale(Year) +
           scale(log(SD)) +
           scale(log(FlockSize)) +
           scale(log(BodyMass)) +
           scale(sin(rad)) + scale(cos(rad)) +
           # scale(Day)+
           scale(Temp) +
           scale(StringencyIndex) +
           (1 | weekday) + (scale(StringencyIndex) | genus) + (1 | Species) + (1 | sp_day_year) + (scale(StringencyIndex) | Country) + (1 | IDLocality) + (1 | sp_loc),
         data = s[Nsp > 9], REML = FALSE
         )
         est_m03b = est_out(m03b, "03b) (1|weekday) + (scale(StringencyIndex)|genus) + (1|Species) + (1|sp_day_year) + (scale(StringencyIndex)|Country) + (1|IDLocality) +(1|sp_loc); >9/species")

        m03c = lmer(scale(log(FID)) ~
           scale(Year) +
           scale(log(SD)) +
           scale(log(FlockSize)) +
           scale(log(BodyMass)) +
           scale(sin(rad)) + scale(cos(rad)) +
           # scale(Day)+
           scale(Temp) +
           scale(StringencyIndex) +
           (scale(StringencyIndex) | Country) + (1 | IDLocality),
         data = s[Nsp > 9], REML = FALSE
         )
         # (1|Year), (1|genus), (1|sp_day_year), (1|sp_loc),  explain nothing - could stay
         est_m03c = est_out(m03c, "03c) (scale(StringencyIndex)|Country) + 1|IDLocality); >9/species")

  # prepare estimates google mobility 
    mg01a=lmer(scale(log(FID))~
      scale(Year)+       
      scale(log(SD))+
      scale(log(FlockSize))+
      scale(log(BodyMass))+
      scale(sin(rad)) + scale(cos(rad)) + 
      #scale(Day)+
      scale(Temp)+
      scale(parks_percent_change_from_baseline)+
      (1|weekday) + (1|genus) +(1|Species) + (1|sp_day_year) + (scale(parks_percent_change_from_baseline)|Country) + (1|IDLocality)+(1|sp_loc),  
      data = ss, REML = FALSE
      ) 
      # (1|Year), (1|genus), (1|sp_day_year), (1|sp_loc),  explain nothing - could stay
    est_mg01a = est_out(mg01a, "01a) (1|weekday) + (1|genus) +(1|Species) + (1|sp_day_year) + (scale(parks_percent_change_from_baseline)|Country) + (1|IDLocality)+(1|sp_loc)")
    
    mg01b <- lmer(scale(log(FID)) ~
     scale(Year) +
     scale(log(SD)) +
     scale(log(FlockSize)) +
     scale(log(BodyMass)) +
     scale(sin(rad)) + scale(cos(rad)) +
     # scale(Day)+
     scale(Temp) +
     scale(parks_percent_change_from_baseline) +
     (1|weekday) + (scale(parks_percent_change_from_baseline)| genus) + (1 | Species) + (1 | sp_day_year) + 
     (scale(parks_percent_change_from_baseline)| Country) + (1 | IDLocality) + (1 | sp_loc),
     data = ss, REML = FALSE
     ) 

     est_mg01b = est_out(mg01b, "01b) (1|weeekday) + (scale(parks_percent_change_from_baseline)|genus)+(1|Species)+(1|sp_day_year) + (scale(parks_percent_change_from_baseline)|Country) + (1|IDLocality) +(1|sp_loc)")
     
    mg01c <- lmer(scale(log(FID)) ~
      scale(Year) +
      scale(log(SD)) +
      scale(log(FlockSize)) +
      scale(log(BodyMass)) +
      scale(sin(rad)) + scale(cos(rad)) +
      # scale(Day)+
      scale(Temp) +
      scale(parks_percent_change_from_baseline) +
      (scale(parks_percent_change_from_baseline) | Country) + (1 | IDLocality),
    data = ss, REML = FALSE
    )
    est_mg01c = est_out(mg01c, "01c) (scale(parks_percent_change_from_baseline)|Country) + (1|IDLocality)")

    mg02a = lmer(scale(log(FID)) ~
      scale(Year) +
      scale(log(SD)) +
      scale(log(FlockSize)) +
      scale(log(BodyMass)) +
      scale(sin(rad)) + scale(cos(rad)) +
      # scale(Day)+
      scale(Temp) +
      scale(parks_percent_change_from_baseline) +
      (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (scale(parks_percent_change_from_baseline) | Country) + (1 | IDLocality) + (1 | sp_loc),
    data = ss[Nsp>4], REML = FALSE
    )
    # (1|Year), (1|genus), (1|sp_day_year), (1|sp_loc),  explain nothing - could stay
    est_mg02a = est_out(mg02a, "02a) (1|weekday) + (1|genus) +(1|Species) + (1|sp_day_year) + (scale(parks_percent_change_from_baseline)|Country) + (1|IDLocality)+(1|sp_loc); >4/specie")

    mg02b <- lmer(scale(log(FID)) ~
      scale(Year) +
      scale(log(SD)) +
      scale(log(FlockSize)) +
      scale(log(BodyMass)) +
      scale(sin(rad)) + scale(cos(rad)) +
      # scale(Day)+
      scale(Temp) +
      scale(parks_percent_change_from_baseline) +
      (1 | weekday) + (scale(parks_percent_change_from_baseline) | genus) + (1 | Species) + (1 | sp_day_year) +
      (scale(parks_percent_change_from_baseline) | Country) + (1 | IDLocality) + (1 | sp_loc),
    data = ss[Nsp>4], REML = FALSE
    )

    est_mg02b = est_out(mg02b, "02b) (1|weeekday) + (scale(parks_percent_change_from_baseline)|genus)+(1|Species)+(1|sp_day_year) + (scale(parks_percent_change_from_baseline)|Country) + (1|IDLocality) +(1|sp_loc); >4/specie")

    mg02c <- lmer(scale(log(FID)) ~
      scale(Year) +
      scale(log(SD)) +
      scale(log(FlockSize)) +
      scale(log(BodyMass)) +
      scale(sin(rad)) + scale(cos(rad)) +
      # scale(Day)+
      scale(Temp) +
      scale(parks_percent_change_from_baseline) +
      (scale(parks_percent_change_from_baseline) | Country) + (1 | IDLocality),
    data = ss[Nsp>4], REML = FALSE
    )
    est_mg02c = est_out(mg02c, "02c) (scale(parks_percent_change_from_baseline)|Country) + (1|IDLocality); >4/specie")

    mg03a = lmer(scale(log(FID)) ~
      scale(Year) +
      scale(log(SD)) +
      scale(log(FlockSize)) +
      scale(log(BodyMass)) +
      scale(sin(rad)) + scale(cos(rad)) +
      # scale(Day)+
      scale(Temp) +
      scale(parks_percent_change_from_baseline) +
      (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (scale(parks_percent_change_from_baseline) | Country) + (1 | IDLocality) + (1 | sp_loc),
    data = ss[Nsp > 9], REML = FALSE
    )
    # (1|Year), (1|genus), (1|sp_day_year), (1|sp_loc),  explain nothing - could stay
    est_mg03a = est_out(mg03a, "03a) (1|weekday) + (1|genus) +(1|Species) + (1|sp_day_year) + (scale(parks_percent_change_from_baseline)|Country) + (1|IDLocality)+(1|sp_loc); >9/specie")

    mg03b <- lmer(scale(log(FID)) ~
      scale(Year) +
      scale(log(SD)) +
      scale(log(FlockSize)) +
      scale(log(BodyMass)) +
      scale(sin(rad)) + scale(cos(rad)) +
      # scale(Day)+
      scale(Temp) +
      scale(parks_percent_change_from_baseline) +
      (1 | weekday) + (scale(parks_percent_change_from_baseline) | genus) + (1 | Species) + (1 | sp_day_year) +
      (scale(parks_percent_change_from_baseline) | Country) + (1 | IDLocality) + (1 | sp_loc),
    data = ss[Nsp > 9], REML = FALSE
    )

    est_mg03b = est_out(mg03b, "03b) (1|weeekday) + (scale(parks_percent_change_from_baseline)|genus)+(1|Species)+(1|sp_day_year) + (scale(parks_percent_change_from_baseline)|Country) + (1|IDLocality) +(1|sp_loc); >9/specie")

    mg03c <- lmer(scale(log(FID)) ~
      scale(Year) +
      scale(log(SD)) +
      scale(log(FlockSize)) +
      scale(log(BodyMass)) +
      scale(sin(rad)) + scale(cos(rad)) +
      # scale(Day)+
      scale(Temp) +
      scale(parks_percent_change_from_baseline) +
      (scale(parks_percent_change_from_baseline) | Country) + (1 | IDLocality),
    data = ss[Nsp > 9], REML = FALSE
    )
    est_mg03c = est_out(mg03c, "03c) (scale(parks_percent_change_from_baseline)|Country) + (1|IDLocality); >9/specie")

  # prepare estimates number of humans
    mh01a=lmer(scale(log(FID))~  
      scale(log(SD))+
      scale(log(FlockSize))+
      scale(log(BodyMass))+
      scale(sin(rad)) + scale(cos(rad)) + 
      #scale(Day)+
      scale(Temp)+
      scale(Human)+
      (1 | Year) + (1|weekday) + (1|genus) +(1|Species) + (1|sp_day_year) + (scale(Human)|Country) + (1|IDLocality)+(1|sp_loc),  
      data = dh, REML = FALSE
      ) 
    est_mh01a = est_out(mh01a, "01a) (1|Year) + (1|weekday) + (1|genus) +(1|Species) + (1|sp_day_year) + (scale(Human)|Country) + (1|IDLocality)+(1|sp_loc)")
    
    mh01b <- lmer(scale(log(FID)) ~
     scale(log(SD)) +
     scale(log(FlockSize)) +
     scale(log(BodyMass)) +
     scale(sin(rad)) + scale(cos(rad)) +
     # scale(Day)+
     scale(Temp) +
     scale(Human) +
     (1|Year) + (1|weekday) + (0+scale(Human)| genus) + (1 | Species) + (0+scale(Human) | sp_day_year) + 
     (scale(Human)| Country) + (1 | IDLocality) + (1 | sp_loc),
     data = dh, REML = FALSE
     ) 
     est_mh01b = est_out(mh01b, "01b) (1|Year) + (1|weeekday) + (0 + scale(Human)|genus)+(1|Species)+(0 + scale(Human)|sp_day_year) + (scale(Human)|Country) + (1|IDLocality) +(1|sp_loc)")
     
    mh01c <- lmer(scale(log(FID)) ~
      scale(log(SD)) +
      scale(log(FlockSize)) +
      scale(log(BodyMass)) +
      scale(sin(rad)) + scale(cos(rad)) +
      # scale(Day)+
      scale(Temp) +
      scale(Human) +
      (1|Year) + (1|genus) + (1|sp_day_year) + (scale(Human) | Country) + (1 | IDLocality) +(1|sp_loc),
    data = dh, REML = FALSE
    )
    est_mh01c = est_out(mh01c, "01c) (1|Year) + (1|genus) + (1|sp_day_year) + (scale(Human)|Country) + (1|IDLocality) +(1|sp_loc)")

    mh02a = lmer(scale(log(FID)) ~
      scale(log(SD)) +
      scale(log(FlockSize)) +
      scale(log(BodyMass)) +
      scale(sin(rad)) + scale(cos(rad)) +
      # scale(Day)+
      scale(Temp) +
      scale(Human) +
      (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (scale(Human) | Country) + (1 | IDLocality) + (1 | sp_loc),
    data = dh[Nsp>4], REML = FALSE
    )
    est_mh02a = est_out(mh02a, "02a) (1|Year) + (1|weekday) + (1|genus) +(1|Species) + (1|sp_day_year) + (scale(Human)|Country) + (1|IDLocality)+(1|sp_loc); >4/specie")

    mh02b <- lmer(scale(log(FID)) ~
      scale(log(SD)) +
      scale(log(FlockSize)) +
      scale(log(BodyMass)) +
      scale(sin(rad)) + scale(cos(rad)) +
      # scale(Day)+
      scale(Temp) +
      scale(Human) +
      (1 | Year) + (1 | weekday) + (scale(Human) | genus) + (1 | Species) + (1 | sp_day_year) +
      (scale(Human) | Country) + (1 | IDLocality) + (1 | sp_loc),
    data = dh[Nsp>4], REML = FALSE
    )

    est_mh02b = est_out(mh02b, "02b) (1|Year) + (1|weeekday) + (scale(Human)|genus)+(1|Species)+(1|sp_day_year) + (scale(Human)|Country) + (1|IDLocality) +(1|sp_loc); >4/specie")

    mh02c <- lmer(scale(log(FID)) ~
      scale(log(SD)) +
      scale(log(FlockSize)) +
      scale(log(BodyMass)) +
      scale(sin(rad)) + scale(cos(rad)) +
      # scale(Day)+
      scale(Temp) +
      scale(Human) +
      (1|Year) + (1|genus) + (1|sp_day_year) + (scale(Human) | Country) + (1 | IDLocality) +(1|sp_loc),
    data = dh[Nsp>4], REML = FALSE
    )

    est_mh02c = est_out(mh02c, "02c) (1|Year) + (1|genus) + (1|sp_day_year) + (scale(Human)|Country) + (1|IDLocality) +(1|sp_loc); >4/specie")

    mh03a = lmer(scale(log(FID)) ~
      scale(log(SD)) +
      scale(log(FlockSize)) +
      scale(log(BodyMass)) +
      scale(sin(rad)) + scale(cos(rad)) +
      # scale(Day)+
      scale(Temp) +
      scale(Human) +
      (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (scale(Human) | Country) + (1 | IDLocality) + (1 | sp_loc),
    data = dh[Nsp > 9], REML = FALSE
    )
    est_mh03a = est_out(mh03a, "03a) (1|Year) + (1|weekday) + (1|genus) +(1|Species) + (1|sp_day_year) + (scale(parks_percent_change_from_baseline)|Country) + (1|IDLocality)+(1|sp_loc); >9/specie")

    mh03b <- lmer(scale(log(FID)) ~
      scale(log(SD)) +
      scale(log(FlockSize)) +
      scale(log(BodyMass)) +
      scale(sin(rad)) + scale(cos(rad)) +
      # scale(Day)+
      scale(Temp) +
      scale(Human) +
      (1 | Year) + (1 | weekday) + (scale(Human) | genus) + (1 | Species) + (1 | sp_day_year) +
      (scale(Human) | Country) + (1 | IDLocality) + (1 | sp_loc),
    data = dh[Nsp > 9], REML = FALSE
    )

    est_mh03b = est_out(mh03b, "03b) (1|Year) + (1|weeekday) + (scale(Human)|genus)+(1|Species)+(1|sp_day_year) + (scale(Human)|Country) + (1|IDLocality) +(1|sp_loc); >9/specie")

    mh03c <- lmer(scale(log(FID)) ~
      scale(log(SD)) +
      scale(log(FlockSize)) +
      scale(log(BodyMass)) +
      scale(sin(rad)) + scale(cos(rad)) +
      # scale(Day)+
      scale(Temp) +
      scale(Human) +
      (1|Year) + (1|genus) + (1|sp_day_year) + (scale(Human) | Country) + (1 | IDLocality) +(1|sp_loc),
    data = dh[Nsp > 9], REML = FALSE
    )
    est_mh03c = est_out(mh03c, "03c) (1|Year) + (1|genus) + (1|sp_day_year) + (scale(Human) | Country) + (1 | IDLocality) +(1|sp_loc); >9/specie")

    # export
      save(file = here::here("Data/dat_Fig_S2_estimates.Rdata"), 
      est_m1a, est_m1b, est_m1c, est_m1d, est_m1e, est_m2a, est_m2b, est_m2c, est_m3a, est_m3b, est_m3c, 
      est_m01a, est_m01b, est_m01c, est_m02a, est_m02b, est_m02c, est_m03a, est_m03b, est_m03c, 
      est_mg01a, est_mg01b, est_mg01c, est_mg02a, est_mg02b, est_mg02c, est_mg03a, est_mg03b, est_mg03c,
      est_mh01a, est_mh01b, est_mh01c, est_mh02a, est_mh02b, est_mh02c, est_mh03a, est_mh03b, est_mh03c) # load(here::here("Data/dat_Fig_S2_estimates.Rdata"))
     
     # prepare plot for Period
     xs = rbind(est_m1a, est_m1b, est_m1c, est_m1d, est_m1e, est_m2a, est_m2b, est_m2c, est_m3a, est_m3b, est_m3c)
     #gsub("scale\\(Covid\\)", "Period", "(scale(Covid)|Country)")
     xs[, model := gsub("scale\\(Covid\\)", "Period", xs$model)]
     xs[, model := gsub("Year", "year", xs$model)]
     xs[, model := gsub("Species", "species", xs$model)]
     xs[, model := gsub("sp_day_year", "species within day & year", xs$model)]
     xs[, model := gsub("IDLocality", "site", xs$model)]
     xs[, model := gsub("sp_loc", "species within site", xs$model)]

     gs2_ =
       ggplot(xs[predictor == "scale(Covid)"], aes(y = model, x = estimate, col = model)) +
       geom_vline(xintercept = 0, col = "grey30", lty = 3) +
       geom_errorbar(aes(xmin = lwr, xmax = upr, col = model), width = 0, position = position_dodge(width = 0.01)) +
       # ggtitle ("Sim based")+
       geom_point(position = position_dodge(width = 0.01)) +
       # scale_colour_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
       # scale_fill_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
       scale_color_viridis(discrete = TRUE, guide = guide_legend(reverse = TRUE)) +
       scale_fill_viridis(discrete = TRUE, guide = guide_legend(reverse = TRUE)) +
       scale_y_discrete(limits = rev) +
       coord_fixed(ratio = 0.05, xlim = c(-0.23, 0.15)) +
       # scale_shape(guide = guide_legend(reverse = TRUE)) +
       # scale_x_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2,2, by = 1), labels = seq(-2,2, by = 1)) +
       labs(y = NULL, x = "Period (before vs during shutdowns)", tag = "a)") + # title = "a)",Effect of Period (before/during shutdown)
       # ylim(c(0,100))+
       # coord_flip()+
       theme_bw() +
       theme(
         legend.position = "none",
         #plot.title.position = "plot",
         plot.subtitle = element_text(size = 7),
         plot.title = element_text(size = 7),
         plot.tag = element_text(size = 7),
         legend.title = element_text(size = 7),
         legend.text = element_text(size = 6),
         ## legend.spacing.y = unit(0.1, 'cm'),
         legend.key.height = unit(0.5, "line"),
         # plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r =0.5, unit =  "pt"),
         panel.grid = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.line = element_line(colour = ax_lines, size = 0.25),
         axis.line.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.ticks.x = element_line(colour = ax_lines, size = 0.25),
         axis.ticks.length = unit(1, "pt"),
         axis.text.x = element_text(colour = "black", size = 6),
         axis.text.y = element_text(colour = "black", size = 7),
         axis.title = element_text(size = 7)
       )
     #gs2
     # ggsave(here::here('Outputs/Figure_Sy.png'),g, width = 30, height =5, units = 'cm')
     
     # prepare plot for Stringency
     xs0 = rbind(est_m01a, est_m01b, est_m01c, est_m02a, est_m02b, est_m02c, est_m03a, est_m03b, est_m03c)
     xs0[, model := gsub("scale\\(StringencyIndex\\)", "Stringency Index", model)]
     xs0[, model := gsub("Year", "year", model)]
     xs0[, model := gsub("Species", "species", model)]
     xs0[, model := gsub("sp_day_year", "species within day & year", model)]
     xs0[, model := gsub("IDLocality", "site", model)]
     xs0[, model := gsub("sp_loc", "species within site", model)]
     g0_ =
       ggplot(xs0[predictor == "scale(StringencyIndex)"], aes(y = model, x = estimate, col = model)) +
       geom_vline(xintercept = 0, col = "grey30", lty = 3) +
       geom_errorbar(aes(xmin = lwr, xmax = upr, col = model), width = 0, position = position_dodge(width = 0.01)) +
       # ggtitle ("Sim based")+
       geom_point(position = position_dodge(width = 0.01)) +
       # scale_colour_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
       # scale_fill_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
       scale_color_viridis(discrete = TRUE, guide = guide_legend(reverse = TRUE)) +
       scale_fill_viridis(discrete = TRUE, guide = guide_legend(reverse = TRUE)) +
       scale_y_discrete(limits = rev) +
       coord_fixed(ratio = 0.05, xlim = c(-0.23, 0.15)) +
       # scale_shape(guide = guide_legend(reverse = TRUE)) +
       # scale_x_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2,2, by = 1), labels = seq(-2,2, by = 1)) +
       labs(y = NULL, x = "Stringency index", tag = "b)") + # title = "b) Effect of ") +
       # ylim(c(0,100))+
       # coord_flip()+
       theme_bw() +
       theme(
         legend.position = "none",
         plot.subtitle = element_text(size = 7),
         plot.title = element_text(size = 7),
         plot.tag = element_text(size = 7),
         legend.title = element_text(size = 7),
         legend.text = element_text(size = 6),
         ## legend.spacing.y = unit(0.1, 'cm'),
         legend.key.height = unit(0.5, "line"),
         # plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r =0.5, unit =  "pt"),
         panel.grid = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.line = element_line(colour = ax_lines, size = 0.25),
         axis.line.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.ticks.x = element_line(colour = ax_lines, size = 0.25),
         axis.ticks.length = unit(1, "pt"),
         axis.text.x = element_text(colour = "black", size = 6),
         axis.text.y = element_text(colour = "black", size = 7),
         axis.title = element_text(size = 7)
       )
     #g0
     # ggsave(here::here('Outputs/Figure_Sz.png'),g0, width = 30, height =5, units = 'cm')
     
     # prepare plot for Google
     xg0 = rbind(est_mg01a, est_mg01b, est_mg01c, est_mg02a, est_mg02b, est_mg02c, est_mg03a, est_mg03b, est_mg03c)
     xg0[, model := gsub("scale\\(parks_percent_change_from_baseline\\)", "Google Mobility", model)]
     xg0[, model := gsub("Year", "year", model)]
     xg0[, model := gsub("Species", "species", model)]
     xg0[, model := gsub("sp_day_year", "species within day & year", model)]
     xg0[, model := gsub("IDLocality", "site", model)]
     xg0[, model := gsub("sp_loc", "species within site", model)]
     gg0_ =
       ggplot(xg0[predictor == "scale(parks_percent_change_from_baseline)"], aes(y = model, x = estimate, col = model)) +
       geom_vline(xintercept = 0, col = "grey30", lty = 3) +
       geom_errorbar(aes(xmin = lwr, xmax = upr, col = model), width = 0, position = position_dodge(width = 0.01)) +
       # ggtitle ("Sim based")+
       geom_point(position = position_dodge(width = 0.01)) +
       # scale_colour_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
       # scale_fill_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
       scale_color_viridis(discrete = TRUE, guide = guide_legend(reverse = TRUE)) +
       scale_fill_viridis(discrete = TRUE, guide = guide_legend(reverse = TRUE)) +
       scale_y_discrete(limits = rev) +
       coord_fixed(ratio = 0.05, xlim = c(-0.23, 0.15)) +
       # scale_shape(guide = guide_legend(reverse = TRUE)) +
       # scale_x_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2,2, by = 1), labels = seq(-2,2, by = 1)) +
       labs(y = NULL, x = "Google Mobility", tag = "c)") + # title = "b) Effect of ") +
       # ylim(c(0,100))+
       # coord_flip()+
       theme_bw() +
       theme(
         legend.position = "none",
         plot.subtitle = element_text(size = 7),
         plot.title = element_text(size = 7),
         plot.tag = element_text(size = 7),
         legend.title = element_text(size = 7),
         legend.text = element_text(size = 6),
         ## legend.spacing.y = unit(0.1, 'cm'),
         legend.key.height = unit(0.5, "line"),
         # plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r =0.5, unit =  "pt"),
         panel.grid = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.line = element_line(colour = ax_lines, size = 0.25),
         axis.line.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.ticks.x = element_line(colour = ax_lines, size = 0.25),
         axis.ticks.length = unit(1, "pt"),
         axis.text.x = element_text(colour = "black", size = 6),
         axis.text.y = element_text(colour = "black", size = 7),
         axis.title = element_text(size = 7)
       )
     #gg0
     # ggsave(here::here('Outputs/Figure_Sz.png'),g0, width = 30, height =5, units = 'cm')

     # prepare plot for number of humans
     xh0 = rbind(est_mh01a, est_mh01b, est_mh01c, est_mh02a, est_mh02b, est_mh02c, est_mh03a, est_mh03b, est_mh03c)
     xh0[, model := gsub("scale\\(Human\\)", "# of humans", model)]
     xh0[, model := gsub("Year", "year", model)]
     xh0[, model := gsub("Species", "species", model)]
     xh0[, model := gsub("sp_day_year", "species within day & year", model)]
     xh0[, model := gsub("IDLocality", "site", model)]
     xh0[, model := gsub("sp_loc", "species within site", model)]
     gh0_ =
       ggplot(xh0[predictor == "scale(Human)"], aes(y = model, x = estimate, col = model)) +
       geom_vline(xintercept = 0, col = "grey30", lty = 3) +
       geom_errorbar(aes(xmin = lwr, xmax = upr, col = model), width = 0, position = position_dodge(width = 0.01)) +
       # ggtitle ("Sim based")+
       geom_point(position = position_dodge(width = 0.01)) +
       # scale_colour_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
       # scale_fill_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
       scale_color_viridis(discrete = TRUE, guide = guide_legend(reverse = TRUE)) +
       scale_fill_viridis(discrete = TRUE, guide = guide_legend(reverse = TRUE)) +
       scale_y_discrete(limits = rev) +
       coord_fixed(ratio = 0.05, xlim = c(-0.23, 0.15)) +
       # scale_shape(guide = guide_legend(reverse = TRUE)) +
       # scale_x_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2,2, by = 1), labels = seq(-2,2, by = 1)) +
       labs(y = NULL, x = "# of humans\n[Standardised effect sizes on\nflight initiation distances]", tag = "d)") + # title = "b) Effect of ") +
       # ylim(c(0,100))+
       # coord_flip()+
       theme_bw() +
       theme(
         legend.position = "none",
         plot.subtitle = element_text(size = 7),
         plot.title = element_text(size = 7),
         plot.tag = element_text(size = 7),
         legend.title = element_text(size = 7),
         legend.text = element_text(size = 6),
         ## legend.spacing.y = unit(0.1, 'cm'),
         legend.key.height = unit(0.5, "line"),
         # plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r =0.5, unit =  "pt"),
         panel.grid = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.line = element_line(colour = ax_lines, size = 0.25),
         axis.line.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.ticks.x = element_line(colour = ax_lines, size = 0.25),
         axis.ticks.length = unit(1, "pt"),
         axis.text.x = element_text(colour = "black", size = 6),
         axis.text.y = element_text(colour = "black", size = 7),
         axis.title = element_text(size = 7)
       )
     #gg0
     # ggsave(here::here('Outputs/Figure_Sz.png'),g0, width = 30, height =5, units = 'cm')

     # combine
     if(save_plot==TRUE){
     ggsave(here::here("Outputs/Fig_S1_rev_v7.png"), rbind(ggplotGrob(gs2_), ggplotGrob(g0_), ggplotGrob(gg0_), ggplotGrob(gh0_)), width = 30, height = 20, units = "cm")
     }
     
     grid::grid.draw(rbind(ggplotGrob(gs2_), ggplotGrob(g0_), ggplotGrob(gg0_), ggplotGrob(gh0_)))

#' <a name="F_S1">
#' **Figure S1 | Comparing estimates from alternative models.**</a> Changes in avian tolerance towards humans in response to (**a**) Period (before vs during the COVID-19 shutdowns) (**b**) stringency of governmental measures, (**c**) Google Mobility and (**d**) number of humans during an escape distance trial. The dots with horizontal lines represent the estimated standardised effect size and their 95% confidence intervals based on the joint posterior distribution of 5,000 simulated values generated by the *sim* function  from the *arm* R-package ([Gelman et al 2016](https://www.cambridge.org/highereducation/books/data-analysis-using-regression-and-multilevel-hierarchical-models/32A29531C7FD730C3A68951A17C9D983#overview)) from the output of the mixed models (for details see [Table S2](#T_S2a)). The name of each effect size highlights the corresponding model in [Table S2a](#T_S2a) for (a), [Table S2b](#T_S2b) for (b), [Table S2c](#T_S2c) for (c) and [Table S2d](#T_S2d) for (d), the random structure of the specific model, if applicable, the condition used to reduce the dataset, and sample size. Depicted are effect sizes based on full (01) and reduced datasets with ≥5 (02) or ≥10 observations per species and period (03). For plots of model assumptions see [below](#MS_TS1a). Note that effect sizes are small and estimates close to zero.
#'
#'   
#'  
#'
#' <a name="T_S2a">
#' **Table S2a | Alternative models on escape distance given Period** </a>
#'
    m1a_ = m_out(name = "Table S2a - 1a", dep = "Escape distance", model = m1a, nsim = 5000)
    m1b_ = m_out(name = "Table S2a - 1b", dep = "Escape distance", model = m1b, nsim = 5000)
    m1c_ = m_out(name = "Table S2a - 1c", dep = "Escape distance", model = m1c, nsim = 5000)
    m1d_ = m_out(name = "Table S2a - 1d", dep = "Escape distance", model = m1d, nsim = 5000)
    m1e_ = m_out(name = "Table S2a - 1e", dep = "Escape distance", model = m1e, nsim = 5000)

    m2a_ = m_out(name = "Table S2a - 2a", dep = "Escape distance", model = m2a, nsim = 5000)
    m2b_ = m_out(name = "Table S2a - 2b", dep = "Escape distance", model = m2b, nsim = 5000)
    m2c_ = m_out(name = "Table S2a - 2c", dep = "Escape distance", model = m2c, nsim = 5000)
    m3a_ = m_out(name = "Table S2a - 3a", dep = "Escape distance", model = m3a, nsim = 5000)
    m3b_ = m_out(name = "Table S2a - 3b", dep = "Escape distance", model = m3b, nsim = 5000)
    m3c_ = m_out(name = "Table S2a - 3c", dep = "Escape distance", model = m3c, nsim = 5000)
     
    out1 = rbind(m1a_, m1b_, m1c_, m1d_, m1e_, m2a_, m2b_, m2c_, m3a_, m3b_, m3c_, fill = TRUE)
    out1[is.na(out1)] = ""
    out1[, effect := gsub("scale\\(Covid\\)", "Period", effect)]
    out1[, effect := gsub("scale\\(Year\\)", "year", effect)]
    out1[, effect := gsub("scale\\(log\\(SD\\)\\)", "starting distance (ln)", effect)]
    out1[, effect := gsub("scale\\(Temp\\)", "temperature", effect)]
    out1[, effect := gsub("scale\\(log\\(FlockSize\\)\\)", "flock size (ln)", effect)]
    out1[, effect := gsub("scale\\(log\\(BodyMass\\)\\)", "body mass (ln)", effect)]
    out1[, effect := gsub("scale\\(sin\\(rad\\)\\)", "time (sine of radians)", effect)]
    out1[, effect := gsub("scale\\(cos\\(rad\\)\\)", "time (cosine of radians)", effect)]
    out1[, effect := gsub("Species", "species", effect)]
    out1[, effect := gsub("Country", "country", effect)]
    out1[, effect := gsub("Year", "year", effect)]
    out1[, effect := gsub("sp_day_year", "species within day & year", effect)]
    out1[, effect := gsub("IDLocality", "site", effect)]
    out1[, effect := gsub("sp_loc", "species within site", effect)]
    out1[type == "random" & grepl("Period", effect, fixed = TRUE), effect := paste("Period (slope) |", gsub(" Period", "", effect))]

    out1$R2_mar=out1$R2_con=NULL
    fwrite(file = here::here("Outputs/Table_S2a.csv"), out1)

    out1$response = out1$error_structure = NULL
    out1[model != "", model := paste0('0',substring(model, 13))]
    setnames(out1, old = c("estimate_r", "lwr_r", "upr_r"), new = c("estimate", "lower", "upper"))
    out1 %>%
       kbl() %>%
       kable_paper("hover", full_width = F) %>%
         scroll_box(width = "100%", height = "650px")
#' Note that (1a) model is the one reported in the main text.
#'
#' <br />
#'  
#' <a name="T_S2b">
#' **Table S2b | Alternative models on escape distance given Stringency index** </a>
    m01a_ = m_out(name = "Table S2b - 1a", dep = "Escape distance", model = m01a, nsim = 5000)
    m01b_ = m_out(name = "Table S2b - 1b", dep = "Escape distance", model = m01b, nsim = 5000)
    m01c_ = m_out(name = "Table S2b - 1c", dep = "Escape distance", model = m01c, nsim = 5000)

    m02a_ = m_out(name = "Table S2b - 2a", dep = "Escape distance", model = m02a, nsim = 5000)
    m02b_ = m_out(name = "Table S2b - 2b", dep = "Escape distance", model = m02b, nsim = 5000)
    m02c_ = m_out(name = "Table S2b - 2c", dep = "Escape distance", model = m02c, nsim = 5000)

    m03a_ = m_out(name = "Table S2b - 3a", dep = "Escape distance", model = m03a, nsim = 5000)
    m03b_ = m_out(name = "Table S2b - 3b", dep = "Escape distance", model = m03b, nsim = 5000)
    m03c_ = m_out(name = "Table S2b - 3c", dep = "Escape distance", model = m03c, nsim = 5000)

    out2 = rbind(m01a_, m01b_, m01c_, m02a_, m02b_, m02c_, m03a_, m03b_, m03c_, fill = TRUE)
    out2[is.na(out2)] = ""
    out2[, effect := gsub("scale\\(StringencyIndex\\)", "Stringency index", effect)]
    out2[, effect := gsub("scale\\(Year\\)", "year", effect)]
    out2[, effect := gsub("scale\\(log\\(SD\\)\\)", "starting distance (ln)", effect)]
    out2[, effect := gsub("scale\\(Temp\\)", "temperature", effect)]
    out2[, effect := gsub("scale\\(log\\(FlockSize\\)\\)", "flock size (ln)", effect)]
    out2[, effect := gsub("scale\\(log\\(BodyMass\\)\\)", "body mass (ln)", effect)]
    out2[, effect := gsub("scale\\(sin\\(rad\\)\\)", "time (sine of radians)", effect)]
    out2[, effect := gsub("scale\\(cos\\(rad\\)\\)", "time (cosine of radians)", effect)]
    out2[, effect := gsub("Species", "species", effect)]
    out2[, effect := gsub("Country", "country", effect)]
    out2[, effect := gsub("Year", "year", effect)]
    out2[, effect := gsub("sp_day_year", "species within day & year", effect)]
    out2[, effect := gsub("IDLocality", "site", effect)]
    out2[, effect := gsub("sp_loc", "species within site", effect)]
    out2[type == "random" & grepl("Stringency index", effect, fixed = TRUE), effect := paste("Stringency index (slope) |", gsub(" Stringency index", "", effect))]
    out2$R2_mar = out2$R2_con = NULL

    fwrite(file = here::here("Outputs/Table_S2b.csv"), out2)

    out2$response = out2$error_structure = NULL
    out2[model != "", model := paste0("0", substring(model, 13))]
    setnames(out2, old = c("estimate_r", "lwr_r", "upr_r"), new = c("estimate", "lower", "upper"))
    out2 %>%
      kbl() %>%
      kable_paper("hover", full_width = F) %>%
        scroll_box(width = "100%", height = "650px")
#' Note that (1a) model is the one reported in the main text.
#'   
#'   
#'
#' <a name="T_S2c">
#' **Table S2c | Alternative models on escape distance given Google Mobility**</a>
#'
    # google
     mg01a_ = m_out(name = "Table S2c - 1a", dep = "Escape distance", model = mg01a, nsim = 5000, R2 = FALSE)
     mg01b_ = m_out(name = "Table S2c - 1b", dep = "Escape distance", model = mg01b, nsim = 5000, R2 = FALSE)
     mg01c_ = m_out(name = "Table S2c - 1c", dep = "Escape distance", model = mg01c, nsim = 5000, R2 = FALSE)

     mg02a_ = m_out(name = "Table S2c - 2a", dep = "Escape distance", model = mg02a, nsim = 5000)
     mg02b_ = m_out(name = "Table S2c - 2b", dep = "Escape distance", model = mg02b, nsim = 5000)
     mg02c_ = m_out(name = "Table S2c - 2c", dep = "Escape distance", model = mg02c, nsim = 5000)

     mg03a_ = m_out(name = "Table S2c - 3a", dep = "Escape distance", model = mg03a, nsim = 5000)
     mg03b_ = m_out(name = "Table S2c - 3b", dep = "Escape distance", model = mg03b, nsim = 5000)
     mg03c_ = m_out(name = "Table S2c - 3c", dep = "Escape distancey", model = mg03c, nsim = 5000)

     out3 = rbind(mg01a_, mg01b_, mg01c_, mg02a_, mg02b_, mg02c_, mg03a_, mg03b_, mg03c_, fill = TRUE)
     out3[is.na(out3)] = ""
      out3[, effect := gsub("scale\\(parks_percent_change_from_baseline\\)", "Google Mobility", effect)]
      out3[, effect := gsub("scale\\(Year\\)", "year", effect)]
      out3[, effect := gsub("scale\\(log\\(SD\\)\\)", "starting distance (ln)", effect)]
      out3[, effect := gsub("scale\\(Temp\\)", "temperature", effect)]
      out3[, effect := gsub("scale\\(log\\(FlockSize\\)\\)", "flock size (ln)", effect)]
      out3[, effect := gsub("scale\\(log\\(BodyMass\\)\\)", "body mass (ln)", effect)]
      out3[, effect := gsub("scale\\(sin\\(rad\\)\\)", "time (sine of radians)", effect)]
      out3[, effect := gsub("scale\\(cos\\(rad\\)\\)", "time (cosine of radians)", effect)]
       out3[, effect := gsub("Species", "species", effect)]
       out3[, effect := gsub("Country", "country", effect)]
       out3[, effect := gsub("Year", "year", effect)]
       out3[, effect := gsub("sp_day_year", "species within day & year", effect)]
       out3[, effect := gsub("IDLocality", "site", effect)]
       out3[, effect := gsub("sp_loc", "species within site", effect)]
       out3[type == "random" & grepl("Google Mobility", effect, fixed = TRUE), effect := paste("Google Mobility (slope) |", gsub(" Google Mobility", "", effect))]
       out3$R2_mar = out3$R2_con = NULL
     fwrite(file = here::here("Outputs/Table_S2c.csv"), out3)
     out3$response = out3$error_structure = NULL
     out3[model != "", model := paste0("0", substring(model, 13))]
     setnames(out3, old = c("estimate_r", "lwr_r", "upr_r"), new = c("estimate", "lower", "upper"))
     out3 %>%
          kbl() %>%
          kable_paper("hover", full_width = F) %>%
            scroll_box(width = "100%", height = "650px")
#' Note that (1a) model is the one reported in the main text.
#'  
#'   
#'
#' <a name="T_S2d">
#' **Table S2d | Alternative models on escape distance given # of humans**</a>
#'
    # humans
     mh01a_ = m_out(name = "Table S2d - 1a", dep = "Escape distance", model = mh01a, nsim = 5000)
     mh01b_ = m_out(name = "Table S2d - 1b", dep = "Escape distance", model = mh01b, nsim = 5000)
     mh01c_ = m_out(name = "Table S2d - 1c", dep = "Escape distance", model = mh01c, nsim = 5000)

     mh02a_ = m_out(name = "Table S2d - 2a", dep = "Escape distance", model = mh02a, nsim = 5000)
     mh02b_ = m_out(name = "Table S2d - 2b", dep = "Escape distance", model = mh02b, nsim = 5000)
     mh02c_ = m_out(name = "Table S2d - 2c", dep = "Escape distance", model = mh02c, nsim = 5000)

     mh03a_ = m_out(name = "Table S2d - 3a", dep = "Escape distance", model = mh03a, nsim = 5000)
     mh03b_ = m_out(name = "Table S2d - 3b", dep = "Escape distance", model = mh03b, nsim = 5000)
     mh03c_ = m_out(name = "Table S2d - 3c", dep = "Escape distancey", model = mh03c, nsim = 5000)

     out4 = rbind(mh01a_, mh01b_, mh01c_, mh02a_, mh02b_, mh02c_, mh03a_, mh03b_, mh03c_, fill = TRUE)
     out4[is.na(out4)] = ""
      out4[, effect := gsub("scale\\(Human\\)", "# of humans", effect)]
      out4[, effect := gsub("scale\\(Year\\)", "year", effect)]
      out4[, effect := gsub("scale\\(log\\(SD\\)\\)", "starting distance (ln)", effect)]
      out4[, effect := gsub("scale\\(Temp\\)", "temperature", effect)]
      out4[, effect := gsub("scale\\(log\\(FlockSize\\)\\)", "flock size (ln)", effect)]
      out4[, effect := gsub("scale\\(log\\(BodyMass\\)\\)", "body mass (ln)", effect)]
      out4[, effect := gsub("scale\\(sin\\(rad\\)\\)", "time (sine of radians)", effect)]
      out4[, effect := gsub("scale\\(cos\\(rad\\)\\)", "time (cosine of radians)", effect)]
       out4[, effect := gsub("Species", "species", effect)]
       out4[, effect := gsub("Country", "country", effect)]
       out4[, effect := gsub("Year", "year", effect)]
       out4[, effect := gsub("sp_day_year", "species within day & year", effect)]
       out4[, effect := gsub("IDLocality", "site", effect)]
       out4[, effect := gsub("sp_loc", "species within site", effect)]
       out4[type == "random" & grepl("# of humans", effect, fixed = TRUE), effect := paste("# of humans (slope) |", gsub(" # of humans", "", effect))]
       out4$R2_mar = out4$R2_con = NULL
     fwrite(file = here::here("Outputs/Table_S2d.csv"), out4)
     out4$response = out4$error_structure = NULL
     out4[model != "", model := paste0("0", substring(model, 13))]
     setnames(out4, old = c("estimate_r", "lwr_r", "upr_r"), new = c("estimate", "lower", "upper"))
     out4 %>%
          kbl() %>%
          kable_paper("hover", full_width = F) %>%
            scroll_box(width = "100%", height = "650px")
#' Note that (1a) model is the one reported in the main text.
#'
#' ***
#' 
#' ### Species-site specific trends
#'
#+ Fig_2_site_human, fig.width=18*.393701, fig.height = 18*.393701
dh[, NspL := .N, by = "sp_loc"]
dhl <- dh[NspL > 9]
dhl[, Country := factor(Country, levels = rev(c("Finland", "Poland", "Czechia", "Hungary")))]
dhl[, sp_C_loc2 := paste(scinam_abb, Country, IDLocality, sep = "\n")]

col3_ <- c("#357EBDFF", "#D43F3AFF", "#46B8DAFF", "#5CB85CFF", "#EEA236FF", "#9632B8FF", "#9632B8FF")[7:1]
col_gssh <- col3_[4:7]

gssh <-
  ggplot(dhl, aes(x = Human, y = FID)) +
  # stat_smooth(method = 'rlm', se = FALSE, col = 'black', lwd = 0.5)+
  geom_point(pch = 21, alpha = 0.7, aes(fill = Country), col = "grey20") +
  stat_smooth(se = FALSE, aes(colour = "Locally weighted\nsmoothing"), lwd = 0.5) + # show_guide=TRUE
  facet_wrap(~sp_C_loc2, ncol = 9) +
  scale_fill_manual(values = col_gssh, guide = guide_legend(reverse = TRUE)) +
  scale_colour_manual(values = col_fit) + #c("grey60")) +#
  scale_x_continuous("Hourly human levels\n[# of humans]", expand = c(0, 0)) +
  scale_y_continuous("Flight initiation distance [m]", expand = c(0, 0), trans = "log10") +
  coord_cartesian(xlim = c(-3, 40)) +
  theme_bw() + theme_MB2 +
  theme(
    legend.position = c(1, 0.0625),
    legend.spacing.y = unit(-0.8, "cm"),
  )
 

ggssh <- ggplotGrob(gssh) # gg$layout$name
ggssh <- gtable_filter_remove(ggssh, name = c("axis-b-2-8", "axis-b-4-8", "axis-b-6-8" , "axis-b-8-7"), trim = FALSE) # paste0("axis-b-", c(2, 4), "-9")

if (save_plot == TRUE) {
  ggsave(here::here("Outputs/Fig_2_width-180mm.png"), ggssh, width = 18, height = 18, unit = "cm", dpi = 600)
}
grid::grid.draw(ggssh)

#' <a name="F_2">
#' **Figure 2 | Variation in avian tolerance toward hourly human numbers during the escape trial across species and sites.**</a> Dots represent single escape distance observations of species at specific sites (e.g. specific park or cemetery) during the escape trial and not corrected for other factors such as starting distance of the observer. Dot colour highlights the country. Yellow lines represent locally weighted smoothing, a non-parametric local regression fitted with the *ggplot* function of the *ggplot2* package ([Wickham 2016](https://ggplot2-book.org)), highlighting heterogenous (and usually unclear – close to zero) within- and between- species trends. Some species lack trend lines because data distribution hindered the smoothing and visualised are only data for species-site combinations with ≥10 escape distance observations and where number (#) of humans was estimated.The y-axes are on the log-scale. Panels are ordered alphabetically according to species names, then country and site identifier. Abbreviated genus names represent Col = Columba, Den = Dendrocopos, Eri = Erithacus, Fri = Fringilla, Gar = Garrulus , Lar = Larus, Lus = Luscinia, Pas = Passer, Str = Streptopelia, Stu = Sturnus, and abbreviated species name megarhync = megarhynchos.
#'
#+ Fig_3_site_google, fig.width=18*.393701, fig.height = (18+18/8)*.393701
ss[, NspL := .N, by = "sp_loc"]
ssl <- ss[NspL > 9]
ssl[, Country := factor(Country, levels = rev(c("Finland", "Poland", "Czechia", "Hungary","Australia")))]
ssl[, sp_C_loc2 := paste(scinam_abb, Country, IDLocality, sep = "\n")]
ssl[, Google_mobility := parks_percent_change_from_baseline]

col3_ <- c("#357EBDFF", "#D43F3AFF", "#46B8DAFF", "#5CB85CFF", "#EEA236FF", "#9632B8FF", "#9632B8FF")[7:1]
col3_ssl <- col3_[3:7]

gssl <-
  ggplot(ssl, aes(x = Google_mobility, y = FID)) +
  # stat_smooth(method = 'rlm', se = FALSE, col = 'black', lwd = 0.5)+
  geom_point(pch = 21, alpha = 0.7, aes(fill = Country), col = "grey20") +
  stat_smooth(se = FALSE, aes(colour = "Locally weighted\nsmoothing"), lwd = 0.5) + # show_guide=TRUE
  facet_wrap(~sp_C_loc2, ncol = 9) + # mayby do 9
  scale_fill_manual(values = col3_ssl, guide = guide_legend(reverse = TRUE)) +
  scale_colour_manual(values = col_fit) + # c("grey60")) +
  scale_x_continuous("Daily human levels\n[Google Mobility]", expand = c(0, 0)) +
  scale_y_continuous("Flight initiation distance [m]", expand = c(0, 0), trans = "log10") +
  theme_bw() + theme_MB2 +
  theme(
    legend.position = c(1, 0.0625),
    legend.spacing.y = unit(-0.95, "cm"),
  )

ggssl <- ggplotGrob(gssl) # gg$layout$name
#ggssl <- gtable_filter_remove(ggssl, name = c("axis-b-2-10", "axis-b-4-10", "axis-b-6-10"), trim = FALSE) # paste0("axis-b-", c(2, 4), "-9")

if (save_plot == TRUE) {
  ggsave(here::here("Outputs/Fig_3_width-180mm.png"), ggssl, width = 18, height = 18+18/8, unit = "cm", dpi = 600)
}
grid::grid.draw(ggssl)

#' <a name="F_3">
#' **Figure 3 | Variation in avian tolerance toward daily human levels (Google Mobility) across species and sites.**</a> Dots represent single escape distance observations of species at specific sites (e.g. specific park or cemetery) and not corrected for other factors such as starting distance of the observer. Dot colour highlights the country. Yellow lines represent locally weighted smoothing, a non-parametric local regression fitted with the *ggplot* function of the *ggplot2* package ([Wickham 2016](https://ggplot2-book.org), highlighting heterogenous (and usually unclear – close to zero) within- and between- species trends. Some species lack trend lines because data distribution hindered the smoothing and visualised are only data for species-site combinations with ≥10 escape distance observations, for which Google Mobility data were available. The y-axes are on the log-scale. Panels are ordered alphabetically according to species names, then country and site identifier. Abbreviation in the species names represent Aca. chrysorrh. = Acanthiza chrysorrhoa, Acr = Acridotheres, Anas platyrhy. = Anas platyrhynchos,  Ant. caruncula. = Anthochaera carunculata, Che = Chenonetta, Col = Columba, Den = Dendrocopos, Fri = Fringilla, Gal = Gallinula, Gra = Grallina, Gym = Gymnorhina, Lar. novaehol. = Larus novaehollandiae, Lic penicilla = Lichenostomus penicillatus, Lus. megarhync. = Luscinia megarhynchos, Man. melanocep. = Manorina melanocephala, Ocy = Ocyphaps, Pas = Passer, Phy. novaeholl. = Phylidonyris novaehollandiae, Por. = Porphyrio, Rhi = Rhipidura, Sti = Stigmatopelia, and Stu = Sturnus.
#'
#+ Fig_4_site_stringency, fig.width=18*.393701, fig.height = (18+18/8)*.393701
s[, NspL := .N, by = "sp_loc"]
sl <- s[NspL > 9]
sl[, Country := factor(Country, levels = rev(c("Finland", "Poland", "Czechia", "Hungary","Australia")))]
sl[, sp_C_loc2 := paste(scinam_abb, Country, IDLocality, sep = "\n")]

col3_ <- c("#357EBDFF", "#D43F3AFF", "#46B8DAFF", "#5CB85CFF", "#EEA236FF", "#9632B8FF", "#9632B8FF")[7:1]
col_sl <- col3_[3:7]

gsl <-
  ggplot(sl, aes(x = StringencyIndex, y = FID)) +
  # stat_smooth(method = 'rlm', se = FALSE, col = 'black', lwd = 0.5)+
  geom_point(pch = 21, alpha = 0.7, aes(fill = Country), col = "grey20") +
  stat_smooth(se = FALSE, aes(colour = "Locally weighted\nsmoothing"), lwd = 0.5) + # show_guide=TRUE
  facet_wrap(~sp_C_loc2, ncol = 9) +
  scale_fill_manual(values = col_sl, guide = guide_legend(reverse = TRUE)) +
  scale_colour_manual(values = col_fit) +
  scale_x_continuous("Weekly human levels\n[Stringency index]", expand = c(0, 0)) +
  scale_y_continuous("Flight initiation distance [m]", expand = c(0, 0), trans = "log10") +
  theme_bw() + theme_MB2 +
  theme(
      legend.position = c(1, 0.0625),
      legend.spacing.y = unit(-0.95, "cm"),
    )

ggsl <- ggplotGrob(gsl) # gg$layout$name
ggsl <- gtable_filter_remove(ggsl, name = c("axis-b-2-9", "axis-b-4-8", "axis-b-6-8", "axis-b-8-8"), trim = FALSE) # paste0("axis-b-", c(2, 4), "-9")

if (save_plot == TRUE) {
  ggsave(here::here("Outputs/Fig_4_width-180mm.png"), ggsl, width = 18, height = 18+18/8, unit = "cm", dpi = 600)
}

grid::grid.draw(ggsl)

#' <a name="F_4">
#' **Figure 4 | Variation in avian tolerance toward weekly human levels (proxied by Stringency index) across species and sites.**</a> Dots represent single escape distance observations of species at specific sites (e.g. specific park or cemetery) and not corrected for other factors such as starting distance of the observer. Dot colour highlights the country. Yellow lines represent locally weighted smoothing, a non-parametric local regression fitted with the *ggplot* function of the *ggplot2* package ([Wickham 2016](https://ggplot2-book.org)), highlighting heterogenous (and usually unclear – close to zero) within- and between- species trends. Some species lack trend lines because data distribution hindered the smoothing and visualised are only data for species-site combinations with ≥10 escape distance observations, for which Stringency index data were available. The y-axes are on the log-scale. Panels are ordered alphabetically according to species names, then country and site identifier Abbreviation in the species names represent Aca. chrysorrh. = Acanthiza chrysorrhoa, Acr = Acridotheres, Anas platyrhy. = Anas platyrhynchos,  Ant. caruncula. = Anthochaera carunculata, Che = Chenonetta, Col = Columba, Den = Dendrocopos, Fri = Fringilla, Gal = Gallinula, Gra = Grallina, Gym = Gymnorhina, Lar. novaehol. = Larus novaehollandiae, Lic penicilla = Lichenostomus penicillatus, Lus. megarhync. = Luscinia megarhynchos, Man. melanocep. = Manorina melanocephala, Ocy = Ocyphaps, Pas = Passer, Phy. novaeholl. = Phylidonyris novaehollandiae, Por. = Porphyrio, Rhi = Rhipidura, Sti = Stigmatopelia, and Stu = Sturnus.
#'
#+ Fig_5_site_period, fig.width=17*.393701, fig.height = 14*.393701
pxx = pp[N_during > 4 & N_before > 4]
dxx = d[paste(IDLocality, Species) %in% paste(pxx$IDLocality, pxx$Species)]
# table(dxx$IDLocality, dxx$Year)
#length(unique(pxx$IDLocality))
#length(unique(pxx$Species))
#sum(pxx$N_during) + sum(pxx$N_before)

dxx[, Country := factor(Country, levels = rev(c("Finland", "Poland", "Czechia", "Hungary", "Australia")))]
dxx[, sp_C_loc2 := paste(gsub("[_]", " ", Species), Country, IDLocality, sep = "\n")]
dxx[, genus := sub("_.*", "", Species)]
dxx[Covid == 0, period := "before COVID-19"]
dxx[Covid == 1, period := "during COVID-19"]

col3_ <- c("#357EBDFF", "#D43F3AFF", "#46B8DAFF", "#5CB85CFF", "#EEA236FF", "#9632B8FF", "#9632B8FF")[7:1]
col3__ <- col3_[3:7]

g_f5 = 
ggplot(dxx, aes(x = as.factor(Year), y = FID, col = Country)) +
  geom_boxplot(outlier.size = 0.5) +
  #geom_rect(data=NULL, aes(xmin = 3.5, xmax = Inf, ymin = -Inf, ymax = Inf), color = "grey60", fill = "grey60")+
  scale_y_continuous("Flight initiation distance [m]", trans = "log10") +
  annotate(geom = "rect",
             xmin = 3.5,
             xmax = +Inf,
             ymin = 0.6,
             ymax = +Inf,
             color = "grey95", fill = "grey95") +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~sp_C_loc2) +
  
  scale_x_discrete("Year", guide = guide_axis(angle = 45)) +
  # scale_color_continuous() +
  scale_colour_manual(values = col3__, guide = guide_legend(reverse = TRUE)) +
  #scale_fill_manual(values = c("white", "lightgrey")) +
  theme_bw() + theme_MB2 +
  theme(legend.position = "right",
        legend.justification = c(0, 0.94),
        legend.margin = margin(0,0,0,0, unit="cm"),
        legend.box.margin = margin(5, 5, 5, 5))

if (save_plot == TRUE) {
  ggsave(here::here("Outputs/Fig_5_width-170mm.png"), g_f5, width = 17, height = 14, units = "cm") # width = 18, height = 16, 
}

g_f5

#' <a name="F_5">
#' **Figure 5 | Between-year variation in avian tolerance toward humans across species and sites.**</a> Panels are ordered alphabetically according to species names, then country and site identifier within each country (e.g. specific park or cemetery). Boxplots outline colour highlights country, background colour indicates Period (white: before the COVID-19 shutdowns; grey: during the COVID-19 shutdowns). Boxplots depict median (horizontal line inside the box), the 25th and 75th percentiles (box) ± 1.5 times the interquartile range or the minimum/maximum value, whichever is smaller (bars), and the outliers (dots). Included are only species–site combinations with ≥5 observations per Period. The y-axes are on the log-scale. Note the lack of consistent shutdowns effects within and between species, sites and countries.
#' 
#' ***
#' 
#' ### Species trends
#+ Fig_S2_species_human, fig.width=0.393701*15.25, fig.height = .393701*6 * 19 / 9
    #fig.width= 0.393701*18, fig.height = 0.393701*5 * 19 / 9
dh[, NspC := .N, by = "sp_country"]
dhc <- dh[NspC > 9]
dhc[, sp2 := gsub(" ", "\n", sp)]
dhc[, Country := factor(Country, levels = rev(c("Finland", "Poland", "Czechia", "Hungary")))]

col3_ <- c("#357EBDFF", "#D43F3AFF", "#46B8DAFF", "#5CB85CFF", "#EEA236FF", "#9632B8FF", "#9632B8FF")[7:1]
col_gh2 <- col3_[4:7]

# two rows labels
gh2 <-
  ggplot(dhc, aes(x = Human, y = FID)) +
  # stat_smooth(method = 'rlm', se = FALSE, col = 'black', lwd = 0.5)+
  geom_point(pch = 21, alpha = 0.7, aes(fill = Country), col = "grey20") +
  stat_smooth(se = FALSE, aes(colour = "Locally weighted\nsmoothing"), lwd = 0.5) + # show_guide=TRUE
  facet_wrap(~sp2, ncol = 6) +
  scale_fill_manual(values = col_gh2, guide = guide_legend(reverse = TRUE)) +
  scale_colour_manual(values = col_fit) +
  scale_x_continuous("Hourly human levels\n[# of humans]", expand = c(0, 0)) +
  scale_y_continuous("Flight initiation distance [m]", expand = c(0, 0), trans = "log10") +
  coord_cartesian(xlim = c(-5, 70)) +
  # annotate("text", x = 1, y = 1, label = c(rep("", 52),"Observation"), hjust = -0.08, size = 1) +
  # labs(title = "Species means per sampling location")+
  # scale_colour_manual(values=c('grey60'))+
  # scale_color_manual(name = 'try', values = c('LOESS smoothed = "grey60"'))+
  theme_MB +
  theme(
    plot.title = element_text(size = 7),
    strip.background = element_blank(),
    # strip.text.x = element_text(size = 4.5, color="grey30",  margin=margin(1,1,1,1,"mm")),
    # panel.spacing = unit(1, "mm"),
    legend.position = c(1, 0.03),
    legend.justification = c(1, 0.4),
    #legend.position = c(1, 0.14),
    #legend.justification = c(1, 0),
    legend.title = element_blank(),
    # legend.spacing.y = unit(-0.78, "cm")
    # legend.spacing.y = unit(0.02, "cm") use if LOESS smooth text as legend
    legend.spacing.y = unit(-0.9, "cm")
  )

gth2 <- ggplotGrob(gh2) # gth2$layout$name
#gthx2 <- gtable_filter_remove(gth2, name = c("axis-b-2-4", "axis-b-4-4", "axis-b-6-4", "axis-b-8-4"), trim = FALSE) # paste0("axis-b-", c(2, 4), "-9")

if (save_plot == TRUE) {
  ggsave(here::here("Outputs/Fig_S2_width-152mm_humans.png"), gth2, width = 15.25, height = 6 * 19 / 9, unit = "cm", dpi = 600)
}
grid::grid.draw(gth2)
#' <a name="F_S2">
#' **Figure S2 | Species-specific avian tolerance towards humans in relation to # of humans during escape trail.**</a> Each dot represents a single escape distance observation (not corrected for other factors such as starting distance of the observer) with the # of humans within 50m radius. Dot colour highlights the country. Yellow lines represent locally weighted smoothing, a non-parametric local regression fitted with the ggplot function of ggplot2 package ([Wickham 2016](https://ggplot2-book.org)), highlighting heterogenous (and usually unclear – close to zero) within- and between- species trends. Some species lack trend lines because data distribution hindered the smoothing and visualised are only data for species with ≥10 escape distance observationshere number (#) of humans was estimated. The y-axes are on the log-scale. 
#'
#+ Fig_S3_species_google, fig.width=0.393701*15.25, fig.height = 0.393701*19
ss[, NspC := .N, by = "sp_country"]
ssc <- ss[NspC > 9]
ssc[, sp2 := gsub(" ", "\n", sp)]
ssc[, Google_mobility := parks_percent_change_from_baseline]
ssc[, Country := factor(Country, levels = rev(c("Finland", "Poland", "Czechia", "Hungary", "Australia")))]

col3_ <- c("#357EBDFF", "#D43F3AFF", "#46B8DAFF", "#5CB85CFF", "#EEA236FF", "#9632B8FF", "#9632B8FF")[7:1]
col3__ = col3_[3:7]
# two rows labels
gt2 <-
  ggplot(ssc, aes(x = Google_mobility, y = FID)) +
  # stat_smooth(method = 'rlm', se = FALSE, col = 'black', lwd = 0.5)+
  geom_point(pch = 21, alpha = 0.7, aes(fill = Country), col = "grey20") +
  stat_smooth(se = FALSE, aes(colour = "Locally weighted\nsmoothing"), lwd = 0.5) + # show_guide=TRUE
  facet_wrap(~sp2, ncol = 6) +
  scale_fill_manual(values = col3__, guide = guide_legend(reverse = TRUE)) +
  scale_colour_manual(values = col_fit) +
  scale_x_continuous("Daily human levels\n[Google Mobility]", expand = c(0, 0)) +
  scale_y_continuous("Flight initiation distance [m]", expand = c(0, 0), trans = "log10") +
  theme_MB +
  theme(
    plot.title = element_text(size = 7),
    strip.background = element_blank(),
    # strip.text.x = element_text(size = 4.5, color="grey30",  margin=margin(1,1,1,1,"mm")),
    # panel.spacing = unit(1, "mm"),
    legend.position = c(1, 0.03),
    legend.justification = c(1, 0.4),
    legend.title = element_blank(),
    # legend.spacing.y = unit(-0.78, "cm")
    # legend.spacing.y = unit(0.02, "cm") use if LOESS smooth text as legend
    legend.spacing.y = unit(-0.9, "cm")
  )

gtg2 <- ggplotGrob(gt2) # gg$layout$name
#gtgx2 <- gtable_filter_remove(gtg2, name = c("axis-b-2-9", "axis-b-5-8"), trim = FALSE) # paste0("axis-b-", c(2, 4), "-9")

if (save_plot == TRUE) {
  ggsave(here::here("Outputs/Fig_S3_width-152mm_Google.png"), gtg2, width = 15.24, height = 19, unit = "cm", dpi = 600)
}
grid::grid.draw(gtg2)
#' <a name="F_S3">
#' **Figure S3 | Species-specific avian tolerance towards humans in relation to relative daily human levels (Google Mobility).**</a> Each dot represents a single escape distance observation (not corrected for other factors such as starting distance of the observer) and a day-specific value of Google Mobility for parks in a given country. Dot colour highlights the country. Yellow lines represent locally weighted smoothing, a non-parametric local regression fitted with the ggplot function of ggplot2 package ([Wickham 2016](https://ggplot2-book.org)), highlighting heterogenous (and usually unclear – close to zero) within- and between- species trends. Some species lack trend lines because data distribution hindered the smoothing and visualised are only data for species with ≥10 escape distance observations, for which Google Mobility data were available. The y-axes are on the log-scale. 
#'
#+ Fig_S4_species_string, fig.width=6, fig.height = 7.5
ss_sp <- s[Nsp > 9]
ss_sp[, sp2 := gsub(" ", "\n", sp)]
ss_sp[, Country := factor(Country, levels = rev(c("Finland", "Poland", "Czechia", "Hungary", "Australia")))]
gs2 <-
  ggplot(ss_sp, aes(x = StringencyIndex, y = FID)) +
  # stat_smooth(method = 'rlm', se = FALSE, col = 'black', lwd = 0.5)+
  geom_point(pch = 21, alpha = 0.7, aes(fill = Country), col = "grey20") +
  stat_smooth(se = FALSE, aes(colour = "Locally weighted\nsmoothing"), lwd = 0.5) + # show_guide=TRUE
  facet_wrap(~sp2, ncol = 6) +
  scale_fill_manual(values = col3__, guide = guide_legend(reverse = TRUE)) +
  # scale_fill_viridis(discrete = TRUE, guide = guide_legend(reverse = FALSE)) +
  scale_x_continuous("Weekly human levels\n[Stringency index of governmental COVID-19 restrictions]", expand = c(0, 0)) +
  scale_y_continuous("Flight initiation distance [m]", expand = c(0, 0), trans = "log10") +
  # annotate("text", x = 1, y = 1, label = c(rep("", 52),"Observation"), hjust = -0.08, size = 1) +
  # labs(title = "Species means per sampling location")+
  scale_colour_manual(values = col_fit) +
  # scale_color_manual(name = 'try', values = c('LOESS smoothed = "grey60"'))+
  theme_MB +
  theme(
    plot.title = element_text(size = 7),
    strip.background = element_blank(),
    # strip.text.x = element_text(size = 4.5, color="grey30",  margin=margin(1,1,1,1,"mm")),
    # panel.spacing = unit(1, "mm"),
    legend.position = c(1, 0.05),
    legend.justification = c(1, 0),
    legend.title = element_blank(),
    # legend.spacing.y = unit(-0.78, "cm")
    # legend.spacing.y = unit(0.02, "cm") use if LOESS smooth text as legend
    legend.spacing.y = unit(-0.9, "cm")
  )

gsg2 <- ggplotGrob(gs2) # gg$layout$name
gsgx2 <- gtable_filter_remove(gsg2, name = paste0("axis-b-", c(2, 4), "-9"), trim = FALSE)

if (save_plot == TRUE) {
  ggsave(here::here("Outputs/Fig_S4_width-152mm.png"), gsgx2, width = 15.24, height = 19, unit = "cm", dpi = 600)
}
grid::grid.draw(gsgx2)
#' <a name="F_S4">
#' **Figure S4 | Species-specific avian tolerance towards humans in relation to weekly human levels (stringency of antipandemic governmental restrictions during COVID-19 shutdowns quantified as a Stringency index).**</a> Each dot represents a single escape distance observation (not corrected for other factors such as starting distance of the observer) and a day-specific value of governmental Stringency index in a given country. Dot colour highlights the country. Yellow lines represent locally weighted smoothing, a non-parametric local regression fitted with the ggplot function of *ggplot2* package ([Wickham 2016](https://ggplot2-book.org)), highlighting heterogenous (and usually unclear – close to zero) within- and between- species trends. Some species lack trend lines because data distribution hindered the smoothing and visualised are only data for species with ≥10 escape distance observations, for which Stringency index data were available. The y-axes are on the log-scale.
#'
#+ Fig_S5_species_compare, fig.width=5.3, fig.height = 7
dxx <- d[paste(IDLocality, Species) %in% paste(pp$IDLocality, pp$Species)]
# length(dxx[, unique(paste(sp_loc, Covid))])
m <- lm(log(FID) ~ log(SD), dxx)
dxx[, resid_FID := resid(m)]
a <- dxx[, .(mean(resid_FID), sd(resid_FID), mean(FID), .N), by = .(Country, IDLocality, genus, Species, sp_loc, Covid)]
setnames(a, old = c("V1", "V2", "V3"), new = c("resid_FID_avg", "SD", "FID_avg"))
a[is.na(SD), SD := 0]

aw <- reshape(a, idvar = c("Country", "IDLocality", "genus", "Species", "sp_loc"), timevar = "Covid", direction = "wide")
aw[, Species := gsub("[_]", " ", Species)]
aw <- merge(aw, t, all.x = TRUE)

aw[, sp2 := gsub(" ", "\n", Species)]
ann_text <- data.frame(
  FID_avg.0 = 8, FID_avg.1 = 10, lab = "Text",
  Species = factor("Aegithalos caudatus", levels = levels(as.factor(aw$Species)))
)
ann_text$sp2 <- gsub(" ", "\n", ann_text$Species)
g3 <-
  ggplot(aw, aes(x = FID_avg.0, y = FID_avg.1)) +
  # geom_errorbar(aes(ymin = FID_avg.1-SD.1, ymax = FID_avg.1+SD.1, col = Country), width = 0) +
  # geom_errorbar(aes(xmin = FID_avg.0-SD.0, xmax = FID_avg.0+SD.0, col = Country), width = 0) +
  # geom_point(pch = 21, alpha = 0.7, aes(col = Country)) +
  geom_point(pch = 21, alpha = 0.7, aes(fill = Country), col = "grey20", alpha = 0.8) + # col = "white") +
  # ggtitle ("Sim based")+
  geom_abline(intercept = 0, slope = 1, lty = 3, col = "grey80") +
  geom_text(data = ann_text, label = "No difference", col = "grey80", angle = 45, size = 2) +
  facet_wrap(~sp2) +
  # geom_phylopic(data = o, aes(image = uid),  color = "grey80", size = o$size) + # ,
  # scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = FALSE))  +
  scale_fill_manual(values = col3__, guide = guide_legend(reverse = TRUE)) +
  scale_x_continuous("Before COVID-19 shutdown - flight initiation distance [m]", expand = c(0, 0), trans = "log10") +
  scale_y_continuous("During COVID-19 shutdown - flight initiation distance [m]", expand = c(0, 0), trans = "log10") +
  labs(title = "Species means per sampling location") +
  theme_MB +
  theme(
    plot.title = element_text(size = 7),
    strip.background = element_blank(),
    # panel.spacing = unit(1, "mm"),
    legend.position = c(0.96, 0.0),
    legend.justification = c(1, 0)
  )
gg3 <- ggplotGrob(g3) # gg2$layout$name
ggx3 <- gtable_filter_remove(gg3,
  name = c(paste0("axis-b-", c(2, 4), "-7"), "axis-b-6-6"),
  trim = FALSE
)
if (save_plot == TRUE) {
  ggsave(here::here("Outputs/Fig_S5.png"), ggx3, width = 13.5, height = 17.5, unit = "cm", dpi = 600) # 11.43cm
}
grid::grid.draw(ggx3)
#' <a name="F_S5">
#' **Figure S5 | Avian tolerance towards humans before and during the COVID-19 shutdowns according to species.**</a> Dots represent means or single escape distance observations of species at specific sites (e.g. specific park or cemetery) with data for both periods (i.e. before and during the shutdowns) and not corrected for other factors such as starting distance of the observer. Dot colour highlights the country. Dotted lines indicate no difference, dots above the lines indicate lower tolerance towards humans (i.e. longer escape distances), dots below the lines indicate lower tolerance before than during the COVID-19 shutdowns. Panels are ordered alphabetically. The axes are on the log-scale.
#'
#' ***
#' 
#' ### Correlations among predictors
#+ Fig_S6_cor_pred, fig.width=9, fig.height = 9
d[, sin_rad := sin(rad)]
d[, cos_rad := cos(rad)]

dp <- d[, c("Covid", "StringencyIndex", "parks_percent_change_from_baseline", "Human", "SD_ln", "flock_ln", "body_ln", "sin_rad", "cos_rad", "Temp", "Day")]
setnames(dp, old = c("Covid", "StringencyIndex", "parks_percent_change_from_baseline", "Human", "SD_ln", "flock_ln", "body_ln", "sin_rad", "cos_rad", "Temp", "Day"), new = c("Period", "Stringency Index", "Google Mobility", "# of humans", "Starting distance\nln(m)", "Flock size\nln()", "Body mass\nln(g)", "Sine\nof radians", "Cosine\nof radians", "Temperature\n°C", "Day"))

# dp <- d[, c("SD_ln", "flock_ln", "body_ln", "sin_rad", "cos_rad", "Temp", "Day")]
# setnames(dp, old = c("SD_ln", "flock_ln", "body_ln", "sin_rad", "cos_rad", "Temp", "Day"), new = c("Starting distance\nln(m)", "Flock size\nln(m)", "Body mass\nln(m)", "Sine\nof radians", "Cosine\nof radians", "Temperature\n°C", "Day"))

#if (save_plot == TRUE) {
#  png(here::here("Outputs/Fig_S2_rev.png"), width = 19, height = 19, units = "cm", bg = "transparent", res = 600)
#  chart.Correlation(dp, histogram = TRUE, pch = 19, alpha = 0.5)
#  mtext("Single observations", side = 3, line = 3)
#  dev.off()
#}
chart.Correlation(dp, histogram = TRUE, pch = 19, alpha = 0.5)
mtext("Single observations", side = 3, line = 3)

#' <a name="F_S6">
#' **Figure S6 | Pairwise correlations among predictors of escape distance used in this study.**</a> On the diagonal: histograms and density lines (red) for each variable. Above diagonal: Pearson’s correlation with stars indicating significance. Below diagonal: the bivariate scatterplots, with each dot representing a single observation and a red line representing smoothed fit. Note that first four predictors were never entered togetheer in a single model with escape distance as dependent variable. Created with *chart.Correlation* function from R-package *PerformanceAnalytics* ([Peterson & Carl 2020](https://CRAN.R-project.org/package=PerformanceAnalytics)).
#'
#' ***
#' 
#' ### Number of humans ~ Google & Stringency
#+ Fig_6, fig.width=16*0.393701, fig.height = 8*0.393701
sh[, Country := factor(Country, levels = (c("Finland", "Czechia", "Hungary")))]
ssh[, Country := factor(Country, levels = (c("Finland", "Czechia", "Hungary")))]

# predictions for humans ~ stringency
lsh <- list()
shf <- sh[Country == "Finland"]
shfi <- lmer(Human ~
  Year +
  StringencyIndex +
  (scale(StringencyIndex) | year_weekday),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = shf, REML = FALSE
)
bsim <- sim(shfi, n.sim = nsim)
v <- apply(bsim@fixef, 2, quantile, prob = c(0.5))
newD <- data.frame(Year = mean(shf$Year), StringencyIndex = seq(min(shf$StringencyIndex), max(shf$StringencyIndex), length.out = 100)) # values to predict for
X <- model.matrix(~ Year + StringencyIndex, data = newD) # exactly the model which was used has to be specified here
newD$pred <- (X %*% v)
predmatrix <- matrix(nrow = nrow(newD), ncol = nsim)
for (j in 1:nsim) predmatrix[, j] <- (X %*% bsim@fixef[j, ])
newD$lwr <- apply(predmatrix, 1, quantile, prob = 0.025)
newD$upr <- apply(predmatrix, 1, quantile, prob = 0.975)
newD$pred <- apply(predmatrix, 1, quantile, prob = 0.5)
newD$Country <- "Finland"
newD$Year <- NULL
lsh[[1]] <- newD

shc <- sh[Country == "Czechia"]
shcz <- lmer(Human ~
  StringencyIndex +
  (scale(StringencyIndex) | weekday),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = shc, REML = FALSE
)
bsim <- sim(shcz, n.sim = nsim)
v <- apply(bsim@fixef, 2, quantile, prob = c(0.5))
newD <- data.frame(StringencyIndex = seq(min(shc$StringencyIndex), max(shc$StringencyIndex), length.out = 100)) # values to predict for
X <- model.matrix(~StringencyIndex, data = newD) # exactly the model which was used has to be specified here
newD$pred <- (X %*% v)
predmatrix <- matrix(nrow = nrow(newD), ncol = nsim)
for (j in 1:nsim) predmatrix[, j] <- (X %*% bsim@fixef[j, ])
newD$lwr <- apply(predmatrix, 1, quantile, prob = 0.025)
newD$upr <- apply(predmatrix, 1, quantile, prob = 0.975)
newD$pred <- apply(predmatrix, 1, quantile, prob = 0.5)
newD$Country <- "Czechia"
lsh[[2]] <- newD

sh_h <- sh[Country == "Hungary"]
shhu <- lmer(Human ~
  Year +
  StringencyIndex +
  (scale(StringencyIndex) | year_weekday),
data = sh_h, REML = FALSE
)

bsim <- sim(shhu, n.sim = nsim)
v <- apply(bsim@fixef, 2, quantile, prob = c(0.5))
newD <- data.frame(Year = mean(sh_h$Year), StringencyIndex = seq(min(sh_h$StringencyIndex), max(sh_h$StringencyIndex), length.out = 100)) # values to predict for
X <- model.matrix(~ Year + StringencyIndex, data = newD) # exactly the model which was used has to be specified here
newD$pred <- (X %*% v)
predmatrix <- matrix(nrow = nrow(newD), ncol = nsim)
for (j in 1:nsim) predmatrix[, j] <- (X %*% bsim@fixef[j, ])
newD$lwr <- apply(predmatrix, 1, quantile, prob = 0.025)
newD$upr <- apply(predmatrix, 1, quantile, prob = 0.975)
newD$pred <- apply(predmatrix, 1, quantile, prob = 0.5)
newD$Country <- "Hungary"
newD$Year <- NULL
lsh[[3]] <- newD

h_s <- data.table(do.call(rbind, lsh))
h_s[, Country := factor(Country, levels = (c("Finland", "Czechia", "Hungary")))]

# predictions for log10(humans+0.01) ~ stringency
lsh_ <- list()

shf <- sh[Country == "Finland"]
shfi_ln <- lmer(log10(Human+0.01) ~
  Year +
  StringencyIndex +
  (scale(StringencyIndex) | year_weekday),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = shf, REML = FALSE
)
bsim <- sim(shfi_ln, n.sim = nsim)
v <- apply(bsim@fixef, 2, quantile, prob = c(0.5))
newD <- data.frame(Year = mean(shf$Year), StringencyIndex = seq(min(shf$StringencyIndex), max(shf$StringencyIndex), length.out = 100)) # values to predict for
X <- model.matrix(~ Year + StringencyIndex, data = newD) # exactly the model which was used has to be specified here
newD$pred <- (X %*% v)
predmatrix <- matrix(nrow = nrow(newD), ncol = nsim)
for (j in 1:nsim) predmatrix[, j] <- (X %*% bsim@fixef[j, ])
newD$lwr <- apply(predmatrix, 1, quantile, prob = 0.025)
newD$upr <- apply(predmatrix, 1, quantile, prob = 0.975)
newD$pred <- apply(predmatrix, 1, quantile, prob = 0.5)
newD$Country <- "Finland"
newD$Year <- NULL
lsh_[[1]] <- newD

shc <- sh[Country == "Czechia"]
shcz_ln <- lmer(log10(Human+0.01) ~
  StringencyIndex +
  (scale(StringencyIndex) | weekday),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = shc, REML = FALSE
)
bsim <- sim(shcz_ln, n.sim = nsim)
v <- apply(bsim@fixef, 2, quantile, prob = c(0.5))
newD <- data.frame(StringencyIndex = seq(min(shc$StringencyIndex), max(shc$StringencyIndex), length.out = 100)) # values to predict for
X <- model.matrix(~StringencyIndex, data = newD) # exactly the model which was used has to be specified here
newD$pred <- (X %*% v)
predmatrix <- matrix(nrow = nrow(newD), ncol = nsim)
for (j in 1:nsim) predmatrix[, j] <- (X %*% bsim@fixef[j, ])
newD$lwr <- apply(predmatrix, 1, quantile, prob = 0.025)
newD$upr <- apply(predmatrix, 1, quantile, prob = 0.975)
newD$pred <- apply(predmatrix, 1, quantile, prob = 0.5)
newD$Country <- "Czechia"
lsh_[[2]] <- newD

sh_h <- sh[Country == "Hungary"]
shhu_ln <- lmer(log10(Human+0.01) ~
  Year +
  StringencyIndex +
  (scale(StringencyIndex) | year_weekday),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = sh_h, REML = FALSE
)

bsim <- sim(shhu_ln, n.sim = nsim)
v <- apply(bsim@fixef, 2, quantile, prob = c(0.5))
newD <- data.frame(Year = mean(sh_h$Year), StringencyIndex = seq(min(sh_h$StringencyIndex), max(sh_h$StringencyIndex), length.out = 100)) # values to predict for
X <- model.matrix(~ Year + StringencyIndex, data = newD) # exactly the model which was used has to be specified here
newD$pred <- (X %*% v)
predmatrix <- matrix(nrow = nrow(newD), ncol = nsim)
for (j in 1:nsim) predmatrix[, j] <- (X %*% bsim@fixef[j, ])
newD$lwr <- apply(predmatrix, 1, quantile, prob = 0.025)
newD$upr <- apply(predmatrix, 1, quantile, prob = 0.975)
newD$pred <- apply(predmatrix, 1, quantile, prob = 0.5)
newD$Country <- "Hungary"
newD$Year <- NULL
lsh_[[3]] <- newD

h_s_ln <- data.table(do.call(rbind, lsh_))
h_s_ln[, Country := factor(Country, levels = (c("Finland", "Czechia", "Hungary")))]
h_s_ln[, pred_o := 10^pred - 0.01]
h_s_ln[, lwr_o := 10^lwr - 0.01]
h_s_ln[, upr_o := 10^upr - 0.01]

# predictions for humans ~ Google
lshg <- list()

sshf <- ssh[Country == "Finland"]
sshfi <- lmer(Human ~
  Year +
  parks_percent_change_from_baseline +
  (scale(parks_percent_change_from_baseline) | year_weekday),
data = sshf, REML = FALSE
)
bsim <- sim(sshfi, n.sim = nsim)
v <- apply(bsim@fixef, 2, quantile, prob = c(0.5))
newD <- data.frame(Year = mean(sshf$Year), parks_percent_change_from_baseline = seq(min(sshf$parks_percent_change_from_baseline), max(sshf$parks_percent_change_from_baseline), length.out = 100)) # values to predict for
X <- model.matrix(~ Year + parks_percent_change_from_baseline, data = newD) # exactly the model which was used has to be specified here
newD$pred <- (X %*% v)
predmatrix <- matrix(nrow = nrow(newD), ncol = nsim)
for (j in 1:nsim) predmatrix[, j] <- (X %*% bsim@fixef[j, ])
newD$lwr <- apply(predmatrix, 1, quantile, prob = 0.025)
newD$upr <- apply(predmatrix, 1, quantile, prob = 0.975)
newD$pred <- apply(predmatrix, 1, quantile, prob = 0.5)
newD$Country <- "Finland"
newD$Year <- NULL
lshg[[1]] <- newD

sshc <- ssh[Country == "Czechia"]
sshcz <- lmer(Human ~
  parks_percent_change_from_baseline +
  (scale(parks_percent_change_from_baseline) | weekday),
# (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
data = sshc, REML = FALSE
)
bsim <- sim(sshcz, n.sim = nsim)
v <- apply(bsim@fixef, 2, quantile, prob = c(0.5))
newD <- data.frame(parks_percent_change_from_baseline = seq(min(sshc$parks_percent_change_from_baseline), max(sshc$parks_percent_change_from_baseline), length.out = 100)) # values to predict for
X <- model.matrix(~parks_percent_change_from_baseline, data = newD) # exactly the model which was used has to be specified here
newD$pred <- (X %*% v)
predmatrix <- matrix(nrow = nrow(newD), ncol = nsim)
for (j in 1:nsim) predmatrix[, j] <- (X %*% bsim@fixef[j, ])
newD$lwr <- apply(predmatrix, 1, quantile, prob = 0.025)
newD$upr <- apply(predmatrix, 1, quantile, prob = 0.975)
newD$pred <- apply(predmatrix, 1, quantile, prob = 0.5)
newD$Country <- "Czechia"
lshg[[2]] <- newD

ssh_h <- ssh[Country == "Hungary"]
sshhu <- lmer(Human ~
  Year +
  parks_percent_change_from_baseline +
  ((parks_percent_change_from_baseline) | year_weekday),
data = ssh_h, REML = FALSE
)

bsim <- sim(sshhu, n.sim = nsim)
v <- apply(bsim@fixef, 2, quantile, prob = c(0.5))
newD <- data.frame(Year = mean(ssh_h$Year), parks_percent_change_from_baseline = seq(min(ssh_h$parks_percent_change_from_baseline), max(ssh_h$parks_percent_change_from_baseline), length.out = 100)) # values to predict for
X <- model.matrix(~ Year + parks_percent_change_from_baseline, data = newD) # exactly the model which was used has to be specified here
newD$pred <- (X %*% v)
predmatrix <- matrix(nrow = nrow(newD), ncol = nsim)
for (j in 1:nsim) predmatrix[, j] <- (X %*% bsim@fixef[j, ])
newD$lwr <- apply(predmatrix, 1, quantile, prob = 0.025)
newD$upr <- apply(predmatrix, 1, quantile, prob = 0.975)
newD$pred <- apply(predmatrix, 1, quantile, prob = 0.5)
newD$Country <- "Hungary"
newD$Year <- NULL
lshg[[3]] <- newD

h_sg <- data.table(do.call(rbind, lshg))
h_sg[, Country := factor(Country, levels = (c("Finland", "Czechia", "Hungary")))]

# predictions for log(humans+0.01) ~ Google
lshg_ <- list()

sshf <- ssh[Country == "Finland"]
sshfi_ln <- lmer(log10(Human+0.01) ~
  Year +
  parks_percent_change_from_baseline +
  (scale(parks_percent_change_from_baseline) | year_weekday),
data = sshf, REML = FALSE
)
bsim <- sim(sshfi_ln, n.sim = nsim)
v <- apply(bsim@fixef, 2, quantile, prob = c(0.5))
newD <- data.frame(Year = mean(sshf$Year), parks_percent_change_from_baseline = seq(min(sshf$StringencyIndex), max(sshf$parks_percent_change_from_baseline), length.out = 100)) # values to predict for
X <- model.matrix(~ Year + parks_percent_change_from_baseline, data = newD) # exactly the model which was used has to be specified here
newD$pred <- (X %*% v)
predmatrix <- matrix(nrow = nrow(newD), ncol = nsim)
for (j in 1:nsim) predmatrix[, j] <- (X %*% bsim@fixef[j, ])
newD$lwr <- apply(predmatrix, 1, quantile, prob = 0.025)
newD$upr <- apply(predmatrix, 1, quantile, prob = 0.975)
newD$pred <- apply(predmatrix, 1, quantile, prob = 0.5)
newD$Country <- "Finland"
newD$Year <- NULL
lshg_[[1]] <- newD

sshc <- ssh[Country == "Czechia"]
sshcz_ln <- lmer(log10(Human+0.01) ~
  parks_percent_change_from_baseline +
  (scale(parks_percent_change_from_baseline) | weekday),
data = sshc, REML = FALSE
)
bsim <- sim(sshcz_ln, n.sim = nsim)
v <- apply(bsim@fixef, 2, quantile, prob = c(0.5))
newD <- data.frame(parks_percent_change_from_baseline = seq(min(sshc$parks_percent_change_from_baseline), max(sshc$parks_percent_change_from_baseline), length.out = 100)) # values to predict for
X <- model.matrix(~parks_percent_change_from_baseline, data = newD) # exactly the model which was used has to be specified here
newD$pred <- (X %*% v)
predmatrix <- matrix(nrow = nrow(newD), ncol = nsim)
for (j in 1:nsim) predmatrix[, j] <- (X %*% bsim@fixef[j, ])
newD$lwr <- apply(predmatrix, 1, quantile, prob = 0.025)
newD$upr <- apply(predmatrix, 1, quantile, prob = 0.975)
newD$pred <- apply(predmatrix, 1, quantile, prob = 0.5)
newD$Country <- "Czechia"
lshg_[[2]] <- newD

ssh_h <- ssh[Country == "Hungary"]
sshhu_ln <- lmer(log10(Human+0.01) ~
  Year +
  parks_percent_change_from_baseline +
  (scale(parks_percent_change_from_baseline) | year_weekday),
data = ssh_h, REML = FALSE
)

bsim <- sim(sshhu_ln, n.sim = nsim)
v <- apply(bsim@fixef, 2, quantile, prob = c(0.5))
newD <- data.frame(Year = mean(ssh_h$Year), parks_percent_change_from_baseline = seq(min(ssh_h$parks_percent_change_from_baseline), max(ssh_h$parks_percent_change_from_baseline), length.out = 100)) # values to predict for
X <- model.matrix(~ Year + parks_percent_change_from_baseline, data = newD) # exactly the model which was used has to be specified here
newD$pred <- (X %*% v)
predmatrix <- matrix(nrow = nrow(newD), ncol = nsim)
for (j in 1:nsim) predmatrix[, j] <- (X %*% bsim@fixef[j, ])
newD$lwr <- apply(predmatrix, 1, quantile, prob = 0.025)
newD$upr <- apply(predmatrix, 1, quantile, prob = 0.975)
newD$pred <- apply(predmatrix, 1, quantile, prob = 0.5)
newD$Country <- "Hungary"
newD$Year <- NULL
lshg_[[3]] <- newD

h_g_ln <- data.table(do.call(rbind, lshg_))
h_g_ln[, Country := factor(Country, levels = (c("Finland", "Czechia", "Hungary")))]
h_g_ln[, pred_o := 10^pred - 0.01]
h_g_ln[, lwr_o := 10^lwr - 0.01]
h_g_ln[, upr_o := 10^upr - 0.01]

# plot
col_h <- c("#357EBDFF", "#D43F3AFF", "#46B8DAFF", "#5CB85CFF", "#EEA236FF", "#9632B8FF", "#9632B8FF")[7:1]
col_h <- col_h[c(4, 5, 7)] #show_col(col_h)
# original stringency
p_hs_o <-
  ggplot(h_s, aes(x = StringencyIndex, y = pred, col = Country)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Country, color = NULL), alpha = .15) +
  # geom_point(aes(y = Human, fill = Country), data = sh, pch = 21, col = "grey20",alpha = 0.5)+
  geom_jitter(aes(y = Human, fill = Country), data = sh, pch = 21, col = "grey20", width = 0.5, height = 0.1, alpha = 0.5) +
  geom_line(lwd = 1) +
  labs(subtitle = "", y = "# of humans", x = "Stringency index") +
  coord_cartesian(xlim = c(25, 75), ylim = c(0, 50)) +
  facet_wrap(~Country) +
  # scale_color_locuszoom()+
  # scale_fill_locuszoom(guide = "none")
  scale_x_continuous(breaks = round(seq(25, 75, by = 25), 1)) +
  # scale_y_continuous(breaks = log(c(0.01, 1, 10, 50)), labels = c(0, 1, 10, 50)) +
  # scale_y_continuous(breaks = round(seq(-100, 175, by = 25), 1)) +
  scale_colour_manual(
    values = col_h, guide = guide_legend(reverse = TRUE, override.aes = list(size = 0)),
    labels = paste(
      "<span style='color:",
      col_h,
      "'>",
      levels(h_s$Country),
      "</span>"
    )
  ) +
  scale_fill_manual(values = col_h, guide = "none") +
  theme_bw() +
  theme(
    legend.text = element_markdown(size = 6),
    legend.position = "none",
    legend.title = element_blank(),
    # legend.spacing.y = unit(0.1, 'cm'),
    legend.key.height = unit(0.5, "line"),
    legend.key.size = unit(0, "line"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-10, 1, -10, -10),
    # legend.position=c(0.5,1.6),
    plot.title = element_text(color = "grey", size = 7),
    plot.subtitle = element_text(color = "grey60", size = 6),
    plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r = 0.5, unit = "pt"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = ax_lines, size = 0.25),
    axis.ticks = element_line(colour = ax_lines, size = 0.25),
    # axis.text.x = element_text()
    axis.ticks.length = unit(1, "pt"),
    axis.text = element_text(, size = 6),
    axis.title = element_text(size = 7),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 6, color = "grey30", margin = margin(1, 1, 1, 1, "mm"))
  )
#p_hs_o
#if (save_plot == TRUE) {
#  ggsave(here::here("Outputs/Fig_HS_temp_rev_widht_70mm_v2_origi.png"), p_hs_o #+ theme(plot.subtitle = element_blank())
#    ,
#    width = 10, height = 4.5, unit = "cm", dpi = 600
#  )
#}

# log stringency
p_hs_ln = ggplot(h_s_ln, aes(x = StringencyIndex, y = pred, col = Country)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Country, color = NULL), alpha = .15) +
  # geom_point(aes(y = Human, fill = Country), data = sh, pch = 21, col = "grey20",alpha = 0.5)+
  geom_jitter(aes(y = log10(Human + 0.01), fill = Country), data = sh, pch = 21, col = "grey20", width = 0.5, height = 0.1, alpha = 0.5) +
  geom_line(lwd = 1) +
  labs(subtitle = "", y = "# of humans", x = "Stringency index") +
  coord_cartesian(xlim = c(25, 75), ylim = c(log10(0.01), log10(50))) +
  facet_wrap(~Country) +
  # scale_color_locuszoom()+
  # scale_fill_locuszoom(guide = "none")
  scale_x_continuous(breaks = round(seq(25, 75, by = 25), 1)) +
  scale_y_continuous(breaks = log10(c(0.01, 1, 10, 50)), labels = c(0, 1, 10, 50)) +
  # scale_y_continuous(breaks = round(seq(-100, 175, by = 25), 1)) +
  scale_colour_manual(
    values = col_h, guide = guide_legend(reverse = TRUE, override.aes = list(size = 0)),
    labels = paste(
      "<span style='color:",
      col_h,
      "'>",
      levels(h_s$Country),
      "</span>"
    )
  ) +
  scale_fill_manual(values = col_h, guide = "none") +
  theme_bw() +
  theme(
    legend.text = element_markdown(size = 6),
    legend.position = "none",
    legend.title = element_blank(),
    # legend.spacing.y = unit(0.1, 'cm'),
    legend.key.height = unit(0.5, "line"),
    legend.key.size = unit(0, "line"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-10, 1, -10, -10),
    # legend.position=c(0.5,1.6),
    plot.title = element_text(color = "grey", size = 7),
    plot.subtitle = element_text(color = "grey60", size = 6),
    plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r = 0.5, unit = "pt"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = ax_lines, size = 0.25),
    axis.ticks = element_line(colour = ax_lines, size = 0.25),
    # axis.text.x = element_text()
    axis.ticks.length = unit(1, "pt"),
    axis.text = element_text(, size = 6),
    axis.title = element_text(size = 7),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 6, color = "grey30", margin = margin(1, 1, 1, 1, "mm"))
  )
#p_hs_ln
#if (save_plot == TRUE) {
#  ggsave(here::here("Outputs/Fig_HS_ln.png"), p_hs_ln #+ theme(plot.subtitle = element_blank())
#    ,
#    width = 10, height = 4.5, unit = "cm", dpi = 600
#  )
#}

# original google
p_hg_o <-
  ggplot(h_sg, aes(x = parks_percent_change_from_baseline, y = pred, col = Country)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Country, color = NULL), alpha = .15) +
  # geom_point(aes(y = Human, fill = Country), data = sh, pch = 21, col = "grey20",alpha = 0.5)+
  geom_jitter(aes(y = Human, fill = Country), data = ssh, pch = 21, col = "grey20", width = 0.5, height = 0.1, alpha = 0.5) +
  geom_line(lwd = 1) +
  labs(subtitle = "Original scale", y = "# of humans", x = "Google Mobility") +
  coord_cartesian(xlim = c(-50, 150), ylim = c(0, 50)) +
  facet_wrap(~Country) +
  # scale_color_locuszoom()+
  # scale_fill_locuszoom(guide = "none")
  #scale_x_continuous(breaks = round(seq(25, 75, by = 25), 1)) +
  # scale_y_continuous(breaks = log(c(0.01, 1, 10, 50)), labels = c(0, 1, 10, 50)) +
  # scale_y_continuous(breaks = round(seq(-100, 175, by = 25), 1)) +
  scale_colour_manual(
    values = col_h, guide = guide_legend(reverse = TRUE, override.aes = list(size = 0)),
    labels = paste(
      "<span style='color:",
      col_h,
      "'>",
      levels(h_s$Country),
      "</span>"
    )
  ) +
  scale_fill_manual(values = col_h, guide = "none") +
  theme_bw() +
  theme(
    legend.text = element_markdown(size = 6),
    legend.position = "none",
    legend.title = element_blank(),
    # legend.spacing.y = unit(0.1, 'cm'),
    legend.key.height = unit(0.5, "line"),
    legend.key.size = unit(0, "line"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-10, 1, -10, -10),
    # legend.position=c(0.5,1.6),
    plot.title = element_text(color = "grey", size = 7),
    plot.subtitle = element_text(color = "grey60", size = 6),
    plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r = 0.5, unit = "pt"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = ax_lines, size = 0.25),
    axis.ticks = element_line(colour = ax_lines, size = 0.25),
    # axis.text.x = element_text()
    axis.ticks.length = unit(1, "pt"),
    axis.text = element_text(, size = 6),
    axis.title = element_text(size = 7),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 6, color = "grey30", margin = margin(1, 1, 1, 1, "mm"))
  )
#p_hg_o
#if (save_plot == TRUE) {
#  ggsave(here::here("Outputs/Fig_HG_rev_widht_70mm_v2_origi.png"), p_hg_o #+ theme(plot.subtitle = element_blank())
#    ,
#    width = 10, height = 4.5, unit = "cm", dpi = 600
#  )
#}

# log google
p_hg_ln <- ggplot(h_g_ln, aes(x = parks_percent_change_from_baseline, y = pred, col = Country)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Country, color = NULL), alpha = .15) +
  # geom_point(aes(y = Human, fill = Country), data = sh, pch = 21, col = "grey20",alpha = 0.5)+
  geom_jitter(aes(y = log10(Human + 0.01), fill = Country), data = ssh, pch = 21, col = "grey20", width = 0.5, height = 0.1, alpha = 0.5) +
  geom_line(lwd = 1) +
  labs(subtitle = "Log-scale", y = "# of humans", x = "Google Mobility") +
  coord_cartesian(xlim = c(-50, 150), ylim = c(log10(0.01), log10(50))) +
  facet_wrap(~Country) +
  # scale_color_locuszoom()+
  # scale_fill_locuszoom(guide = "none")
  #scale_x_continuous(breaks = round(seq(25, 75, by = 25), 1)) +
  scale_y_continuous(breaks = log10(c(0.01, 1, 10, 50)), labels = c(0, 1, 10, 50)) +
  # scale_y_continuous(breaks = round(seq(-100, 175, by = 25), 1)) +
  scale_colour_manual(
    values = col_h, guide = guide_legend(reverse = TRUE, override.aes = list(size = 0)),
    labels = paste(
      "<span style='color:",
      col_h,
      "'>",
      levels(h_s$Country),
      "</span>"
    )
  ) +
  scale_fill_manual(values = col_h, guide = "none") +
  theme_bw() +
  theme(
    legend.text = element_markdown(size = 6),
    legend.position = "none",
    legend.title = element_blank(),
    # legend.spacing.y = unit(0.1, 'cm'),
    legend.key.height = unit(0.5, "line"),
    legend.key.size = unit(0, "line"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-10, 1, -10, -10),
    # legend.position=c(0.5,1.6),
    plot.title = element_text(color = "grey", size = 7),
    plot.subtitle = element_text(color = "grey60", size = 6),
    plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r = 0.5, unit = "pt"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = ax_lines, size = 0.25),
    axis.ticks = element_line(colour = ax_lines, size = 0.25),
    # axis.text.x = element_text()
    axis.ticks.length = unit(1, "pt"),
    axis.text = element_text(, size = 6),
    axis.title = element_text(size = 7),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 6, color = "grey30", margin = margin(1, 1, 1, 1, "mm"))
  )
#p_hg_ln
#if (save_plot == TRUE) {
#  ggsave(here::here("Outputs/Fig_GS_ln.png"), p_hs_ln #+ theme(plot.subtitle = element_blank())
#    ,
#    width = 10, height = 4.5, unit = "cm", dpi = 600
#  )
#}

# combine & export
sshh <- ggarrange(
  p_hg_o + rremove("grid")+  rremove("ylab") + xlab('') + theme(plot.margin = margin(r = 4)), 
  p_hs_o + rremove("grid") + rremove("ylab") + xlab('')+ theme(plot.margin = margin(r = 4)) , 
  p_hg_ln + rremove("grid") + rremove("ylab")+ theme(plot.margin = margin(r = 4)) , 
  p_hs_ln + rremove("grid") +  rremove("ylab")+ theme(plot.margin = margin(r = 4)) ,
  ncol = 2, nrow = 2
  )

fig_6 = annotate_figure(sshh + rremove("grid"),
  left = textGrob("# of humans", rot = 90, gp = gpar(cex = 0.6))
)

if (save_plot == TRUE) {
ggsave(file = here::here("Outputs/Fig_6_width-160mm.png"), fig_6, width = 16, height = 8, units = "cm",bg = "white")
}

fig_6

#' <a name="F_6">
#' **Figure 6 | Human numbers within 50 meters during the escape trial in association with daily human presence in parks (Google Mobility, left) and stringency of antipandemic governmental restrictions (Stringency index, right).**</a> Dots represent individual data points (**top** panels on original scale, **bottom** on log-scale), jittered to increase visibility. Lines with shaded areas represent predictions with 95%CIs from mixed effect models controlled for the year (in case of Finland and Hungary) and non-independence of data points by including weekday within the year as a random intercept and Stringency index or Google Mobility as a random slope (Table [S3a](#T_S3a), [S3b](#T_S3b), [S3c](#T_S3c), [S3d](#T_S3d)). Human numbers during the escape trials shutdowns were missing for Australia (all year) and Poland (during shutdowns). Note the weak and country-specific associations.
#' 
#' ***
#' 
#' <a name="T_S3a">
#' **Table S3a | Number of humans in relation to Stringency index**</a>
llb = list()

shf <- sh[Country == "Finland"]
shfi_ <- lmer(scale(Human) ~
  scale(Year) +
  scale(StringencyIndex) +
  (scale(StringencyIndex) | year_weekday),
data = shf, REML = FALSE
)
llb[[1]] = m_out(name = "Table S3a - FI", dep = "# of humans", model = shfi_, nsim = 5000)

shc <- sh[Country == "Czechia"]
shcz_ <- lmer(scale(Human) ~
  scale(StringencyIndex) +
  (scale(StringencyIndex) | weekday),
data = shc, REML = FALSE
)
llb[[2]] = m_out(name = "Table S3a - CZ", dep = "# of humans", model = shcz_, nsim = 5000)
llb[[2]]$R2_mar <- llb[[2]]$R2_con <- NULL

sh_h <- sh[Country == "Hungary"]
shhu_ <- lmer(scale(Human) ~
  scale(Year) +
  scale(StringencyIndex) +
  (scale(StringencyIndex) | year_weekday),
data = sh_h, REML = FALSE
)
llb[[3]] = m_out(name = "Table S3a - HU", dep = "# of humans", model = shhu_, nsim = 5000)
llb[[3]]$R2_mar = llb[[3]]$R2_con = NULL

out_t3b = data.table(do.call(rbind, llb))
out_t3b[is.na(out_t3b)] <- ""
out_t3b$R2_mar = out_t3b$R2_con = NULL
out_t3b[, effect := gsub("scale\\(Year\\)", "year", effect)]
out_t3b[, effect := gsub("scale\\(StringencyIndex\\)", "Stringency index", effect)]
out_t3b[, effect := gsub("year_weekday", "weekday within year", effect)]
out_t3b[type == "random" & grepl("Stringency index", effect, fixed = TRUE), effect := paste("Stringency index (slope) |", gsub(" Stringency index", "", effect))]
fwrite(file = here::here("Outputs/Table_S3a_rev.csv"), out_t3b)

out_t3b$error_structure = out_t3b$response = NULL
out_t3b[model != "", model := c(
  "Finland",
  "Czechia",
  "Hungary"
)]
setnames(out_t3b, old = c("estimate_r", "lwr_r", "upr_r"), new = c("estimate", "lower", "upper"))
out_t3b %>%
  kbl() %>%
  kable_paper("hover", full_width = F) %>%
    scroll_box(width = "100%")
#'
#' ***
#'
#' <a name="T_S3b">
#' **Table S3b | Number of humans (ln-transformed) in relation to Stringency index**</a>
llc = list()

shf <- sh[Country == "Finland"]
shfi_ln_ <- lmer(scale(log(Human+0.01)) ~
  scale(Year) +
  scale(StringencyIndex) +
  (scale(StringencyIndex) | year_weekday),
data = shf, REML = FALSE
)
llc[[1]] = m_out(name = "Table S3b - FI", dep = "# of humans (ln)", model = shfi_ln_, nsim = 5000)
llc[[1]]$R2_mar <- llc[[1]]$R2_con <- NULL

shc <- sh[Country == "Czechia"]
shcz_ln_ <- lmer(scale(log(Human+0.01))~
  scale(StringencyIndex) +
  (scale(StringencyIndex) | weekday),
data = shc, REML = FALSE
)
llc[[2]] = m_out(name = "Table S3b - CZ", dep = "# of humans (ln)", model = shcz_ln_, nsim = 5000)
llc[[2]]$R2_mar <- llc[[2]]$R2_con <- NULL

sh_h <- sh[Country == "Hungary"]
shhu_ln_ <- lmer(scale(log(Human+0.01)) ~
  scale(Year) +
  scale(StringencyIndex) +
  (scale(StringencyIndex) | year_weekday),
data = sh_h, REML = FALSE
)
llc[[3]] = m_out(name = "Table S3b - HU", dep = "# of humans (ln)", model = shhu_ln_, nsim = 5000)

out_t3c = data.table(do.call(rbind, llc))
out_t3c[is.na(out_t3c)] <- ""
out_t3c$R2_mar = out_t3c$R2_con = NULL
out_t3c[, effect := gsub("scale\\(Year\\)", "year", effect)]
out_t3c[, effect := gsub("scale\\(StringencyIndex\\)", "Stringency index", effect)]
out_t3c[, effect := gsub("year_weekday", "weekday within year", effect)]
out_t3c[type == "random" & grepl("Stringency index", effect, fixed = TRUE), effect := paste("Stringency index (slope) |", gsub(" Stringency index", "", effect))]
fwrite(file = here::here("Outputs/Table_S3b_rev.csv"), out_t3c)

out_t3c$error_structure = out_t3b$response = NULL
out_t3c[model != "", model := c(
  "Finland",
  "Czechia",
  "Hungary"
)]
setnames(out_t3c, old = c("estimate_r", "lwr_r", "upr_r"), new = c("estimate", "lower", "upper"))
out_t3c %>%
  kbl() %>%
  kable_paper("hover", full_width = F) %>%
  scroll_box(width = "100%")
#'
#' ***
#' 
#' <a name="T_S3c">
#' **Table S3c | Number of humans in relation to Google Mobility**</a>
lld = list()

sshf <- ssh[Country == "Finland"]
shfi_g <- lmer(scale(Human) ~
  scale(Year) +
  scale(parks_percent_change_from_baseline) +
  (scale(parks_percent_change_from_baseline) | year_weekday),
data = sshf, REML = FALSE
)
lld[[1]] = m_out(name = "Table S3c - FI", dep = "# of humans", model = shfi_g, nsim = 5000, R2 = FALSE)

sshc <- ssh[Country == "Czechia"]
shcz_g <- lmer(scale(Human) ~
  scale(parks_percent_change_from_baseline) +
  (scale(parks_percent_change_from_baseline) | weekday),
data = sshc, REML = FALSE
)
lld[[2]] = m_out(name = "Table S3c - CZ", dep = "# of humans", model = shcz_g, nsim = 5000, , R2 = FALSE)

ssh_h <- ssh[Country == "Hungary"]
shhu_g <- lmer(scale(Human) ~
  scale(Year) +
  scale(parks_percent_change_from_baseline) +
  (scale(parks_percent_change_from_baseline) | year_weekday),
data = ssh_h, REML = FALSE
)
lld[[3]] = m_out(name = "Table S3c - HU", dep = "# of humans", model = shhu_g, nsim = 5000, R2 = FALSE)

out_t3d = data.table(do.call(rbind, lld))
out_t3d[is.na(out_t3d)] <- ""
out_t3d[, effect := gsub("scale\\(Year\\)", "year", effect)]
out_t3d[, effect := gsub("scale\\(parks_percent_change_from_baseline\\)", "Google Mobility", effect)]
out_t3d[, effect := gsub("year_weekday", "weekday within year", effect)]
out_t3d[type == "random" & grepl("Google Mobility", effect, fixed = TRUE), effect := paste("Google Mobility (slope) |", gsub(" Google Mobility", "", effect))]
fwrite(file = here::here("Outputs/Table_S3c_rev.csv"), out_t3d)

out_t3d$error_structure = out_t3d$response = NULL
out_t3d[model != "", model := c(
  "Finland",
  "Czechia",
  "Hungary"
)]
setnames(out_t3d, old = c("estimate_r", "lwr_r", "upr_r"), new = c("estimate", "lower", "upper"))
out_t3d %>%
  kbl() %>%
  kable_paper("hover", full_width = F) %>%
    scroll_box(width = "100%")
#'
#' ***
#'
#' <a name="T_S3d">
#' **Table S3d | Number of humans (ln-transformed) in relation to Google Mobility**</a>
lle = list()

sshf <- ssh[Country == "Finland"]
shfi_g_ln <- lmer(scale(log(Human + 0.01)) ~
  scale(Year) +
  scale(parks_percent_change_from_baseline) +
  (scale(parks_percent_change_from_baseline) | year_weekday),
data = sshf, REML = FALSE
)
lle[[1]] = m_out(name = "Table S3d - FI", dep = "# of humans (ln)", model = shfi_g_ln, nsim = 5000, R2 = FALSE)

sshc <- ssh[Country == "Czechia"]
shcz_g_ln <- lmer(scale(log(Human + 0.01)) ~
  scale(parks_percent_change_from_baseline) +
  (scale(parks_percent_change_from_baseline) | weekday),
data = sshc, REML = FALSE
)
lle[[2]] = m_out(name = "Table S3d - CZ", dep = "# of humans (ln)", model = shcz_g_ln, nsim = 5000, R2 = FALSE)

ssh_h <- ssh[Country == "Hungary"]
shhu_g_ln <- lmer(scale(log(Human + 0.01)) ~
  scale(Year) +
  scale(parks_percent_change_from_baseline) +
  (scale(parks_percent_change_from_baseline) | year_weekday),
data = ssh_h, REML = FALSE
)
lle[[3]] = m_out(name = "Table S3d - HU", dep = "# of humans (ln)", model = shhu_g_ln, nsim = 5000, R2 = FALSE)

out_t3e = data.table(do.call(rbind, lle))
out_t3e[is.na(out_t3e)] <- ""
out_t3e[, effect := gsub("scale\\(Year\\)", "year", effect)]
out_t3e[, effect := gsub("scale\\(parks_percent_change_from_baseline\\)", "Google Mobility", effect)]
out_t3e[, effect := gsub("year_weekday", "weekday within year", effect)]
out_t3e[type == "random" & grepl("Google Mobility", effect, fixed = TRUE), effect := paste("Google Mobility (slope) |", gsub(" Google Mobility", "", effect))]
fwrite(file = here::here("Outputs/Table_S3d_rev.csv"), out_t3e)

out_t3e$error_structure = out_t3e$response = NULL
out_t3e[model != "", model := c(
  "Finland",
  "Czechia",
  "Hungary"
)]
setnames(out_t3e, old = c("estimate_r", "lwr_r", "upr_r"), new = c("estimate", "lower", "upper"))
out_t3e %>%
  kbl() %>%
  kable_paper("hover", full_width = F) %>%
    scroll_box(width = "100%")
#' ***
#'
#' ### Google Mobility ~ Stringency index
#+ Fig_7_gsfig_pred, fig.width=7/2.5, fig.height = 6/2.5
# Predictions
l <- list()
sscz <- ss[Country == "Czechia"]
cz <- lmer(
  parks_percent_change_from_baseline ~
    StringencyIndex +
    (scale(StringencyIndex) | weekday),
  # (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
  data = sscz, REML = FALSE
)
bsim <- sim(cz, n.sim = nsim)
v <- apply(bsim@fixef, 2, quantile, prob = c(0.5))
newD <- data.frame(StringencyIndex = seq(min(sscz$StringencyIndex), max(sscz$StringencyIndex), length.out = 100)) # values to predict for
X <- model.matrix(~StringencyIndex, data = newD) # exactly the model which was used has to be specified here
newD$pred <- (X %*% v)
predmatrix <- matrix(nrow = nrow(newD), ncol = nsim)
for (j in 1:nsim) predmatrix[, j] <- (X %*% bsim@fixef[j, ])
newD$lwr <- apply(predmatrix, 1, quantile, prob = 0.025)
newD$upr <- apply(predmatrix, 1, quantile, prob = 0.975)
newD$pred <- apply(predmatrix, 1, quantile, prob = 0.5)
newD$Country <- "Czechia"
l[[1]] <- newD

s[, year_weekday := paste(Year, weekday)]
sf <- ss[Country == "Finland"]
fi <- lmer(
  parks_percent_change_from_baseline ~
    Year +
    StringencyIndex +
    (scale(StringencyIndex) | year_weekday),
  # (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
  data = sf, REML = FALSE
)
bsim <- sim(fi, n.sim = nsim)
v <- apply(bsim@fixef, 2, quantile, prob = c(0.5))
newD <- data.frame(Year = mean(sf$Year), StringencyIndex = seq(min(sf$StringencyIndex), max(sf$StringencyIndex), length.out = 100)) # values to predict for
X <- model.matrix(~ Year + StringencyIndex, data = newD) # exactly the model which was used has to be specified here
newD$pred <- (X %*% v)
predmatrix <- matrix(nrow = nrow(newD), ncol = nsim)
for (j in 1:nsim) predmatrix[, j] <- (X %*% bsim@fixef[j, ])
newD$lwr <- apply(predmatrix, 1, quantile, prob = 0.025)
newD$upr <- apply(predmatrix, 1, quantile, prob = 0.975)
newD$pred <- apply(predmatrix, 1, quantile, prob = 0.5)
newD$Country <- "Finland"
newD$Year <- NULL
l[[2]] <- newD

s_h <- ss[Country == "Hungary"]
hu <- lmer(
  parks_percent_change_from_baseline ~
    Year +
    StringencyIndex +
    (scale(StringencyIndex) | year_weekday),
  # (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
  data = s_h, REML = FALSE
)

bsim <- sim(hu, n.sim = nsim)
v <- apply(bsim@fixef, 2, quantile, prob = c(0.5))
newD <- data.frame(Year = mean(s_h$Year), StringencyIndex = seq(min(s_h$StringencyIndex), max(s_h$StringencyIndex), length.out = 100)) # values to predict for
X <- model.matrix(~ Year + StringencyIndex, data = newD) # exactly the model which was used has to be specified here
newD$pred <- (X %*% v)
predmatrix <- matrix(nrow = nrow(newD), ncol = nsim)
for (j in 1:nsim) predmatrix[, j] <- (X %*% bsim@fixef[j, ])
newD$lwr <- apply(predmatrix, 1, quantile, prob = 0.025)
newD$upr <- apply(predmatrix, 1, quantile, prob = 0.975)
newD$pred <- apply(predmatrix, 1, quantile, prob = 0.5)
newD$Country <- "Hungary"
newD$Year <- NULL
l[[3]] <- newD

sp <- ss[Country == "Poland"]
pl <- lmer(
  parks_percent_change_from_baseline ~
    Year +
    StringencyIndex +
    (scale(StringencyIndex) | year_weekday),
  # (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
  data = sp, REML = FALSE
)
bsim <- sim(pl, n.sim = nsim)
v <- apply(bsim@fixef, 2, quantile, prob = c(0.5))
newD <- data.frame(Year = mean(sp$Year), StringencyIndex = seq(min(sp$StringencyIndex), max(sp$StringencyIndex), length.out = 100)) # values to predict for
X <- model.matrix(~ Year + StringencyIndex, data = newD) # exactly the model which was used has to be specified here
newD$pred <- (X %*% v)
predmatrix <- matrix(nrow = nrow(newD), ncol = nsim)
for (j in 1:nsim) predmatrix[, j] <- (X %*% bsim@fixef[j, ])
newD$lwr <- apply(predmatrix, 1, quantile, prob = 0.025)
newD$upr <- apply(predmatrix, 1, quantile, prob = 0.975)
newD$pred <- apply(predmatrix, 1, quantile, prob = 0.5)
newD$Country <- "Poland"
newD$Year <- NULL
l[[4]] <- newD

sa <- ss[Country == "Australia"]
au <- lmer(
  parks_percent_change_from_baseline ~
    Year +
    StringencyIndex +
    (scale(StringencyIndex) | year_weekday),
  # (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
  data = sa, REML = FALSE
)

bsim <- sim(au, n.sim = nsim)
v <- apply(bsim@fixef, 2, quantile, prob = c(0.5))
newD <- data.frame(Year = mean(sa$Year), StringencyIndex = seq(min(sa$StringencyIndex), max(sa$StringencyIndex), length.out = 100)) # values to predict for
X <- model.matrix(~ Year + StringencyIndex, data = newD) # exactly the model which was used has to be specified here
newD$pred <- (X %*% v)
predmatrix <- matrix(nrow = nrow(newD), ncol = nsim)
for (j in 1:nsim) predmatrix[, j] <- (X %*% bsim@fixef[j, ])
newD$lwr <- apply(predmatrix, 1, quantile, prob = 0.025)
newD$upr <- apply(predmatrix, 1, quantile, prob = 0.975)
newD$pred <- apply(predmatrix, 1, quantile, prob = 0.5)
newD$Country <- "Australia"
newD$Year <- NULL
l[[5]] <- newD

# Figure G_S
g_s <- data.table(do.call(rbind, l))
g_s[, Country := factor(Country, levels = rev(c("Finland", "Poland", "Czechia", "Hungary", "Australia")))]

col3_ <- c("#357EBDFF", "#D43F3AFF", "#46B8DAFF", "#5CB85CFF", "#EEA236FF", "#9632B8FF", "#9632B8FF")[7:1]
col3__ <- col3_[3:7]
p <-
  ggplot(g_s, aes(x = StringencyIndex, y = pred, col = Country)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Country, color = NULL), alpha = .15) +
  geom_jitter(aes(y = parks_percent_change_from_baseline, fill = Country), data = s, pch = 21, col = "grey20", width = 0.7, height = 3, alpha = 0.5) +
  geom_line(lwd = 1) +
  labs(subtitle = "Mixed model per country predicitons", y = "Google Mobility\n[% change from baseline]", x = "Stringency index") +
  # scale_color_locuszoom()+
  # scale_fill_locuszoom(guide = "none")
  scale_x_continuous(breaks = round(seq(25, 75, by = 25), 1)) +
  scale_y_continuous(breaks = round(seq(-100, 200, by = 50), 1)) +
  # scale_y_continuous(breaks = round(seq(-100, 175, by = 25), 1)) +
  scale_colour_manual(
    values = col3__, guide = guide_legend(reverse = TRUE, override.aes = list(size = 0)),
    labels = paste(
      "<span style='color:",
      col3__,
      "'>",
      levels(g_s$Country),
      "</span>"
    )
  ) +
  scale_fill_manual(values = col3__, guide = "none") +
  theme_bw() +
  theme(
    legend.text = element_markdown(size = 6),
    # legend.position = "right",
    legend.title = element_blank(),
    # legend.spacing.y = unit(0.1, 'cm'),
    legend.key.height = unit(0.5, "line"),
    legend.key.size = unit(0, "line"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-10, 1, -10, -10),
    # legend.position=c(0.5,1.6),
    plot.title = element_text(color = "grey", size = 7),
    plot.subtitle = element_text(color = "grey60", size = 6),
    plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r = 0.5, unit = "pt"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = ax_lines, size = 0.25),
    axis.ticks = element_line(colour = ax_lines, size = 0.25),
    # axis.text.x = element_text()
    axis.ticks.length = unit(1, "pt"),
    axis.text = element_text(, size = 6),
    axis.title = element_text(size = 7)
  )

if (save_plot == TRUE) {
  ggsave(here::here("Outputs/Fig_7_width-70mm.png"), p + theme(plot.subtitle = element_blank()), width = 7, height = 6, unit = "cm", dpi = 600)
}

p

#' <a name="F_7">
#' **Figure 7 | Association between daily human presence in parks (Google Mobility) and the stringency of antipandemic governmental restrictions (Stringency index).**</a> Lines with shaded areas represent predictions with 95%CIs from country-specific mixed effect models controlled for the year and non-independence of data points by including weekday within the year as random intercept and Stringency index as a random slope (Table [S3e](#T_S3e)). Dots represent raw data, jittered to increase visibility, for days within which we collected escape distances in each city. Colour indicates country. Note the generally negative but weak association between Google Mobility and Stringency index.
#'
#' <a name="T_S3e">
#' **Table S3e | Google Mobility in relation to Stringency index**</a>
ll <- list()

ssf <- ss[Country == "Finland"]
fi <- lmer(
  scale(parks_percent_change_from_baseline) ~
    scale(Year) +
    scale(StringencyIndex) +
    (scale(StringencyIndex) | year_weekday),
  data = ssf, REML = FALSE
)
ll[[1]] <- m_out(name = "Table S3e - FI", dep = "Google Mobility", model = fi, nsim = 5000, R2 = FALSE)

ssp <- ss[Country == "Poland"]
pl <- lmer(
  scale(parks_percent_change_from_baseline) ~
    scale(Year) +
    scale(StringencyIndex) +
    (scale(StringencyIndex) | year_weekday),
  data = ssp, REML = FALSE
)
ll[[2]] <- m_out(name = "Table S3e - PL", dep = "Google Mobility", model = pl, nsim = 5000, R2 = FALSE)

sscz <- ss[Country == "Czechia"]
cz <- lmer(
  scale(parks_percent_change_from_baseline) ~
    scale(StringencyIndex) +
    (scale(StringencyIndex) | weekday),
  data = sscz, REML = FALSE
)
ll[[3]] <- m_out(name = "Table S3e - CZ", dep = "Google Mobility", model = cz, nsim = 5000, R2 = FALSE)

ss_h <- ss[Country == "Hungary"]
hu <- lmer(
  scale(parks_percent_change_from_baseline) ~
    scale(Year) +
    scale(StringencyIndex) +
    (scale(StringencyIndex) | year_weekday),
  data = ss_h, REML = FALSE
)
ll[[4]] <- m_out(name = "Table S3e - HU", dep = "Google Mobility", model = hu, nsim = 5000, R2 = FALSE)

ssa <- ss[Country == "Australia"]
au <- lmer(
  scale(parks_percent_change_from_baseline) ~
    scale(Year) +
    scale(StringencyIndex) +
    (scale(StringencyIndex) | year_weekday),
  # (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
  data = ssa, REML = FALSE
)
ll[[5]] <- m_out(name = "Table S3e - AU", dep = "Google Mobility", model = au, nsim = 5000, R2 = FALSE)

out_g_s <- data.table(do.call(rbind, ll))
out_g_s[is.na(out_g_s)] <- ""
out_g_s$R2_mar <- out_g_s$R2_con <- NULL
out_g_s[, effect := gsub("scale\\(Year\\)", "year", effect)]
out_g_s[, effect := gsub("scale\\(StringencyIndex\\)", "Stringency index", effect)]
out_g_s[, effect := gsub("year_weekday", "weekday within year", effect)]
out_g_s[type == "random" & grepl("Stringency index", effect, fixed = TRUE), effect := paste("Stringency index (slope) |", gsub(" Stringency index", "", effect))]
fwrite(file = here::here("Outputs/Table_S3e_rev.csv"), out_g_s)

out_g_s$error_structure <- out_g_s$response <- NULL
out_g_s[model != "", model := c(
  "Finland",
  "Poland",
  "Czechia",
  "Hungary",
  "Australia"
)]
setnames(out_g_s, old = c("estimate_r", "lwr_r", "upr_r"), new = c("estimate", "lower", "upper"))
out_g_s %>%
  kbl() %>%
  kable_paper("hover", full_width = F) %>%
  scroll_box(width = "100%")
#'
#' ### Exploration of human presence
#+ Fig_8, fig.width=16*0.393701, fig.height = 11*0.39370
g_ <- fread(here::here("Data/google_mobility.txt")) # fwrite(d, here::here('Data/data.txt'), sep ='\t')
g_[, Year := as.integer(substring(date, nchar(date) - 3, nchar(date)))]
g_[nchar(date) == 9, date := paste0("0", date)]
g_[, date_ := as.Date(date, format = "%d.%m.%Y")]
g_[, Day := yday(date_)]
setnames(g_, old = "country_region", new = "Country")
g_[, weekday := weekdays(date_)]
g_[, Country := factor(Country, levels = (c("Finland", "Poland", "Czechia", "Hungary", "Australia")))]

ann_text_gg <- data.frame(
  parks_percent_change_from_baseline = 15, Year = 2021, lab = "Text",
  Country = factor("Finland", levels = (c("Finland", "Poland", "Czechia", "Hungary", "Australia")))
)

g0 <- ggplot(g_, aes(x = parks_percent_change_from_baseline, fill = factor(Year))) +
  geom_histogram(position = "dodge", color = NA) +
  facet_wrap(~Country, nrow = 5) +
  # scale_y_continuous(trans = 'log')+
  scale_fill_manual(values = c("orange", "skyblue", "black"), guide = "none") +
  geom_vline(xintercept = 0, lty = 3, col = "#991616") +
  geom_text(data = ann_text_gg, aes(y = 100), label = "Baseline", col = "#991616", size = 1.75, hjust = 0) +
  labs(subtitle = "\nDistributions", x = "Google Mobility\n[% change in human presence]", y = "\n# of days") +
  theme_bw()+theme_MB2 +
  theme(
    plot.subtitle = element_text(size = 7, color="grey30"),
    axis.title.y = element_text(vjust=3.5),
    strip.text.x = element_text(size = 6, color = "grey30", margin = margin(1, 1, 1, 1, "mm")),
    # panel.spacing = unit(1, "mm"),
    legend.position = "none" # c(1, 0.01),
    # legend.spacing.y = unit(-0.78, "cm")
    # legend.spacing.y = unit(0.02, "cm") use if LOESS smooth text as legend
    # legend.spacing.y = unit(-0.9, "cm"),
    #axis.text.x = element_text(colour = "grey30", size = 6),
    #axis.text.y = element_text(colour = "grey30", size = 6)
  )

g1 <- ggplot(g_, aes(x = Day, y = parks_percent_change_from_baseline, col = factor(Year))) +
  geom_line() +
  facet_wrap(~Country, nrow = 5) +
  # scale_y_continuous(trans = 'log')+
  coord_cartesian(ylim = c(-100, 300)) +
  scale_color_manual(values = c("orange", "skyblue", "black"), guide = "none") +
  labs(subtitle = "\nRaw data", x = "Day\n ", y = "Google Mobility\n[% change in human presence]") +
  theme_bw()+theme_MB2 +
  theme(
    plot.subtitle = element_text(size = 7, color="grey30"),
    axis.title.y = element_text(vjust=3.5),
    strip.text.x = element_text(size = 6, color = "grey30", margin = margin(1, 1, 1, 1, "mm")),
    legend.position = "none" # c(1, 0.01),
    )


g2 <- ggplot(g_, aes(x = Day, y = parks_percent_change_from_baseline, col = factor(Year))) +
  stat_smooth(se = FALSE) +
  facet_wrap(~Country, nrow = 5) +
  # scale_y_continuous(trans = 'log')+
  coord_cartesian(ylim = c(-100, 300)) +
  scale_color_manual(values = c("orange", "skyblue", "black"), labels = c("2020", "2021", "2022 (post-COVID-19)"), name = "Year") +
  labs(subtitle = "Locally estimated\nscatterplot smoothing", x = "Day\n ") +
  theme_bw() +  
  theme(
    title = element_text(size=8, colour="grey30"),
    plot.subtitle = element_text(size = 7, color="grey30"),
  
    strip.text.x = element_text(size = 6, color = "grey30", margin = margin(1, 1, 1, 1, "mm")),
    strip.background = element_blank(), 

    legend.text = element_text(color = "grey30", size = 6),
    legend.title=element_text(size=6), 
    legend.justification = c(0, 1),
    legend.margin = margin(0,0,0,0, unit="cm"),
    legend.box.margin = margin(5, 5, 5, 5),

    #legend.spacing.y = unit(-0.9, "cm"),
    legend.background = element_blank(),
    legend.key = element_rect(colour = NA, fill = NA),
    legend.key.height= unit(0.5,"line"),
    legend.key.width = unit(0.25, "cm"),
  
    axis.title = element_text(size = 7),
    axis.text = element_text(size=6),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "grey70", linewidth = 0.1),
    axis.ticks.length = unit(1, "pt"),

    panel.spacing = unit(0, "mm"),
    panel.background=element_blank(),
    panel.border = element_rect(colour="grey70", linewidth=0.1, fill = NA),
    panel.grid = element_blank()           
  )

g1$labels$x <- g2$labels$x <- " \n "

#if(save_plot==TRUE){
#png(here::here("Outputs/Fig_8_width-160mm.png"),  width = 16, height = 11, unit = "cm", res = 600)
#(g0|g1|g2)
#grid::grid.draw(grid::textGrob('Day in the breading season\n ', y = 0.055, x = 0.6, gp = gpar(fontsize =7)))
#dev.off()
#}

g0|g1|g2; grid::grid.draw(grid::textGrob('Day in the breading season\n ', y = 0.055, x = 0.6, gp = gpar(fontsize =7)))

#' <a name="F_8">
#' **Figure 8 | Changes in human levels in parks within and between years and countries.**</a> **Left panels** represent distributions (histograms) of daily human levels (Google Mobility), **middle panels** the raw daily data across years, and **right panels** locally estimated scatterplot smoothing. Dotted vertical line in “Distribution” panels indicates baseline value of human levels, separating negative values that represent decreased human presence and positive values that indicate increased human presence when compared with the country and weekday specific baseline human levels estimated as the median from 3 January - 6 February 2020 (see also Methods; for weekday-specific patterns see Fig. [S7](#F_S7)). The zeros on the x-axes of middle and right panels represent the beginning of the breading season. Note that Google Mobility data were unavailable for the years before the COVID-19 pandemic (i.e. before 2020) but the year 2022 was without shutdowns in the studied countries.
#' 
#'  
#+ Fig_S7_gm_week_year, fig.width=18*.393701, fig.height=13*.393701
g_w <- fread(here::here("Data/google_mobility.txt")) # fwrite(d, here::here('Data/data.txt'), sep ='\t')
g_w[, Year := as.integer(substring(date, nchar(date) - 3, nchar(date)))]
g_w[nchar(date) == 9, date := paste0("0", date)]
g_w[, date_ := as.Date(date, format = "%d.%m.%Y")]
g_w[, Day := yday(date_)]
setnames(g_w, old = "country_region", new = "Country")
g_w[, weekday := weekdays(date_)]
g_w[, Country := factor(Country, levels = (c("Finland", "Poland", "Czechia", "Hungary", "Australia")))]
g_w[, weekday := factor(weekday, levels = (c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")))]
g_s5 =
ggplot(g_w, aes(x = Day, y = parks_percent_change_from_baseline, col = factor(Year))) +
    geom_line() +
    facet_grid(rows = vars(Country), cols = vars(weekday)) +
    # scale_y_continuous(trans = 'log')+
    scale_color_manual(values = c("orange", "skyblue", "black"), labels = c("2020", "2021", "2022 (post-COVID-19)"), name = "Year") +
    labs(x = "Day in the breading season", y = "Google Mobility\n[% change in human presence]")+
    theme_bw() +  
  theme(
    title = element_text(size=8, colour="grey30"),
    plot.subtitle = element_text(size = 7, color="grey30"),
  
    strip.text.x = element_text(size = 6, color = "grey30", margin = margin(1, 1, 1, 1, "mm")),
    strip.text.y = element_text(size = 6, color = "grey30", margin = margin(1, 1, 1, 1, "mm"), angle = 90),
    strip.background = element_blank(), 

    legend.text = element_text(color = "grey30", size = 6),
    legend.title=element_text(size=6), 
    legend.justification = c(0, 1),
    legend.margin = margin(0,0,0,0, unit="cm"),
    legend.box.margin = margin(5, 5, 5, 5),

    #legend.spacing.y = unit(-0.9, "cm"),
    legend.background = element_blank(),
    legend.key = element_rect(colour = NA, fill = NA),
    legend.key.height= unit(0.5,"line"),
    legend.key.width = unit(0.25, "cm"),
  
    axis.title = element_text(size = 7),
    axis.text = element_text(size=6),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "grey70", linewidth = 0.1),
    axis.ticks.length = unit(1, "pt"),

    panel.spacing = unit(0, "mm"),
    panel.background=element_blank(),
    panel.border = element_rect(colour="grey70", linewidth=0.1, fill = NA),
    panel.grid = element_blank()           
  )

if(save_plot==TRUE){
ggsave(here::here("Outputs/Fig_S7.png"), g_s5, width = 18, height = 13, unit = "cm", dpi = 600)
}
g_s5
#' <a name="F_S7">
#' **Figure S7 | Changes in daily human presence (Google Mobility) in parks across weekdays and years.**</a> Depicted are raw data connected by lines. Note that Google Mobility data were not freely available for years before the COVID-19 pandemic (i.e. before 2020) but 2022 was a year without shutdowns in the studied countries.
#' 
#' ***
#'
#' ***
#'
#' ### Testing for phylo-signal in residuals
#' Example code. Don't run, as it takes too long. To make the Table S4, load the saved mcmc outputs instead.
#+ prep_phylo, eval = FALSE, results = 'hide', warning=FALSE
# read trees 
trees = read.tree("Data/trees.tre")

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
  (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) + (scale(StringencyIndex) | Country) + (1 | IDLocality) + (1 | sp_loc),
data = s, REML = FALSE
)
s[, res := resid(m01a)]
# google
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

# # of humans
mhs <- lmer(scale(log(FID)) ~
  scale(log(SD)) +
  scale(log(FlockSize)) +
  scale(log(BodyMass)) +
  scale(sin(rad)) + scale(cos(rad)) +
  # scale(Day)+
  scale(Temp) +
  scale(Human) +
  (1 | Year) + (1 | weekday) + (1 | genus) + (1 | Species) + (1 | sp_day_year) +
  (scale(Human) | Country) + (1 | IDLocality) + (1 | sp_loc),
data = dh, REML = FALSE
)
dh[, res := resid(mhs)]

# MCC tree and covariance matrix
namesd <- d %>% distinct(Species)
row.names(namesd) <- namesd$Species
namecheck_fidd <- name.check(trees[[1]], namesd)
trees_fidd <- lapply(trees, drop.tip, namecheck_fidd$tree_not_data)
class(trees_fidd) <- "multiPhylo"
length(trees_fidd[[1]]$tip.label)
tree_fidd <- maxCladeCred(trees_fidd)
tree_fidd$node.label <- NULL
inv.phylo_d <- inverseA(tree_fidd, nodes = "ALL", scale = TRUE)
phyloresd <- solve(inv.phylo_d$Ainv)
rownames(phyloresd) <- rownames(inv.phylo_d$Ainv)

namess <- s %>% distinct(Species)
row.names(namess) <- namess$Species
namecheck_fids <- name.check(trees[[1]], namess)
trees_fids <- lapply(trees, drop.tip, namecheck_fids$tree_not_data)
class(trees_fids) <- "multiPhylo"
length(trees_fids[[1]]$tip.label)
tree_fids <- maxCladeCred(trees_fids)
tree_fids$node.label <- NULL
inv.phylo_s <- inverseA(tree_fids, nodes = "ALL", scale = TRUE)
phyloress <- solve(inv.phylo_s$Ainv)
rownames(phyloress) <- rownames(inv.phylo_s$Ainv)

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

namesdh <- dh %>% distinct(Species)
row.names(namesdh) <- namesdh$Species
namecheck_fiddh <- name.check(trees[[1]], namesdh)
trees_fiddh <- lapply(trees, drop.tip, namecheck_fiddh$tree_not_data)
class(trees_fiddh) <- "multiPhylo"
length(trees_fiddh[[1]]$tip.label)
tree_fiddh <- maxCladeCred(trees_fiddh)
tree_fiddh$node.label <- NULL
inv.phylo_dh <- inverseA(tree_fiddh, nodes = "ALL", scale = TRUE)
phyloresdh <- solve(inv.phylo_dh$Ainv)
rownames(phyloresdh) <- rownames(inv.phylo_dh$Ainv)

# brms

# period - no - need to run (takes time), if you load the above DAT_brms.Rdata
# without phylo
priors <- get_prior(res_period ~ 0 + Intercept + (1 | Species), data = d) #
mP_no = brm(
  form = res_period ~ 0 + Intercept + (1 | Species),
  data = d,
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
plot(mP_no, ask = FALSE)
pp_check(mP_no, ndraws = 100)
mcmc_plot(mP_no, type = "acf")
summary(mP_no)
# with phylo
priors <- get_prior(res_period ~ 0 + Intercept + (1 | Species) + (1 | gr(scinam, cov = A)), data = d, data2 = list(A = phyloresd))
mP_yes = brm(
  form = res_period ~ 0 + Intercept + (1 | Species) + (1 | gr(scinam, cov = A)),
  data = d,
  data2 = list(A = phyloresd),
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
plot(mP_yes, ask = FALSE)
pp_check(mP_yes, ndraws = 100)
# pp_check(m_yes, type = "scatter_avg_grouped", group = "Species") +  geom_abline(intercept = 0, slope = 1 , color = "red", lty = 2)
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
priors <- get_prior(res ~ 0 + Intercept + (1 | Species), data = s) #
m_no = brm(
  form = res ~ 0 + Intercept + (1 | Species),
  data = s,
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
plot(m_no, ask = FALSE)
pp_check(m_no, ndraws = 100)
mcmc_plot(m_no, type = "acf")
summary(m_no)
# with phylo
priors <- get_prior(res ~ 0 + Intercept + (1 | Species) + (1 | gr(scinam, cov = A)), data = s, data2 = list(A = phyloress))
m_yes = brm(
  form = res ~ 0 + Intercept + (1 | Species) + (1 | gr(scinam, cov = A)),
  data = s,
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
plot(m_yes, ask = FALSE)
pp_check(m_yes, ndraws = 100)
# pp_check(m_yes, type = "scatter_avg_grouped", group = "Species") +  geom_abline(intercept = 0, slope = 1 , color = "red", lty = 2)
mcmc_plot(m_yes, type = "acf")
summary(m_yes)

# lambda
hypothesis(m_yes, "sd_scinam__Intercept^2 / (sd_scinam__Intercept^2 + sd_Species__Intercept^2) = 0", class = NULL) # lambda
hypothesis(m_yes, "sd_scinam__Intercept^2 / (sd_scinam__Intercept^2 + sd_Species__Intercept^2 + sigma^2) = 0", class = NULL)

v_sc <- (VarCorr(m_yes, summary = FALSE)$scinam$sd)^2
v_sp <- (VarCorr(m_yes, summary = FALSE)$Species$sd)^2
v_r <- (VarCorr(m_yes, summary = FALSE)$residual$sd)^2
summary(as.mcmc(v_sc / (v_sc + v_sp + v_r)))
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

v_sc <- (VarCorr(mG_yes, summary = FALSE)$scinam$sd)^2
v_sp <- (VarCorr(mG_yes, summary = FALSE)$Species$sd)^2
v_r <- (VarCorr(mG_yes, summary = FALSE)$residual$sd)^2
summary(as.mcmc(v_sc / (v_sc + v_sp + v_r)))
Mode(as.mcmc(v_sc / (v_sc + v_sp + v_r)))

# compare stringency models
# 1. LOOic
loo(mG_yes)
loo(mG_no)

# 2. Bayes factor
bayes_factor(mG_no, mG_yes)

# 3. Poterior probability
post_prob(mG_yes, mG_no)

# # of humans  - no - need to run (takes time), if you load the above DAT_brms.Rdata
# without phylo
priors <- get_prior(res ~ 0 + Intercept + (1 | Species), data = dh) #
mH_no = brm(
  form = res ~ 0 + Intercept + (1 | Species),
  data = dh,
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
plot(mH_no, ask = FALSE)
pp_check(mH_no, ndraws = 100)
mcmc_plot(mH_no, type = "acf")
summary(mH_no)

# with phylo
priors <- get_prior(res ~ 0 + Intercept + (1 | Species) + (1 | gr(scinam, cov = A)), data = dh, data2 = list(A = phyloresdh))
mH_yes = brm(
  form = res ~ 0 + Intercept + (1 | Species) + (1 | gr(scinam, cov = A)),
  data = dh,
  data2 = list(A = phyloresdh),
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
plot(mH_yes, ask = FALSE)
pp_check(mH_yes, ndraws = 100)
# pp_check(m_yes, type = "scatter_avg_grouped", group = "Species") +  geom_abline(intercept = 0, slope = 1 , color = "red", lty = 2)
mcmc_plot(mH_yes, type = "acf")
summary(mH_yes)

# lambda
hypothesis(mH_yes, "sd_scinam__Intercept^2 / (sd_scinam__Intercept^2 + sd_Species__Intercept^2) = 0", class = NULL) # lambda
hypothesis(mH_yes, "sd_scinam__Intercept^2 / (sd_scinam__Intercept^2 + sd_Species__Intercept^2 + sigma^2) = 0", class = NULL)

v_sc <- (VarCorr(mH_yes, summary = FALSE)$scinam$sd)^2
v_sp <- (VarCorr(mH_yes, summary = FALSE)$Species$sd)^2
v_r <- (VarCorr(mH_yes, summary = FALSE)$residual$sd)^2
summary(as.mcmc(v_sc / (v_sc + v_sp + v_r)))
Mode(as.mcmc(v_sc / (v_sc + v_sp + v_r)))

# compare stringency models
# 1. LOOic
loo(mH_yes)
loo(mH_no)

# 2. Bayes factor
bayes_factor(mH_no, mH_yes)

# 3. Poterior probability
post_prob(mH_yes, mH_no)

# export models
save(file = 'Data/dat_brms_P.Rdata', mP_no, mP_yes)
save(file = 'Data/dat_brms_S.Rdata', m_no, m_yes)
save(file = 'Data/dat_brms_G.Rdata', mG_no, mG_yes)
save(file = 'Data/dat_brms_H.Rdata', mH_no, mH_yes)

load(here::here("Data/dat_brms_P.Rdata"))
load(here::here("Data/dat_brms_S.Rdata"))
load(here::here("Data/dat_brms_G.Rdata"))
load(here::here("Data/dat_brms_H.Rdata"))

v_sc <- (VarCorr(mP_yes, summary = FALSE)$scinam$sd)^2
v_sp <- (VarCorr(mP_yes, summary = FALSE)$Species$sd)^2
v_r <- (VarCorr(mP_yes, summary = FALSE)$residual$sd)^2

s_v_sc <- (VarCorr(m_yes, summary = FALSE)$scinam$sd)^2
s_v_sp <- (VarCorr(m_yes, summary = FALSE)$Species$sd)^2
s_v_r <- (VarCorr(m_yes, summary = FALSE)$residual$sd)^2

g_v_sc <- (VarCorr(mG_yes, summary = FALSE)$scinam$sd)^2
g_v_sp <- (VarCorr(mG_yes, summary = FALSE)$Species$sd)^2
g_v_r <- (VarCorr(mG_yes, summary = FALSE)$residual$sd)^2

h_v_sc <- (VarCorr(mH_yes, summary = FALSE)$scinam$sd)^2
h_v_sp <- (VarCorr(mH_yes, summary = FALSE)$Species$sd)^2
h_v_r <- (VarCorr(mH_yes, summary = FALSE)$residual$sd)^2
# summary(as.mcmc(v_sc / (v_sc + v_sp + v_r)))

ts5 <- data.table(
  "Residuals\nfrom" = c(
    "Table S2a - model 01a (Period)",
    "Table S2b - model 01a (Stringency index)",
    "Table S2c - model 01a (Google Mobility)",
    "Table S2d - model 01a (# of humans)"
  ),
  "Variance explained by phylogeny\n(95%CI)" = c(
    paste0(paste(round(quantile(as.mcmc(v_sc / (v_sc + v_sp + v_r)), probs = c(0.025, 0.975)) * 100, 2), collapse = "-"), "%"),
    paste0(paste(round(quantile(as.mcmc(s_v_sc / (s_v_sc + s_v_sp + s_v_r)), probs = c(0.025, 0.975)) * 100, 2), collapse = "-"), "%"),
    paste0(paste(round(quantile(as.mcmc(g_v_sc / (g_v_sc + g_v_sp + g_v_r)), probs = c(0.025, 0.975)) * 100, 2), collapse = "-"), "%"),
    paste0(paste(round(quantile(as.mcmc(h_v_sc / (h_v_sc + h_v_sp + h_v_r)), probs = c(0.025, 0.975)) * 100, 2), collapse = "-"), "%")
  ),
  "Bayes\nfactor" = c(
    round(bayes_factor(mP_no, mP_yes)$bf),
    round(bayes_factor(m_no, m_yes)$bf),
    round(bayes_factor(mG_no, mG_yes)$bf),
    round(bayes_factor(mH_no, mH_yes)$bf)
  ),
  "Probability of non-phylogenetic model" = c(
    round(post_prob(mP_yes, mP_no)[2], 2),
    round(post_prob(m_yes, m_no)[2], 2),
    round(post_prob(mG_yes, mG_no)[2], 2),
    round(post_prob(mH_yes, mH_no)[2], 2)
  )
)
#save(file = "Data/T_S4.Rdata", ts5)
fwrite(file = here::here("Outputs/Table_S7.csv"), ts5)

#+t_s7, echo=FALSE,results='hide', warning=FALSE, message=FALSE
#' <a name="T_S4">
#' **Table S4 | Variance in residuals explained by phylogeny and comparison of models on residulas without and with control for phylogeny.**</a>
ts5 = fread(here::here('Outputs/Table_S7.csv'))
ts5$R2_mar <- ts5$R2_con <- NULL
ts5 %>%
  kbl(align=c('l', 'r', 'r','r')) %>%
  kable_paper("hover", full_width = F)
#'
#' Testing whether residuals of the original a-models from Table [S1a](#T_S1a), [S1b](#T_S1b), [S1c](#T_S1c), and [S1d](#T_S1d) are confounded by phylogeny. **Variance explained by phylogeny (95%CI)** represents percentage of variance explained by phylogeny in an intercept only model fitted to residuals of the original models in STAN (Stan Development Team 2022) using ‘brm’ function from ‘brms’ R-package (Bürkner 2017, Bürkner 2018) with phylogeny and species as random effects. **Bayes factor** in favour of model without phylogeny and **probability of non-phylogenetic model** in comparison to a model with phylogeny. The Gelman-Rubin diagnostics was 1 for all models, indicating model convergence (Brooks & Gelman 1998). Note that for all cases, the model without phylogeny fits residuals better than a model with phylogeny, which justifies our use of simple original models, not controlled for phylogeny (Table [S1a](#T_S1a), [S1b](#T_S1b), [S1c](#T_S1c), and [S1d](#T_S1d)).
#'
#'*** 
#' 
#' ### Model assumptions  {#MS_}
#+ mode_ass_png, eval = FALSE, results = 'hide', warning=FALSE
# generates pngs with model assumptions that are later loaded within the html 
# Table S1a
m_ass(name = "Table S1a - full a", mo = mhs, dat = dh, fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Human"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1a - full b", mo = mhx, dat = dh, fixed = c("FlockSize", "BodyMass", "rad", "rad", "Temp", "Human"), trans = c("log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1a - CZ a", mo = chs, dat = dh[Country == "Czechia"], fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Human"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1a - CZ b", mo = chx, dat = dh[Country == "Czechia"], fixed = c("Year", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Human"), trans = c("", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1a - FI a", mo = fhs, dat = dh[Country == "Finland"], fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Human"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1a - FI b", mo = fhx, dat = dh[Country == "Finland"], fixed = c("Year", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Human"), trans = c("", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1a - HU a", mo = hhs, dat = dh[Country == "Hungary"], fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Human"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1a - HU b", mo = hhx, dat = dh[Country == "Hungary"], fixed = c("Year", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Human"), trans = c("", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1a - PL a", mo = phs, dat = dh[Country == "Poland"], fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Human"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1a - PL b", mo = phx, dat = dh[Country == "Poland"], fixed = c("FlockSize", "BodyMass", "rad", "rad", "Temp", "Human"), trans = c("log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))

# Table S1b
m_ass(name = "Table S1b - full a", mo = mgs, dat = ss, fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "GoogleMobility"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1b - full b", mo = mgx, dat = ss, fixed = c("Year", "FlockSize", "BodyMass", "rad", "rad", "Temp", "GoogleMobility"), trans = c("", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1b - FI a", mo = fgs, dat = ss[Country == "Finland"], fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "GoogleMobility"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1b - FI b", mo = fgx, dat = ss[Country == "Finland"], fixed = c("Year", "FlockSize", "BodyMass", "rad", "rad", "Temp", "GoogleMobility"), trans = c("", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1b - PL a", mo = pgs, dat = ss[Country == "Poland"], fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "GoogleMobility"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1b - PL b", mo = pgx, dat = ss[Country == "Poland"], fixed = c("Year", "FlockSize", "BodyMass", "rad", "rad", "Temp", "GoogleMobility"), trans = c("", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1b - CZ a", mo = cgs, dat = ss[Country == "Czechia"], fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "GoogleMobility"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1b - CZ b", mo = cgx, dat = ss[Country == "Czechia"], fixed = c("FlockSize", "BodyMass", "rad", "rad", "Temp", "GoogleMobility"), trans = c("log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1b - HU a", mo = hgs, dat = ss[Country == "Hungary"], fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "GoogleMobility"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1b - HU b", mo = hgx, dat = ss[Country == "Hungary"], fixed = c("Year", "FlockSize", "BodyMass", "rad", "rad", "Temp", "GoogleMobility"), trans = c("", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1b - AU a", mo = ags, dat = ss[Country == "Australia"], fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "GoogleMobility"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1b - AU b", mo = agx, dat = ss[Country == "Australia"], fixed = c("Year", "FlockSize", "BodyMass", "rad", "rad", "Temp", "GoogleMobility"), trans = c("", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))

# Table S1c
m_ass(name = "Table S1c - full a", mo = mss, dat = s, fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "StringencyIndex"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1c - full b", mo = msx, dat = s, fixed = c("Year", "FlockSize", "BodyMass", "rad", "rad", "Temp", "StringencyIndex"), trans = c("", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1c - FI a", mo = fss, dat = s[Country == "Finland"], fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "StringencyIndex"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1c - FI b", mo = fsx, dat = s[Country == "Finland"], fixed = c("Year", "FlockSize", "BodyMass", "rad", "rad", "Temp", "StringencyIndex"), trans = c("", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1c - PL a", mo = pss, dat = s[Country == "Poland"], fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "StringencyIndex"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1c - PL b", mo = psx, dat = s[Country == "Poland"], fixed = c("Year", "FlockSize", "BodyMass", "rad", "rad", "Temp", "StringencyIndex"), trans = c("","log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1c - CZ a", mo = css, dat = s[Country == "Czechia"], fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "StringencyIndex"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1c - CZ b", mo = csx, dat = s[Country == "Czechia"], fixed = c("FlockSize", "BodyMass", "rad", "rad", "Temp", "StringencyIndex"), trans = c("log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1c - HU a", mo = hss, dat = s[Country == "Hungary"], fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "StringencyIndex"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1c - HU b", mo = hsx, dat = s[Country == "Hungary"], fixed = c("Year", "FlockSize", "BodyMass", "rad", "rad", "Temp", "StringencyIndex"), trans = c("", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1c - AU a", mo = ass, dat = s[Country == "Australia"], fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "StringencyIndex"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1c - AU b", mo = asx, dat = s[Country == "Australia"], fixed = c("Year", "FlockSize", "BodyMass", "rad", "rad", "Temp", "StringencyIndex"), trans = c("", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
# Table S1d
m_ass(name = "Table S1d - full a", mo = ms, dat = d, fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Covid"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1d - full b", mo = mx, dat = d, fixed = c("FlockSize", "BodyMass", "rad", "rad", "Temp", "Covid"), trans = c("log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1d - CZ a", mo = cs, dat = d[Country == "Czechia"], fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Covid"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1d - CZ b", mo = cx, dat = d[Country == "Czechia"], fixed = c("FlockSize", "BodyMass", "rad", "rad", "Temp", "Covid"), trans = c("log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1d - FI a", mo = fs, dat = d[Country == "Finland"], fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Covid"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1d - FI b", mo = fx, dat = d[Country == "Finland"], fixed = c("FlockSize", "BodyMass", "rad", "rad", "Temp", "Covid"), trans = c("log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1d - HU a", mo = hs, dat = d[Country == "Hungary"], fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Covid"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1d - HU b", mo = hx, dat = d[Country == "Hungary"], fixed = c("FlockSize", "BodyMass", "rad", "rad", "Temp", "Covid"), trans = c("log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1d - PL a", mo = ps, dat = d[Country == "Poland"], fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Covid"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1d - PL b", mo = px, dat = d[Country == "Poland"], fixed = c("FlockSize", "BodyMass", "rad", "rad", "Temp", "Covid"), trans = c("log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1d - AU a", mo = as, dat = d[Country == "Australia"], fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Covid"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S1d - AU b", mo = ax, dat = d[Country == "Australia"], fixed = c("FlockSize", "BodyMass", "rad", "rad", "Temp", "Covid"), trans = c("log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))

# Table S2a
m_ass(name = "Table S2a - 1a", mo = m1a, dat = d, fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Covid"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2a - 1b", mo = m1b, dat = d, fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Covid"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2a - 1c", mo = m1c, dat = d, fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Covid"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2a - 1d", mo = m1d, dat = d, fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Covid"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2a - 1e", mo = m1d, dat = d, fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Covid"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))

dx4 <- dd[N_during > 4 & N_before > 4]
m_ass(name = "Table S2a - 2a", mo = m2a, dat = d[Species %in% dx4$Species], fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Covid"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2a - 2b", mo = m2b, dat = d[Species %in% dx4$Species], fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Covid"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2a - 2c", mo = m2c, dat = d[Species %in% dx4$Species], fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Covid"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))

dx9 = dd[N_during > 9 & N_before > 9]
m_ass(name = "Table S2a - 3a", mo = m3a, dat = d[Species %in% dx9$Species], fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Covid"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2a - 3b", mo = m3b, dat = d[Species %in% dx9$Species], fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Covid"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2a - 3c", mo = m3c, dat = d[Species %in% dx9$Species], fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Covid"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))

# Table S2b
m_ass(name = "Table S2b - 1a", mo = m01a, dat = s, fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "StringencyIndex"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2b - 1b", mo = m01b, dat = s, fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "StringencyIndex"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2b - 1c", mo = m01c, dat = s, fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "StringencyIndex"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))

m_ass(name = "Table S2b - 2a", mo = m02a, dat = s[Nsp > 4], fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "StringencyIndex"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2b - 2b", mo = m02b, dat = s[Nsp > 4], fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "StringencyIndex"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2b - 2c", mo = m02c, dat = s[Nsp > 4], fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "StringencyIndex"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))

m_ass(name = "Table S2b - 3a", mo = m03a, dat = s[Nsp > 9], fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "StringencyIndex"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2b - 3b", mo = m03b, dat = s[Nsp > 9], fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "StringencyIndex"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2b - 3c", mo = m03c, dat = s[Nsp > 9], fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "StringencyIndex"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))

# Table S2c
m_ass(name = "Table S2c - 1a", mo = mg01a, dat = ss, fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "parks_percent_change_from_baseline"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2c - 1b", mo = mg01b, dat = ss, fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "parks_percent_change_from_baseline"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2c - 1c", mo = mg01c, dat = ss, fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "parks_percent_change_from_baseline"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))

m_ass(name = "Table S2c - 2a", mo = mg02a, dat = ss[Nsp > 4], fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "parks_percent_change_from_baseline"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2c - 2b", mo = mg02b, dat = ss[Nsp > 4], fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "parks_percent_change_from_baseline"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2c - 2c", mo = mg02c, dat = ss[Nsp > 4], fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "parks_percent_change_from_baseline"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))

m_ass(name = "Table S2c - 3a", mo = mg03a, dat = ss[Nsp > 9], fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "parks_percent_change_from_baseline"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2c - 3b", mo = mg03b, dat = ss[Nsp > 9], fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "parks_percent_change_from_baseline"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2c - 3c", mo = mg03c, dat = ss[Nsp > 9], fixed = c("Year", "SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "parks_percent_change_from_baseline"), trans = c("", "log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))

# Table S2d
m_ass(name = "Table S2d - 1a", mo = mh01a, dat = dh, fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Human"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2d - 1b", mo = mh01b, dat = dh, fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Human"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2d - 1c", mo = mh01c, dat = dh, fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Human"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))

m_ass(name = "Table S2d - 2a", mo = mh02a, dat = dh[Nsp > 4], fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Human"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2d - 2b", mo = mh02b, dat = dh[Nsp > 4], fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Human"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2d - 2c", mo = mh02c, dat = dh[Nsp > 4], fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Human"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))

m_ass(name = "Table S2d - 3a", mo = mh03a, dat = dh[Nsp > 9], fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Human"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2d - 3b", mo = mh03b, dat = dh[Nsp > 9], fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Human"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S2d - 3c", mo = mh03c, dat = dh[Nsp > 9], fixed = c("SD", "FlockSize", "BodyMass", "rad", "rad", "Temp", "Human"), trans = c("log", "log", "log", "sin", "cos", "", ""), outdir = here::here("Outputs/modelAss/"))

# Table S3a
m_ass(name = "Table S3a - FI", mo = fi, dat = ssf, fixed = c("Year", "StringencyIndex"), trans = c("", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S3a - PL", mo = pl, dat = ssp, fixed = c("Year", "StringencyIndex"), trans = c("", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S3a - CZ", mo = cz, dat = sscz, fixed = c("StringencyIndex"), trans = c(""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S3a - HU", mo = hu, dat = ss_h, fixed = c("Year", "StringencyIndex"), trans = c("", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S3a - AU", mo = au, dat = ssa, fixed = c("Year", "StringencyIndex"), trans = c("", ""), outdir = here::here("Outputs/modelAss/"))

# Table S3b
m_ass(name = "Table S3b - FI", mo = shfi_, dat = shf, fixed = c("Year", "StringencyIndex"), trans = c("", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S3b - CZ", mo = shcz_, dat = shc, fixed = c("Year", "StringencyIndex"), trans = c("", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S3b - HU", mo = shhu_, dat = sh_h, fixed = c("Year", "StringencyIndex"), trans = c("", ""), outdir = here::here("Outputs/modelAss/"))

# Table S3c
m_ass(name = "Table S3c - FI", mo = shfi_ln_, dat = shf, fixed = c("Year", "StringencyIndex"), trans = c("", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S3c - CZ", mo = shcz_ln_, dat = shc, fixed = c("Year", "StringencyIndex"), trans = c("", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S3c - HU", mo = shhu_ln_, dat = sh_h, fixed = c("Year", "StringencyIndex"), trans = c("", ""), outdir = here::here("Outputs/modelAss/"))

# Table S3d
m_ass(name = "Table S3d - FI", mo = shfi_g, dat = sshf, fixed = c("Year", "parks_percent_change_from_baseline"), trans = c("", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S3d - CZ", mo = shcz_g, dat = sshc, fixed = c("Year", "parks_percent_change_from_baseline"), trans = c("", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S3d - HU", mo = shhu_g, dat = ssh_h, fixed = c("Year", "parks_percent_change_from_baseline"), trans = c("", ""), outdir = here::here("Outputs/modelAss/"))

# Table S3e
m_ass(name = "Table S3e - FI", mo = shfi_g_ln, dat = sshf, fixed = c("Year", "parks_percent_change_from_baseline"), trans = c("", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S3e - CZ", mo = shcz_g_ln, dat = sshc, fixed = c("Year", "parks_percent_change_from_baseline"), trans = c("", ""), outdir = here::here("Outputs/modelAss/"))
m_ass(name = "Table S3e - HU", mo = shhu_g_ln, dat = ssh_h, fixed = c("Year", "parks_percent_change_from_baseline"), trans = c("", ""), outdir = here::here("Outputs/modelAss/"))

#' Titles indicate specific models from specific Tables and highlight the model formula in lmer syntax   
#' <br />
#'  
#' #### for Table S1a | Escape distance in relations to # of human numbers during an escape distance trial {#MS_TS1a}
knitr::include_graphics(here::here("Outputs/modelAss/Table S1a - full a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S1a - full b.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S1a - FI a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S1a - FI b.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S1a - PL a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S1a - PL b.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S1a - CZ a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S1a - CZ b.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S1a - HU a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S1a - HU b.png"))

#' <br />
#'  
#' #### for Table S1b | Escape distance in relations to Google Mobility {#MS_TS1b}
knitr::include_graphics(here::here("Outputs/modelAss/Table S1b - full a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S1b - full b.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S1b - FI a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S1b - FI b.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S1b - PL a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S1b - PL b.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S1b - CZ a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S1b - CZ b.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S1b - HU a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S1b - HU b.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S1b - AU a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S1b - AU b.png"))
#' <br />
#'  
#' #### for Table S1c | Escape distance in relations to Stringency index {#MS_TS1c}
knitr::include_graphics(here::here("Outputs/modelAss/Table S1c - full a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S1c - full b.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S1c - FI a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S1c - FI b.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S1c - PL a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S1c - PL b.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S1c - CZ a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S1c - CZ b.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S1c - HU a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S1c - HU b.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S1c - AU a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S1c - AU b.png"))
#' <br />
#'
#' #### for Table S1d | Escape distance in relations to Period {#MS_TS1d}
knitr::include_graphics(here::here("Outputs/modelAss/Table S1d - full a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S1d - full b.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S1d - FI a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S1d - FI b.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S1d - PL a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S1d - PL b.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S1d - CZ a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S1d - CZ b.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S1d - HU a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S1d - HU b.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S1d - AU a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S1d - AU b.png"))
#' <br />
#'  
#' #### for Table S2a | Alternative models on escape distance given Period  
knitr::include_graphics(here::here("Outputs/modelAss/Table S2a - 1a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2a - 1b.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2a - 1c.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2a - 1d.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2a - 1e.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S2a - 2a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2a - 2b.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2a - 2c.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S2a - 3a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2a - 3b.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2a - 3c.png"))
#' <br />
#'  
#' #### for Table S2b | Alternative models on escape distance given Stringency index  
knitr::include_graphics(here::here("Outputs/modelAss/Table S2b - 1a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2b - 1b.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2b - 1c.png"))


knitr::include_graphics(here::here("Outputs/modelAss/Table S2b - 2a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2b - 2b.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2b - 2c.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S2b - 3a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2b - 3b.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2b - 3c.png"))
#' <br />
#'   
#' #### for Table S2c | Alternative models on escape distance given Google Mobility  
knitr::include_graphics(here::here("Outputs/modelAss/Table S2c - 1a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2c - 1b.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2c - 1c.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S2c - 2a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2c - 2b.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2c - 2c.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S2c - 3a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2c - 3b.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2c - 3c.png"))
#' <br />  
#'    
#' #### for Table S2d | Alternative models on escape distance given # of humans  
knitr::include_graphics(here::here("Outputs/modelAss/Table S2d - 1a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2d - 1b.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2d - 1c.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S2d - 2a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2d - 2b.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2d - 2c.png"))

knitr::include_graphics(here::here("Outputs/modelAss/Table S2d - 3a.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2d - 3b.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S2d - 3c.png"))
#' <br />
#'   
#' #### for Table S3a  | Google Mobility in relation to Stringency index
knitr::include_graphics(here::here("Outputs/modelAss/Table S3a - FI.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S3a - PL.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S3a - CZ.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S3a - HU.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S3a - AU.png"))
#' <br />
#'   
#' #### for Table S3b  | Number of humans in relation to Stringency index
knitr::include_graphics(here::here("Outputs/modelAss/Table S3b - FI.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S3b - CZ.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S3b - HU.png"))
#' <br />
#'   
#' #### for Table S3c | Number of humans (ln-transformed) in relation to Stringency index  
knitr::include_graphics(here::here("Outputs/modelAss/Table S3c - FI.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S3c - CZ.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S3c - HU.png"))
#' <br />
#'   
#' #### for Table S3d  | Number of humans in relation to Google Mobility
knitr::include_graphics(here::here("Outputs/modelAss/Table S3d - FI.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S3d - CZ.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S3d - HU.png"))
#' <br />
#'   
#' #### for Table S3e  | Number of humans (ln-transformed) in relation to Google Mobility  
knitr::include_graphics(here::here("Outputs/modelAss/Table S3e - FI.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S3e - CZ.png"))
knitr::include_graphics(here::here("Outputs/modelAss/Table S3e - HU.png"))
#' <br />
#'  
#' *** 
#' 
#' ### References
#' - Bulla, M., Blumstein, D.T., Benedetti, Y., Floigl, K., Jokimäki, J., Kaisanlahti-Jokimäki, M.-L., Markó, G., Morelli, F., Siretckaia, A., Szakony, S., Weston, M.A., Zeid, F.A., Tryjanowski, P., Albrecht, T. & Mikula, P. (2022). Supporting information for 'Urban birds' flight responses were unaffected by the COVID-19 shutdowns'. Open Science Framework https://doi.org/10.17605/OSF.IO/WUZH7.
#' - Bürkner PC. (2018). Advanced Bayesian multilevel modeling with the R package brms. R Journal 10:395–411. https://doi.org/10.32614/RJ-2018-017
#' - Bürkner PC. (2017). brms: An R package for Bayesian multilevel models using Stan. J Stat Softw 80:1–28. https://doi.org/10.18637/jss.v080.i01
#' - Gelman, A., Su, Y.-S., Yajima, M., Hill, J., Pittau, M., Kerman, J., Zheng, T., & Vincent, D. (2016). Data Analysis using Regression and Multilevel/Hierarchical Models. In CRAN Repository (1.8-6.; pp. 1–53). https://www.cambridge.org/highereducation/books/data-analysis-using-regression-and-multilevel-hierarchical-models/32A29531C7FD730C3A68951A17C9D983#overview
#' - Hale T., Angrist N., Goldszmidt R., Kira B., Petherick A., Phillips T., Webster S., Cameron-Blake E., Hallas L., Majumdar S., Tatlow H. (2021). A global panel database of pandemic policies (Oxford COVID-19 Government Response Tracker). Nature Human Behaviour 2021 5:4 5:529–538. https://doi.org/10.1038/s41562-021-01079-8.
#' - Peterson, B. G., & Carl, P. (2020). PerformanceAnalytics: Econometric Tools for Performance and Risk Analysis. R package version 2.0.4. https://CRAN.R-project.org/package=PerformanceAnalytics
#' - Stan Development Team (2022). Stan Modeling Language Users Guide and Reference Manual, Version 2.28. https://mc-stan.org/users/documentation/
#' - Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag, New York. https://ggplot2-book.org
#'
#'***
#'
#' ### Session info
pander(sessionInfo()) # pander(sessionInfo(), compact = FALSE)