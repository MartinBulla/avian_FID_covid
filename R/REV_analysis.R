#' ---
#' title: "Results - revision"
#' author: "Martin Bulla"
#' date: "`r Sys.time()`"
#' output:
#'     html_document:
#'         toc: true
#'         toc_float: true
#'         toc_depth: 5
#'         code_folding: hide
#'         link-citations: yes
#' ---

#+ r setup, include=FALSE
knitr::opts_chunk$set(message = FALSE, warning = FALSE, cache = FALSE)

#' ###### Code to load tools and data
#+ start, echo = T, results = 'hide', warning=FALSE
 # packages
    require(arm)
    require(data.table)
    require(effects)
    require(foreach)
    require(here)
    require(ggimage)
    require(ggplot2)
    require(ggpubr)
    require(ggsci)
    require(ggtext)
    require(grid)
    require(gtable)
    require(MASS)
    require(multcomp)
    require(optimx)
    require(performance)  
    require(PerformanceAnalytics)
    require(png)
    require(RColorBrewer)
    require(rmeta)
    require(rphylopic)
    require(scales)
    require(viridis)
 # constants
    round_ = 3 # number of decimal places to round model coefficients
    nsim = 5000 # number of simulations to extract estimates and 95%CrI
    ax_lines = "grey60" # defines color of the axis lines
    #colors <- c("#999999", "#E69F00", "#56B4E9") #viridis(3)
    set.seed(42)
    #width_ = .7 # spacing between error bars
    #col_ = c(brewer.pal(n =12, name = "Paired"), 'grey30','grey80')
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
        round_ = 3, nsim = 5000, aic = FALSE, save_sim = 'Data/model_sim/', back_tran = FALSE, perc_ = 1){
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

         ri=data.table(type='random %',effect=pred, estimate_r=round(100*q050/sum(q050)), lwr_r=round(100*q025/sum(q025)), upr_r=round(100*q975/sum(q975)))
           
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
        if(type == "mixed" & nrow(x[type=='random %' & estimate_r =='0%'])==0){
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
        png(paste(outdir,name, ".png", sep=""), width=width_,height=height_,units="in",res=600) # width = 6
        par(mfrow=c(4, n_col),tcl = -0.08, cex = 0.5, cex.main = 0.9,#ceiling(n/n_col),n_col)
            oma = c(1,1,2,1),
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
              #i = "lat_abs"
            j=sub("\\).*", "", sub(".*\\(", "",i)) 
            scatter[length(scatter)+1]=j
          }
          x = data.frame(scatter=unique(scatter)[2:length(unique(scatter))],
                          log_ = grepl("log",rownames(summary(mo)$coef)[2:length(unique(scatter))]), stringsAsFactors = FALSE)
          for (i in 1:length(fixed)){
              jj =fixed[i]
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
       
       mtext(stringr::str_wrap(paste(paste0(name," model: "), slot(mo,"call")[1],'(',slot(mo,"call")[2],sep=''), width = ceiling(nchar(paste(slot(mo,"call")[1],'(',slot(mo,"call")[2],sep=''))/2)+10), side = 3, line = 0, cex=0.5,outer = TRUE, col = 'darkblue') #ceiling(nchar(paste(slot(mo,"call")[1],'(',slot(mo,"call")[2],sep=''))/2)
       if(PNG==TRUE){dev.off()}
      }
    
  
 # data
    o  =  fread(here::here('Data/phylopic.txt'))
    setnames(o, old = c('Name', 'Code'), new = c('genus2', 'uid'))

    t = fread(here::here('Data/taxonomy.txt'))
    
    g = fread(here::here('Data/google_mobility.txt')) #fwrite(d, here::here('Data/data.txt'), sep ='\t')
    g[, Year := as.integer(substring(date, nchar(date)-3, nchar(date)))]
    g[nchar(date)==9, date:=paste0('0',date)]
    g[, date_ :=as.Date(date, format = '%d.%m.%Y')]
    g[, Day :=yday(date_)]
    g[country_region!='Ausstralia', Day := Day-92 +1] # 1 April = start of breeding season (1st day) = 92 day of the year 
    g[country_region=='Ausstralia', Day := Day-228 +1] # 15 Augusst = start of breeding season (1st day) = 228 day of the year 
    setnames(g, old = 'country_region', new ='Country')

    d = fread(here::here('Data/data_corrected.txt')) #fwrite(d, here::here('Data/data.txt'), sep ='\t')
    # add data and weekdays
    x = fread(here::here("Data/date.txt"))[, .(IDObs, Date_corr)]
    x[, date_:=as.Date(Date_corr, "%d.%m.%Y" )]
    x[, weekday := weekdays(date_)]
    d =  merge(d, x[, .(IDObs, date_, weekday)], by = "IDObs")
    #d[, date_ := as.Date(Day, origin = paste0(Year, '-01-01'))]

    # adjust correct assignment of season (Year) for Australia
    d[Country == 'Australia' & Year == 2020 & Covid == 0, Year:=2019]
    d[Country == 'Australia' & Year == 2021 & Day>139, Year:=2020]
    d[Country == 'Australia' & Year == 2022 & Day>139, Year:=2021]

    d = d[order(Year, IDLocality, Day, Hour)]
    d[Country %in% c("Czech_Republic", "Czech Republic"), Country := "Czechia"]
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
    d[, weekday := weekdays(date_)]
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
    # add google mobility
    s = merge(s, g[,.(Country,  date_, parks_percent_change_from_baseline)], all.x = TRUE)
    s[, year_ := as.character(Year)]  


    ss = s[!is.na(parks_percent_change_from_baseline)]
    ss[, country_year := paste(Country, Year)] #table(paste(s$Country, s$Year))   
    ss[parks_percent_change_from_baseline<0, google := 'before_zero']
    ss[parks_percent_change_from_baseline>0, google := 'after_zero']
    ss[, sp_country_google:= paste(sp_country, google)]

    g[, weekday := weekdays(date_)]

#' #### Title: Urban birds' flight responses were largely unaffected by the COVID-19 shutdowns

#' ## disstribution of FID ~ year
#+ fid_yr, fig.width=8, fig.height = 6
px = pp[N_during > 4 & N_before > 4]
dxx = d[paste(IDLocality, Species) %in% paste(px$IDLocality, px$Species)]
# table(dxx$IDLocality, dxx$Year)
length(unique(px$IDLocality))
length(unique(px$Species))
sum(px$N_during) + sum(px$N_before)

dxx[, Country := factor(Country, levels = rev(c("Finland", "Poland", "Czechia", "Hungary", "Australia")))]
dxx[, sp_C_loc2 := paste(gsub("[_]", " ", Species), Country, IDLocality, sep = "\n")]
dxx[, genus := sub("_.*", "", Species)]
dxx[Covid == 0, period := "before COVID-19"]
dxx[Covid == 1, period := "during COVID-19"]

col3_ <- c("#357EBDFF", "#D43F3AFF", "#46B8DAFF", "#5CB85CFF", "#EEA236FF", "#9632B8FF", "#9632B8FF")[7:1]
col3__ <- col3_[3:7]
g =
  ggplot(dxx, aes(x = as.factor(Year), y = FID, col = Country, fill = period)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~sp_C_loc2) +
  scale_y_continuous("Flight initiation distance [m]", trans = "log10") +
  scale_x_discrete("Year", guide = guide_axis(angle = 45)) +
  #scale_color_continuous() +
  scale_colour_manual(values = col3__, guide = guide_legend(reverse = TRUE))+
  scale_fill_manual(values = c("white", "lightgrey")) +
  theme_MB +
  theme(
    plot.title = element_text(size = 7),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 5, color = "grey30", margin = margin(1, 1, 1, 1, "mm")),
    # panel.spacing = unit(1, "mm"),
    legend.position = "none", # c(1, 0.01),
    legend.justification = c(1, 0),
    legend.title = element_blank(),
    # legend.spacing.y = unit(-0.78, "cm")
    # legend.spacing.y = unit(0.02, "cm") use if LOESS smooth text as legend
    legend.spacing.y = unit(-0.9, "cm"),
    axis.text.x = element_text(colour = "grey30", size = 6),
    axis.text.y = element_text(colour = "grey30", size = 6)
  )

ggsave(here::here("Outputs/Fig_2_rev_v2.png"), g, width = 18, height = 16, units = "cm")

#' **Figure 2 | Temporal and variation in the flight initiation distance across species.** Each heading denotes the scientific name of the species, country and unique site ID within each country (city). Boxplots outline colour highlights country, fill colour indicates Period (white: pre-COVID-19, grey: during COVID-19). Boxplots depict median (horizontal line inside the box), the 25th and 75th percentiles (box), the 25th and 75th percentiles ±1.5 times the interquartile range or the minimum/maximum value, whichever is smaller (bars), and the outliers (dots). Included are only species-site combinations with ≥5 observations per period. Note the log-scale in y-axis and the lack of consistent shutdown effects within and between species as well as within and between the countries.

#' ## prediction for FID ~ Period
# predictions
# full
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

# CZ - singular fits only due to genera estimated as zero (removing it changes no results)
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

# FI
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

# HU - singular fits only due to sp_loc estimated as zero (removing it changes no results)
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

# AU - singular fits only due to Year and random slope estimated as zero (removing those changes no results)
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

# PL 
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

#+ est_1, fig.width=3, fig.height = 2.5
load(here::here("Data/dat_est_rev.Rdata"))
o[predictor %in% c("scale(Covid)"), predictor := "Period"]
oo <- o[predictor %in% c("Period")]
oo[, N:=as.numeric(sub('.*N = ', '', model))]
# add meta-analytical mean
  oo_s = oo[control_for_starting_distance == 'yes']
  met = summary(meta.summaries(d = oo_s$estimate, se = oo_s$sd, method = "fixed", weights = oo_s$N))$summci
  oo_met = data.table(predictor = "Period", estimate = met[2], lwr = met[1], upr = met[3], sd = NA, model = NA, control_for_starting_distance = "yes", Country = "Combined\n(metanalytical)", N = NA)

  oo_sx = oo[control_for_starting_distance == "no"]
  metx = summary(meta.summaries(d = oo_sx$estimate, se = oo_sx$sd, method = "fixed", weights = oo_sx$N))$summci
  oo_metx = data.table(predictor = "Period", estimate = metx[2], lwr = metx[1], upr = metx[3], sd = NA, model = NA, control_for_starting_distance = "no", Country = "Combined\n(metanalytical)", N = NA)
  
  oo = rbind(oo, oo_met, oo_metx)
    
oo[, Country := factor(Country, levels = rev(c("Finland", "Poland", "Czechia", "Hungary", "Australia", "Combined\n(metanalytical)", "All\n(mixed model)")))]

width_ <- .5 # spacing between error bars

#col_ <- c(brewer.pal(n = 12, name = "Paired"), "grey30", "grey80")
#Tol_bright <- c("#EE6677", "#228833", "#4477AA", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB")
#Tol_muted <- c("#88CCEE", "#44AA99", "#117733", "#332288", "#DDCC77", "#999933", "#CC6677", "#882255", "#AA4499", "#DDDDDD")
#Tol_light <- c("#BBCC33", "#AAAA00", "#77AADD", "#EE8866", "#EEDD88", "#FFAABB", "#99DDFF", "#44BB99", "#DDDDDD")

# From Color Universal Design (CUD): https://jfly.uni-koeln.de/color/
#Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
#col_ = Okabe_Ito[7:1]
# JAMA and LocusZoom modified order
#col_ =  c("#374E55FF", "#374E55FF", "#DF8F44FF", "#79AF97FF", "#00A1D5FF", "#B24745FF",  "#80796BFF") #"#6A6599FF",
#col_ <- c("#357EBDFF", "#9632B8FF", "#46B8DAFF", "#5CB85CFF", "#EEA236FF", "#D43F3AFF", "#D43F3AFF")[7:1] # "#D43F3AFF", "#B8B8B8FF"
col_ = c("#357EBDFF", "#D43F3AFF", "#46B8DAFF", "#5CB85CFF", "#EEA236FF", "#9632B8FF", "#9632B8FF")[7:1] # "#D43F3AFF", "#B8B8B8FF"
#show_col(col_)

#p = 
ggplot(oo, aes(x = estimate, y = Country, col = Country, shape = control_for_starting_distance)) +
    geom_vline(xintercept = 0, color = "grey", linetype = "dotted") +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0, position = ggstance::position_dodgev(width_)) +
    # geom_point(position = ggstance::position_dodgev(.6)) +
    geom_point(position = position_dodge(width = width_), bg = "white", size = 1.1) +
    # scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
    # scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) +
    # geom_text( aes(x = n_pos,label = N), vjust = 0, size = 1.75, position = ggstance::position_dodgev(width_))+ # 3 positions for 3 bars
    # annotate("text", x=log10(3), y=85, label= "Used", col = "grey30", size = 2.5)+

    scale_shape_manual(name = "Controlled for\nstarting distance", guide = guide_legend(reverse = TRUE), values = c(21, 19)) +
    #scale_color_jama(guide = "none")+ #, palette = 'light'
    scale_color_manual(guide = "none", values = col_) + #guide_legend(reverse = TRUE)
    scale_x_continuous(breaks = round(seq(-0.6, 1.2, by = 0.3), 1)) +
    ylab("") +
    xlab("Standardized effect size of period\n[on flight initiation distance]") +
    # coord_cartesian(xlim = c(-.15, .15)) +
    # scale_x_continuous(breaks = round(seq(-.15, .15, by = 0.05),2)) +
    theme_bw() +
    theme(
        legend.position = "right",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        # legend.spacing.y = unit(0.1, 'cm'),
        legend.key.height = unit(0.5, "line"),
        legend.margin = margin(0, 0, 0, 0),
        # legend.position=c(0.5,1.6),
        plot.title = element_text(color = "grey", size = 7),
        plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r = 0.5, unit = "pt"),
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

#ggsave("Outputs/Fig_rev_width_CustomLocusZoom_v2.png", width = 8, height = 6, unit = "cm", dpi = 600)
#ggsave("Outputs/Fig_rev_width_CustomJAMAv1.png", width = 8, height = 6, unit = "cm", dpi = 600)
#ggsave("Outputs/Fig_rev_width_Okabe_v2.png", width = 8, height = 6, unit = "cm", dpi = 600)
#ggsave("Outputs/Fig_rev_width_UChicago_v3.png", width = 8, height = 6, unit = "cm", dpi = 600)
#ggsave("Outputs/Fig_rev_width_d3_v2.png", width = 8, height = 6, unit = "cm", dpi = 600)
#ggsave("Outputs/Fig_rev_width_nejm_v2.png", width = 8, height = 6, unit = "cm", dpi = 600)
#ggsave("Outputs/Fig_rev_width_jama_v2.png", width = 8, height = 6, unit = "cm", dpi = 600)
#ggsave("Outputs/Fig_rev_width_jco_v2.png", width = 8, height = 6, unit = "cm", dpi = 600)
#ggsave("Outputs/Fig_rev_width_npg_v2.png", width = 8, height = 6, unit = "cm", dpi = 600)

# show used colors
#gg <- ggplot_build(p)
#col_ = unique(gg$data[[3]]["colour"])$colour
#show_col(col_)

#' **Fig. 1 | The effect size with 95%CIs of Period (0 – before, 1 – during shutdowns) on flight initiation distance (ln-transformed).**  The general model on all data was controlled for starting distance (ln-transformed; filled circles) or not (empty circles), flock size (ln-transformed), temperature (also a proxy for a day within the breeding season: r Pearson = 0.48; Fig. S1) and time of day. To account for circular properties of time, time was transformed into radians (2 × time × π/24) and fitted as sine and cosine of radians (Bulla et al. 2016). All continuous variables were standardised by subtracting the mean and dividing by the standard deviation. The multicollinearity was small as correlations between predictors were  weak (r<XX, Fig. S1). To account for the non-independence of data points (Schielzeth & Forstmeier 2009; Barr et al. 2013), we attempted to fit random intercepts of year, weekday, genus, species, species at a given day and year, country, site, and species within a site, while fitting Period as random slope within site. Note that fitting Period as random slope at other random intercepts produces similar results (see Fig Sxx). We used this approach with a full dataset with all observations (n = 6369), as well as with conservative datasets, one with at least five observations per species and Period (i.e. at least five observations before and five during the COVID-19 shutdowns; n = 5260), the other with at least 10 observations per species and each Period (n = 5106). In other words, conservative datasets included species sampled in the both Periods. Albeit random structure of the full model accounts for the potential country specific effect (country accounted for xx% of variance, Table xx), depicted are also estimates from country specific models, with same predictors and random structure as the full model (excluding the random intercept of country and in case of Poland also site, as specific sites were not noted in Poland). Using meta.summaries function from rmeta R-package we used the country estimates, their standard deviation, and sample size per country to estimate the meta-analytical mean, which reflects the results based on the full mixed effect model. Also, albeit the response of humans to shutdowns were likely country specific and presence of humans in parks might have increased in some countries, decreased in others or did not change, the escape distances do not reflect such changes. In other words, the birds are either inflexible or the change in human behavior due to shutdowns was not strong enought, which might have been the case - see Fig. YY and ZZ- decide what to show in this figure (I think 2020-2022 changes in Google Mobility for each country including the fits and refering to supplementary figure)
#' 
#+ line, fig.width=4, fig.height = 6
g0 = ggplot(g, aes(x = parks_percent_change_from_baseline, fill = factor(Year))) +
  geom_histogram(position = "dodge") +
  # scale_y_continuous(trans = 'log')+
  scale_fill_manual(values = c("orange", "skyblue", "black"), guide = 'none') +
  geom_vline(xintercept = 0, lty = 3, col = "red") +
  labs(subtitle = "Distribution") +
  facet_wrap(~Country, nrow = 5)

g1 = ggplot(g, aes(x = Day, y = parks_percent_change_from_baseline, col = factor(Year))) +
  geom_line() +
  facet_wrap(~Country, nrow = 5) +
  labs(subtitle = "Raw data") +
  # scale_y_continuous(trans = 'log')+
  coord_cartesian(ylim = c(-100, 300))+
  scale_color_manual(values = c("orange", "skyblue", "black"), guide = 'none')

g2 = ggplot(g, aes(x = Day, y = parks_percent_change_from_baseline, col = factor(Year))) +
     stat_smooth() +
     facet_wrap(~Country, nrow = 5) +
     labs(subtitle = "Loess") +
     # scale_y_continuous(trans = 'log')+
     coord_cartesian(ylim = c(-100, 300))+
     scale_color_manual(values = c("orange", "skyblue", "black"))+
      theme(
         axis.text.y = element_blank(),
         axis.title.y = element_blank()
       )
  ggarrange(
    g0, g1, g2,
    ncol = 3, widths = c(1, 1, 1.1)
  )
#'
#' **Fig. ZZ | Changes in human presence (Google Mobility) in parks across year and between years.** Left plots represent the raw data, right plots LOESS smoothed curves. Google Mobility is absent for years before COVID-19. Nevertheless, 2022 was a year without shutdowns in the studied countries. Assuming that human activities levels might have been similar to pre-COVID-19 years (which may not be the case), human acctivity in the shutdown years (2020 and 2021) decreased. However, such decrease might be irrelevant for birds as the day-to-day variatioin in human presence seems larger than the general decrease in activity.
#'
#' ##  Google Mobility vs Stringency
#+ gsfig, fig.width=4.5, fig.height = 3.5
col2_ = c("#357EBDFF", "#D43F3AFF", "#46B8DAFF", "#5CB85CFF", "#EEA236FF", "#9632B8FF", "#9632B8FF")[7:1]
col2__ = col2_[3:7]
ggplot(s, aes(x = StringencyIndex, y = parks_percent_change_from_baseline, col = Country)) +
    stat_smooth(method = "lm") +
    stat_cor(method = "pearson", size = 2) +
    geom_point()+
    scale_color_manual(values = col2__) +
    labs(subtitle = "simple lm & Pearson's R")

#+ gsfig_mod, fig.width=3, fig.height = 2.5
# Predictions 
l = list()
 sc = s[Country == "Czechia"]
 cz <- lmer(parks_percent_change_from_baseline ~
    StringencyIndex + 
   (scale(StringencyIndex)|weekday),
 # (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
 data = sc, REML = FALSE
 )
 bsim <- sim(cz, n.sim = nsim)
 v <- apply(bsim@fixef, 2, quantile, prob = c(0.5))
 newD <- data.frame(StringencyIndex = seq(min(sc$StringencyIndex), max(sc$StringencyIndex), length.out = 100)) # values to predict for
 X <- model.matrix(~StringencyIndex, data = newD) # exactly the model which was used has to be specified here
 newD$pred <- (X %*% v)
 predmatrix <- matrix(nrow = nrow(newD), ncol = nsim)
 for (j in 1:nsim) predmatrix[, j] <- (X %*% bsim@fixef[j, ])
 newD$lwr <- apply(predmatrix, 1, quantile, prob = 0.025)
 newD$upr <- apply(predmatrix, 1, quantile, prob = 0.975)
 newD$pred <- apply(predmatrix, 1, quantile, prob = 0.5)
 newD$Country = 'Czechia'
 l[[1]] = newD
 
 s[, year_weekday :=paste(Year, weekday)]
 sf = s[Country == "Finland"]
 fi <- lmer(parks_percent_change_from_baseline ~
    Year+
    StringencyIndex + 
   (scale(StringencyIndex)|year_weekday),
 # (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
 data = sf, REML = FALSE
 )

 bsim <- sim(fi, n.sim = nsim)
 v <- apply(bsim@fixef, 2, quantile, prob = c(0.5))
 newD <- data.frame(Year = mean(sf$Year), StringencyIndex = seq(min(sf$StringencyIndex), max(sf$StringencyIndex), length.out = 100)) # values to predict for
 X <- model.matrix(~Year + StringencyIndex, data = newD) # exactly the model which was used has to be specified here
 newD$pred <- (X %*% v)
 predmatrix <- matrix(nrow = nrow(newD), ncol = nsim)
 for (j in 1:nsim) predmatrix[, j] <- (X %*% bsim@fixef[j, ])
 newD$lwr <- apply(predmatrix, 1, quantile, prob = 0.025)
 newD$upr <- apply(predmatrix, 1, quantile, prob = 0.975)
 newD$pred <- apply(predmatrix, 1, quantile, prob = 0.5)
 newD$Country = 'Finland'
 newD$Year = NULL
 l[[2]] = newD


 sh <- s[Country == "Hungary"]
 hu <- lmer(parks_percent_change_from_baseline ~
     Year +
     StringencyIndex +
     (scale(StringencyIndex) | year_weekday),
 # (1 | Year) + (1 | weekday) + (1|genus) + (1 | Species) + (1 | sp_day_year) + (1 | IDLocality),
 data = sh, REML = FALSE
 )
 bsim <- sim(hu, n.sim = nsim)
 v <- apply(bsim@fixef, 2, quantile, prob = c(0.5))
 newD <- data.frame(Year = mean(sh$Year), StringencyIndex = seq(min(sh$StringencyIndex), max(sh$StringencyIndex), length.out = 100)) # values to predict for
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
 
 sp <- s[Country == "Poland"]
 pl <- lmer(parks_percent_change_from_baseline ~
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

 sa <- s[Country == "Australia"]
 au <- lmer(parks_percent_change_from_baseline ~
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

sp = data.table(do.call(rbind,l))
sp[, Country := factor(Country, levels = rev(c("Finland", "Poland", "Czechia", "Hungary", "Australia")))]

#+ gsfig_pred, fig.width=4, fig.height = 3.5
col3_ = c("#357EBDFF", "#D43F3AFF", "#46B8DAFF", "#5CB85CFF", "#EEA236FF", "#9632B8FF", "#9632B8FF")[7:1]
col3__ = col3_[3:7]
#p = 
ggplot(sp, aes(x = StringencyIndex, y = pred, col = Country)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Country, color = NULL), alpha = .15) +
  geom_jitter(aes(y = parks_percent_change_from_baseline, fill = Country), data = s, pch = 21, col = 'grey20', width = 0.7, height = 3, alpha = 0.5) +
  geom_line(lwd = 1) +
  labs(subtitle = "Mixed model per country predicitons", y = "Google Mobility\n[% change from baseline]", x = "Stringency index") +
   # scale_color_locuszoom()+
   # scale_fill_locuszoom(guide = "none")
  scale_x_continuous(breaks = round(seq(25, 75, by = 25), 1)) +
  scale_y_continuous(breaks = round(seq(-100, 200, by = 50), 1)) +
  #scale_y_continuous(breaks = round(seq(-100, 175, by = 25), 1)) +
  scale_colour_manual(values = col3__, guide = guide_legend(reverse = TRUE, override.aes = list(size = 0)), 
            labels = paste("<span style='color:",
                                   col3__,
                                   "'>",
                                   levels(sp$Country),
                                   "</span>")
            ) +
  scale_fill_manual(values = col3__, guide = "none") +
  theme_bw() +
  theme(
    legend.text = element_markdown(size = 6),
    #legend.position = "right",
    legend.title = element_blank(),
    # legend.spacing.y = unit(0.1, 'cm'),
    legend.key.height = unit(0.5, "line"),
    legend.key.size = unit(0, "line"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin=margin(-10,1,-10,-10),
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

# for export
p <-
  ggplot(sp, aes(x = StringencyIndex, y = pred, col = Country)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Country, color = NULL), alpha = .15) +
  geom_jitter(aes(y = parks_percent_change_from_baseline, fill = Country), data = s, pch = 21, col = "grey20", width = 0.7, height = 3, alpha = 0.5) +
  geom_line(lwd = 1) +
  labs( y = "Google Mobility\n[% change from baseline]", x = "Stringency index") +
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
      levels(sp$Country),
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

ggsave("Outputs/Fig_G-S_rev_width_CustomLocusZoom_year-weekday_v2.png", p, width = 8, height = 6, unit = "cm", dpi = 600)

#' **Fig. 2 | Relationship betweeen Google Mobility and Stringency Index.** Lines with shaded areas represent predicted relationships from country-specific mixed effect models controlled for the year and non-independence of data points by including weekday within the year as random intercept and Stringency Index as a random slope. Dots represent raw data, jittered to increase visibility. Colors indicate country. The predictions reveal generally negative and week relationship between Gooble Mobility (human activity) in parks and seveareness of governmental restrictions. 

#+ gsfig_day, fig.width=5, fig.height = 8
ggplot(s, aes(x = StringencyIndex, y = parks_percent_change_from_baseline, col = weekday)) +
    facet_wrap(~Country, nrow = 5, scales = "free_y") +
    geom_jitter(width = 2, height = 1, pch = 21, alpha = 0.5) +
    stat_smooth(method = "lm", se = FALSE) +
    stat_cor(method = "pearson", size = 2)+
    labs(subtitle = "simple lm & Pearson's R")

#+ gm_week_year, fig.width=8, fig.height = 6
  g[, weekday := factor(weekday, levels = (c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")))]
  
  ggplot(g, aes(x = Day, y = parks_percent_change_from_baseline, col = factor(Year))) +
         geom_line() +
         facet_grid(rows = vars(Country), cols = vars(weekday)) +
         # scale_y_continuous(trans = 'log')+
         scale_color_manual(values = c("orange", "skyblue", "black"))
#'
#+ str_weekday, fig.width=2.5, fig.height = 3
s[, weekday := factor(weekday, levels = (c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")))]
ggplot(s, aes(y = StringencyIndex, x = weekday)) +
  geom_boxplot()
#'
#+ str_weekday_2, fig.width=2.5, fig.height = 3.5
ggplot(s, aes(y = StringencyIndex, x = weekday, fill = year_)) +
    geom_boxplot()
    