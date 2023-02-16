#' ---
#' title: "METHODS: FID ~ google mobiliity"
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
knitr::opts_chunk$set(message = FALSE, warning = FALSE, cache = TRUE)

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
    require(grid)
    require(gtable)
    require(MASS)
    require(multcomp)
    require(optimx)
    require(performance)  
    require(PerformanceAnalytics)
    require(png)
    require(RColorBrewer)
    require(rphylopic)
    require(viridis)
 # constants
    round_ = 3 # number of decimal places to round model coefficients
    nsim = 5000 # number of simulations to extract estimates and 95%CrI
    ax_lines = "grey60" # defines color of the axis lines
    #colors <- c("#999999", "#E69F00", "#56B4E9") #viridis(3)
    set.seed(42)
    width_ = .7 # spacing between error bars
    col_ = c(brewer.pal(n =12, name = "Paired"), 'grey30','grey80')
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
    g[country_region == "Czechia",  country_region := 'Czech Republic']
    g[, Year := as.integer(substring(date, nchar(date)-3, nchar(date)))]
    g[nchar(date)==9, date:=paste0('0',date)]
    g[, date_ :=as.POSIXct(date, format = '%d.%m.%y')]
    g[, Day :=yday(date_)]
    g[country_region!='Ausstralia', Day := Day-92 +1] # 1 April = start of breeding season (1st day) = 92 day of the year 
    g[country_region=='Ausstralia', Day := Day-228 +1] # 15 Augusst = start of breeding season (1st day) = 228 day of the year 
    setnames(g, old = 'country_region', new ='Country')

    d = fread(here::here('Data/data_corrected.txt')) #fwrite(d, here::here('Data/data.txt'), sep ='\t')
  
    # adjust correct assignment of season (Year) for Australia
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
    # add google mobility
    s = merge(s, g[,.(Country,  Year, date_,Day, parks_percent_change_from_baseline)], all.x = TRUE)
    s[, year_ := as.character(Year)]  
    s[, weekday := weekdays(date_)]

    ss = s[!is.na(parks_percent_change_from_baseline)]
    ss[, country_year := paste(Country, Year)] #table(paste(s$Country, s$Year))   
#' Google Mobility Index uses February baseline values for each city and day of the week and reports % changes for each day of the week. In other words, a 10% increase on Monday may mean differnt human mobility than 10% increase on Sunday and this may differ between cities. We can control for this in the models but it is worth keeping in mind.
#'<b>  

#' ## Distributions
#+ hist, fig.width=4, fig.height = 6
     ggplot(g, aes(x = parks_percent_change_from_baseline, fill = factor(Year))) +
       geom_histogram(position = "dodge") +
       #scale_y_continuous(trans = 'log')+
       scale_fill_manual(values = c('orange', 'skyblue', 'black')) +
       geom_vline(xintercept = 0, lty = 3, col = 'red')+
       facet_wrap(~Country, nrow = 5) 
#+ box, fig.width=6, fig.height = 3.5 
    ggplot(g, aes(x = Country, y = parks_percent_change_from_baseline, fill = factor(Year))) +
       geom_boxplot()+
       #scale_y_continuous(trans = 'log')+
       scale_fill_manual(values = c('orange', 'skyblue', 'black'))
#+ line, fig.width=4, fig.height = 6
       ggplot(g, aes(x = Day, y = parks_percent_change_from_baseline, col = factor(Year))) +
       geom_line()+
        facet_wrap(~Country, nrow = 5) +
       #scale_y_continuous(trans = 'log')+
       scale_color_manual(values = c('orange', 'skyblue', 'black'))
#+ line_2, fig.width=4, fig.height = 6
      ggplot(g, aes(x = Day, y = parks_percent_change_from_baseline, col = factor(Year))) +
       stat_smooth()+
        facet_wrap(~Country, nrow = 5) +
        labs(subtitle = 'Loess' )+
       #scale_y_continuous(trans = 'log')+
       scale_color_manual(values = c('orange', 'skyblue', 'black'))
#+ line_3, fig.width=4, fig.height = 7
      ggplot(g, aes(x = Day, y = parks_percent_change_from_baseline, col = factor(Year))) +
       stat_smooth(method = 'lm')+
        facet_wrap(~Country, nrow = 5) +
         labs(subtitle = 'linear fit' )+
       #scale_y_continuous(trans = 'log')+
       scale_color_manual(values = c('orange', 'skyblue', 'black'))

#' **!! Although histograms show slight shift (perhaps decline) in the distribution of human mobility in 2020 (21) in compariison to post covid 2022 and boxplots perhaps confirm it, the day to day variattion in human activity is far greater than the covid vs. post-covid differrences.  Although we cannot be sure because about pre-covid human mobility, the available Google Mobility data indicate that human mobility within the parks might have not changed much during COVID and hence that the relevance of our study might be compromised. **  
#' <b>  

#' ##  Google Mobility vs stringency
#+ gsfig, fig.width=5, fig.height = 3.5
   ggplot(s, aes(x = StringencyIndex, y = parks_percent_change_from_baseline, col = Country)) + 
      stat_smooth(method = 'lm', se = FALSE)+
      stat_cor(method="pearson", size =2)+
      geom_point() 

#' ## FID ~ Google Mobility
#' We are not sure whether the test between FID and Google Mobility is meaningful because it  looks at whether FID is a truly plastic trait that changes according to daily changes in human mobility. Nevertheless,  I have tested for that, in general and for each country separately. Also, as raw data indicated that there might be a quadratic effect, I spedified also quadratic models and models that use only negative Google Mobility index and only positive Google Mobility index values.
#' <b>  

#' ### Quick and dirty exploration with ggplot
#+ wd, fig.width=10, fig.height = 3
    ss[, Nsp := .N, by ='sp']  
    ss[, weekday:=factor(weekday, levels = c('Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday'))]
    ggplot(ss, aes(x = parks_percent_change_from_baseline, y = FID, col = Country, lty = year_, pch = year_)) +     
        stat_smooth(method='rlm', lwd =0.5) +
        facet_wrap(~weekday, nrow =1)
   
#+ sp, fig.width=6, fig.height = 4
    ggplot(ss[Nsp>5], aes(y = FID, x = parks_percent_change_from_baseline, groups = sp, col = Country)) + 
        stat_smooth(method = 'lm', se = FALSE) +
        labs(subtitle  = "regresions for specis with >5 data points")
        #theme(legend.position = 'none')
      # geom_point(size =0.5, pch = 1) + 
        #facet_wrap(~sp, scales ='free_y')  
        
#' ## Model outputs
# PREDICTIONS
  # full model
     ss[, country_weekday:=paste(Country,weekday)]
     ms=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
        (scale(parks_percent_change_from_baseline)|genus)+ (1|Species)+(1|weekday) + (1|sp_day_year) + (scale(parks_percent_change_from_baseline)|Country) + (1|IDLocality) +(1|sp_loc),
        data = ss, REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  

    est_ms = est_out(ms, 'ALL: (scale(google_mob)|genus)+(1|Species)+(1|weekday) +(1|sp_day_year) + (scale(google_mob)|Country) + (1|IDLocality) +(1|sp_loc)')    
    est_ms[, control_for_starting_distance:='yes']

    mx=lmer(scale(log(FID))~ 
        scale(Year)+ 
        #scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
        (scale(parks_percent_change_from_baseline)|genus)+ (1|Species)+(1|weekday) + (1|sp_day_year) + (scale(parks_percent_change_from_baseline)|Country) + (1|IDLocality) +(1|sp_loc),
        data = ss, REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  

    est_mx = est_out(mx, 'ALL: (scale(google_mob)|genus)+(1|Species)+(1|weekday) +(1|sp_day_year) + (scale(google_mob)|Country) + (1|IDLocality) +(1|sp_loc)')    
    est_mx[, control_for_starting_distance:='no']
  
    msp=lmer(scale(log(FID))~  # singular due to weekday, else all ok
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(poly(parks_percent_change_from_baseline,2))+ 
        (poly(parks_percent_change_from_baseline,2)|genus)+ (1|Species)+(1|weekday)+(1|sp_day_year) + (poly(parks_percent_change_from_baseline,2)|Country) + (1|IDLocality) +(1|sp_loc),
        data = ss, REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
    est_msp = est_out(msp, 'ALL: (poly(google_mob,2)|genus)+(1|Species)+(1|weekday) +(1|sp_day_year) + (poly(google_mob,2)|Country) + (1|IDLocality) +(1|sp_loc)') 
    est_msp[, control_for_starting_distance:='yes']

    mxp=lmer(scale(log(FID))~  # singular due to weekday, else all ok
        scale(Year)+ 
       # scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(poly(parks_percent_change_from_baseline,2))+ 
        (poly(parks_percent_change_from_baseline,2)|genus)+ (1|Species)+(1|weekday)+(1|sp_day_year) + (poly(parks_percent_change_from_baseline,2)|Country) + (1|IDLocality) +(1|sp_loc),
        data = ss, REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
    est_mxp = est_out(mxp, 'ALL: (poly(google_mob,2)|genus)+(1|Species)+(1|weekday) +(1|sp_day_year) + (poly(google_mob,2)|Country) + (1|IDLocality) +(1|sp_loc)') 
    est_mxp[, control_for_starting_distance:='no']

  # PL
    ps=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
        (1|weekday)+(1|Species)+(1|sp_day_year),
        data = ss[Country=='Poland'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
    est_ps = est_out(ps, 'Poland: (1|weekday)+(1|Species)+(1|sp_day_year)')    
    est_ps[, control_for_starting_distance:='yes']
    px=lmer(scale(log(FID))~ 
        scale(Year)+ 
        #scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
        (1|weekday)+(1|Species)+(1|sp_day_year),
        data = ss[Country=='Poland'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
    est_px = est_out(px, 'Poland: (1|weekday)+(1|Species)+(1|sp_day_year)')    
    est_px[, control_for_starting_distance:='no']

    psp=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(poly(parks_percent_change_from_baseline,2))+ 
        (1|weekday)+(1|Species)+(1|sp_day_year),
        data = ss[Country=='Poland'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
    est_psp = est_out(psp, 'Poland - quadratic: (1|weekday)+(1|Species)+(1|sp_day_year)')    
    est_psp[, control_for_starting_distance:='yes']   
    
    pxp=lmer(scale(log(FID))~ 
        scale(Year)+ 
        #scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(poly(parks_percent_change_from_baseline,2))+ 
        (1|weekday)+(1|Species)+(1|sp_day_year),
        data = ss[Country=='Poland'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  

    est_pxp = est_out(pxp, 'Poland - quadratic: (1|weekday)+(1|Species)+(1|sp_day_year)')    
    est_pxp[, control_for_starting_distance:='no']
  # HU
   hs=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
       (1|weekday) + (1|Species)+(1|sp_day_year) + (1|IDLocality),
        data = ss[Country=='Hungary'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
    est_hs = est_out(hs, 'Hungary: (1|weekday)+(1|Species)+(1|sp_day_year)+ (1|IDLocality)')    
    est_hs[, control_for_starting_distance:='yes']
    hx=lmer(scale(log(FID))~ 
        scale(Year)+ 
       # scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
       (1|weekday) + (1|Species)+(1|sp_day_year) + (1|IDLocality),
        data = ss[Country=='Hungary'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )      
    est_hx = est_out(hx, 'Hungary: (1|weekday)+(1|Species)+(1|sp_day_year)+ (1|IDLocality)')    
    est_hx[, control_for_starting_distance:='no']

    hsp=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(poly(parks_percent_change_from_baseline,2))+ 
       (1|weekday) + (1|Species)+(1|sp_day_year) + (1|IDLocality),
        data = ss[Country=='Hungary'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
    est_hsp = est_out(hsp, 'Hungary - quadratic: (1|weekday)+(1|Species)+(1|sp_day_year)+ (1|IDLocality)')    
    est_hsp[, control_for_starting_distance:='yes']
    
    hxp=lmer(scale(log(FID))~ 
        scale(Year)+ 
       # scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(poly(parks_percent_change_from_baseline,2))+  
       (1|weekday) + (1|Species)+(1|sp_day_year) + (1|IDLocality),
        data = ss[Country=='Hungary'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )      
    est_hxp = est_out(hxp, 'Hungary - quadratic: (1|weekday)+(1|Species)+(1|sp_day_year)+ (1|IDLocality)')    
    est_hxp[, control_for_starting_distance:='no']
  # CZ
   cs=lmer(scale(log(FID))~ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
       (1|weekday) + (1|Species)+(1|sp_day_year) + (1|IDLocality),
        data = ss[Country=='Czech Republic'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
    est_cs = est_out(cs, 'Czechia: (1|weekday)+(1|Species)+(1|sp_day_year)+ (1|IDLocality)')    
    est_cs[, control_for_starting_distance:='yes']
    cx=lmer(scale(log(FID))~ 
       # scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
       (1|weekday) + (1|Species)+(1|sp_day_year) + (1|IDLocality),
        data = ss[Country=='Czech Republic'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )      
    est_cx = est_out(cx, 'Czechia: (1|weekday)+(1|Species)+(1|sp_day_year)+ (1|IDLocality)')    
    est_cx[, control_for_starting_distance:='no']

    csp=lmer(scale(log(FID))~ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(poly(parks_percent_change_from_baseline,2))+       
        (1|weekday) + (1|Species)+(1|sp_day_year) + (1|IDLocality),
        data = ss[Country=='Czech Republic'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
    est_csp = est_out(csp, 'Czechia - quadratic: (1|weekday)+(1|Species)+(1|sp_day_year)+ (1|IDLocality)')    
    est_csp[, control_for_starting_distance:='yes']
    
    cxp=lmer(scale(log(FID))~ 
       # scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(poly(parks_percent_change_from_baseline,2))+ 
       (1|weekday) + (1|Species)+(1|sp_day_year) + (1|IDLocality),
        data = ss[Country=='Czech Republic'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )      
    est_cxp = est_out(cxp, 'Czechia - quadratic: (1|weekday)+(1|Species)+(1|sp_day_year)+ (1|IDLocality)')    
    est_cxp[, control_for_starting_distance:='no']

  # FI
   fs=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
       (1|weekday) + (1|Species)+(1|sp_day_year) + (1|IDLocality),
        data = ss[Country=='Finland'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
    est_fs = est_out(fs, 'Finland: (1|weekday)+(1|Species)+(1|sp_day_year)+ (1|IDLocality)')    
    est_fs[, control_for_starting_distance:='yes']
    fx=lmer(scale(log(FID))~ 
        scale(Year)+ 
       # scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
       (1|weekday) + (1|Species)+(1|sp_day_year) + (1|IDLocality),
        data = ss[Country=='Finland'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )      
    est_fx = est_out(fx, 'Finland: (1|weekday)+(1|Species)+(1|sp_day_year)+ (1|IDLocality)')    
    est_fx[, control_for_starting_distance:='no']

    fsp=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(poly(parks_percent_change_from_baseline,2))+ 
       (1|weekday) + (1|Species)+(1|sp_day_year) + (1|IDLocality),
        data = ss[Country=='Finland'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
    est_fsp = est_out(fsp, 'Finland - quadratic: (1|weekday)+(1|Species)+(1|sp_day_year)+ (1|IDLocality)')    
    est_fsp[, control_for_starting_distance:='yes']
    
    fxp=lmer(scale(log(FID))~ 
        scale(Year)+ 
       # scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(poly(parks_percent_change_from_baseline,2))+  
       (1|weekday) + (1|Species)+(1|sp_day_year) + (1|IDLocality),
        data = ss[Country=='Finland'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )      
    est_fxp = est_out(fxp, 'Finland - quadratic: (1|weekday)+(1|Species)+(1|sp_day_year)+ (1|IDLocality)')    
    est_fxp[, control_for_starting_distance:='no']



  # AU
   as=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
       (1|weekday) + (1|Species)+(1|sp_day_year) + (1|IDLocality),
        data = ss[Country=='Australia'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
    est_as = est_out(as, 'Australia: (1|weekday)+(1|Species)+(1|sp_day_year)+ (1|IDLocality)')    
    est_as[, control_for_starting_distance:='yes']
    ax=lmer(scale(log(FID))~ 
        scale(Year)+ 
       # scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
       (1|weekday) + (1|Species)+(1|sp_day_year) + (1|IDLocality),
        data = ss[Country=='Australia'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )      
    est_ax = est_out(ax, 'Australia: (1|weekday)+(1|Species)+(1|sp_day_year)+ (1|IDLocality)')    
    est_ax[, control_for_starting_distance:='no']

    asp=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(poly(parks_percent_change_from_baseline,2))+ 
       (1|weekday) + (1|Species)+(1|sp_day_year) + (1|IDLocality),
        data = ss[Country=='Australia'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
    est_asp = est_out(asp, 'Australia - quadratic: (1|weekday)+(1|Species)+(1|sp_day_year)+ (1|IDLocality)')    
    est_asp[, control_for_starting_distance:='yes']
    
    axp=lmer(scale(log(FID))~ 
        scale(Year)+ 
       # scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(poly(parks_percent_change_from_baseline,2))+  
       (1|weekday) + (1|Species)+(1|sp_day_year) + (1|IDLocality),
        data = ss[Country=='Australia'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )      
    est_axp = est_out(axp, 'Australia - quadratic: (1|weekday)+(1|Species)+(1|sp_day_year)+ (1|IDLocality)')    
    est_axp[, control_for_starting_distance:='no']

  # combine
    o = rbind(est_ms, est_mx, est_msp, est_mxp,
                    est_as, est_ax, est_asp, est_axp,
                    est_cs, est_cx, est_csp, est_cxp,
                    est_hs, est_hx, est_hsp, est_hxp,
                    est_ps, est_px, est_psp, est_pxp,
                    est_fs, est_fx, est_fsp, est_fxp)
    save(o, file = here::here('Data/dat_est_rev.Rdata'))
  
#' ### compare linear and quadratic model with AIC
    AIC(hs,hsp) 
    AIC(cs,csp) 
    AIC(ps,psp) 
    AIC(as,asp) 
    AIC(fs,fsp) 
#' Quadratic  (indicated with 'p' in the model name) is never  better than linear.  

#' ### PLOT estimates
#+ est_1, fig.width=10, fig.height = 5
  load(here::here('Data/dat_est_rev.Rdata'))
  o[predictor%in%c('scale(poly(parks_percent_change_from_baseline, 2))1','scale(parks_percent_change_from_baseline)'), predictor:='google mobility\n(linear)']
  o[predictor%in%'scale(poly(parks_percent_change_from_baseline, 2))2', predictor:='google mobility\n(quadratic)']
  oo = o[predictor %in% c('google mobility\n(linear)', 'google mobility\n(quadratic)')]
  oo[, predictor := factor(predictor, levels = rev(c('google mobility\n(linear)', 'google mobility\n(quadratic)')))]
 ggplot(oo, aes(x = estimate, y = predictor, color = model, shape = control_for_starting_distance)) +
    geom_vline(xintercept = 0, color = "grey", linetype = "dotted") +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0, position = ggstance::position_dodgev(width_)) +
    #geom_point(position = ggstance::position_dodgev(.6)) +
    geom_point(position = position_dodge(width = width_), bg = 'white', size = 1.1) +
    #scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
    #scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
    #geom_text( aes(x = n_pos,label = N), vjust = 0, size = 1.75, position = ggstance::position_dodgev(width_))+ # 3 positions for 3 bars
    # annotate("text", x=log10(3), y=85, label= "Used", col = "grey30", size = 2.5)+
  
    scale_shape_manual(name = 'Controlled for starting distance',  guide = guide_legend(reverse = TRUE), values = c(23,21))+
    scale_color_manual(name = "Model", guide = guide_legend(reverse = TRUE), values = col_)+
    scale_x_continuous(breaks = round(seq(-0.6, 0.9, by = 0.1),1))+
    ylab("") +
    xlab("Standardized effect size\n[on flight initiation distance]") +
    #coord_cartesian(xlim = c(-.15, .15)) +
   #scale_x_continuous(breaks = round(seq(-.15, .15, by = 0.05),2)) +
    theme_bw() +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        # legend.spacing.y = unit(0.1, 'cm'),
        legend.key.height = unit(0.5, "line"),
        legend.margin=margin(0,0,0,0),
        #legend.position=c(0.5,1.6),
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
  
#' ### PLOT estimates for linear models only
#+ est_2, fig.width=10, fig.height = 3  
  load(here::here('Data/dat_est_rev.Rdata'))
  o[predictor%in%c('poly(parks_percent_change_from_baseline, 2)1','scale(parks_percent_change_from_baseline)'), predictor:='google mobility\n(linear)']
  oo = o[predictor %in% c('google mobility\n(linear)') & !model%in%o[predictor%in%'poly(parks_percent_change_from_baseline, 2)2', unique(model)]]
  ggplot(oo, aes(x = estimate, y = model, shape = control_for_starting_distance)) +
    geom_vline(xintercept = 0, color = "grey", linetype = "dotted") +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0, position = ggstance::position_dodgev(width_)) +
    #geom_point(position = ggstance::position_dodgev(.6)) +
    geom_point(position = position_dodge(width = width_), bg = 'white', size = 1.1) +
    #scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
    #scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
    #geom_text( aes(x = n_pos,label = N), vjust = 0, size = 1.75, position = ggstance::position_dodgev(width_))+ # 3 positions for 3 bars
    # annotate("text", x=log10(3), y=85, label= "Used", col = "grey30", size = 2.5)+
  
    scale_shape_manual(name = 'Controlled for starting distance',  guide = guide_legend(reverse = TRUE), values = c(23,21))+
   # scale_color_manual(name = "Model", guide = guide_legend(reverse = TRUE), values = col_)+
    ylab("") +
    xlab("Standardized effect size\n[on flight initiation distance]") +
    #coord_cartesian(xlim = c(-.15, .15)) +
   #scale_x_continuous(breaks = round(seq(-.15, .15, by = 0.05),2)) +
    theme_bw() +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        # legend.spacing.y = unit(0.1, 'cm'),
        legend.key.height = unit(0.5, "line"),
        legend.margin=margin(0,0,0,0),
        #legend.position=c(0.5,1.6),
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
  
#' ### SHOW that quadratic models missfit the data
#+ est_3, fig.width=4, fig.height = 3
    # generate predictions for each country 
    nsim_= 5000
    ss[, sin_ :=sin(rad)]
    ss[, cos_ :=cos(rad)]
    ss[, poly_lin := poly(parks_percent_change_from_baseline,2)[,1]]
    ss[, poly_qua := poly(parks_percent_change_from_baseline,2)[,2]]
    # HU, FI
        p1 = foreach(j = c( 'Finland','Hungary'), .combine = rbind) %do% {
        # j = 'Finland'
        ssp = ss[Country == j]
        psp=lmer(FID_ln~ 
        Year+ 
        SD_ln + 
        flock_ln +
        body_ln +
        sin_ + cos_  +
        Temp+ 
        poly_lin + poly_qua + 
            (1|weekday)+(1|Species)+(1|sp_day_year) + (1|IDLocality),
            data = ssp, REML = FALSE,  
            control = lmerControl( 
                optimizer ='optimx', optCtrl=list(method='nlminb')) 
            )   
        bsim = sim(psp, n.sim=nsim_)  
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5)) # coefficients
        newD=data.frame(
            Year = mean(ssp$Year),
            SD_ln = mean(ssp$SD_ln),
            flock_ln = mean(ssp$flock_ln),
            body_ln = mean(ssp$body_ln),
            rad = mean(ssp$rad),
            Temp = mean(ssp$Temp),
            parks_percent_change_from_baseline = seq(min(ssp$parks_percent_change_from_baseline), max(ssp$parks_percent_change_from_baseline), length.out = 100)) # values to predict for
        newD$sin_ = sin(newD$rad)
        newD$cos_ = cos(newD$rad)
        newD$poly_lin = poly(newD$parks_percent_change_from_baseline,2)[,1]    
        newD$poly_qua = poly(newD$parks_percent_change_from_baseline,2)[,2]
        
        X <- model.matrix(~ 
        Year+ 
        SD_ln + 
        flock_ln +
        body_ln +
        sin_ + cos_  +
        Temp+ 
        poly_lin + poly_qua 
        , data=newD) # exactly the model which was used has to be specified here 
        
        # calculate predicted values and creditability intervals
        newD$pred <-X%*%v #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim) 
                for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
                newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
        newD$Country = j
        return(newD)
        }
    # Poland 
    ssp = ss[Country=='Poland']
    psp=lmer(FID_ln~ 
        Year+ 
        SD_ln + 
        flock_ln +
        body_ln +
        sin_ + cos_  +
        Temp+ 
        poly_lin + poly_qua + 
            (1|weekday)+(1|Species)+(1|sp_day_year),
            data = ssp, REML = FALSE,  
            control = lmerControl( 
                optimizer ='optimx', optCtrl=list(method='nlminb')) 
            )   

        bsim = sim(psp, n.sim=nsim_)  
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5)) # coefficients
        newD=data.frame(
            Year = mean(ssp$Year),
            SD_ln = mean(ssp$SD_ln),
            flock_ln = mean(ssp$flock_ln),
            body_ln = mean(ssp$body_ln),
            rad = mean(ssp$rad),
            Temp = mean(ssp$Temp),
            parks_percent_change_from_baseline = seq(min(ssp$parks_percent_change_from_baseline), max(ssp$parks_percent_change_from_baseline), length.out = 100)) # values to predict for
        newD$sin_ = sin(newD$rad)
        newD$cos_ = cos(newD$rad)
        newD$poly_lin = poly(newD$parks_percent_change_from_baseline,2)[,1]    
        newD$poly_qua = poly(newD$parks_percent_change_from_baseline,2)[,2]
        
        X <- model.matrix(~ 
            Year+ 
        SD_ln + 
        flock_ln +
        body_ln +
        sin_ + cos_  +
        Temp+ 
        poly_lin + poly_qua 
            , data=newD) # exactly the model which was used has to be specified here 
        
        # calculate predicted values and creditability intervals
        newD$pred <-X%*%v #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim) 
                for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
                newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
        newD$Country = 'Poland'
        p2 = newD

    # CZ 
    ssp = ss[Country=='Czech Republic']
    psp=lmer(FID_ln~ 
        SD_ln + 
        flock_ln +
        body_ln +
        sin_ + cos_  +
        Temp+ 
        poly_lin + poly_qua + 
            (1|weekday)+(1|Species)+(1|sp_day_year)+ (1|IDLocality),,
            data = ssp, REML = FALSE,  
            control = lmerControl( 
                optimizer ='optimx', optCtrl=list(method='nlminb')) 
            )   

        bsim = sim(psp, n.sim=nsim_)  
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5)) # coefficients
        newD=data.frame(
            SD_ln = mean(ssp$SD_ln),
            flock_ln = mean(ssp$flock_ln),
            body_ln = mean(ssp$body_ln),
            rad = mean(ssp$rad),
            Temp = mean(ssp$Temp),
            parks_percent_change_from_baseline = seq(min(ssp$parks_percent_change_from_baseline), max(ssp$parks_percent_change_from_baseline), length.out = 100)) # values to predict for
        newD$sin_ = sin(newD$rad)
        newD$cos_ = cos(newD$rad)
        newD$poly_lin = poly(newD$parks_percent_change_from_baseline,2)[,1]    
        newD$poly_qua = poly(newD$parks_percent_change_from_baseline,2)[,2]
        
        X <- model.matrix(~ 
        SD_ln + 
        flock_ln +
        body_ln +
        sin_ + cos_  +
        Temp+ 
        poly_lin + poly_qua 
            , data=newD) # exactly the model which was used has to be specified here 
        
        # calculate predicted values and creditability intervals
        newD$pred <-X%*%v #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim) 
                for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
                newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
        newD$Country = 'Czech Republic'
        newD$Year = 2021
        p3= newD

    # AU
        ssp = ss[Country == 'Australia']
        psp=lmer(FID_ln~ 
        Year+ 
        SD_ln + 
        flock_ln +
        body_ln +
        sin_ + cos_  +
        Temp+ 
        poly_lin + poly_qua + 
        (1|Species)+(1|sp_day_year) + (1|IDLocality),
            data = ssp, REML = FALSE,  
            control = lmerControl( 
                optimizer ='optimx', optCtrl=list(method='nlminb')) 
            )   
        bsim = sim(psp, n.sim=nsim_)  
        v = apply(bsim@fixef, 2, quantile, prob=c(0.5)) # coefficients
        newD=data.frame(
            Year = mean(ssp$Year),
            SD_ln = mean(ssp$SD_ln),
            flock_ln = mean(ssp$flock_ln),
            body_ln = mean(ssp$body_ln),
            rad = mean(ssp$rad),
            Temp = mean(ssp$Temp),
            parks_percent_change_from_baseline = seq(min(ssp$parks_percent_change_from_baseline), max(ssp$parks_percent_change_from_baseline), length.out = 100)) # values to predict for
        newD$sin_ = sin(newD$rad)
        newD$cos_ = cos(newD$rad)
        newD$poly_lin = poly(newD$parks_percent_change_from_baseline,2)[,1]    
        newD$poly_qua = poly(newD$parks_percent_change_from_baseline,2)[,2]
        
        X <- model.matrix(~ 
        Year+ 
        SD_ln + 
        flock_ln +
        body_ln +
        sin_ + cos_  +
        Temp+ 
        poly_lin + poly_qua 
        , data=newD) # exactly the model which was used has to be specified here 
        
        # calculate predicted values and creditability intervals
        newD$pred <-X%*%v #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim) 
                for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
                newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
        newD$Country = 'Australia'
        p4 = newD
    
    # combine
        p = data.table(rbind(p1,p2,p3,p4))
        p [, FID:=exp(pred)]
        p [, FID_lwr:=exp(lwr)]
        p [, FID_upr:=exp(upr)]
        p[, Google_mobility:=parks_percent_change_from_baseline]  

  # plot
    ggplot(p, aes(x = Google_mobility, y = pred, col = Country, fill = Country)) + 
      geom_ribbon(aes(ymin=lwr, ymax=upr, x=Google_mobility, fill = Country),  alpha = 0.2, show.legend = NA, colour = NA) +
      geom_line(aes(Google_mobility, y = pred, col = Country)) +
      #scale_fill_viridis_c(option = "plasma", name = "°C\nwhen\npredated") +
      #scale_x_continuous(expand = c(0, 0), lim = c(30,60),  name = "Mid-day temperature [°C]") +
      scale_y_continuous(  name = "ln(Flight initiation distance)") +
      #coord_cartesian(xlim = c(30,61), clip = 'off') + 
      #labs(tag = "(c)") +
      theme_MB +
      theme(  
            legend.text=element_text(size=5),
            legend.title=element_text(size=6, hjust = 0.5)
            )  
  
# ' ### Global quadratic model
#+ est_2a, fig.width=3, fig.height = 3
 # generate predictions foo full model
   nsim_= 5000
   ss[, sin_ :=sin(rad)]
   ss[, cos_ :=cos(rad)]
   ss[, poly_lin := poly(parks_percent_change_from_baseline,2)[,1]]
   ss[, poly_qua := poly(parks_percent_change_from_baseline,2)[,2]]
   ssp =ss
   psp=lmer(FID_ln~ 
       Year+ 
       SD_ln + 
       flock_ln +
       body_ln +
       sin_ + cos_  +
       Temp+ 
       poly_lin + poly_qua + 
       (poly(parks_percent_change_from_baseline,2)|genus)+ (1|Species)+(1|weekday)+(1|sp_day_year) + (poly(parks_percent_change_from_baseline,2)|Country) + (1|IDLocality) +(1|sp_loc),
        data = ssp, REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )   
     bsim = sim(psp, n.sim=nsim_)  
     v = apply(bsim@fixef, 2, quantile, prob=c(0.5)) # coefficients
     newD=data.frame(
        Year = mean(ssp$Year),
        SD_ln = mean(ssp$SD_ln),
        flock_ln = mean(ssp$flock_ln),
        body_ln = mean(ssp$body_ln),
        rad = mean(ssp$rad),
        Temp = mean(ssp$Temp),
        parks_percent_change_from_baseline = seq(min(ssp$parks_percent_change_from_baseline), max(ssp$parks_percent_change_from_baseline), length.out = 100)) # values to predict for
     newD$sin_ = sin(newD$rad)
     newD$cos_ = cos(newD$rad)
     newD$poly_lin = poly(newD$parks_percent_change_from_baseline,2)[,1]    
     newD$poly_qua = poly(newD$parks_percent_change_from_baseline,2)[,2]
    
     X <- model.matrix(~ 
       Year+ 
       SD_ln + 
       flock_ln +
       body_ln +
       sin_ + cos_  +
       Temp+ 
       poly_lin + poly_qua 
      , data=newD) # exactly the model which was used has to be specified here 
    
     # calculate predicted values and creditability intervals
      newD$pred <-X%*%v #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
      predmatrix <- matrix(nrow=nrow(newD), ncol=nsim) 
			for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
			newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
			newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
      p= data.table(newD)
      p [, FID:=exp(pred)]
      p [, FID_lwr:=exp(lwr)]
      p [, FID_upr:=exp(upr)]
      p[, Google_mobility:=parks_percent_change_from_baseline]  
   # plot
    ggplot(p, aes(x = Google_mobility, y = pred)) + 
      geom_ribbon(aes(ymin=lwr, ymax=upr, x=Google_mobility),  alpha = 0.2, show.legend = NA, colour = NA) +
      geom_line(aes(Google_mobility, y = pred)) +
      #scale_fill_viridis_c(option = "plasma", name = "°C\nwhen\npredated") +
      #scale_x_continuous(expand = c(0, 0), lim = c(30,60),  name = "Mid-day temperature [°C]") +
      scale_y_continuous(  name = "ln(Flight initiation distance)") +
      #coord_cartesian(xlim = c(30,61), clip = 'off') + 
      #labs(tag = "(c)") +
      theme_MB
#+ est_2b, fig.width=3, fig.height = 3
    ggplot(p, aes(x = Google_mobility, y = FID)) + 
      geom_ribbon(aes(ymin=FID_lwr, ymax=FID_upr, x=Google_mobility),  alpha = 0.2, show.legend = NA, colour = NA) +
      geom_line(aes(Google_mobility, y = FID)) +
      #scale_fill_viridis_c(option = "plasma", name = "°C\nwhen\npredated") +
      #scale_x_continuous(expand = c(0, 0), lim = c(30,60),  name = "Mid-day temperature [°C]") +
      scale_y_continuous(  name = "Flight initiation distance [m]") +
      #coord_cartesian(xlim = c(30,61), clip = 'off') + 
      #labs(tag = "(c)") +
      theme_MB

#' ### SPECIES-specific regressions
#+ est_4, fig.width=4, fig.height = 3.5 
 # predictions 
  ssc <- ss[, sp_country_N := .N, by = sp_country]
  ssc = ssc[sp_country_N>=15]
  ssc[, log_FID :=log(FID)]
  ssc[, log_SD :=log(SD)]
  nsim = 1000
  vv = list()
  u = foreach(i = unique(ssc$sp_country), .combine = rbind) %do% {
    # i = "Luscinia_megarhynchos\nHungary" "Gallinula_tenebrosa\nAustralia" "Columba_livia\nPoland" #ssc$sp_country[1]
    ssci = ssc[sp_country == i]
    if(ssci$Country[1] =='Poland' | ssci$sp_country[1]%in%c( "Gallinula_tenebrosa\nAustralia","Luscinia_megarhynchos\nHungary")){
      mi <- lmer(log_FID ~
        log_SD +
        parks_percent_change_from_baseline +
        (1 | weekday), data = ssci, REML = FALSE)
    }else{
      mi = lmer(log_FID ~ 
            log_SD+ 
            parks_percent_change_from_baseline +
            (1 | weekday) + (1 | IDLocality), data = ssci, REML = FALSE)
    }        
    bsim = sim(mi, n.sim = nsim)
    v = apply(bsim@fixef, 2, quantile, prob = c(0.5)) # coefficients
    newD = data.frame(
      log_SD = mean(ssci$log_SD),
      parks_percent_change_from_baseline = seq(min(ssci$parks_percent_change_from_baseline), max(ssci$parks_percent_change_from_baseline), length.out = 100)
    ) # values to predict for
    X <- model.matrix(~log_SD +parks_percent_change_from_baseline, data = newD) # exactly the model which was used has to be specified here

    # calculate predicted values and creditability intervals
    newD$pred <- (X %*% v)
    predmatrix <- matrix(nrow = nrow(newD), ncol = nsim)
    for (j in 1:nsim) {
      predmatrix[, j] <- (X %*% bsim@fixef[j, ])
    }
    #predmatrix[predmatrix < 0] <- 0
    newD$pred <- as.numeric(newD$pred[, 1])
    newD$lwr <- apply(predmatrix, 1, quantile, prob = 0.025)
    newD$upr <- apply(predmatrix, 1, quantile, prob = 0.975)
    newD$sp_country = as.character(i)
    newD$Species= as.character(ssci$sp[1])
    newD$Country = as.character(ssci$Country[1])
    #print(i)
    v = data.frame(v)  
    names(v) = 'estimate'
    v$predictor =c('intercept','log_sd','google_mobiliity')
    v = data.table(v)
    v$sp_country = as.character(i)
    v$Species= as.character(ssci$sp[1])
    v$Country = as.character(ssci$Country[1])
    vv[[i]] = data.table(v)
    return(newD)
    }
  u = data.table(u)
  vvv =  do.call(rbind, vv) 
  u[, FID:=exp(pred)] 
  ggplot(u, aes(x = parks_percent_change_from_baseline, y = FID, col = Country, group = sp_country)) +
    geom_line() +
    labs(x = 'Google mobility index\n[% change from baseline]', y = 'Escape distance [m]', title = 'For species with N ≥ 15') +
    scale_x_continuous(breaks = seq(-80, 140, by =20)) + 
    theme_bw() +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        # legend.spacing.y = unit(0.1, 'cm'),
        legend.key.height = unit(0.5, "line"),
        legend.margin = margin(0, 0, 0, 1),
        # legend.position=c(0.5,1.6),
        plot.title = element_text(color ="grey40", size = 8),
        #plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r = 0.5, unit = "pt"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = ax_lines, size = 0.25),
        #axis.line.y = element_blank(),
        #axis.ticks.y = element_blank(),
        axis.ticks = element_line(colour = ax_lines, size = 0.25),
        # axis.text.x = element_text()
        axis.ticks.length = unit(1, "pt"),
        # axis.text = element_text(, size = 6),
        axis.text = element_text(colour = "black", size = 7),
        axis.title = element_text(size = 8)
      )

#' ### SPECIES- specific' regressions for - and +  google index
#+ est_5, fig.width=5, fig.height = 3  
 # predictions 
  ss[parks_percent_change_from_baseline<0, google := 'before_zero']
  ss[parks_percent_change_from_baseline>0, google := 'after_zero']
  ss[, sp_country_google:= paste(sp_country, google)]
  ssg <- ss[, sp_country_google_N := .N, by = sp_country_google]
  ssg = ssg[sp_country_google_N >= 15]
  ssg[, log_FID := log(FID)]
  ssg[, log_SD := log(SD)]
  nsim = 1000
  ba = list()
  uu = foreach(i = unique(ssg$sp_country_google), .combine = rbind) %do% {
    # i = "Luscinia_megarhynchos\nHungary" "Gallinula_tenebrosa\nAustralia" "Columba_livia\nPoland" #ssc$sp_country[1]
    ssgi = ssg[sp_country_google == i]
    if(ssgi$Country[1] =='Poland' | ssgi$sp_country[1]%in%c( "Gallinula_tenebrosa\nAustralia","Luscinia_megarhynchos\nHungary")){
      mi <- lmer(log_FID ~
        log_SD +
        parks_percent_change_from_baseline +
        (1 | weekday), data = ssgi, REML = FALSE)
    }else{
      mi = lmer(log_FID ~ 
            log_SD+ 
            parks_percent_change_from_baseline +
            (1 | weekday) + (1 | IDLocality), data = ssgi, REML = FALSE)
    }        
    bsim = sim(mi, n.sim = nsim)
    v = apply(bsim@fixef, 2, quantile, prob = c(0.5)) # coefficients
    newD = data.frame(
      log_SD = mean(ssgi$log_SD),
      parks_percent_change_from_baseline = seq(min(ssgi$parks_percent_change_from_baseline), max(ssgi$parks_percent_change_from_baseline), length.out = 100)
    ) # values to predict for
    X <- model.matrix(~log_SD +parks_percent_change_from_baseline, data = newD) # exactly the model which was used has to be specified here

    # calculate predicted values and creditability intervals
    newD$pred <- (X %*% v)
    predmatrix <- matrix(nrow = nrow(newD), ncol = nsim)
    for (j in 1:nsim) {
      predmatrix[, j] <- (X %*% bsim@fixef[j, ])
    }
    #predmatrix[predmatrix < 0] <- 0
    newD$pred <- as.numeric(newD$pred[, 1])
    newD$lwr <- apply(predmatrix, 1, quantile, prob = 0.025)
    newD$upr <- apply(predmatrix, 1, quantile, prob = 0.975)
    newD$sp_country_google = as.character(i)
    newD$sp_country = as.character(ssgi$sp_country[1])
    newD$Species = as.character(ssgi$sp[1])
    newD$Country = as.character(ssgi$Country[1])
    newD$google = as.character(ssgi$google[1])
    #print(i)
    v = data.frame(v)  
    names(v) = 'estimate'
    v$predictor =c('intercept','log_sd','google_mobiliity')
    v = data.table(v)
    v$sp_country_google = as.character(i)
    v$sp_country = as.character(i)
    v$Species= as.character(ssci$sp[1])
    v$Country = as.character(ssci$Country[1])
    v$google = as.character(ssgi$google[1])
    ba[[i]] = data.table(v)
    #l[[i]] = newD
    return(data.table(newD))
    }
  bav =  do.call(rbind, ba) 
  uu[, FID:=exp(pred)]

 # plotting 
  ggplot(uu, aes(x = parks_percent_change_from_baseline, y = FID, col = Country, group = sp_country_google)) +
    geom_line() +
    labs(x = 'Google mobility index\n[% change from baseline]', y = 'Escape distance [m]', title = 'For species with N ≥ 15 per period') +
    scale_x_continuous(breaks = seq(-80, 140, by =20)) + 
    theme_bw() +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        # legend.spacing.y = unit(0.1, 'cm'),
        legend.key.height = unit(0.5, "line"),
        legend.margin = margin(0, 0, 0, 1),
        # legend.position=c(0.5,1.6),
        plot.title = element_text(color ="grey40", size = 8),
        #plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r = 0.5, unit = "pt"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = ax_lines, size = 0.25),
        #axis.line.y = element_blank(),
        #axis.ticks.y = element_blank(),
        axis.ticks = element_line(colour = ax_lines, size = 0.25),
        # axis.text.x = element_text()
        axis.ticks.length = unit(1, "pt"),
        # axis.text = element_text(, size = 6),
        axis.text = element_text(colour = "black", size = 7),
        axis.title = element_text(size = 8)
      )
  
  #' ### Global mode estimates for - & + index 
#+ pred_g, fig.width=10, fig.height = 4  
 # model before
    ssb = ss[google == 'before_zero']
    #summary(factor(ssb$Country))
    mbs=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
        (scale(parks_percent_change_from_baseline)|genus)+ (1|Species)+(1|weekday) + (1|sp_day_year) + (1|Country) + (1|IDLocality) +(1|sp_loc),
        data = ssb, REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  

    est_mbs = est_out(mbs, 'Before zero: (scale(google_mob)|genus)+(1|Species)+(1|weekday) +(1|sp_day_year) + (1|Country) + (1|IDLocality) +(1|sp_loc)')    
    est_mbs[, control_for_starting_distance:='yes']

    mbx=lmer(scale(log(FID))~ 
        scale(Year)+ 
        #scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
        (scale(parks_percent_change_from_baseline)|genus)+ (1|Species)+(1|weekday) + (1|sp_day_year) + (1|Country) + (1|IDLocality) +(1|sp_loc),
        data = ssb, REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  

    est_mbx = est_out(mbx, 'Before zero: (scale(google_mob)|genus)+(1|Species)+(1|weekday) +(1|sp_day_year) + (1|Country) + (1|IDLocality) +(1|sp_loc)')    
    est_mbx[, control_for_starting_distance:='no']
  
    mbsp=lmer(scale(log(FID))~  # singular due to weekday, else all ok
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(poly(parks_percent_change_from_baseline,2))+ 
        (scale(poly(parks_percent_change_from_baseline,2))|genus)+ (1|Species)+(1|weekday)+(1|sp_day_year) + (1|Country) + (1|IDLocality) +(1|sp_loc),
        data = ssb, REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
    est_mbsp = est_out(mbsp, 'Before zero: (poly(google_mob,2)|genus)+(1|Species)+(1|weekday) +(1|sp_day_year) + (1|Country) + (1|IDLocality) +(1|sp_loc)') 
    est_mbsp[, control_for_starting_distance:='yes']

    mbxp=lmer(scale(log(FID))~  # singular due to weekday, else all ok
        scale(Year)+ 
        #scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(poly(parks_percent_change_from_baseline,2))+ 
        (scale(poly(parks_percent_change_from_baseline,2))|genus)+ (1|Species)+(1|weekday)+(1|sp_day_year) + (1|Country) + (1|IDLocality) +(1|sp_loc),
        data = ssb, REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
    est_mbxp = est_out(mbxp, 'Before zero: (poly(google_mob,2)|genus)+(1|Species)+(1|weekday) +(1|sp_day_year) + (1|Country) + (1|IDLocality) +(1|sp_loc)') 
    est_mbxp[, control_for_starting_distance:='no']

 # model after
    ssa = ss[google == 'after_zero']
    mas=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
        (scale(parks_percent_change_from_baseline)|genus)+ (1|Species)+(1|weekday) + (1|sp_day_year) + (scale(parks_percent_change_from_baseline)|Country) + (1|IDLocality) +(1|sp_loc),
        data = ssa, REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  

    est_mas = est_out(mas, 'After zero: (scale(google_mob)|genus)+(1|Species)+(1|weekday) +(1|sp_day_year) + (scale(google_mob)|Country) + (1|IDLocality) +(1|sp_loc)')    
    est_mas[, control_for_starting_distance:='yes']

    max=lmer(scale(log(FID))~ 
        scale(Year)+ 
        #scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
        (scale(parks_percent_change_from_baseline)|genus)+ (1|Species)+(1|weekday) + (1|sp_day_year) + (scale(parks_percent_change_from_baseline)|Country) + (1|IDLocality) +(1|sp_loc),
        data = ssa, REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  

    est_max = est_out(max, 'After zero: (scale(google_mob)|genus)+(1|Species)+(1|weekday) +(1|sp_day_year) + (scale(google_mob)|Country) + (1|IDLocality) +(1|sp_loc)')    
    est_max[, control_for_starting_distance:='no']
  
    masp=lmer(scale(log(FID))~  # singular due to weekday, else all ok
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(poly(parks_percent_change_from_baseline,2))+ 
        (scale(poly(parks_percent_change_from_baseline,2))|genus)+ (1|Species)+(1|weekday)+(1|sp_day_year) + (scale(poly(parks_percent_change_from_baseline,2))|Country) + (1|IDLocality) +(1|sp_loc),
        data = ssa, REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
    est_masp = est_out(masp, 'After zero: (poly(google_mob,2)|genus)+(1|Species)+(1|weekday) +(1|sp_day_year) + (poly(google_mob,2)|Country) + (1|IDLocality) +(1|sp_loc)') 
    est_masp[, control_for_starting_distance:='yes']

    maxp=lmer(scale(log(FID))~  # singular due to weekday, else all ok
        scale(Year)+ 
       # scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(poly(parks_percent_change_from_baseline,2))+ 
        (1|genus)+ (1|Species)+(1|weekday)+(1|sp_day_year) + (1|Country) + (1|IDLocality) +(1|sp_loc),
        data = ssa, REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
    est_maxp = est_out(maxp, 'After zero: (poly(google_mob,2)|genus)+(1|Species)+(1|weekday) +(1|sp_day_year) + (poly(google_mob,2)|Country) + (1|IDLocality) +(1|sp_loc)') 
    est_maxp[, control_for_starting_distance:='no']
 # combine
    oab = rbind(est_mbs, est_mbx, est_mbsp, est_mbxp,
                   est_mas, est_max, est_masp, est_maxp)
    save(oab, file = here::here('Data/dat_est_bef-aft_rev.Rdata'))
 # plot
  load(here::here('Data/dat_est_bef-aft_rev.Rdata'))
  oab[predictor%in%c('scale(poly(parks_percent_change_from_baseline, 2))1','scale(parks_percent_change_from_baseline)'), predictor:='google mobility\n(linear)']
  oab[predictor%in%'scale(poly(parks_percent_change_from_baseline, 2))2', predictor:='google mobility\n(quadratic)']
  ooab = oab[predictor %in% c('google mobility\n(linear)', 'google mobility\n(quadratic)')]
  ooab[, predictor := factor(predictor, levels = rev(c('google mobility\n(linear)', 'google mobility\n(quadratic)')))]
  ggplot(ooab, aes(x = estimate, y = predictor, color = model, shape = control_for_starting_distance)) +
    geom_vline(xintercept = 0, color = "grey", linetype = "dotted") +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0, position = ggstance::position_dodgev(width_)) +
    #geom_point(position = ggstance::position_dodgev(.6)) +
    geom_point(position = position_dodge(width = width_), bg = 'white', size = 1.1) +
    #scale_color_viridis(discrete=TRUE, begin=0, end = 0.5)  +
    #scale_fill_viridis(discrete=TRUE, begin=0, end = 0.5) + 
    #geom_text( aes(x = n_pos,label = N), vjust = 0, size = 1.75, position = ggstance::position_dodgev(width_))+ # 3 positions for 3 bars
    # annotate("text", x=log10(3), y=85, label= "Used", col = "grey30", size = 2.5)+
  
    scale_shape_manual(name = 'Controlled for starting distance',  guide = guide_legend(reverse = TRUE), values = c(23,21))+
    scale_color_manual(name = "Model", guide = guide_legend(reverse = TRUE), values = col_)+
    scale_x_continuous(breaks = round(seq(-0.6, 0.9, by = 0.1),1))+
    ylab("") +
    xlab("Standardized effect size\n[on flight initiation distance]") +
    #coord_cartesian(xlim = c(-.15, .15)) +
   #scale_x_continuous(breaks = round(seq(-.15, .15, by = 0.05),2)) +
    theme_bw() +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        # legend.spacing.y = unit(0.1, 'cm'),
        legend.key.height = unit(0.5, "line"),
        legend.margin=margin(0,0,0,0),
        #legend.position=c(0.5,1.6),
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
  
#' ### Figure 3 - google alternative
#+ fig3_g, fig.width=12, fig.height = 12  
  ss[, NspC := .N, by ='sp_country']
  ssc = ss[NspC>9]
  ssc[, sp2 := gsub(" ", "\n", sp)]
  ssc[, Google_mobility:=parks_percent_change_from_baseline]
 
 # two rows labels
    g2 = 
    ggplot(ssc, aes(x = Google_mobility, y = FID,)) +
      #stat_smooth(method = 'rlm', se = FALSE, col = 'black', lwd = 0.5)+
      geom_point(pch = 21, alpha = 0.7, aes(fill = Country), col = 'white') + 
      stat_smooth(se = FALSE, aes(col = Country), lwd = 0.5)+ # show_guide=TRUE
      facet_wrap(~sp2, ncol = 6) +
      scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = FALSE))  +
      scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = FALSE))  +
      scale_x_continuous("Google mobility index", expand = c(0, 0)) +
      scale_y_continuous("Flight initiation distance [m]", expand = c(0, 0), trans = 'log10') +

      #annotate("text", x = 1, y = 1, label = c(rep("", 52),"Observation"), hjust = -0.08, size = 1) +
      #labs(title = "Species means per sampling location")+
      #scale_colour_manual(values=c('grey60'))+
      #scale_color_manual(name = 'try', values = c('LOESS smoothed = "grey60"'))+
      theme_MB  +
      theme(
          plot.title = element_text(size=7),
          strip.background = element_blank(),
          #strip.text.x = element_text(size = 4.5, color="grey30",  margin=margin(1,1,1,1,"mm")),
          #panel.spacing = unit(1, "mm"),
          legend.position = c(1, 0.01),
          legend.justification = c(1, 0),
          legend.title = element_blank(),
          #legend.spacing.y = unit(-0.78, "cm")
          #legend.spacing.y = unit(0.02, "cm") use if LOESS smooth text as legend
          legend.spacing.y = unit(-0.9, "cm")
          ) 

    gg2 <- ggplotGrob(g2) #gg$layout$name
    ggx2 <- gtable_filter_remove(gg2, name = c('axis-b-2-9','axis-b-5-8'), trim = FALSE)# paste0("axis-b-", c(2, 4), "-9")
    grid.draw(ggx2)
