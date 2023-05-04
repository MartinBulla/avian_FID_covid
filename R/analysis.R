# Load tools & data
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
    require(performance)  
    require(PerformanceAnalytics)
    require(png)
    require(rphylopic)
    require(viridis)
  # constants
    round_ = 3 # number of decimal places to round model coefficients
    nsim = 5000 # number of simulations to extract estimates and 95%CrI
    ax_lines = "grey60" # defines color of the axis lines
    #colors <- c("#999999", "#E69F00", "#56B4E9") #viridis(3)
    set.seed(42)
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
    o  =  fread('Data/phylopic.txt')
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



# Figure S1
   d[, sin_rad:=sin(rad)]
   d[, cos_rad:=cos(rad)]
   d[, Period:=Covid]

   dp = d[,c('Period', 'StringencyIndex','SD_ln', 'flock_ln', 'body_ln', 'sin_rad', 'cos_rad','Temp', 'Day')]
   setnames(dp,old = c('StringencyIndex','SD_ln', 'flock_ln', 'body_ln', 'sin_rad', 'cos_rad','Temp', 'Day'), new = c('Stringency\nindex','Starting distance\nln(m)', 'Flock size\nln (m)', 'Body mass\nln(m)', 'Sine\n of radians', 'Cosine\nof radians','Temperature\n°C', 'Day'))
   
   png("Outputs/Fig_S1.png", width =19, height = 19, units = "cm", bg = "transparent", res = 600)
   chart.Correlation(dp, histogram=TRUE, pch=19, alpha = 0.5)
   mtext("Single observations", side=3, line=3)
   dev.off()
# Figure Sx - correlation betwen stringency and google mobility
   g=
   ggplot(s, aes(x = StringencyIndex, y = parks_percent_change_from_baseline, col = Country)) + 
      stat_smooth(method = 'lm', se = FALSE)+
      stat_cor(method="pearson", size =2)+
      geom_point() 
   ggsave(here::here('Outputs/Fig_Sx.png'),g, width = 14, height =12, units = 'cm')

  ggplot(s, aes(x = parks_percent_change_from_baseline)) + 
      geom_histogram() 

  g2=
   ggplot(s, aes(x = parks_percent_change_from_baseline, col = Country)) + 
      geom_histogram() + 
      facet_wrap(~Country, nrow = 5) +
      theme(legend.position = 'none')
    ggsave(here::here('Outputs/Fig_Sx_hist.png'),g2, width = 10, height =14, units = 'cm')   
# Figure Sy - fid and google mobility
    s[, country_year := paste(Country, Year)] #table(paste(s$Country, s$Year))   
    s[, year_ := as.character(Year)] #table(paste(s$Country, s$Year))   
    ss = s[!is.na(parks_percent_change_from_baseline), Nsp := .N, by ='sp']  
    ss[, weekday:=factor(weekday, levels = c('Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday'))]

    gy = 
    ggplot(s, aes(y = FID, x = parks_percent_change_from_baseline, col = Country)) + 
      stat_smooth()+
      #stat_cor(method="pearson", size =2) +
      scale_y_continuous(trans = 'log')
    ggsave(here::here('Outputs/Fig_Sy_ln.png'),gy, width = 10, height =14, units = 'cm')   
  
    gy2 = 
    ggplot(s, aes(y = FID, x = parks_percent_change_from_baseline, col = Country)) + 
      stat_smooth()
      ggsave(here::here('Outputs/Fig_Sy.png'),gy2, width = 10, height =14, units = 'cm')   

   
    ggplot(s, aes(x = parks_percent_change_from_baseline, col = year_)) + 
      geom_histogram() + 
      facet_wrap(~Country, nrow = 5) 



  ggplot(s, aes(y = FID, x = parks_percent_change_from_baseline)) + 
        stat_smooth() +
      # geom_point(size =0.5, pch = 1) + 
        facet_wrap(~sp, scales ='free_y')

  ggplot(ss, aes(x = Day, y = parks_percent_change_from_baseline, col = Country, lty = year_, pch = year_)) + 
        geom_point()+
        stat_smooth(method='rlm', lwd =0.5) +
        facet_wrap(~weekday, nrow =2)
  ggsave(here::here('Outputs/Fig_Sz_weekdays.png'),width = 20, height =12, units = 'cm')   

  ggplot(ss, aes(y = FID, x = parks_percent_change_from_baseline, col = Country, lty = year_, pch = year_)) + 
        stat_smooth(method="gam",formula = y ~ s(x, bs = "cs", k=3))+
        #stat_smooth(method='rlm', lwd =0.5) +
        facet_wrap(~weekday, nrow =2) + 
        scale_y_continuous(trans = 'log')

  ggsave(here::here('Outputs/Fig_Sz_weekdays_FID.png'),width = 20, height =12, units = 'cm')   

  ggplot(ss, aes(y = FID, x = parks_percent_change_from_baseline, lty = year_)) + 
        stat_smooth(method="gam",formula = y ~ s(x, bs = "cs", k=5))+
        #stat_smooth(method='rlm', lwd =0.5) +
        facet_wrap(~weekday, nrow =7, strip.position="right") + 
        scale_y_continuous(trans = 'log')

  ggsave(here::here('Outputs/Fig_Sz_weekdays_FID_global.png'),width = 10, height =20, units = 'cm')   

  ggplot(s, aes(x = Day, y = parks_percent_change_from_baseline, col = Country, lty = year_)) + 
        stat_smooth(method="gam",formula = y ~ s(x, bs = "cs", k=6))

  ggplot(s, aes(x = Day, y = parks_percent_change_from_baseline, col = Country, lty = year_)) + 
        stat_smooth(method = 'lm') 
  ggsave(here::here('Outputs/Fig_Sz.png'),width = 12, height =10, units = 'cm')   

  ggplot(s, aes(x = Day, y = FID, col = Country, lty = year_)) + 
        stat_smooth(method = 'lm')  +
        scale_y_continuous(trans = 'log')
  ggsave(here::here('Outputs/Fig_Sz_FID.png'),width = 12, height =10, units = 'cm')   

  gy2 = 
    ggplot(ss[Nsp>5], aes(y = FID, x = parks_percent_change_from_baseline, groups = sp, col = Country)) + 
        stat_smooth(method = 'lm', se = FALSE) +
        labs(subtitle  = "regresions for specis with >5 data points")
        #theme(legend.position = 'none')
      # geom_point(size =0.5, pch = 1) + 
        #facet_wrap(~sp, scales ='free_y')
        ggsave(here::here('Outputs/Fig_Sy_sp_country.png'),gy2, width = 10, height =14, units = 'cm')   

    ggplot(s, aes(x = parks_percent_change_from_baseline, col = year_)) + 
        geom_histogram(position = "dodge") + 
        facet_wrap(~IDLocality) 

      gy3 = 
      ggplot(s, aes(y = FID, x = parks_percent_change_from_baseline, col = country_year)) + 
        stat_smooth()
        ggsave(here::here('Outputs/Fig_Sy.png'),gy2, width = 10, height =14, units = 'cm')   

      
# Figure Sz - hour distritributiion of observations
   ggplot(d, aes(x = Hour)) + geom_histogram()
      
# Figure S4 - year trend for >4 observations/site before and during covid
   px = pp[N_during>4 & N_before>4]
   dxx = d[paste(IDLocality, Species) %in% paste(px$IDLocality, px$Species)]
    #table(dxx$IDLocality, dxx$Year)
    length(unique(px$IDLocality))  
    length(unique(px$Species)) 
    sum(px$N_during) + sum(px$N_before)

   dxx[, sp_C_loc2 := paste(gsub('[_]', ' ', Species), Country, IDLocality, sep ='\n')]
   dxx[, genus := sub("_.*", "", Species)]
   dxx[Covid == 0 , period := 'before COVID-19']
   dxx[Covid == 1 , period := 'during COVID-19']
   g = 
   ggplot(dxx, aes(x = as.factor(Year), y = FID, col = Year, fill = period)) + 
    geom_boxplot() + facet_wrap(~sp_C_loc2) + 
    scale_y_continuous("Flight initiation distance [m]", trans = 'log10') + 
    scale_x_discrete("Year", guide = guide_axis(angle = 45)) +
    scale_color_continuous()+
    scale_fill_manual(values = c('white', 'lightgrey'))+
    theme_MB  +
    theme(
          plot.title = element_text(size=7),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 5, color="grey30",  margin=margin(1,1,1,1,"mm")),
          #panel.spacing = unit(1, "mm"),
          legend.position = "none",#c(1, 0.01),
          legend.justification = c(1, 0),
          legend.title = element_blank(),
          #legend.spacing.y = unit(-0.78, "cm")
          #legend.spacing.y = unit(0.02, "cm") use if LOESS smooth text as legend
          legend.spacing.y = unit(-0.9, "cm"),
          axis.text.x = element_text(colour="grey30", size = 6),
          axis.text.y=element_text(colour="grey30", size = 6)
          ) 

   ggsave(here::here('Outputs/Fig_S4_rev.png'),g, width = 18, height =16, units = 'cm')

# not used - Figure Sub - year trend for >0 observations/site before and during covid
   px = pp[N_during>0 & N_before>0]
   dxx = d[paste(IDLocality, Species) %in% paste(px$IDLocality, px$Species)]
    table(dxx$IDLocality, dxx$Year)

   dxx[, sp_C_loc2 := paste(gsub('[_]', ' ', Species), Country, IDLocality, sep ='\n')]
   dxx[, genus := sub("_.*", "", Species)]
   g = 
   ggplot(dxx, aes(x = as.factor(Year), y = FID, col = Year)) + 
    geom_boxplot() + facet_wrap(~sp_C_loc2, ncol = 15) + 
    scale_y_continuous("Flight initiation distance [m]", trans = 'log10') + 
    scale_x_discrete("Year", guide = guide_axis(angle = 45)) +
    scale_color_continuous()+
    theme_MB  +
    theme(
          plot.title = element_text(size=7),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 5, color="grey30",  margin=margin(1,1,1,1,"mm")),
          #panel.spacing = unit(1, "mm"),
          legend.position = "none",#c(1, 0.01),
          legend.justification = c(1, 0),
          legend.title = element_blank(),
          #legend.spacing.y = unit(-0.78, "cm")
          #legend.spacing.y = unit(0.02, "cm") use if LOESS smooth text as legend
          legend.spacing.y = unit(-0.9, "cm"),
          axis.text.x = element_text(colour="grey30", size = 6),
          axis.text.y=element_text(colour="grey30", size = 6)
          ) 

   #ggsave(here::here('Outputs/Fig_Su_atLeast1.png'),g, width = 24, height =30, units = 'cm')

   #ggplot(d, aes(x = Day, y = FID)) +  geom_smooth(aes(col = as.factor(Year))) + scale_y_continuous(trans = 'log10') +  scale_color_viridis(discrete=TRUE)
   #ggplot(d, aes(x = Day, y = FID)) +  geom_smooth(aes(col = as.factor(Year))) + scale_y_continuous(trans = 'log10') +  scale_color_viridis(discrete=TRUE) + facet_wrap(~Country, ncol = 5)

# Figure 1, S2, Tables S1 & S2
  
  # prepare estimates Period
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
  # prepare estimates Stringency 
     m01a=lmer(scale(log(FID))~ 
        scale(Year)+ 
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
        # (1|Year) explains nothing - could stay 
     est_m01a = est_out(m01a, '01a) (scale(StringencyIndex)|genus)+(1|Species)+(1|sp_day_year) + (1|Country) + (scale(StringencyIndex)|IDLocality) +(1|sp_loc)')
     m01b=lmer(scale(log(FID))~
           scale(Year)+       
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
          scale(Year)+       
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
         scale(Year)+ 
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
          scale(Year)+ 
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
          scale(Year)+ 
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
        scale(Year)+ 
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
          scale(Year)+ 
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
          scale(Year)+ 
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
        (scale(parks_percent_change_from_baseline)|genus)+(1|Species)+(1|sp_day_year) + (1|Country) + (scale(parks_percent_change_from_baseline)|IDLocality) +(1|sp_loc),
        data = s, REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
        # (1|Year) explains nothing - could stay 
     est_mg01a = est_out(mg01a, 'g01a) (scale(parks_percent_change_from_baseline)|genus)+(1|Species)+(1|sp_day_year) + (1|Country) + (scale(parks_percent_change_from_baseline)|IDLocality) +(1|sp_loc)')
     mg01b=lmer(scale(log(FID))~
           scale(Year)+       
          scale(log(SD))+
          scale(log(FlockSize))+
          scale(log(BodyMass))+
          scale(sin(rad)) + scale(cos(rad)) + 
          #scale(Day)+
          scale(Temp)+
          scale(parks_percent_change_from_baseline)+
          (1|genus) +(1|Species) + (1|sp_day_year) + (1|Country) + (scale(parks_percent_change_from_baseline)|IDLocality)+(1|sp_loc),  
          data = s, REML = FALSE, 
          control = lmerControl(
              optimizer ='optimx', optCtrl=list(method='nlminb'))
          ) 
          # (1|Year), (1|genus), (1|sp_day_year), (1|sp_loc),  explain nothing - could stay
     est_mg01b = est_out(mg01b, 'g01b) (1|genus) +(1|Species) + (1|sp_day_year) + (1|Country) + (scale(parks_percent_change_from_baseline)|IDLocality)+(1|sp_loc)')  
     mg01c=lmer(scale(log(FID))~
          scale(Year)+       
          scale(log(SD))+
          scale(log(FlockSize))+
          scale(log(BodyMass))+
          scale(sin(rad)) + scale(cos(rad)) + 
          #scale(Day)+
          scale(Temp)+
          scale(parks_percent_change_from_baseline)+
          (1|Country) + (scale(parks_percent_change_from_baseline)|IDLocality),  
          data = s, REML = FALSE, 
          control = lmerControl(
              optimizer ='optimx', optCtrl=list(method='nlminb'))
          ) 
          # (1|Year), (1|genus), (1|sp_day_year), (1|sp_loc),  explain nothing - could stay
     est_mg01c = est_out(mg01c, 'g01c) (1|Country) + (scale(parks_percent_change_from_baseline)|IDLocality)')  

     mg02a=lmer(scale(log(FID))~ 
         scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
        (scale(parks_percent_change_from_baseline)|genus)+(1|Species)+(1|sp_day_year) + 
        (1|Country) + (scale(parks_percent_change_from_baseline)|IDLocality) +(1|sp_loc),   
        data = s[Nsp>4], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
        # (1|Year) explains nothing - could stay 
     est_mg02a = est_out(mg02a, 'g02a) (scale(parks_percent_change_from_baseline)|genus)+(1|Species)+(1|sp_day_year) + (1|Country) + (scale(parks_percent_change_from_baseline)|IDLocality) +(1|sp_loc);>4/species')   
     mg02b=lmer(scale(log(FID))~
          scale(Year)+ 
          scale(log(SD))+
          scale(log(FlockSize))+
          scale(log(BodyMass))+
          scale(sin(rad)) + scale(cos(rad)) + 
          #scale(Day)+
          scale(Temp)+
          scale(parks_percent_change_from_baseline)+
          (1|genus) +(1|Species) + (1|sp_day_year) + (1|Country) + (scale(parks_percent_change_from_baseline)|IDLocality)+(1|sp_loc),  
          data = s[Nsp>4], REML = FALSE, 
          control = lmerControl(
              optimizer ='optimx', optCtrl=list(method='nlminb'))
          ) 
     est_mg02b = est_out(mg02b, 'g02b)  (1|genus) +(1|Species) + (1|sp_day_year) + (1|Country) + (scale(parks_percent_change_from_baseline)|IDLocality)+(1|sp_loc); >4/species') 
     mg02c=lmer(scale(log(FID))~
          scale(Year)+ 
          scale(log(SD))+
          scale(log(FlockSize))+
          scale(log(BodyMass))+
          scale(sin(rad)) + scale(cos(rad)) + 
          #scale(Day)+
          scale(Temp)+
          scale(parks_percent_change_from_baseline)+
          (1|Country) + (scale(parks_percent_change_from_baseline)|IDLocality),  
          data = s[Nsp>4], REML = FALSE, 
          control = lmerControl(
              optimizer ='optimx', optCtrl=list(method='nlminb'))
          ) 
     est_mg02c = est_out(mg02c, 'g02c)  (1|Country) + (scale(parks_percent_change_from_baseline)|IDLocality); >4/species') 
     
     mg03a=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
        (scale(parks_percent_change_from_baseline)|genus)+(1|Species)+(1|sp_day_year) + 
        (1|Country) + (scale(parks_percent_change_from_baseline)|IDLocality) +(1|sp_loc),   
        data = s[Nsp>9], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
        # (1|Year) explains nothing - could stay 
     est_mg03a = est_out(mg03a, 'g03a) (scale(parks_percent_change_from_baseline)|genus)+(1|Species)+(1|sp_day_year) + (1|Country) + (scale(parks_percent_change_from_baseline)|IDLocality) +(1|sp_loc);>9/species') 
     mg03b=lmer(scale(log(FID))~
          scale(Year)+ 
          scale(log(SD))+
          scale(log(FlockSize))+
          scale(log(BodyMass))+
          scale(sin(rad)) + scale(cos(rad)) + 
          #scale(Day)+
          scale(Temp)+
          scale(parks_percent_change_from_baseline)+
          (1|genus) +(1|Species) + (1|sp_day_year) + (1|Country) + (scale(parks_percent_change_from_baseline)|IDLocality)+(1|sp_loc),  
          data = s[Nsp>9], REML = FALSE, 
          control = lmerControl(
              optimizer ='optimx', optCtrl=list(method='nlminb'))
          ) 
     est_mg03b = est_out(mg03b, 'g03b)  (1|genus) +(1|Species) + (1|sp_day_year) + (1|Country) + (scale(parks_percent_change_from_baseline)|IDLocality)+(1|sp_loc); >9/species') 
     mg03c=lmer(scale(log(FID))~
          scale(Year)+ 
          scale(log(SD))+
          scale(log(FlockSize))+
          scale(log(BodyMass))+
          scale(sin(rad)) + scale(cos(rad)) + 
          #scale(Day)+
          scale(Temp)+
          scale(parks_percent_change_from_baseline)+
          (1|Country) + (scale(parks_percent_change_from_baseline)|IDLocality),  
          data = s[Nsp>9], REML = FALSE, 
          control = lmerControl(
              optimizer ='optimx', optCtrl=list(method='nlminb'))
          ) 
     est_mg03c = est_out(mg03c, 'g03c)  (1|Country) + (scale(parks_percent_change_from_baseline)|IDLocality); >9/species') 
 
  # export dataset with model residuals for test of phylo-signal (key models are Table S1 1d and S2 1c)
     d[,res := resid(m1d)]
     d_ = d[,.(Species,res)]
     s[, res := resid(m01c)]
     s_ = s[,.(Species,res)]
     save(file = 'Data/DAT_res.Rdata', d_, s_)

     d[,res_m1a := resid(m1a)]
     d[,res_m1b := resid(m1b)]
     d[,res_m1c := resid(m1c)]
     d[,res_m1d := resid(m1d)]
     
     dx = dd[N_during>4 & N_before >4]
     d5 = d[Species %in% dx$Species]
     d5[,res_m2a := resid(m2a)]
     d5[,res_m2b := resid(m2b)]

     dx = dd[N_during>9 & N_before >9]
     d10 = d[Species %in% dx$Species]
     d10[,res_m3a := resid(m3a)]
     d10[,res_m3b := resid(m3b)]
     

     s[, res_m1a := resid(m01a)]
     s[, res_m1b := resid(m01b)]
     s[, res_m1c := resid(m01c)]

     s5= s[Nsp>4]
     s5[, res_m2a := resid(m02a)]
     s5[, res_m2b := resid(m02b)]
     s5[, res_m2c := resid(m02c)]

     s10= s[Nsp>9]
     s10[, res_m3a := resid(m03a)]
     s10[, res_m3b := resid(m03b)]
     s10[, res_m3c := resid(m03c)]

     d_ = d[,.(Species,res_m1a, res_m1b, res_m1c, res_m1d)]
     d5 = d5[.,(Species, res_m2a, res_m2b)]
     d10 = d10[.,(Species, res_m3a, res_m3c)]

     s_ = s[,.(Species,res_m1a, res_m1b, res_m1c)]
     s5 = s5[,.(Species,res_m2a, res_m2b,res_m2c)]
     s10 = s10[,.(Species,res_m3a, res_m3b, res_m3c)]

     #save(file = 'Data/DAT_res.Rdata', d_, s_, d5, d10, s5, s10)
  
  # Figure 1
    xc = rbind(est_m1d, est_m2b,est_m3b,est_m01c, est_m02c,est_m03c)
    xc = xc[predictor %in% c('scale(Covid)', 'scale(StringencyIndex)')]
    xc[predictor %in% 'scale(Covid)', predictor := 'a) Period - before vs during COVID-19 shutdown)']
    xc[predictor %in% 'scale(StringencyIndex)', predictor := 'b) Stringency of governmental COVID-19\n    restrictions']
    xc[, N:=c('N = 6369; all data', 'N = 5260; ≥5 observations/species/period', 'N = 5106; ≥10 observations/species/period',
              'N = 3676; all data', 'N = 3573; ≥5 observations/species', 'N = 3425; ≥10 observations/species')]
    xc[, N := factor(N, levels = c('N = 6369; all data', 'N = 5260; ≥5 observations/species/period', 'N = 5106; ≥10 observations/species/period',
              'N = 3676; all data', 'N = 3573; ≥5 observations/species', 'N = 3425; ≥10 observations/species'))]
    
    col_ = c(rep(viridis(1, alpha = 1, begin = 0.3, end = 0.4, direction = 1, option = "D"),3),
              rep(viridis(1, alpha = 1, begin = 0.8, end = 0.8, direction = 1, option = "D"),3))
    #col_ = c(viridis(3, alpha = 1, begin = 0.25, end = 0.4, direction = 1, option = "D"),viridis(3, alpha = 1, begin = 0.75, end = 0.9, direction = 1, option = "D"))
      #show_col(viridis(6, alpha = 1, begin = 0, end = 1, direction = 1, option = "D"))
    g = 
    ggplot(xc, aes(y = N, x = estimate, col = N)) +
          geom_vline(xintercept = 0, col = "grey30", lty =3)+
          geom_errorbar(aes(xmin = lwr, xmax = upr, col = N), width = 0, position = position_dodge(width = 0.01) ) +
          #ggtitle ("Sim based")+
          geom_point(position = position_dodge(width = 0.01)) +
          #scale_colour_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          #scale_fill_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          #scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
          scale_color_manual(values = col_,guide = guide_legend(reverse = TRUE))  +
          #scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE)) + 
          scale_fill_manual(values = col_,guide = guide_legend(reverse = TRUE)) + 
          scale_y_discrete(limits=rev, position = "right")+
          #coord_fixed(ratio = 0.05)+ #,xlim = c(-0.23, 0.15)
          #scale_shape(guide = guide_legend(reverse = TRUE)) + 
          #scale_x_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2,2, by = 1), labels = seq(-2,2, by = 1)) +
          labs(y = NULL ,x = "Standardized effect sizes ") + #title = "a)",Effect of Period (before/during shutdown)
          #ylim(c(0,100))+
          #coord_flip()+
          facet_wrap(~predictor, nrow = 2, scales = 'free_y')+
          theme_MB +
          theme( legend.position ="none",
                plot.title = element_text(size=7),
                plot.tag = element_text(size=7),
                legend.title=element_text(size=7), 
                legend.text=element_text(size=6),
                ##legend.spacing.y = unit(0.1, 'cm'), 
                legend.key.height= unit(0.5,"line"),
                #plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r =0.5, unit =  "pt"),
                panel.grid = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                strip.background = element_blank(),
                strip.text = element_text(hjust = 0),
                axis.line = element_line(colour = ax_lines, size = 0.25),
                axis.line.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x= element_line( colour = ax_lines, size = 0.25),
                axis.ticks.length = unit(1, "pt"),
                axis.text.x = element_text(colour="grey30", size = 6),
                axis.text.y=element_text(colour="grey30", size = 6),
                axis.title=element_text(size=7)
                )
    g
    ggsave(here::here('Outputs/Fig_1_width-92mm.png'),g, width = 9.2, height =6, units = 'cm')   

  # Figure 1 - revision
    xc = rbind(est_m1d, est_m2b,est_m3b,est_m01c, est_m02c,est_m03c,est_mg01c, est_mg02c,est_mg03c)
    xc = xc[predictor %in% c('scale(Covid)', 'scale(StringencyIndex)', 'scale(parks_percent_change_from_baseline)')]
    xc[predictor %in% 'scale(Covid)', predictor := 'a) Period - before vs during COVID-19 shutdown)']
    xc[predictor %in% 'scale(StringencyIndex)', predictor := 'b) Stringency of governmental COVID-19\n    restrictions']
    xc[predictor %in% 'scale(parks_percent_change_from_baseline)', predictor := 'c) Google mobilitty for parks given baseline']
    xc[, N:=c('N = 6369; all data', 'N = 5260; ≥5 observations/species/period', 'N = 5106; ≥10 observations/species/period',
              'N = 3676; all data', 'N = 3573; ≥5 observations/species', 'N = 3425; ≥10 observations/species',
              'N = 3644; all data', 'N = 3545; ≥5 observations/species', 'N = 3399; ≥10 observations/species'
              )]
    xc[, N := factor(N, levels = c(
              'N = 6369; all data', 'N = 5260; ≥5 observations/species/period', 'N = 5106; ≥10 observations/species/period',
              'N = 3676; all data', 'N = 3573; ≥5 observations/species', 'N = 3425; ≥10 observations/species',
               'N = 3644; all data', 'N = 3545; ≥5 observations/species', 'N = 3399; ≥10 observations/species'))]

    col_ = c(rep(viridis(1, alpha = 1, begin = 0.3, end = 0.3, direction = 1, option = "D"),3),
              rep(viridis(1, alpha = 1, begin = 0.55, end = 0.55, direction = 1, option = "D"),3),
              rep(viridis(1, alpha = 1, begin = 0.8, end = 0.8, direction = 1, option = "D"),3))
    #col_ = c(viridis(3, alpha = 1, begin = 0.25, end = 0.4, direction = 1, option = "D"),viridis(3, alpha = 1, begin = 0.75, end = 0.9, direction = 1, option = "D"))
      #show_col(viridis(6, alpha = 1, begin = 0, end = 1, direction = 1, option = "D"))
      #show_col(viridis(1, alpha = 1, begin = 0.8, end = 0.8, direction = 1, option = "D"))
    g = 
    ggplot(xc, aes(y = N, x = estimate, col = N)) +
          geom_vline(xintercept = 0, col = "grey30", lty =3)+
          geom_errorbar(aes(xmin = lwr, xmax = upr, col = N), width = 0, position = position_dodge(width = 0.01) ) +
          #ggtitle ("Sim based")+
          geom_point(position = position_dodge(width = 0.01)) +
          #scale_colour_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          #scale_fill_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          #scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
          scale_color_manual(values = col_,guide = guide_legend(reverse = TRUE))  +
          #scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE)) + 
          scale_fill_manual(values = col_,guide = guide_legend(reverse = TRUE)) + 
          scale_y_discrete(limits=rev, position = "right")+
          #coord_fixed(ratio = 0.05)+ #,xlim = c(-0.23, 0.15)
          #scale_shape(guide = guide_legend(reverse = TRUE)) + 
          #scale_x_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2,2, by = 1), labels = seq(-2,2, by = 1)) +
          labs(y = NULL ,x = "Standardized effect sizes ") + #title = "a)",Effect of Period (before/during shutdown)
          #ylim(c(0,100))+
          #coord_flip()+
          facet_wrap(~predictor, nrow = 3, scales = 'free_y')+
          theme_MB +
          theme( legend.position ="none",
                plot.title = element_text(size=7),
                plot.tag = element_text(size=7),
                legend.title=element_text(size=7), 
                legend.text=element_text(size=6),
                ##legend.spacing.y = unit(0.1, 'cm'), 
                legend.key.height= unit(0.5,"line"),
                #plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r =0.5, unit =  "pt"),
                panel.grid = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                strip.background = element_blank(),
                strip.text = element_text(hjust = 0),
                axis.line = element_line(colour = ax_lines, size = 0.25),
                axis.line.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.ticks.x= element_line( colour = ax_lines, size = 0.25),
                axis.ticks.length = unit(1, "pt"),
                axis.text.x = element_text(colour="grey30", size = 6),
                axis.text.y=element_text(colour="grey30", size = 6),
                axis.title=element_text(size=7)
                )
    g
    ggsave(here::here('Outputs/Fig_1_width-92mm_rev.png'),g, width = 9.2, height =9, units = 'cm')   

   # Figure S2  
    # prepare plot for Period
      xs = rbind(est_m1a,est_m1b,est_m1c,est_m1d, est_m2a, est_m2b,est_m3a, est_m3b)
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
          coord_fixed(ratio = 0.05,xlim = c(-0.23, 0.15))+
          #scale_shape(guide = guide_legend(reverse = TRUE)) + 
          #scale_x_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2,2, by = 1), labels = seq(-2,2, by = 1)) +
          labs(y = NULL ,x = "Period (before/during shutdown)\n ",  tag = 'a)') + #title = "a)",Effect of Period (before/during shutdown)
          #ylim(c(0,100))+
          #coord_flip()+
          theme_bw() +
          theme( legend.position ="none",
                plot.title = element_text(size=7),
                plot.tag = element_text(size=7),
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
      #ggsave(here::here('Outputs/Figure_Sy.png'),g, width = 30, height =5, units = 'cm')    
    # prepare plot for Stringency
      xs0 = rbind(est_m01a,est_m01b,est_m01c, est_m02a, est_m02b,est_m02c, est_m03a, est_m03b, est_m03c)
      g0=
      ggplot(xs0[predictor == 'scale(StringencyIndex)'], aes(y = model, x = estimate, col = model)) +
          geom_vline(xintercept = 0, col = "grey30", lty =3)+
          geom_errorbar(aes(xmin = lwr, xmax = upr, col = model), width = 0, position = position_dodge(width = 0.01) ) +
          #ggtitle ("Sim based")+
          geom_point(position = position_dodge(width = 0.01)) +
          #scale_colour_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          #scale_fill_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
          scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE)) + 
          scale_y_discrete(limits=rev)+
          coord_fixed(ratio = 0.05, xlim = c(-0.23, 0.15))+
          #scale_shape(guide = guide_legend(reverse = TRUE)) + 
          #scale_x_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2,2, by = 1), labels = seq(-2,2, by = 1)) +
          labs(y = NULL ,x = "Shutdown stringency\n[standardized effect sizes]",  tag = 'b)')+#title = "b) Effect of ") +
          #ylim(c(0,100))+
          #coord_flip()+
          theme_bw() +
          theme(legend.position ="none",
                plot.title = element_text(size=7),
                plot.tag = element_text(size=7),
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
      g0
      #ggsave(here::here('Outputs/Figure_Sz.png'),g0, width = 30, height =5, units = 'cm')
    # prepare plot for Google 
      xs0 = rbind(est_mg01a,est_mg01b,est_mg01c, est_mg02a, est_mg02b,est_mg02c, est_mg03a, est_mg03b, est_mg03c)
      gg0=
      ggplot(xs0[predictor == 'scale(parks_percent_change_from_baseline)'], aes(y = model, x = estimate, col = model)) +
          geom_vline(xintercept = 0, col = "grey30", lty =3)+
          geom_errorbar(aes(xmin = lwr, xmax = upr, col = model), width = 0, position = position_dodge(width = 0.01) ) +
          #ggtitle ("Sim based")+
          geom_point(position = position_dodge(width = 0.01)) +
          #scale_colour_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          #scale_fill_brewer(type = 'qual', palette = 'Paired',guide = guide_legend(reverse = TRUE))+
          scale_color_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
          scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE)) + 
          scale_y_discrete(limits=rev)+
          coord_fixed(ratio = 0.05, xlim = c(-0.23, 0.15))+
          #scale_shape(guide = guide_legend(reverse = TRUE)) + 
          #scale_x_continuous(limits = c(-2, 2), expand = c(0, 0), breaks = seq(-2,2, by = 1), labels = seq(-2,2, by = 1)) +
          labs(y = NULL ,x = "Shutdown's google mobility in parks\n[standardized effect sizes]",  tag = 'b)')+#title = "b) Effect of ") +
          #ylim(c(0,100))+
          #coord_flip()+
          theme_bw() +
          theme(legend.position ="none",
                plot.title = element_text(size=7),
                plot.tag = element_text(size=7),
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
      gg0
      #ggsave(here::here('Outputs/Figure_Sz.png'),g0, width = 30, height =5, units = 'cm')
    # combine
     ggsave(here::here('Outputs/Fig_S2_rev.png'),rbind(ggplotGrob(g),ggplotGrob(g0),ggplotGrob(gg0)), width = 30, height =15, units = 'cm') 

  # TABLES
    m1a_ = m_out(name = "Table S1 - 1a", dep = 'Period',model = m1a, nsim = 5000)
    m1b_ = m_out(name = "Table S1 - 1b", dep = 'Period',model = m1b, nsim = 5000)
    m1c_ = m_out(name = "Table S1 - 1c", dep = 'Period',model = m1c, nsim = 5000)
    m1d_ = m_out(name = "Table S1 - 1d", dep = 'Period',model = m1d, nsim = 5000)

    m2a_ = m_out(name = "Table S1 - 2a", dep = 'Period',model = m2a, nsim = 5000)
    m2b_ = m_out(name = "Table S1 - 2b", dep = 'Period',model = m2b, nsim = 5000)
    m3a_ = m_out(name = "Table S1 - 3a", dep = 'Period',model = m3a, nsim = 5000)
    m3b_ = m_out(name = "Table S1 - 3b", dep = 'Period',model = m3b, nsim = 5000)

    
    m01a_ = m_out(name = "Table S2 - 1a", dep = 'Stringency Index',model = m01a, nsim = 5000)
    m01b_ = m_out(name = "Table S2 - 1b", dep = 'Stringency Index',model = m01b, nsim = 5000)
    m01c_ = m_out(name = "Table S2 - 1c", dep = 'Stringency Index',model = m01c, nsim = 5000)
   
    m02a_ = m_out(name = "Table S2 - 2a", dep = 'Stringency Index',model = m02a, nsim = 5000)
    m02b_ = m_out(name = "Table S2 - 2b", dep = 'Stringency Index',model = m02b, nsim = 5000)
    m02c_ = m_out(name = "Table S2 - 2c", dep = 'Stringency Index',model = m02c, nsim = 5000)
    
    m03a_ = m_out(name = "Table S2 - 3a", dep = 'Stringency Index',model = m03a, nsim = 5000)
    m03b_ = m_out(name = "Table S2 - 3b", dep = 'Stringency Index',model = m03b, nsim = 5000)
    m03c_ = m_out(name = "Table S2 - 3c", dep = 'Stringency Index',model = m03c, nsim = 5000)

    mg01a_ = m_out(name = "Table S3 - 1a", dep = 'Stringency Index',model = mg01c, nsim = 5000)
    

    out1 = rbind(m1a_, m1b_, m1c_, m1d_, m2a_, m2b_, m3a_, m3b_,fill = TRUE)
    out1[is.na(out1)] = ""
    fwrite(file = "./Outputs/Table_S1.csv", out1)

    out2 = rbind(m01a_, m01b_, m01c_, m02a_, m02b_, m02c_, m03a_, m03b_, m03c_,fill = TRUE)
    out2[is.na(out2)] = ""
    fwrite(file = "./Outputs/Table_S2.csv", out2)

  # modelAss
    m_ass(name = "Table S1 - 1a", mo = m1a, fixed = c('SD', 'FlockSize', 'BodyMass', 'rad','rad','Temp','Covid'), trans = c("log","log","log","sin","cos","","") , categ = 'Covid', outdir = 'Outputs/modelAss/')
    m_ass(name = "Table S1 - 1b", mo = m1b, fixed = c('SD', 'FlockSize', 'BodyMass', 'rad','rad','Temp','Covid'), trans = c("log","log","log","sin","cos","","") , categ = 'Covid', outdir = 'Outputs/modelAss/')
    m_ass(name = "Table S1 - 1c", mo = m1c, fixed = c('SD', 'FlockSize', 'BodyMass', 'rad','rad','Temp','Covid'), trans = c("log","log","log","sin","cos","","") , categ = 'Covid', outdir = 'Outputs/modelAss/')
    m_ass(name = "Table S1 - 1d", mo = m1d, fixed = c('SD', 'FlockSize', 'BodyMass', 'rad','rad','Temp','Covid'), trans = c("log","log","log","sin","cos","","") , categ = 'Covid', outdir = 'Outputs/modelAss/')

    m_ass(name = "Table S1 - 2a", mo = m2a, dat = d[Species %in% dd[N_during>4 & N_before >4, Species]], fixed = c('SD', 'FlockSize', 'BodyMass', 'rad','rad','Temp','Covid'), trans = c("log","log","log","sin","cos","","") , categ = 'Covid', outdir = 'Outputs/modelAss/')
    m_ass(name = "Table S1 - 2b", mo = m2b, dat = d[Species %in% dd[N_during>4 & N_before >4, Species]],fixed = c('SD', 'FlockSize', 'BodyMass', 'rad','rad','Temp','Covid'), trans = c("log","log","log","sin","cos","","") , categ = 'Covid', outdir = 'Outputs/modelAss/')

    m_ass(name = "Table S1 - 3a", mo = m3a, dat = d[Species %in% dd[N_during>9 & N_before >9, Species]], fixed = c('SD', 'FlockSize', 'BodyMass', 'rad','rad','Temp','Covid'), trans = c("log","log","log","sin","cos","","") , categ = 'Covid', outdir = 'Outputs/modelAss/')
    m_ass(name = "Table S1 - 3b", mo = m3b, dat = d[Species %in% dd[N_during>9 & N_before >9, Species]], fixed = c('SD', 'FlockSize', 'BodyMass', 'rad','rad','Temp','Covid'), trans = c("log","log","log","sin","cos","","") , categ = 'Covid', outdir = 'Outputs/modelAss/')

    m_ass(name = "Table S2 - 1a", mo = m01a, dat = s, fixed = c('Year','SD', 'FlockSize', 'BodyMass', 'rad','rad','Temp','StringencyIndex'), trans = c("","log","log","log","sin","cos","","") , categ = 'Year', outdir = 'Outputs/modelAss/')
    m_ass(name = "Table S2 - 1b", mo = m01b, dat = s, fixed = c('Year','SD', 'FlockSize', 'BodyMass', 'rad','rad','Temp','StringencyIndex'), trans = c("","log","log","log","sin","cos","","") , categ = 'Year', outdir = 'Outputs/modelAss/')
    m_ass(name = "Table S2 - 1c", mo = m01c, dat = s, fixed = c('Year','SD', 'FlockSize', 'BodyMass', 'rad','rad','Temp','StringencyIndex'), trans = c("","log","log","log","sin","cos","","") , categ = 'Year', outdir = 'Outputs/modelAss/')
    
    m_ass(name = "Table S2 - 2a", mo = m02a, dat = s[Nsp>4],fixed = c('Year','SD', 'FlockSize', 'BodyMass', 'rad','rad','Temp','StringencyIndex'), trans = c("","log","log","log","sin","cos","","") , categ = 'Year', outdir = 'Outputs/modelAss/')
    m_ass(name = "Table S2 - 2b", mo = m02b, dat = s[Nsp>4], fixed = c('Year','SD', 'FlockSize', 'BodyMass', 'rad','rad','Temp','StringencyIndex'), trans = c("","log","log","log","sin","cos","","") , categ = 'Year', outdir = 'Outputs/modelAss/')
    m_ass(name = "Table S2 - 2c", mo = m02c, dat = s[Nsp>4], fixed = c('Year','SD', 'FlockSize', 'BodyMass', 'rad','rad','Temp','StringencyIndex'), trans = c("","log","log","log","sin","cos","","") , categ = 'Year', outdir = 'Outputs/modelAss/')
    m_ass(name = "Table S2 - 3a", mo = m03a, dat = s[Nsp>9], fixed = c('Year','SD', 'FlockSize', 'BodyMass', 'rad','rad','Temp','StringencyIndex'), trans = c("","log","log","log","sin","cos","","") , categ = 'Year', outdir = 'Outputs/modelAss/')
    m_ass(name = "Table S2 - 3b", mo = m03b, dat = s[Nsp>9], fixed = c('Year','SD', 'FlockSize', 'BodyMass', 'rad','rad','Temp','StringencyIndex'), trans = c("","log","log","log","sin","cos","","") , categ = 'Year', outdir = 'Outputs/modelAss/')
    m_ass(name = "Table S2 - 3c", mo = m03c, dat = s[Nsp>9], fixed = c('Year','SD', 'FlockSize', 'BodyMass', 'rad','rad','Temp','StringencyIndex'), trans = c("","log","log","log","sin","cos","","") , categ = 'Year', outdir = 'Outputs/modelAss/')
# Figure 2, S3, S5
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
  
    # Fig 2 & left panel of S5- plot from files
        anas = rasterGrob(change_col('#CCCCCC',readPNG('Data/Pics/Anas.png')))
        columba = rasterGrob(change_col('#CCCCCC',readPNG('Data/Pics/Columba.png')))
        Dendrocopos = rasterGrob(change_col('#CCCCCC',readPNG('Data/Pics/Dendrocopos.png')))
        Larus = rasterGrob(change_col('#CCCCCC',readPNG('Data/Pics/Larus_flip.png')))
        Picus = rasterGrob(change_col('#CCCCCC',readPNG('Data/Pics/Picus.png')))
        Motacilla = rasterGrob(change_col('#CCCCCC',readPNG('Data/Pics/Motacilla.png')))
        Erithacus = rasterGrob(change_col('#CCCCCC',readPNG('Data/Pics/Erithacus.png')))
        Phoenicurus = rasterGrob(change_col('#CCCCCC',readPNG('Data/Pics/Phoenicurus.png')))
        Turdus = rasterGrob(change_col('#CCCCCC',readPNG('Data/Pics/Turdus.png')))
        Sylvia = rasterGrob(change_col('#CCCCCC',readPNG('Data/Pics/Sylvia.png')))
        Parus = rasterGrob(change_col('#CCCCCC',readPNG('Data/Pics/Parus_flip.png')))
        Sitta = rasterGrob(change_col('#CCCCCC',readPNG('Data/Pics/Sitta.png')))
        Pica = rasterGrob(change_col('#CCCCCC',readPNG('Data/Pics/Pica.png')))
        Garrulus = rasterGrob(change_col('#CCCCCC',readPNG('Data/Pics/Garrulus.png')))
        Corvus = rasterGrob(change_col('#CCCCCC',readPNG('Data/Pics/Corvus.png')))
        Sturnus = rasterGrob(change_col('#CCCCCC',readPNG('Data/Pics/Sturnus.png')))
        Passer = rasterGrob(change_col('#CCCCCC',readPNG('Data/Pics/Passer_flip.png')))
        Fringilla = rasterGrob(change_col('#CCCCCC',readPNG('Data/Pics/Fringilla.png')))
        other = rasterGrob(change_col('#CCCCCC',readPNG('Data/Pics/other_flip.png')))

        ann_text <- data.frame(FID_avg.0 = 8, FID_avg.1 = 10,lab = "Text",
                       genus2 = factor('Anas',levels = c('Anas', 'Larus', 'Columba', 'Dendrocopos','Picus', 'Motacilla','Erithacus','Phoenicurus','Turdus', 'Sylvia','Parus','Sitta','Pica','Garrulus','Corvus','Sturnus','Passer','Fringilla','other')))
        ann_text2 <- data.frame(FID_avg.0 = 6, FID_avg.1 = 3,lab = "Text",
                       genus2 = factor('Larus',levels = c('Anas', 'Larus', 'Columba', 'Dendrocopos','Picus', 'Motacilla','Erithacus','Phoenicurus','Turdus', 'Sylvia','Parus','Sitta','Pica','Garrulus','Corvus','Sturnus','Passer','Fringilla','other')))

        aw2 = data.frame(FID_avg.0 = c(11.25,11.25), FID_avg.1 = c(3.5,5.8), genus2 = factor('Larus',levels = c('Anas', 'Larus', 'Columba', 'Dendrocopos','Picus', 'Motacilla','Erithacus','Phoenicurus','Turdus', 'Sylvia','Parus','Sitta','Pica','Garrulus','Corvus','Sturnus','Passer','Fringilla','other')))

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
          annotation_custom2(Phoenicurus, data=o[genus2 == 'Phoenicurus'], xmin = 0.05, xmax =0.35, ymax = 2.7)+
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
          geom_line(data = aw2, col = "grey80", lwd = 0.25)+
          geom_text(data = ann_text,label = "No difference", col = "grey80",angle = 45, size = 2) + 
          geom_text(data = ann_text2,label = "Species mean / site", col = "grey60", size = 2,) +
          facet_wrap(~genus2) +
          #geom_phylopic(data = o, aes(image = uid),  color = "grey80", size = o$size) + # ,
          scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = FALSE))  +
          scale_x_continuous("Before COVID-19 shutdown - flight initiation distance [m]", expand = c(0, 0), trans = 'log10') +
          scale_y_continuous("During COVID-19 shutdown - flight initiation distance [m]", expand = c(0, 0), trans = 'log10') +
          #labs(title = "Species means per sampling location")+
          theme_MB  +
          theme(
                  plot.title = element_text(size=7),
                  strip.background = element_blank(),
                  #panel.spacing = unit(1, "mm"),
                  legend.position = c(1, 0.025),
                  legend.justification = c(1, -0.05)
                  )  
        gg <- ggplotGrob(g) #gg$layout$name
        ggx <- gtable_filter_remove(gg, name = paste0("axis-b-", c(2, 4), "-4"),
                                         trim = FALSE)
        #grid.draw(ggx)
        #ggsave('Outputs/Fig_2_width-114mm.png',ggx, width=4.5,height=4.5,dpi=600) # 11.43cm # with label on top
        ggsave('Outputs/Fig_2_width-122mm.png',ggx, width=4.8,height=4.5,dpi=600) # 12.2cm # with label inside
    
    # Fig S5 right panel
        g2 =     
        ggplot(aw, aes(x = resid_FID_avg.0, y = resid_FID_avg.1)) + 
          geom_point(pch = 21, alpha = 0.7, aes(fill = Country), col = 'white') + 
          geom_abline(intercept = 0, slope = 1, lty =3, col = "grey80")+
          facet_wrap(~genus2) +
          #geom_phylopic(data = o, aes(image = uid),  color = "grey80", size = o$size) + 
          scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = TRUE))  +
          scale_x_continuous("Before COVID-19 shutdown - residual escape distance", expand = c(0, 0)) +
          scale_y_continuous("During COVID-19 shutdown - residual escape distance", expand = c(0, 0)) +
          #labs(title = "Species means per sampling location")+
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
    # Fig S5 combine
        grid.draw(cbind(cbind(ggx,ggx2 ,  size = "last")))
        ggsave('Outputs/Fig_S5.png',cbind(ggx,ggx2 , size = "last"), width=4.8*2,height=4.5,dpi=600)
    
    # Fig S5 legend
        # correlation between mean fid and residual fid
            ggplot(a, aes(x = resid_FID_avg)) + geom_histogram()
            ggplot(a, aes(x=FID_avg, y = resid_FID_avg)) + 
                stat_smooth() + 
                geom_point() +
                stat_cor(aes(label = ..r.label..),  label.x = 0.3, size = 2) +
                scale_x_continuous(trans = 'log')
            ggplot(dxx, aes(x=FID, y = resid_FID)) + 
                stat_smooth() + 
                geom_point() +
                stat_cor(aes(label = ..r.label..),  label.x = 0.3, size = 2) +
                scale_x_continuous(trans = 'log')    
            cor(log(dxx$FID), dxx$resid_FID)
            cor(log(a$FID_avg), a$resid_FID_avg)
            ggplot(aw, aes(x = N.0-N.1)) +geom_histogram()
            nrow(aw[abs(N.0-N.1)>2])
            nrow(aw[!abs(N.0-N.1)>2])
     
    # not used - plot from web
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
    
    # Fig S3
      aw[, sp2 := gsub(" ", "\n", Species)]
      ann_text <- data.frame(FID_avg.0 = 8, FID_avg.1 = 10,lab = "Text",
                       Species = factor('Aegithalos caudatus',levels = levels(as.factor(aw$Species))))
      ann_text$sp2 = gsub(" ", "\n", ann_text$Species)
      g3 = 
        ggplot(aw, aes(x = FID_avg.0, y = FID_avg.1)) + 
          #geom_errorbar(aes(ymin = FID_avg.1-SD.1, ymax = FID_avg.1+SD.1, col = Country), width = 0) +
          #geom_errorbar(aes(xmin = FID_avg.0-SD.0, xmax = FID_avg.0+SD.0, col = Country), width = 0) +
          #geom_point(pch = 21, alpha = 0.7, aes(col = Country)) + 
          geom_point(pch = 21, alpha = 0.7, aes(fill = Country), col = 'white') + 
            #ggtitle ("Sim based")+
          geom_abline(intercept = 0, slope = 1, lty =3, col = "grey80")+
          geom_text(data = ann_text,label = "No difference", col = "grey80",angle = 45, size = 2) + 
          facet_wrap(~sp2) +
          #geom_phylopic(data = o, aes(image = uid),  color = "grey80", size = o$size) + # ,
          scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = FALSE))  +
          scale_x_continuous("Before COVID-19 shutdown - flight initiation distance [m]", expand = c(0, 0), trans = 'log10') +
          scale_y_continuous("During COVID-19 shutdown - flight initiation distance [m]", expand = c(0, 0), trans = 'log10') +
          labs(title = "Species means per sampling location")+
          theme_MB  +
          theme(
                  plot.title = element_text(size=7),
                  strip.background = element_blank(),
                  #panel.spacing = unit(1, "mm"),
                  legend.position = c(0.96, 0.0),
                  legend.justification = c(1, 0)
                  )  
        gg3 <- ggplotGrob(g3) #gg2$layout$name
        ggx3 <- gtable_filter_remove(gg3, name = c(paste0("axis-b-", c(2, 4), "-7"), "axis-b-6-6"),
                                         trim = FALSE)
        #grid.draw(ggx2)
        ggsave('Outputs/Fig_S3_species.png',ggx3, width=13.5,height=17.5,unit = 'cm', dpi=600) # 11.43cm
      # not used S3 version
        g = 
        ggplot(aw, aes(x = FID_avg.0, y = FID_avg.1)) + 
          #geom_errorbar(aes(ymin = FID_avg.1-SD.1, ymax = FID_avg.1+SD.1, col = Country), width = 0) +
          #geom_errorbar(aes(xmin = FID_avg.0-SD.0, xmax = FID_avg.0+SD.0, col = Country), width = 0) +
          #geom_point(pch = 21, alpha = 0.7, aes(col = Country)) + 
          geom_point(pch = 21, alpha = 0.7, aes(fill = Country), col = 'white') + 
            #ggtitle ("Sim based")+
          geom_abline(intercept = 0, slope = 1, lty =3, col = "grey80")+
          geom_text(data = ann_text,label = "No difference", col = "grey80",angle = 45, size = 2) + 
          facet_wrap(~Species) +
          #geom_phylopic(data = o, aes(image = uid),  color = "grey80", size = o$size) + # ,
          scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = FALSE))  +
          scale_x_continuous("Before COVID-19 shutdown - flight initiation distance [m]", expand = c(0, 0), trans = 'log10') +
          scale_y_continuous("During COVID-19 shutdown - flight initiation distance [m]", expand = c(0, 0), trans = 'log10') +
          labs(title = "Species means per sampling location")+
          theme_MB  +
          theme(
                  plot.title = element_text(size=7),
                  strip.background = element_blank(),
                  #panel.spacing = unit(1, "mm"),
                  legend.position = c(1, 0.025),
                  legend.justification = c(1, 0)
                  )  
        gg <- ggplotGrob(g) #gg$layout$name
        ggx <- gtable_filter_remove(gg, name = paste0("axis-b-", c(2, 4), "-4"),
                                         trim = FALSE)
        grid.draw(ggx)
        ggsave('Outputs/Fig_S3_names1row.png',gg, width=4.5*1.75,height=4.5*1.75,dpi=600) # 11.43cm
# Figure 3
  ss = s[Nsp>9]
  ss[, sp2 := gsub(" ", "\n", sp)]
  # not used - one row labels
    g = 
    ggplot(ss, aes(x = StringencyIndex, y = FID)) +
      stat_smooth(se = FALSE, aes(colour = 'Locally weighted\nsmoothing'), lwd = 0.5)+ # show_guide=TRUE
      #stat_smooth(method = 'rlm', se = FALSE, col = 'black', lwd = 0.5)+
      geom_point(pch = 21, alpha = 0.7, aes(fill = Country), col = 'white') + 
      facet_wrap(~sp, ncol = 6) +
      scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = FALSE))  +
      scale_x_continuous("Stringency index of governmental COVID-19 restrictions", expand = c(0, 0)) +
      scale_y_continuous("Flight initiation distance [m]", expand = c(0, 0), trans = 'log10') +
      #annotate("text", x = 1, y = 1, label = c(rep("", 52),"Observation"), hjust = -0.08, size = 1) +
      #labs(title = "Species means per sampling location")+
      scale_colour_manual(values=c('grey60'))+
      #scale_color_manual(name = 'try', values = c('LOESS smoothed = "grey60"'))+
      theme_MB  +
      theme(
          plot.title = element_text(size=7),
          strip.background = element_blank(),
          #strip.text.x = element_text(size = 4.5, color="grey30",  margin=margin(1,1,1,1,"mm")),
          #panel.spacing = unit(1, "mm"),
          legend.position = c(0.98, 0.02),
          legend.justification = c(1, 0),
          legend.title = element_blank(),
          #legend.spacing.y = unit(-0.78, "cm")
          #legend.spacing.y = unit(0.02, "cm") use if LOESS smooth text as legend
          legend.spacing.y = unit(-0.9, "cm")
          ) 

      gg <- ggplotGrob(g) #gg$layout$name
      ggx <- gtable_filter_remove(gg, name = paste0("axis-b-", c(2, 4), "-9"), trim = FALSE)
      #grid.draw(ggx)
      ggsave('Outputs/Fig_3_width-152mm_v3.png',ggx, width=19,height=22.2, unit = 'cm',dpi=600)
  # two rows labels
    g2 = 
    ggplot(ss, aes(x = StringencyIndex, y = FID)) +
      stat_smooth(se = FALSE, aes(colour = 'Locally weighted\nsmoothing'), lwd = 0.5)+ # show_guide=TRUE
      #stat_smooth(method = 'rlm', se = FALSE, col = 'black', lwd = 0.5)+
      geom_point(pch = 21, alpha = 0.7, aes(fill = Country), col = 'white') + 
      facet_wrap(~sp2, ncol = 6) +
      scale_fill_viridis(discrete=TRUE,guide = guide_legend(reverse = FALSE))  +
      scale_x_continuous("Stringency index of governmental COVID-19 restrictions", expand = c(0, 0)) +
      scale_y_continuous("Flight initiation distance [m]", expand = c(0, 0), trans = 'log10') +
      #annotate("text", x = 1, y = 1, label = c(rep("", 52),"Observation"), hjust = -0.08, size = 1) +
      #labs(title = "Species means per sampling location")+
      scale_colour_manual(values=c('grey60'))+
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
    ggx2 <- gtable_filter_remove(gg2, name = paste0("axis-b-", c(2, 4), "-9"), trim = FALSE)
    #grid.draw(ggx)
    ggsave('Outputs/Fig_3_width-152mm_2-row_v1.png',ggx2, width=15.24,height=19, unit = 'cm',dpi=600)



# Testing random structure for reviews
   m=lmer(scale(log(FID))~
          scale(Year)+ 
          scale(log(SD))+
          scale(log(FlockSize))+
          scale(log(BodyMass))+
          scale(sin(rad)) + scale(cos(rad)) + 
          #scale(Day)+
          scale(Temp)+
          scale(parks_percent_change_from_baseline)+
          (scale(parks_percent_change_from_baseline)|Country) + (1|IDLocality),  
          data = s, REML = FALSE, 
          control = lmerControl(
              optimizer ='optimx', optCtrl=list(method='nlminb'))
          ) 
    est_m_test = est_out(m, 'test') 

    m=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
        (scale(parks_percent_change_from_baseline)|genus)+(1|Species)+(1|sp_day_year) + (scale(parks_percent_change_from_baseline)|Country) + (1|IDLocality) +(1|sp_loc),
        data = s, REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
        # (1|Year) explains nothing - could stay 
     est_test2 = est_out(m, 'g01a) (scale(parks_percent_change_from_baseline)|genus)+(1|Species)+(1|sp_day_year) + (1|Country) + (scale(parks_percent_change_from_baseline)|IDLocality) +(1|sp_loc)')

# Testing mobility for reviews
  m=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
        (scale(parks_percent_change_from_baseline)|genus)+ (1|Species)+(1|sp_day_year) + (scale(parks_percent_change_from_baseline)|Country) + (1|IDLocality) +(1|sp_loc),
        data = s, REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  

    est_test3 = est_out(m, '(scale(parks_percent_change_from_baseline)|genus)+(1|Species)+(1|sp_day_year) + (scale(parks_percent_change_from_baseline)|Country) + (1|IDLocality) +(1|sp_loc)')    

    mp=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        poly(parks_percent_change_from_baseline,2)+ 
        (poly(parks_percent_change_from_baseline,2)|genus)+ (1|Species)+(1|sp_day_year) + (poly(parks_percent_change_from_baseline,2)|Country) + (1|IDLocality) +(1|sp_loc),
        data = s[!is.na(parks_percent_change_from_baseline)], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
      summary(mp)
      est_out(mp, '')    

  m=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
       (1|genus) + (1|sp_country)+(1|Species)+(1|sp_day_year) + (scale(parks_percent_change_from_baseline)|Country) + (1|IDLocality) +(1|sp_loc),
        data = s, REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  

    est_test3 = est_out(m, '(1|genus) + (scale(parks_percent_change_from_baseline)|sp_country)+(1|Species)+(1|sp_day_year) + (scale(parks_percent_change_from_baseline)|Country) + (1|IDLocality) +(1|sp_loc)')       

 mp=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        poly(parks_percent_change_from_baseline,2)+ 
      (1|Species)+(1|sp_day_year),
        data = s[Country=='Poland'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
   est_out(mp, 'test)')  

mh=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
       (1|Species)+(1|sp_day_year) + (1|IDLocality),
        data = s[Country=='Hungary'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
summary(mh)
est_out(mh, 'test)')  
mh=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        poly(parks_percent_change_from_baseline,2)+ 
       (1|Species)+(1|sp_day_year) + (1|IDLocality),
        data = s[Country=='Hungary'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
summary(mh)
est_out(mh, 'test)')  


mc=lmer(scale(log(FID))~ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
       (1|Species)+(1|sp_day_year) + (1|IDLocality),
        data = s[Country=='Czech Republic'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
summary(mc)
est_out(mc, 'test)')  
mcp=lmer(scale(log(FID))~ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        poly(parks_percent_change_from_baseline,2)+ 
       (1|Species)+(1|sp_day_year) + (1|IDLocality),
        data = s[Country=='Czech Republic'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
summary(mcp)
est_out(mcp, 'test)')  


ma=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
       (1|Species)+(1|sp_day_year) + (1|IDLocality),
        data = s[Country=='Australia'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
summary(ma)
est_out(ma, 'test)')  
map=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        poly(parks_percent_change_from_baseline,2)+ 
       (1|Species)+(1|sp_day_year) + (1|IDLocality),
        data = s[Country=='Australia'], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
summary(map)
est_out(map, 'test)')  

mf=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        scale(parks_percent_change_from_baseline)+ 
       (1|sp_day_year) + (1|IDLocality),
        data = s[Country=='Finland' & !is.na(parks_percent_change_from_baseline)], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
summary(mf)
est_out(mf, 'test)')  
mfp=lmer(scale(log(FID))~ 
        scale(Year)+ 
        scale(log(SD))+ 
        scale(log(FlockSize))+ 
        scale(log(BodyMass))+ 
        scale(sin(rad)) + scale(cos(rad)) +  
        #scale(Day)+ 
        scale(Temp)+ 
        poly(parks_percent_change_from_baseline,2)+ 
      (1|sp_day_year) + (1|IDLocality),
        data = s[Country=='Finland' & !is.na(parks_percent_change_from_baseline)], REML = FALSE,  
        control = lmerControl( 
            optimizer ='optimx', optCtrl=list(method='nlminb')) 
        )  
summary(mfp)
est_out(mfp, 'test)')  
# END