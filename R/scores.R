##requires hydroGOF package

#obs: matrix or data frame containing the observations
     #stations in columns and days in rows
#sim: matrix or data frame containing the estimates
     #stations in columns and days in rows. This matrix 
     #has to have the same NAs of obs
#alts: numeric vector with the altitude ranges that you
      #want to analyse
#dates: vector of dates (as.Date()) of the same length
       #as the nrow() of obs and sim

scores <- function(obs, sim, alts, dates){
  
  #ordering the stations
  sts <- sts[match(colnames(obs), sts$ID),]
  sim <- sim[,match(colnames(obs), colnames(sim))]
  
  
  ####################
  #10 quantiles stats
  ####################   
  dec <- matrix(NA, ncol = 2, nrow = 10)
  onum <- as.numeric(as.matrix(obs))
  snum <- as.numeric(as.matrix(sim))
  #remove nas
  ww <- which(is.na(onum))
  onum <- onum[-c(ww)]; snum <- snum[-c(ww)]
  #remove zeros
  w <- which(onum > 0  & snum > 0)  
  onum <- onum[w]; snum <- snum[w]
  #compute quantiles
  q <- quantile(onum, 
                probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), 
                na.rm = T)
  q <- c(0, q)
  #breaks
  brks <- cut(onum, q)
  
  ###################
  print('Computing statistics by quantiles')
  errors <- data.frame(matrix(NA,nrow=10,ncol=7),stringsAsFactors=F)
  names(errors) <- c('%OBS','%PRED','MAE','ME','RM','RSD','Range'); rownames(errors) <- paste('d',1:10,sep='')
  dir.create('./val', showWarnings = F)
  
  no <- sum(onum)
  ns <- sum(snum)
  for(i in 1:10){
    o <- onum[which(levels(brks)[i] == brks)]
    s <- snum[which(levels(brks)[i] == brks)]
    errors[i,1] <- round(sum(o) * 100 / no, 1)
    errors[i,2] <- round(sum(s) * 100 / ns, 1)
    errors[i,3] <- round(mae(s,o),2)
    errors[i,4] <- round(me(s,o),2)
    errors[i,5] <- round(mean(s, na.rm = T) / mean(o, na.rm = T), 2)
    errors[i,6] <- round(sd(s, na.rm = T) / sd(o, na.rm = T), 2)
    if(i == 1){ errors[i,7] <- paste('>0 to', round(max(o), 2))} else{
      errors[i,7] <- paste(round(min(o),2),'to',round(max(o),2))}
  }
  
  write.table(errors ,'./val/stats_by_quantiles.txt', quote = F, row.names = F, sep='\t')
  
  ##
  pdf('./val/histogram_obs_vs_pred.pdf',height=4,width=9,pointsize=10,useDingbats=F)
  lo <- log10(onum)
  lsi <-log10(snum)
  hist(log10(onum), main = '', xlab = '', freq = F, xlim = c(min(lo, lsi), max(lo, lsi)))
  par(new = T)
  hist(log10(snum), lty = 2, xaxt = 'n',
       yaxt = 'n', main = '', xlab = 'logPcp (mm)', ylab = '',
       freq = F,
       xlim = c(min(lo, lsi), max(lo, lsi)))
  legend('topleft', lty = c(1,2), legend = c('Observed','Predicted'), 
         bty = 'n')
  dev.off()
  
  ##
  ####################
  #altitude stats
  ####################   
  print('Computing statistics by altitudes')
  #breaks
  brks <- cut(sts$ALT, alts)
  errors <- data.frame(matrix(NA,nrow=length(levels(brks)),ncol=7),stringsAsFactors=F)
  names(errors) <- c('%OBS','%PRED','MAE','ME','RM','RSD','N')
  
  no <- sum(onum)
  ns <- sum(snum)
  for(i in 1:length(levels(brks))){
    o <- as.numeric(as.matrix(obs[,which(levels(brks)[i] == brks)]))
    wna <- which(is.na(o))
    if(length(wna) > 0) o <- o[-c(wna)]
    s <- as.numeric(as.matrix(sim[,which(levels(brks)[i] == brks)]))
    wna <- which(is.na(s))
    if(length(wna) > 0) s <- s[-c(wna)]
    errors[i,1] <- round(sum(o) * 100 / no, 1)
    errors[i,2] <- round(sum(s) * 100 / ns, 1)
    errors[i,3] <- round(mae(s,o),2)
    errors[i,4] <- round(me(s,o),2)
    errors[i,5] <- round(mean(s, na.rm = T) / mean(o, na.rm = T), 2)
    errors[i,6] <- round(sd(s, na.rm = T) / sd(o, na.rm = T), 2)
    errors[i,7] <- length(which(levels(brks)[i] == brks))
  }
  write.table(errors ,'./val/stats_by_altitudes.txt', quote = F, row.names = F, sep='\t')
  
  ###
  ####################
  #number of zeros
  ####################   
  print('Computing statistics of wet/dry days')
  pdf('./val/number of zeros.pdf',height = 6,width = 6)
  d <- density(colSums(obs == 0, na.rm = T))
  plot(d, main = '', 
       xlab = 'Number of zeros',xlim=c(0,max(d$x)), col ='darkgrey')
  lines(density(colSums(sim == 0, na.rm = T)), lty = 2)
  legend('topright',col=c('darkgrey','black'),lty=c(1,2),
         legend = c('Observed','Predicted'),bty='n')
  dev.off()
  
  print('Computing statistics by stations and days')
  pdf('./val/medias_medians_95_estaciones_y_dias.pdf', width=10, height=6, pointsize = 12)
  par(mfrow=c(2,3),mar=c(4,3.9,1,2))
  Lab.palette <- colorRampPalette(c("white", "white", "blue", "orange", "red4"), space = "Lab")
  #means all
  om1=colMeans(obs,na.rm=T)
  pm1=colMeans(sim,na.rm=T)
  smoothScatter(om1,pm1,ylim=c(min(om1,pm1),max(om1,pm1)),xlim=c(min(om1,pm1),max(om1,pm1)),
                xlab='Observed (mm)',ylab='Predicted (mm)',colramp=Lab.palette,nbin=200)
  abline(sd(pm1),1,lty=2,col=colors()[220])
  abline(-sd(pm1),1,lty=2,col=colors()[220])
  legend('bottomright',bty='n',legend=c(paste('Pearson',round(cor(om1,pm1,use="pairwise.complete.obs"),3))))
  abline(0,1) 
  
  #medians wet
  ome=rep(NA,nrow(obs)); pme=rep(NA,nrow(obs))
  o95=rep(NA,nrow(obs)); p95=rep(NA,nrow(obs))
  for(i in 1:ncol(obs)){
    w=which(obs[,i]>0 & sim[,i]>0)
    if(length(w)>10){
      ome[i] = median(obs[w,i])
      pme[i] = median(sim[w,i])
      o95[i] = quantile(obs[w,i],probs=0.95)
      p95[i] = quantile(sim[w,i],probs=0.95)
    } 
  }
  
  smoothScatter(ome,pme,ylim=c(min(ome,pme,na.rm=T),max(ome,pme,na.rm=T)),
                xlim=c(min(ome,pme,na.rm=T),max(ome,pme,na.rm=T)),xlab='Observed (mm)',
                ylab='Predicted (mm)',colramp=Lab.palette,nbin=200)
  abline(sd(pme,na.rm=T),1,lty=2,col=colors()[220])
  abline(-sd(pme,na.rm=T),1,lty=2,col=colors()[220])
  legend('bottomright',bty='n',legend=c(paste('Pearson',round(cor(ome,pme,use="pairwise.complete.obs"),3))))
  abline(0,1) 
  
  #95th wet
  smoothScatter(o95,p95,ylim=c(min(o95,p95,na.rm=T),max(o95,p95,na.rm=T)),
                xlim=c(min(o95,p95,na.rm=T),max(o95,p95,na.rm=T)),xlab='Observed (mm)',
                ylab='Predicted (mm)',colramp=Lab.palette,nbin=200)
  abline(sd(p95,na.rm=T),1,lty=2,col=colors()[220])
  abline(-sd(p95,na.rm=T),1,lty=2,col=colors()[220])
  legend('bottomright',bty='n',legend=c(paste('Pearson',round(cor(o95,p95,use="pairwise.complete.obs"),3))))
  abline(0,1) 
  
  #by days
  ####################################################
  #means all
  om1=rowMeans(obs,na.rm=T)
  pm1=rowMeans(sim,na.rm=T)
  smoothScatter(om1,pm1,ylim=c(min(om1,pm1),max(om1,pm1)),xlim=c(min(om1,pm1),max(om1,pm1)),
                xlab='Observed (mm)',ylab='Predicted (mm)',colramp=Lab.palette,nbin=200)
  abline(sd(pm1),1,lty=2,col=colors()[220])
  abline(-sd(pm1),1,lty=2,col=colors()[220])
  legend('bottomright',bty='n',legend=c(paste('Pearson',round(cor(om1,pm1,use="pairwise.complete.obs"),3))))
  abline(0,1)
  
  #medians wet
  ome=rep(NA,nrow(obs)); pme=rep(NA,nrow(obs))
  o95=rep(NA,nrow(obs)); p95=rep(NA,nrow(obs))
  
  for(i in 1:nrow(obs)){
    w=which(obs[i,]>0 & sim[i,]>0)
    if(length(w)>10){
      ome[i] = median(as.numeric(obs[i,w]))
      pme[i] = median(as.numeric(sim[i,w]))
      o95[i] = quantile(as.numeric(obs[i,w]),probs=0.95)
      p95[i] = quantile(as.numeric(sim[i,w]),probs=0.95)
    } 
  }
  
  smoothScatter(ome,pme,ylim=c(min(ome,pme,na.rm=T),max(ome,pme,na.rm=T)),
                xlim=c(min(ome,pme,na.rm=T),max(ome,pme,na.rm=T)),xlab='Observed (mm)',
                ylab='Predicted (mm)',colramp=Lab.palette,nbin=200)
  abline(sd(pme,na.rm=T),1,lty=2,col=colors()[220])
  abline(-sd(pme,na.rm=T),1,lty=2,col=colors()[220])
  legend('bottomright',bty='n',legend=c(paste('Pearson',round(cor(ome,pme,use="pairwise.complete.obs"),3))))
  abline(0,1) 
  
  #95th wet
  smoothScatter(o95,p95,ylim=c(min(o95,p95,na.rm=T),max(o95,p95,na.rm=T)),
                xlim=c(min(o95,p95,na.rm=T),max(o95,p95,na.rm=T)),xlab='Observed (mm)',
                ylab='Predicted (mm)',colramp=Lab.palette,nbin=200)
  abline(sd(p95,na.rm=T),1,lty=2,col=colors()[220])
  abline(-sd(p95,na.rm=T),1,lty=2,col=colors()[220])
  legend('bottomright',bty='n',legend=c(paste('Pearson',round(cor(o95,p95,use="pairwise.complete.obs"),3))))
  abline(0,1) 
  
  dev.off()
  
  ##monthly validation
  print('Computing statistics by months')
  mon <- substr(dates,1,7)    
  suma_mon <- function(x, mon){return(aggregate(x,by=list(mon),FUN='sum')[,2])}
  omsum <- apply(obs, 2, suma_mon, mon)
  pmsum <- apply(sim, 2, suma_mon, mon)
  
  meses <- substr(unique(mon),6,7)
  errorsMon <- data.frame(matrix(NA,ncol=6,nrow=12))
  names(errorsMon) <- c('MAE','ME','RM','RSD','Pearson','Range')
  
  for(i in 1:12){
    w <- which(i==as.numeric(meses))
    k=pmsum[w,]; k=k[-c(which(is.na(k)))]
    ko=omsum[w,]; ko=ko[-c(which(is.na(ko)))]
    k[which(k==0)]=1
    ko[which(ko==0)]=1
    
    errorsMon[i,1] <- round(mae(k,ko)/10,2)
    errorsMon[i,2] <- round(me(k,ko)/10,2)
    errorsMon[i,3] <- round(mean(k)/mean(ko),2)
    errorsMon[i,4] <- round(sd(k)/sd(ko),2)
    errorsMon[i,5] <- round(cor(k,ko),2)
    errorsMon[i,6] <- paste(round(min(as.numeric(as.matrix(omsum[w,])),na.rm=T)/10,1),'to',round(max(as.numeric(as.matrix(omsum[w,])),na.rm=T)/10,1))
  }
  write.table(errorsMon,'./val/monthly_stats.txt',sep='\t',quote=F,row.names=F)
  
  #tabla de contingencia
  ####################################
  print('Computing contingency table')
  
  onum <- as.numeric(as.matrix(obs))
  pnum <- as.numeric(as.matrix(sim))
  contTable <- matrix(NA,ncol=2,nrow=2)
  colnames(contTable) <- c('Obs=0','Obs>0')
  rownames(contTable) <- c('Pred=0','Pred>0')
  #Obs=0
  kk = pnum[which(onum==0)]
  contTable[1,1]=round(length(which(kk==0))*100/length(which(onum==0)),2)
  contTable[2,1]=round(length(which(kk>0))*100/length(which(onum==0)),2)
  
  #Obs>0
  kk = pnum[which(onum>0)]
  contTable[1,2]=round(length(which(kk==0))*100/length(which(onum>0)),2)
  contTable[2,2]=round(length(which(kk>0))*100/length(which(onum>0)),2)
  
  write.table(contTable,'./val/contingency_table.txt',sep='\t',row.names=T,quote=F) 
  
}