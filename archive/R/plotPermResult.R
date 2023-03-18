plotPermResult = function(x,method=c('growth','counts'),index,calendar='BP',...)
{
  x = x[[index]]

  if (method=='growth') {x = x$obs.growthrate;ylab='Annual Growth Rate'}
  if (method == 'counts'){x = x$obs.cnts;ylab='Site Counts'}
  
  
  # Boom and Bust Handling ####
  booms <- which(x$obs>x$hi)
  busts <- which(x$obs<x$lo)
  baseline <- rep(NA,nrow(x))
  colpts = rep('grey',nrow(x))
  colpts[booms] = 'red'
  colpts[busts] = 'blue'
  
  boomPlot <- baseline
  if (length(booms)>0){ boomPlot[booms]=x$obs[booms] }
  bustPlot <- baseline
  if (length(busts)>0){ bustPlot[busts]=x$obs[busts] }           
  
  boomBlocks <- vector("list")
  counter <- 0
  state <- "off"
  for (i in 1:length(boomPlot)){
    if (!is.na(boomPlot[i])&state=="off"){
      counter <- counter+1
      boomBlocks <- c(boomBlocks,vector("list",1))
      boomBlocks[[counter]] <- vector("list",2)
      boomBlocks[[counter]][[1]] <- boomPlot[i]
      boomBlocks[[counter]][[2]] <- x[i,"BP"]
      state <- "on"
    }
    if (state=="on"){
      if (!is.na(boomPlot[i])){
        boomBlocks[[counter]][[1]] <- c(boomBlocks[[counter]][[1]],boomPlot[i])
        boomBlocks[[counter]][[2]] <- c(boomBlocks[[counter]][[2]],x[i,"BP"])
      }
      if (is.na(boomPlot[i])){
        state <- "off"
      }
    }   
  }
  
  bustBlocks <- vector("list")
  counter <- 0
  state <- "off"
  for (i in 1:length(bustPlot)){
    if (!is.na(bustPlot[i])&state=="off"){
      counter <- counter+1
      bustBlocks <- c(bustBlocks,vector("list",1))
      bustBlocks[[counter]] <- vector("list",2)
      bustBlocks[[counter]][[1]] <- bustPlot[i]
      bustBlocks[[counter]][[2]] <- x[i,"BP"]
      state <- "on"
    }
    if (state=="on"){
      if (!is.na(bustPlot[i])){
        bustBlocks[[counter]][[1]] <- c(bustBlocks[[counter]][[1]],bustPlot[i])
        bustBlocks[[counter]][[2]] <- c(bustBlocks[[counter]][[2]],x[i,"BP"])
      }
      if (is.na(bustPlot[i])){
        state <- "off"
      }
    }   
  }
  
  
  
   
  
  # Actual Plot ####
  if(calendar=='BP')
  {
    trange=x$BP;
    xlim=rev(rev(trange))
    xlab='BP'
  }
  if(calendar=='BCAD')
  {
    trange=rcarbon::BPtoBCAD(x$BP);
    xlim=range(trange)
    xlab = 'BCE/CE'
    xlab = ifelse(all(trange<0),'BCE',xlab)
    xlab = ifelse(all(trange>0),'CE',xlab)
  }
  
  # Base Plot
  plot(0,0,type='n',xlim=xlim,ylim=range(c(x$obs,x$lo,x$hi),na.rm=T),xlab=xlab,ylab=ylab,axes=F,...)
  #Envelope
  polygon(c(trange,rev(trange)),c(x$lo,rev(x$hi)),col='lightgrey',border=NA)
  
  if (length(booms)>0){
    for (i in 1:length(boomBlocks)){
      bbb = unique(boomBlocks[[i]][[2]])
      index = which(x$BP%in%bbb)
      if (calendar=='BCAD'){bbb=BPtoBCAD(bbb)}
      polygon(c(bbb,rev(bbb)),c(x$obs[index],rev(x$hi[index])),border=NA,col=rgb(0.8,0.36,0.36,0.5))
    }  
  }
    
  if (length(busts)>0){
    for (i in 1:length(bustBlocks)){
      bbb = unique(bustBlocks[[i]][[2]])
      index = which(x$BP%in%bbb)
      if (calendar=='BCAD'){bbb=BPtoBCAD(bbb)}
      polygon(c(bbb,rev(bbb)),c(x$obs[index],rev(x$lo[index])),border=NA,col=rgb(0.25,0.41,0.88,0.5))
    }  
  }
  lines(trange,x$obs,lwd=1.5)
  if(method=='growth'){abline(h=0,lty=2,col='darkgrey')}
  
  # Add axes
  axis(2)
  tt=axTicks(1)
  if(any(tt==0)){tt[which(tt==0)]=1}
  axis(1,at=tt,labels=abs(tt))
  box()
}