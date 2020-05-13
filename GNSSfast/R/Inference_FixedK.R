
###################################
######## Inference procedure for a fixed K
Seg_funct_totK <-function(Data,var.est.month,K,graphK,lmin,lyear,threshold,tol){
  var.est.t=var.est.month[as.numeric(Data$month)]

  period=periodic_estimation_tot_init(Data,var.est.t,lyear)
  auxiliar_data <- Data
  auxiliar_data$signal=Data$signal-period$predict
  segmentation=SegMonthlyVarianceK(auxiliar_data,K,graphK,lmin,var.est.t)

  maxIter = 100
  Diff    = 2*tol
  Iter=0

  while ((Diff  > tol) & (Iter < maxIter))
  {
    Iter = Iter +1
    auxiliar_data$signal=Data$signal-segmentation$mean.est.t
    periodi=periodic_estimation_tot(auxiliar_data,var.est.t,lyear)

    auxiliar_data$signal=Data$signal-periodi$predict
    segmi=SegMonthlyVarianceK(auxiliar_data,K,graphK,lmin,var.est.t)

    if (Iter == 2)
    {
      t2 = c(period$predict,segmentation$mean.est.t)
    }
    if (Iter == 3)
    {
      t1 = c(period$predict,segmentation$mean.est.t)
      t0 = c(periodi$predict,segmi$mean.est.t)
      tp0 = (t2-t1)/sum((t1==t2)+((t2-t1)^2)) + (t0-t1)/sum((t1==t0)+((t0-t1)^2))
      tp0 = t1 + tp0 / sum((tp0==0) + tp0^2)
    }
    if (Iter > 3)
    {
      t2 = t1
      t1 = t0
      t0 = c(periodi$predict,segmi$mean.est.t)
      tp1 = tp0
      tp0 = (t2-t1)/sum((t1==t2)+((t2-t1)^2)) + (t0-t1)/sum((t1==t0)+((t0-t1)^2))
      tp0 = t1 + tp0 / sum((tp0==0) + tp0^2)
      Diff = sum((tp0-tp1)^2)
    }

    period=periodi
    segmentation=segmi
  }


  segmK=c()
  segmK$Tmu   = segmentation$Tmu
  segmK$SSwg  = segmentation$SSwg
  segmK$LogLg = segmentation$LogLg
  segmK$f     = period$predict
  segmK$coeff = period$coeff
  return(segmK)

}


Seg_funct_selbK <-function(Data,var.est.month,K,graphK,lmin,lyear,threshold,tol){
  var.est.t=var.est.month[as.numeric(Data$month)]

  period=periodic_estimation_selb_init(Data,var.est.t,lyear,threshold)
  auxiliar_data <- Data
  auxiliar_data$signal=Data$signal-period$predict
  segmentation=SegMonthlyVarianceK(auxiliar_data,K,graphK,lmin,var.est.t)

  maxIter = 100
  Diff  = 2*tol
  Iter  = 0

  while ((Diff  > tol) & (Iter < maxIter))
  {
    Iter = Iter +1
    auxiliar_data$signal=Data$signal-segmentation$mean.est.t
    periodi=periodic_estimation_selb(auxiliar_data,var.est.t,lyear,threshold)

    auxiliar_data$signal=Data$signal-periodi$predict
    segmi=SegMonthlyVarianceK(auxiliar_data,K,graphK,lmin,var.est.t)

    if (Iter == 2)
    {
      t2 = c(period$predict,segmentation$mean.est.t)
    }
    if (Iter == 3)
    {
      t1 = c(period$predict,segmentation$mean.est.t)
      t0 = c(periodi$predict,segmi$mean.est.t)
      tp0 = (t2-t1)/sum((t1==t2)+((t2-t1)^2)) + (t0-t1)/sum((t1==t0)+((t0-t1)^2))
      tp0 = t1 + tp0 / sum((tp0==0) + tp0^2)
    }
    if (Iter > 3)
    {
      t2 = t1
      t1 = t0
      t0 = c(periodi$predict,segmi$mean.est.t)
      tp1 = tp0
      tp0 = (t2-t1)/sum((t1==t2)+((t2-t1)^2)) + (t0-t1)/sum((t1==t0)+((t0-t1)^2))
      tp0 = t1 + tp0 / sum((tp0==0) + tp0^2)
      Diff = sum((tp0-tp1)^2)
    }

    period=periodi
    segmentation=segmi

  }


  segmK=c()
  segmK$Tmu   = segmentation$Tmu
  segmK$SSwg  = segmentation$SSwg
  segmK$LogLg = segmentation$LogLg
  segmK$f     = period$predict
  segmK$coeff = period$coeff
  return(segmK)

}


###################################
######## Robust estimation of the variances
RobEstiMonthlyVariance <- function(Y){
  Kmonth <- length(unique(Y$month))
  z <- stats::aggregate(signal ~ month + year, data = Y, diff)
  sigma.est <- sapply(1:Kmonth,function(i) {
    e <- subset(z,z$month==levels(z$month)[i])
    ee <- unlist(e$signal)
    robustbase::Qn(ee, constant = 1 / (sqrt(2) * stats::qnorm(5/8)))/sqrt(2)

  })
  return(sigma.est)
}




###################################
######## Functions for segmentation
SegMonthlyVarianceK=function(Data,K,graphK,lmin,var.est.t){
  result=list()
  Res.gfpop=c()
  Res.gfpop=gfpop(data = Data$signal, mygraph = graphK, type = "mean",weights=1/var.est.t)


  Tmu=c()
  rupt  = matrix(Inf,ncol = 2 , nrow= K)
  bp    = Res.gfpop$changepoints
  rupt[,2]  = bp
  bp        = bp +1
  rupt[,1]  = c(1, bp[1:K-1])
  Tmu=data.frame(rupt,Res.gfpop$parameters)
  colnames(Tmu) = c("begin","end","mean")
  mean.est.t  = rep(Tmu$mean,diff(c(0,Tmu$end)))


  result$Tmu=Tmu
  result$mean.est.t = mean.est.t
  result$SSwg=Res.gfpop$globalCost
  result$LogLg=sum(log(diff(c(0,Res.gfpop$changepoints))[diff(c(0,Res.gfpop$changepoints))>0]))
  return(result)

}




###################################
######## Functions for functional
periodic_estimation_tot=function(Data,var.est.t,lyear){
  DataF=Data
  DataF$t=c(as.numeric(DataF$date-DataF$date[1]))/86400
  for (i in 1:4){
    cosX=cos(i*DataF$t*(2*pi)/lyear)
    sinX=sin(i*DataF$t*(2*pi)/lyear)
    DataF=cbind(DataF,cosX,sinX)
    colnames(DataF)[(dim(DataF)[2]-1):(dim(DataF)[2])]=c(paste0('cos',i),paste0('sin',i))
  }
  reg=stats::lm(signal~-1+cos1+sin1+cos2+sin2+cos3+sin3+cos4+sin4,weights=1/var.est.t,data=DataF)
  coeff=base::summary(reg)$coefficients[,1]
  result=list()
  result$predict=stats::predict(reg,DataF)
  result$coeff=coeff
  return(result)
}


periodic_estimation_tot_init=function(Data,var.est.t,lyear){
  DataF=Data
  DataF$t=c(as.numeric(DataF$date-DataF$date[1]))/86400
  for (i in 1:4){
    cosX=cos(i*DataF$t*(2*pi)/lyear)
    sinX=sin(i*DataF$t*(2*pi)/lyear)
    DataF=cbind(DataF,cosX,sinX)
    colnames(DataF)[(dim(DataF)[2]-1):(dim(DataF)[2])]=c(paste0('cos',i),paste0('sin',i))
  }
  reg=stats::lm(signal~-1+cos1+sin1+cos2+sin2+cos3+sin3+cos4+sin4,data=DataF)
  coeff=base::summary(reg)$coefficients[,1]
  result=list()
  result$predict=stats::predict(reg,DataF)
  result$coeff=coeff
  return(result)
}


periodic_estimation_selb=function(Data,var.est.t,lyear,threshold=0.001){
  DataF=Data
  DataF$t=c(as.numeric(DataF$date-DataF$date[1]))/86400
  num.col=dim(DataF)[2]
  for (i in 1:4){
    cosX=cos(i*DataF$t*(2*pi)/lyear)
    sinX=sin(i*DataF$t*(2*pi)/lyear)
    DataF=cbind(DataF,cosX,sinX)
    colnames(DataF)[(dim(DataF)[2]-1):(dim(DataF)[2])]=c(paste0('cos',i),paste0('sin',i))
  }
  reg=stats::lm(signal~-1+cos1+sin1+cos2+sin2+cos3+sin3+cos4+sin4,weights=1/var.est.t,data=DataF)
  res.coeff=summary(reg)$coefficients
  names.coeff=rownames(res.coeff)
  #Selection of significant coefficients
  rg=which(res.coeff[,4]<threshold)

  if (length(rg)>=1){
    names.Selected=names.coeff[rg]
    n.Selected=length(names.Selected)
    DataFF=DataF[,c(1:num.col,which(colnames(DataF) %in%  names.Selected ))]
    if (n.Selected >=2){
      a=names.Selected[1]
      for (i in 1:(n.Selected-1)){
        a=paste(a,names.Selected[i+1],sep="+")
      }
    } else {
      a=names.Selected[1]
    }
    reg.Selected=c()
    request=paste(paste0("reg.Selected=stats::lm(signal~-1+",a,",weights=1/var.est.t,data=DataFF)"),sep="")
    eval(parse(text=request))
    pred=reg.Selected$fitted.values
    coeff=reg.Selected$coefficients
  } else {
    pred=rep(0,length(DataF$signal))
    coeff=0
    names(coeff)="no selected coeff"
  }
  result=list()
  result$predict=pred
  result$coeff=coeff
  return(result)
}



periodic_estimation_selb_init=function(Data,var.est.t,lyear,threshold=0.001){
  DataF=Data
  DataF$t=c(as.numeric(DataF$date-DataF$date[1]))/86400
  num.col=dim(DataF)[2]
  for (i in 1:4){
    cosX=cos(i*DataF$t*(2*pi)/lyear)
    sinX=sin(i*DataF$t*(2*pi)/lyear)
    DataF=cbind(DataF,cosX,sinX)
    colnames(DataF)[(dim(DataF)[2]-1):(dim(DataF)[2])]=c(paste0('cos',i),paste0('sin',i))
  }


  reg=stats::lm(signal~-1+cos1+sin1+cos2+sin2+cos3+sin3+cos4+sin4,data=DataF)
  res.coeff=summary(reg)$coefficients
  names.coeff=rownames(res.coeff)
  #Selection of significant coefficients
  rg=which(res.coeff[,4]<threshold)

  if (length(rg)>=1){
    names.Selected=names.coeff[rg]
    n.Selected=length(names.Selected)
    DataFF=DataF[,c(1:num.col,which(colnames(DataF) %in%  names.Selected ))]
    if (n.Selected >=2){
      a=names.Selected[1]
      for (i in 1:(n.Selected-1)){
        a=paste(a,names.Selected[i+1],sep="+")
      }
    } else {
      a=names.Selected[1]
    }
    reg.Selected=c()
    request=paste(paste0("reg.Selected=stats::lm(signal~-1+",a,",data=DataFF)"),sep="")
    eval(parse(text=request))
    pred=reg.Selected$fitted.values
    coeff=reg.Selected$coefficients
  } else {
    pred=rep(0,length(DataF$signal))
    coeff=0
    names(coeff)="no selected coeff"
  }
  result=list()
  result$predict=pred
  result$coeff=coeff
  return(result)
}

GraphBuilding<- function(lmin,Kmax){
  myGraph=list()
  if (lmin==1){
    myGraph[[1]] <- graph(
      Edge(0, 0, "null")
    )
    for (k in 2:Kmax){
      myGraph.add<-graph(
        Edge(k-2, k-1,"std"),
        Edge(k-1, k-1, "null")
      )
      myGraph[[k]] <- rbind(myGraph[[k-1]],myGraph.add)
    }
  }

  if (lmin>1){
    myGraph=list()
    #k=1
    subgraph1 <- c()
    subgraph1 <-graph(Edge(0, paste0("wait",0,".",1,sep=""), "null"))
    subgraph2 <- graph()
    if(lmin>2){
      for (l in 2:(lmin-1)){
        subgraph2 <- rbind(subgraph2,graph(Edge(paste0("wait",0,".",l-1,sep=""), paste0("wait",0,".",l,sep=""))))
      }
    }
    subgraph3 <-graph()
    subgraph3 <- graph(Edge(paste0("wait",0,".",lmin-1,sep=""), paste0("wait",0,".",lmin-1,sep=""), "null"))
    myGraph[[1]] <-rbind(subgraph1,subgraph2,subgraph3)

    #for each k
    for (k in 2:Kmax){
      subgraph1 <- graph()
      subgraph1 <- graph(Edge(myGraph[[k-1]]$state1[dim(myGraph[[k-1]])[1]], k-1, "std"),Edge(k-1, paste0("wait",k-1,".",1,sep=""), "null"))
      subgraph2 <- graph()
      if(lmin>2){
        for (l in 2:(lmin-1)){
          subgraph2 <- rbind(subgraph2,graph(Edge(paste0("wait",k-1,".",l-1,sep=""), paste0("wait",k-1,".",l,sep=""))))
        }
      }
      subgraph3 <-graph()
      subgraph3 <- graph(Edge(paste0("wait",k-1,".",lmin-1,sep=""), paste0("wait",k-1,".",lmin-1,sep=""), "null"))

      myGraph[[k]] <- rbind(myGraph[[k-1]],subgraph1,subgraph2,subgraph3)
      myGraph[[k-1]]<-rbind(myGraph[[k-1]],StartEnd(start = 0, end = myGraph[[k-1]]$state1[dim(myGraph[[k-1]])[1]]))
    }

    myGraph[[Kmax]]<-rbind(myGraph[[Kmax]],StartEnd(start = 0, end = myGraph[[Kmax]]$state1[dim(myGraph[[Kmax]])[1]]))
  }
  return(myGraph)
}
