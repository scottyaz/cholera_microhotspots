## utility functions for tau

reload_source <- function(){
    library(readr)
    library(sp)
    library(sf)
    library(tmap)
    library(dplyr)
    library(tidyr)
    library(IDSpatialStats)
    library(RColorBrewer)
    library(ggplot2)

    source("source/utils.R")
}

## loads point data from Kalemie
load_kalemie_data <- function(){
    kal <- read_csv("data/cleaned_simplified_points_kalemie_jit.csv")
    return(kal)
}

## loads Chad data
load_ndj_data <- function(){
    ndj <- read_csv("data/cleaned_simplified_points_ndjamena_jit.csv")
    return(ndj)
}

##' function for main analaysis
##' @param r.mins vector of minimum distances for tau windows
##' @param r.maxs vector of maximum distances for tau windows
##' @param d.mins vector of min time windows for tau 
##' @param d.maxs vector of max time windows for tau 
##' @param point_mat matrix of x,y and time
run_main_analyis <- function(r.mins=c( 0, seq(5,950,10)),
                          r.maxs=c(40, seq(55,1000,10)),
                          d.mins=c(0,5,10,15,20,25,15,1,0),
                          d.maxs=c(5,10,15,20,25,30,30,5,2),
                          point_mat,
                          filename){
  
  if(!is.matrix(ndj)) {ndj <- as.matrix(ndj)}
  
  set.seed(31255)
  r.mids <- (r.maxs + r.mins)/2
  d.mids <- (d.maxs+d.mins)/2
  
  rc <- list()
  
  for(i in seq_along(d.mins)){
    rc[[i]] <- get_tau_and_cis(day_min=d.mins[i],
                                    day_max=d.maxs[i],
                                    my.mat=point_mat,
                                    r.max=r.maxs,
                                    r.min=r.mins,
                                    n_boots=500,
                                    infectious.process=TRUE
                                    )
    print(i)
  }
  
  saveRDS(rc,file = filename)
}



get_riskzone <- function(dat=ndj,
                         d.mins=c(0,5,10,15,20,25,15,1,0),
                         d.maxs=c(5,10,15,20,25,30,30,5,2),
                         filename="GeneratedData/ndj_riskzone.rds",
                         n_boots=500){
  
  set.seed(45352)
  
  r.mins <- c( 0, seq(5,950,10))
  r.maxs <- c(40, seq(55,1000,10))
  r.mids <- (r.maxs + r.mins)/2
  
  d.mids <- (d.maxs+d.mins)/2
  
  rc <- list()
  for(i in seq_along(d.mins)){
    rc[[i]] <- get_tau_and_cis(day_min=d.mins[i],
                               day_max=d.maxs[i],
                               my.mat=dat,
                               r.max=r.maxs,
                               r.min=r.mins,
                               n_boots=n_boots,
                               infectious.process=TRUE,
                               full_boots=TRUE
    )
    print(i)
  }
  
  cat(paste0("saving outbput to ",filename,"\n"))
  saveRDS(rc,file = filename)
  
  return(rc)
}


##' gets tau estimates and confidence intervals using 
##' IDDSpatial Stats implementation of tau
##' @param day_min 
##' @param day_max 
##' @param my.mat 
##' @param r.max 
##' @param r.min 
##' @param n_boots number of bootstrap iterations
##' @param infectious.process flag if we only want to compute tau considering related cases as those that occur **after** an *index* case
##' @param full_boots flag to return full bootstrap replicate matrix (TRUE) or relavant quantiles for confidence intervals
get_tau_and_cis <- function(day_min,
                            day_max,
                            my.mat,
                            r.max,
                            r.min,
                            n_boots,
                            infectious.process=TRUE,
                            full_boots=FALSE){
  
  
  ## make new relation function 
  if(infectious.process){
    new.gen.func <- function(r1,r2,
                             within.range=c(day_min,day_max)){
      is.in.generation(r1,r2,
                       within.range,
                       infectious=TRUE)
    }
  } else {
    
    new.gen.func <- function(r1,r2,
                             within.range=c(day_min,day_max)){
      is.in.generation(r1,r2,
                       within.range,
                       infectious=FALSE)
    }
  }
  
  
  ## get tau estimate
  tau <- get.tau(posmat=my.mat,
                 fun=new.gen.func,
                 r = r.max,
                 r.low=r.min,
                 comparison.type="independent")
  
  if(!full_boots){
    ## get cis rather than full boot
    boot <- get.tau.ci(
      posmat=my.mat,
      fun=new.gen.func,
      r = r.max,
      r.low=r.min,
      boot.iter=n_boots,
      comparison.type="independent",
      ci.low=0.025,
      ci.high=0.975
    )
  } else {
    ## get full bootstrap estimates matrix
    boot <- get.tau.bootstrap(
      posmat=my.mat,
      fun=new.gen.func,
      r = r.max,
      r.low=r.min,
      boot.iter=n_boots,
      comparison.type="independent")
    
  }
  
  return(list(tau,boot))
}

## the tau function needs a 'relation function' to
## determine how pairs are related
##' @param row1 - case line 1 (expecting 3rd element to be time in vector)
##' @param row2 - case line 2
##' @param within.range - days that case should occur within (lower and upper)
##' @param infectious - if infectous then we only look at time in one directoin
##' @return
##' @author asa
is.in.generation <- function(row1,
                             row2,
                             within.range=c(10,20),
                             infectious=TRUE){
  if(infectious){
    
    if((row1[3]-row2[3]) >= within.range[1] & (row1[3]-row2[3]) < within.range[2]){rc=1}
    else{rc=2}
    
  } else {
    if(abs(row1[3]-row2[3]) >= within.range[1] & abs(row1[3]-row2[3]) < within.range[2]){rc=1}
    else {rc=2}
  }
  
  return(rc)
}

## makes 4-panel figure with taus for
## each location and two time subsets
make_figure_1 <- function(ndj_main,kal_main,r.mids){
    palette(brewer.pal(3,"Dark2"))
    par(mfrow=c(2,2),mar=c(2,1,1,1),oma=c(2,3,2,3))
    simple_tau_plot(ndj_main[[1]],r.mids,col=1,add=FALSE,ylim=c(.5,150),xlim=c(20,450),xlab="",ylab="")
    grid()
    text(28,0.5,"days 0-4")
    text(430,100,"(A)")

    simple_tau_plot(kal_main[[1]],r.mids,col=1,add=FALSE,ylim=c(0.5,150),xlim=c(20,450),yaxt="n",xlab="",ylab="",yaxt="n")
    axis(4)
    grid()
    text(27,0.5,"days 0-4")
    text(430,100,"(C)")

    simple_tau_plot(ndj_main[[8]],r.mids,col=2,add=FALSE,ylim=c(.5,150),xlim=c(20,450),xlab="",ylab="")
    grid()
    text(27,0.5,"days 1-4")
    text(430,100,"(B)")
    simple_tau_plot(kal_main[[8]],r.mids,col=2,add=FALSE,ylim=c(.5,150),xlim=c(20,450),xlab="",ylab="",yaxt="n")
    grid()
    text(27,0.5,"days 1-4")
    text(430,100,"(D)")
    axis(4)

    mtext("N'Djamena",outer=T,side=3,line=-1,at=.25)
    mtext("Kalemie",outer=T,side=3,line=-1,at=.75)

    mtext("Distance (meters)",side=1,outer=T)
    mtext("Relative Risk of Cholera",side=2,outer=T,line=1)
    mtext("Relative Risk of Cholera",side=4,outer=T,line=1)

}

#' Makes figure 2
#' @param kalout_tighttime output from run_main_analyis
#' @param ndjout_tighttime output from run_main_analyis
#'
#' @return ggplot figure
make_figure_2 <- function(kalout_tighttime,ndjout_tighttime){
  
  
  melt_tau <- function(tau_list,r.mids,other_tag){
    rc <- lapply(1:length(tau_list),function(x){
      cbind(tau=tau_list[[x]][[1]],ci_l=tau_list[[x]][[2]][1,],ci_u=tau_list[[x]][[2]][2,],r.mids,other=other_tag[x])
    })
    return(rc %>% do.call('rbind',.))
  }
  
  r.mins <- c( 0, seq(5,950,10))
  r.maxs <- c(40, seq(55,1000,10))
  r.mids <- (r.maxs + r.mins)/2
  
  melted_tau_kal <- melt_tau(kalout_tighttime,r.mids,1:19) %>%
    data.frame %>%
    rename(days=other,distance=r.mids) %>% mutate(location='kalemie')
  
  melted_tau_ndj <- melt_tau(ndjout_tighttime,r.mids,1:19) %>%
    data.frame %>%
    rename(days=other,distance=r.mids) %>% mutate(location='ndjamena')
  
  melted_tau <- rbind(melted_tau_ndj,melted_tau_kal)
  
  key_dists <- c(20,50,150,250,350,450)
  melted_tau %>% filter(distance %in% key_dists, days<=15) %>%
    left_join(data.frame(panel=c('A','B','C','D','E','F'),distance=c(20,50,150,250,350,450),x=rep(1,6),y=rep(100,6)))%>%
    mutate(.,distance=ordered(distance,labels=paste0(unique(distance)," meters"))) %>%
    ggplot(aes(days)) +
    geom_ribbon(aes(ymin=ci_l,ymax=ci_u,group=location,fill=location),alpha=.2)+
    geom_line(aes(y=tau,group=location,color=location)) +
    geom_point(aes(y=tau,group=location,color=location)) +
    facet_wrap(~distance) +
    scale_y_log10(breaks=c(.1,.5,1,2,10,100)) +
    geom_hline(yintercept = 1,lty=2,color="grey40") +
    scale_x_continuous(breaks=seq(0,22,by=2)) +
    theme(legend.position=c(.90,.93)) +
    geom_text(aes(label=panel,x=x,y=y),size=7,alpha=.1) +
    xlab("days from primary case presentation") -> time_plot
  
  return(time_plot + theme_minimal())
}


#' makes a simple plot of tau and bootstrap CI
#' @param dat 
#' @param xmids 
#' @param col 
#' @param add 
#' @param add_abline 
#' @param ... other params to send to plot
#' @return
simple_tau_plot <- function(dat,xmids,col=1,add=FALSE,add_abline=TRUE,...){

    if(add){
        lines(xmids,dat[[1]],col=col)
    } else {
        plot(xmids,dat[[1]],type="l",col=col,log="xy",...)
    }

    polygon(c(xmids,rev(xmids)),c(dat[[2]][1,],rev(dat[[2]][2,])),col=AddAlpha(col,.3),border=FALSE)

    if(add_abline){
        ## add line where the CI stays below 1 for at least 2-steps
        try(abline(v=jitter(r.mids[which(dat[[2]][1,]<1)[3]]),col=col,lty=3,lwd=2))
    }
    abline(h=1,col="grey",lty=2)
}


#' Gets risk zone
#'
#' @param bootmat 
#' @param risk_thresh 
#' @param r.mids 
#'
#' @return
get_risk_zone_dist <- function(bootmat,risk_thresh=1.1,r.mids){

    ## check to make sure r.mids is correct dimension
    if(length(r.mids) != ncol(bootmat)) stop("r.mids and bootmat dimension mismatch")

    apply(bootmat,1,function(x) {
        under_thresh <- which(x<risk_thresh)
        r.mids[under_thresh[which(diff(under_thresh)==1)[1]]]
        })

}

#' Title
#'
#' @param tau_out 
#' @param r.mids 
#'
#' @return
tidy_tau_out <- function(tau_out,r.mids){
    bind_rows(
        data.frame(distance=r.mids,type='median',value=tau_out[[1]]),
        data.frame(distance=r.mids,type='ci_l',value=tau_out[[2]][1,]),
        data.frame(distance=r.mids,type='ci_h',value=tau_out[[2]][2,]))

}


#' Title
#'
#' @param timelist 
#' @param timelist_k 
#'
#' @return
#' @export
#'
#' @examples
make_inhibition_by_time_fig <- function(timelist=timelist,
                                        timelist_k=timelist_k){
  
  par(mfrow=c(2,2),mar=c(.5,.5,.5,.5),oma=c(3,3,1,3),mgp=c(1,.3,0),tck=-.02)
  
  plot(days,sapply(timelist,function(x) x[1,])[1,],
       type="l",log="xy",xaxt="n",ylim=c(.05,150),col=1)
  polygon(c(days,rev(days)),
          c(sapply(timelist,function(x) x[1,])[2,],
            rev(sapply(timelist,function(x) x[1,])[3,])),
          col=AddAlpha(1,.3),
          border=FALSE)
  lines(days,sapply(timelist_k,function(x) x[1,])[1,],col="orange")
  polygon(c(days,rev(days)),
          c(sapply(timelist_k,function(x) x[1,])[2,],
            rev(sapply(timelist_k,function(x) x[1,])[3,])),col=AddAlpha('orange',.3),border=FALSE)
  abline(h=1,lty=2)
  text(18,140,"20 meters")
  axis(3)
  
  plot(days,sapply(timelist,function(x) x[3,])[1,],type="l",log="xy",xaxt="n",ylim=c(.05,150),yaxt="n",col=1)
  polygon(c(days,rev(days)),
          c(sapply(timelist,function(x) x[3,])[2,],
            rev(sapply(timelist,function(x) x[3,])[3,])),col=AddAlpha(1,.3),border=FALSE)
  lines(days,sapply(timelist_k,function(x) x[3,])[1,],col="orange")
  polygon(c(days,rev(days)),
          c(sapply(timelist_k,function(x) x[3,])[2,],
            rev(sapply(timelist_k,function(x) x[3,])[3,])),
          col=AddAlpha('orange',.3),border=FALSE)
  axis(4)
  axis(3)
  abline(h=1,lty=2)
  text(18,140,"40 meters")
  
  plot(days,sapply(timelist,function(x) x[5,])[1,],type="l",log="xy",ylim=c(.05,150),col=1)
  polygon(c(days,rev(days)),
          c(sapply(timelist,function(x) x[5,])[2,],
            rev(sapply(timelist,function(x) x[5,])[3,])),col=AddAlpha(1,.3),border=FALSE)
  lines(days,sapply(timelist_k,function(x) x[5,])[1,],col="orange")
  polygon(c(days,rev(days)),
          c(sapply(timelist_k,function(x) x[5,])[2,],
            rev(sapply(timelist_k,function(x) x[5,])[3,])),col=AddAlpha('orange',.3),border=FALSE)
  abline(h=1,lty=2)
  text(18,140,"100 meters")
  
  plot(days,sapply(timelist,function(x) x[10,])[1,],type="l",log="xy",yaxt="n",ylim=c(.05,150),col=1)
  polygon(c(days,rev(days)),
          c(sapply(timelist,function(x) x[10,])[2,],
            rev(sapply(timelist,function(x) x[10,])[3,])),col=AddAlpha(1,.3),border=FALSE)
  lines(days,sapply(timelist_k,function(x) x[10,])[1,],col="orange")
  polygon(c(days,rev(days)),
          c(sapply(timelist_k,function(x) x[10,])[2,],
            rev(sapply(timelist_k,function(x) x[10,])[3,])),col=AddAlpha('orange',.3),border=FALSE)
  abline(h=1,lty=2)
  axis(4)
  text(18,140,"200 meters")
  mtext("time from 'primary' case",outer=T,side=1,line=1.5)
  mtext(expression(tau),outer=T,side=2,line=1.5)
  mtext(expression(tau),outer=T,side=4,line=1.5)
}

## gets the distance where tau cross null
## returns index
get_threshold_distance <- function(tau_list){
  ci_lows <- tau_list[[2]][1,]
  which_lows <- which(ci_lows<1)
  valid_lows <- which(c(diff(which_lows),100)==1)
  if(length(valid_lows)>0){
    rc <- which_lows[valid_lows[1]]
  } else {
    rc <- NA
  }
  
  return(rc)
}



###################################
## Generic Functions             ##
## from personal .Rprofile       ##
## many functions from elsewhere ##
###################################

##' Adds alpha to a set of colors
##' @param COLORS
##' @param ALPHA
##' @return
AddAlpha <- function(COLORS, ALPHA){
  if(missing(ALPHA)) stop("provide a value for alpha between 0 and 1")
  RGB <- col2rgb(COLORS, alpha=TRUE)
  RGB[4,] <- round(RGB[4,]*ALPHA)
  NEW.COLORS <- rgb(RGB[1,], RGB[2,], RGB[3,], RGB[4,], maxColorValue = 255)
  return(NEW.COLORS)
}


##' wrapper for sending output from function to device
##' (note: this is used by to pdf to create pdfs)
##' @param expr
##' @param dev
##' @param filename
##' @param ...
##' @param verbose
##' @return
##' @author http://nicercode.github.io/blog/2013-07-09-figure-functions/
to.dev <- function(expr, dev, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf("Creating %s\n", filename))
  dev(filename, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

## to make a pdf
to.pdf <- function(expr, filename, ...)
  to.dev(expr, pdf, filename, ...)

## to make a png
to.png <- function(expr, filename, ...)
  to.dev(expr, png, filename,...)


