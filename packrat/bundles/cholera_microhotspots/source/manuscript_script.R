## script file for manuscript
source("source/utils.R")
reload_source()

## bring in main data
kal <- load_kalemie_data()
ndj <- load_ndj_data()

###################
## MAIN ANALYSES ##
###################

## some key numbers for the manuscript
## these are the min and max (and midpoints) for the tau-distance windows
r.mins <- c( 0, seq(5,950,10))
r.maxs <- c(40, seq(55,1000,10))
r.mids <- (r.maxs + r.mins)/2

## time windows
d.mins<- c(0,5,10,15,20,25,15,1,0)
d.maxs<- c(5,10,15,20,25,30,30,5,2)
d.mids <- (d.maxs+d.mins)/2

## these functions run the primary analyses and save the outputs in the <filename>
## this will take a while so we have included the outputs in the repository if needed
run_main_analyis(point_mat=kal %>% as.matrix,filename="generated_data/kal_main.rds")
run_main_analyis(point_mat=ndj %>% as.matrix,filename="generated_data/ndj_main.rds")

## reading outputs
kal_main <- readRDS("generated_data/kal_main.rds")
ndj_main <- readRDS("generated_data/ndj_main.rds")

## looking at a bit more fine scale at time dimension
## similar runs for for generating data for figure 2

## note these were run on cluster to save time
run_main_analyis(point_mat=kal %>% as.matrix,
                 d.mins=c(0,1:18),
                 d.maxs=c(2:20),
                 filename="generated_data/kalout_tighttime.rds")

run_main_analyis(point_mat=ndj %>% as.matrix,
                 d.mins=c(0,1:18),
                 d.maxs=c(2:20),
                 filename="generated_data/ndjout_tighttime.rds")

kalout_tighttime <- readRDS("generated_data/kalout_tighttime.rds")
ndjout_tighttime <- readRDS("generated_data/ndjout_tighttime.rds")


## ------- ##
## Figures ##
## ------- ##

## make figure 1
make_figure_1(ndj_main=ndj_main,kal_main=kal_main,r.mids=r.mids)

## make figure 2
make_figure_2(kalout_tighttime = kalout_tighttime,ndjout_tighttime = ndjout_tighttime)

## ----------------------------------- ##
## some key numbers for the manuscript ##
## ----------------------------------- ##

##################
## first 5-days ##
##################

r.mids[which(ndj_main[[1]][[2]][1,]<1)][2] ## 330
r.mids[which(ndj_main[[8]][[2]][1,]<1)][2] ## 310

r.mids[which(kal_main[[1]][[2]][1,]<1)][2] ## 220
r.mids[which(kal_main[[8]][[2]][1,]<1)][2] ## 90

## but it first crosses at 100
r.mids[which(kal_main[[1]][[2]][1,]<1)][1] ## 210
r.mids[which(kal_main[[8]][[2]][1,]<1)][1] ## 80

## those living within 40-meters of a case
ndj_main[[1]][[1]][1]
ndj_main[[1]][[2]][,1]

kal_main[[1]][[1]][1]
kal_main[[1]][[2]][,1]

## those at 100-meters from a case
ndj_main[[1]][[1]][9]
ndj_main[[1]][[2]][,9]

kal_main[[1]][[1]][9]
kal_main[[1]][[2]][,9]

#################
## day 0 and 1 ##
#################

## those living within 20-meters of a case
ndj_main[[9]][[1]][1]
ndj_main[[9]][[2]][,1]

kal_main[[9]][[1]][1]
kal_main[[9]][[2]][,1]

## extent within the first days
r.mids[which(ndj_main[[9]][[2]][1,]<1)][2] ## 340
r.mids[which(kal_main[[9]][[2]][1,]<1)][2] ## 80


## at 100 meteres
ndj_main[[9]][[1]][9]
ndj_main[[9]][[2]][,9]

kal_main[[9]][[1]][9]
kal_main[[9]][[2]][,9]

####################
## exluding day 0 ##
####################

r.mids[which(ndj_main[[8]][[2]][1,]<1)][2] ## 310
r.mids[which(kal_main[[8]][[2]][1,]<1)][2] ## 90

## those living within 40-meters of a case
ndj_main[[8]][[1]][1]
ndj_main[[8]][[2]][,1]

kal_main[[8]][[1]][1]
kal_main[[8]][[2]][,1]

## those at 100-meters from a case
kal_main[[8]][[1]][9]
kal_main[[8]][[2]][,9]

ndj_main[[8]][[1]][9]
ndj_main[[8]][[2]][,9]


#########################
## longer time periods ##
#########################
## 2-4 weeks before
r.mids[which(ndj_main[[7]][[2]][2,]>1)]
r.mids[which(kal_main[[7]][[2]][2,]>1)]

## those living within 20-meters of a case
1/ndj_main[[7]][[1]][1]
1/ndj_main[[7]][[2]][,1]

1/kal_main[[7]][[1]][1]
1/kal_main[[7]][[2]][,1]

## those at 100-meters from a case
kal_main[[1]][[1]][9]
kal_main[[1]][[2]][,9]

ndj_main[[8]][[1]][9]
ndj_main[[8]][[2]][,9]


######################################################
## Now lets consider zones of increased risk to be  ##
## the max. extent where the RR is greater than 1.1 ##
######################################################

ndj_risk <- get_riskzone(
    dat=ndj %>% as.matrix,
    d.mins=c(0,1,0),
    d.maxs = c(5,5,1),
    filename="generated_data/ndj_riskzone.rds",
    n_boots=1000)

ndj_risk <- readRDS(file="generated_data/ndj_riskzone.rds")
quantile(get_risk_zone_dist(ndj_risk[[1]][[2]],risk_thresh=1.2,r.mids=r.mids),probs=c(.025,.5,.975),na.rm=T)
quantile(get_risk_zone_dist(ndj_risk[[2]][[2]],risk_thresh=1.2,r.mids=r.mids),probs=c(.025,.5,.975),na.rm=T)

kal_risk <- get_riskzone(
    dat=kal %>% as.matrix,
    d.mins=c(0,1,0),
    d.maxs = c(5,5,1),
    filename="generated_data/kal_riskzone.rds",
    n_boots=1000)

kal_risk <- readRDS(file="generated_data/kal_riskzone.rds")
quantile(get_risk_zone_dist(kal_risk[[1]][[2]],risk_thresh=1.2,r.mids=r.mids),probs=c(.025,.5,.975))
quantile(get_risk_zone_dist(kal_risk[[2]][[2]],risk_thresh=1.2,r.mids=r.mids),probs=c(.025,.5,.975))


## now only for day 1
ndj_risk_day1 <- get_riskzone(
    dat=ndj %>% as.matrix,
    d.mins=c(0),
    d.maxs = c(1),
    filename="generated_data/ndj_riskzone_day1.rds",
    n_boots=1000)

kal_risk_day1 <- get_riskzone(
    dat=kal %>% as.matrix,
    d.mins=c(0),
    d.maxs = c(1),
    filename="generated_data/kal_riskzone_day1.rds",
    n_boots=1000)

##########################
## supplemental figures ##
##########################


## ----------------------------------------------------- ##
## Estimating tau at 3 different points in each epidemic ##
## ----------------------------------------------------- ##

## now run analyses

##t1
runme_inf_kal(d.mins=c(0,1),
              d.maxs=c(5,5),
              kal=kal  %>% filter(day <= 200) %>% as.matrix ,
              filename="GeneratedData/kalout_50mWin_10mspace_days0_to_200.rds")

runme_inf_ndj(d.mins=c(0,1),
              d.maxs=c(5,5),
              ndj=ndj  %>% filter(time <= (120-71)) %>% as.matrix ,
              filename="GeneratedData/ndjout_50mWin_10mspace_days0_to_120.rds")

##t2
runme_inf_kal(d.mins=c(0,1),
              d.maxs=c(5,5),
              kal=kal %>% filter(day>200 & day<=300) %>% as.matrix,
              filename="GeneratedData/kalout_50mWin_10mspace_days200_to_300.rds")

runme_inf_ndj(d.mins=c(0,1),
              d.maxs=c(5,5),
              ndj=ndj  %>% filter(time > (120-71) & time <= (150-71)) %>% as.matrix ,
              filename="GeneratedData/ndjout_50mWin_10mspace_days120_to_150.rds")

#t3
runme_inf_kal(d.mins=c(0,1),
              d.maxs=c(5,5),
              kal=kal %>% filter(day>300) %>% as.matrix,
              filename="GeneratedData/kalout_50mWin_10mspace_days300_381.rds")


runme_inf_ndj(d.mins=c(0,1),
              d.maxs=c(5,5),
              ndj=ndj  %>% filter(time > (150-71) & time <=(231-71)) %>% as.matrix ,
              filename="GeneratedData/ndjout_50mWin_10mspace_days150_to_231.rds")


kal_t1 <- readRDS("GeneratedData/kalout_50mWin_10mspace_days0_to_200.rds")
kal_t2 <- readRDS("GeneratedData/kalout_50mWin_10mspace_days200_to_300.rds")
kal_t3 <- readRDS("GeneratedData/kalout_50mWin_10mspace_days300_381.rds")

ndj_t1 <- readRDS("GeneratedData/ndjout_50mWin_10mspace_days0_to_120.rds")
ndj_t2 <- readRDS("GeneratedData/ndjout_50mWin_10mspace_days120_to_150.rds")
ndj_t3 <- readRDS("GeneratedData/ndjout_50mWin_10mspace_days150_to_231.rds")


kal_time_periods <- bind_rows(
    tidy_tau_out(kal_main[[1]],r.mids=r.mids) %>% mutate(time_period="all"),
    tidy_tau_out(kal_t1[[1]],r.mids=r.mids) %>% mutate(time_period="t1"),
    tidy_tau_out(kal_t2[[1]],r.mids=r.mids) %>% mutate(time_period="t2"),
    tidy_tau_out(kal_t3[[1]],r.mids=r.mids) %>% mutate(time_period="t3"))


ndj_time_periods <- bind_rows(
    tidy_tau_out(ndj_main[[1]],r.mids=r.mids) %>% mutate(time_period="all"),
    tidy_tau_out(ndj_t1[[1]],r.mids=r.mids) %>% mutate(time_period="t1"),
    tidy_tau_out(ndj_t2[[1]],r.mids=r.mids) %>% mutate(time_period="t2"),
    tidy_tau_out(ndj_t3[[1]],r.mids=r.mids) %>% mutate(time_period="t3"))


kal_tps <- kal_time_periods %>% ggplot() +
geom_line(data=kal_time_periods %>% filter(type=="median"),aes(x=distance,y=value,color=time_period)) +
geom_ribbon(data=kal_time_periods %>% spread(type,value),aes(x=distance,ymin=ci_l,ymax=ci_h,fill=time_period),alpha=.2) +
scale_y_log10() +
coord_cartesian(xlim=c(0, 500)) + theme_minimal() + ylab("relative cholera risk (tau)")


to.pdf(print(kal_tps + theme(legend.position = "bottom")),filename="figures/kal_time_periods_full.pdf",width=12)


ndj_tps <- ndj_time_periods %>% ggplot() +
geom_line(data=ndj_time_periods %>% filter(type=="median"),aes(x=distance,y=value,color=time_period)) +
geom_ribbon(data=ndj_time_periods %>% spread(type,value),aes(x=distance,ymin=ci_l,ymax=ci_h,fill=time_period),alpha=.2) +
scale_y_log10() +
coord_cartesian(xlim=c(0, 500)) + theme_minimal() + ylab("relative cholera risk (tau)")


to.pdf(print(ndj_tps + theme(legend.position = "bottom")),filename="figures/ndj_time_periods_full.pdf",width=12,height=8)

to.pdf(multiplot(ndj_tps + theme(legend.position = "bottom"),kal_tps + theme(legend.position = "bottom")),filename="figures/time_periods_full.pdf",width=12,height=14)


simple_tau_plot(kal_t1[[1]],
                r.mids,col=1,add=FALSE,,xlab="",ylab="")
    grid()
    text(28,0.5,"days 0-4")
    text(430,100,"(A)")


runme_inf_ndj(d.mins=c(0,1),
              d.maxs=c(5,5),
              filename="GeneratedData/ndjout_50mWin_10mspace_days0_to_120.rds")


