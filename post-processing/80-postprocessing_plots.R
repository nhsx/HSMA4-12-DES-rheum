# Author: Martina Fonseca, NHSE
# Date: 27/09/2022
# Contact: martina.fonseca@nhs.net
#
# R post-processing script for figures used in End-Report (available upon request)
# for the simpy model results from
# The Use Case (V3, an illustrative fictitious trust): Baseline and What-if Scenarios, namely:
# - KPIs
# - Behaviour over simulation time of KPIs (waiting list size, median waiting time, resources occupied)
# The Convergence over number of replications (1-50)
#
# NOTE: Raw simpy model outputs (appointment and patient logs) not provided in repo due to file size
# Script can be adapted to show outputs from new user runs
#
#
# More commented version TBA
#
###############################################################################

setwd("C:/Users/martina.fonseca_nhsx/Documents/HSMA/au-HSMA-DOP/DES-rheum") # amend as appropriate


library(tidyverse)


#### Bundle high-level ####

batchf = 'bundleV3/'
summaryf = 'summary/'
scenarios <- c("out_central_rep30",
               "out_central_rep30_pifu20",
               "out_central_rep30_pifu50",
               "out_central_rep30_pifu30",
               "out_central_rep30_pifu10",
               "out_central_rep30_pifu20_interpifu20",
               "out_central_rep30_ag15",
               "out_central_rep30_pifu20_ag15")

scenariosid <- c("O","A","B1","B2","B3","C","D","E")
scenariosshorthand <- c("Baseline",
                        "20% on PIFU pathway",
                        "50% on PIFU pathway",
                        "30% on PIFU pathway",
                        "10% on PIFU pathway",
                        "20% on PIFU pathway, appointment rate only decreases 16%",
                        "15% of first appointments avoided with A&G",
                        "20% on PIFU pathway, 15% of first appointments avoided with A&G")

scenariosorder <- scenariosshorthand[c(1,3,4,2,5,6,7,8)]


#### Baseline visualisations - observation ####

#subf = 'out_sand'
#subf = 'out_central_15_pifu50'
batchf = 'bundleV3/'
subf = 'out_central_rep30'
mon_app <- read.csv(paste0("./",batchf,subf,"/batch_mon_appointments.csv"))
mon_audit <- read.csv(paste0("./",batchf,subf,"/batch_mon_audit.csv"))

mon_audit_red <- mon_audit %>% select(time,priority.1.patients.waiting,priority.3.patients.waiting,resources.occupied,rep)
mon_audit_red$scenario = "baseline"

### histogram

t_max = mon_app$start_q %>% max()

t_window = 365

mon_app <- mon_app %>% mutate(end_q=start_q+q_time)

ggplot(data=mon_app %>% filter(end_q > t_max - t_window, priority==3), aes(x=q_time,fill=priority ))+
  geom_histogram()+
  facet_wrap(~rep)

ggplot(data=mon_app %>% filter(end_q > t_max - t_window), aes(x=q_time,fill=factor(priority)))+
  geom_histogram()+
  facet_grid(priority~rep)

ggplot(data=mon_app %>% filter(end_q > t_max - t_window), aes(x=q_time,fill=factor(priority)))+
  geom_histogram()

mon_app$priority %>% table()


mystep = 365/4
mon_app <- mon_app %>% mutate(step = end_q %/% mystep * mystep)



log_attributes_l_runsum_l <- mon_app %>%
  filter(priority == 3) %>%
  group_by(rep,step) %>%
  summarise(
    quant_n = c(names(quantile(q_time,c(0.25,0.5,0.75,0.9),na.rm=T)),"mean"),
            quant=c(quantile(q_time,c(0.25,0.5,0.75,0.9),na.rm=T),mean(q_time,na.rm=T)))


log_attributes_l_runsum_l_batch <- log_attributes_l_runsum_l %>%
  group_by(step,quant_n) %>% 
  #mutate(prehandoverNAperc = factor(prehandoverNAperc),
  #       JCToffsiteincrement = factor(JCToffsiteincrement)) %>% 
  summarise(Kp_Tqueuetime=mean(quant,na.rm=T)) %>% 
  ungroup()


ggplot(data=log_attributes_l_runsum_l) +
  geom_boxplot(aes(x=factor(step) , y = quant)) +
  facet_wrap(.~quant_n,nrow=1)+
  labs(x="Time interval", y="Distribution across replications")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggplot(data=log_attributes_l_runsum_l %>% filter(quant_n=="50%")) +
  geom_boxplot(aes(x=factor(round(step/365,2)) , y = quant/7,colour=quant_n)) +
  #facet_wrap(.~quant_n,nrow=1)+
  labs(x="Quarterly time window start (years)", y="RTT median waiting time (weeks)",
       title="",
       footnote = "Distribution across replications")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
                       legend.position="none")+
  theme(text = element_text(size = 16))
  
ggsave(paste0(batchf,summaryf,"baseline_RTTmed_dt.png"),
  unit="cm",
       dpi=300,
       width=20,
       height=20)

ggsave(paste0(batchf,summaryf,"baseline_RTTmed_dt.svg"),
       unit="cm",
       dpi=300,
       width=20,
       height=20)


ggplot(data=mon_audit %>% mutate(a="a")) +
  geom_boxplot(aes(x=factor(round(time/365,1)) , y = priority.3.patients.waiting,colour=a)) +
  labs(x="Audit timepoint (years)", y="RTT waiting list size (no.)",
       footnote = "Distribution across replications")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        legend.position="none")+
  theme(text = element_text(size = 16))+
  scale_x_discrete(breaks=c(5,6,7,8,50))

ggsave(paste0(batchf,summaryf,"/baseline_RTT_WLmean_dt.png"),
       unit="cm",
       dpi=300,
       width=20,
       height=20)

ggsave(paste0(batchf,summaryf,"/baseline_RTT_WLmean_dt.svg"),
       unit="cm",
       dpi=300,
       width=20,
       height=20)



#### Baseline visualisations - observation as 8 years ####

#subf = 'out_sand'
#subf = 'out_central_15_pifu50'
batchf = 'bundleV3/'
subf = 'out_central_rep30_obs8y'
mon_app <- read.csv(paste0("./",batchf,subf,"/batch_mon_appointments.csv"))
mon_audit <- read.csv(paste0("./",batchf,subf,"/batch_mon_audit.csv"))

mon_audit_red <- mon_audit %>% select(time,priority.1.patients.waiting,priority.3.patients.waiting,resources.occupied,rep)
mon_audit_red$scenario = "baseline"

mon_app <- mon_app %>% mutate(end_q=start_q+q_time)
mystep = 365/4
mon_app <- mon_app %>% mutate(step = end_q %/% mystep * mystep)



log_attributes_l_runsum_l <- mon_app %>%
  #filter(priority == 3) %>%
  group_by(type,rep,step) %>%
  summarise(
    quant_n = c(names(quantile(q_time,c(0.25,0.5,0.75,0.9),na.rm=T)),"mean","n"),
    quant=c(quantile(q_time,c(0.25,0.5,0.75,0.9),na.rm=T),mean(q_time,na.rm=T),n()))


log_attributes_l_runsum_l <- log_attributes_l_runsum_l %>% mutate(phase = ifelse(step<5*365,"Warm-up","Observation"))

ggplot(data=log_attributes_l_runsum_l %>% filter(quant_n=="50%")) +
  geom_boxplot(aes(x=factor(round(step/365,2)) , y = quant/7,colour=phase)) +
  facet_wrap(.~type,ncol=1)+
  labs(x="Quarterly time window start (years)", y="Median waiting time (weeks)",
       title="",
       footnote = "Distribution across replications")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        legend.position="bottom")+
  theme(text = element_text(size = 16))

ggsave(paste0(batchf,summaryf,"baseline8_WTmed_type_dt.png"),
       unit="cm",
       dpi=300,
       width=20,
       height=20)

ggsave(paste0(batchf,summaryf,"baseline8_WTmed_type_dt.svg"),
       unit="cm",
       dpi=300,
       width=20,
       height=20)


ggplot(data=log_attributes_l_runsum_l %>% filter(quant_n=="n")) +
  geom_boxplot(aes(x=factor(round(step/365,2)) , y = quant,colour=phase)) +
  facet_wrap(.~type,ncol=1,scales="free_y")+
  labs(x="Quarterly time window start (years)", y="Number seen (patients)",
       title="",
       footnote = "Distribution across replications")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        legend.position="bottom")+
  theme(text = element_text(size = 16))

ggsave(paste0(batchf,summaryf,"baseline8_nseen_type_dt.png"),
       unit="cm",
       dpi=300,
       width=20,
       height=20)

ggsave(paste0(batchf,summaryf,"baseline8_nseen_type_dt.svg"),
       unit="cm",
       dpi=300,
       width=20,
       height=20)


mon_audit_l <- mon_audit %>% pivot_longer(cols=c("priority.1.patients.waiting",
                                                 "priority.3.patients.waiting",
                                                 "resources.occupied"),
                                                  names_to = "quant_n",
                                                  values_to  ="quant")


mon_audit_l <- mon_audit_l %>% mutate(quant_n = ifelse(quant_n == "priority.3.patients.waiting","RTT waiting-list-size",
                                                       ifelse(quant_n == "priority.1.patients.waiting","Follow-up waiting-list-size",
                                                              quant_n)),
                                      )

mon_audit_l <- mon_audit_l %>% mutate(phase = ifelse(time<5*365,"Warm-up","Observation"))


ggplot(data=mon_audit_l) +
  geom_boxplot(aes(x=factor(round(time/365,1)) , y = quant,colour=phase)) +
  labs(x="Audit timepoint (years)", y="",
       footnote = "Distribution across replications")+
  theme_minimal()+
  facet_wrap(.~quant_n,ncol=1,
             scales="free_y")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        legend.position="none")+
  theme(text = element_text(size = 16))+
  scale_x_discrete(breaks=c(0,3,5,6,7,8,50))

ggsave(paste0(batchf,summaryf,"/baseline8_audit_dt.png"),
       unit="cm",
       dpi=300,
       width=20,
       height=20)

ggsave(paste0(batchf,summaryf,"/baseline8_audit_dt.svg"),
       unit="cm",
       dpi=300,
       width=20,
       height=20)



# first to follow-up ratio

data_ffup <- log_attributes_l_runsum_l %>%
  ungroup() %>%
  group_by(rep,step,phase) %>%
  filter(quant_n=="n") %>%
  #summarise(f = sum(quant))
  #mutate(n=as.double(n)) %>%
  summarise(firsttofollowup = sum(quant*(type=="Traditional"))/sum(quant*(type!="Traditional")))

ggplot(data=data_ffup) +
  geom_boxplot(aes(x=factor(round(step/365,2)) , y = firsttofollowup,colour=phase)) +
  labs(x="Quarterly time window start (years)", y="Follow-up ratio (to first)",
       title="",
       footnote = "Distribution across replications")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        legend.position="bottom")+
  theme(text = element_text(size = 16))


ggsave(paste0(batchf,summaryf,"/baseline8_followupfirstratio_dt.svg"),
       unit="cm",
       dpi=300,
       width=20,
       height=20)

ggsave(paste0(batchf,summaryf,"/baseline8_followupfirstratio_dt.png"),
       unit="cm",
       dpi=300,
       width=20,
       height=20)


#### Time trend evolution and comparison ####


mon_app <- data.frame()
mon_audit <- data.frame()
i=1
for (subf in scenarios){
  
  scenarioid <- scenariosid[i]
  scenarioname <- scenariosshorthand[i]
  mon_app_ <- read.csv(paste0("./",batchf,subf,"/batch_mon_appointments.csv"))
  mon_audit_<- read.csv(paste0("./",batchf,subf,"/batch_mon_audit.csv"))
  mon_app_ <- mon_app_ %>% mutate(scenario = subf,scenarioid=scenarioid,scenarioname=scenarioname)
  mon_audit_ <- mon_audit_ %>% select(time,priority.3.patients.waiting,resources.occupied,rep)
  mon_audit_ <- mon_audit_ %>% mutate(scenario = subf,scenarioid=scenarioid,scenarioname=scenarioname)
  mon_audit <- mon_audit %>% bind_rows(mon_audit_)
  mon_app <- mon_app %>% bind_rows(mon_app_)
  i=i+1

}


mystep = 365/4
mon_app <- mon_app %>% mutate(end_q=start_q+q_time)
mon_app <- mon_app %>% mutate(step = end_q %/% mystep * mystep)

log_attributes_l_runsum_l <- mon_app %>%
  filter(priority == 3) %>%
  group_by(scenario,scenarioname,scenarioid,rep,step) %>%
  summarise(
    quant_n = c(names(quantile(q_time,c(0.25,0.5,0.75,0.9),na.rm=T)),"mean"),
    quant=c(quantile(q_time,c(0.25,0.5,0.75,0.9),na.rm=T),mean(q_time,na.rm=T)))


log_attributes_l_runsum_l$scenarioname <- factor(log_attributes_l_runsum_l$scenarioname,
                                                 levels=scenariosorder)


log_attributes_l_runsum_l$scenarioname_ <- unlist(lapply(strwrap(log_attributes_l_runsum_l$scenarioname,
                                                                 width=20, simplify=FALSE), paste, 
                                                         collapse="\n"))
scenariosorder_ <- unlist(lapply(strwrap(scenariosorder,width=20, simplify=FALSE), paste, collapse="\n"))
log_attributes_l_runsum_l$scenarioname_ <- factor(log_attributes_l_runsum_l$scenarioname_,
                                                 levels=scenariosorder_)

ggplot(data=log_attributes_l_runsum_l %>% filter(quant_n=="50%")) +
  geom_boxplot(aes(x=factor(round(step/365,2)) , y = quant/7,colour=scenarioname)) +
  facet_wrap(.~scenarioname_,nrow=2)+
  labs(x="Time (years)", y="RTT median waiting time (weeks)",
       footnote = "Distribution across replications")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        legend.position="none")+
  theme(text = element_text(size = 16))

ggsave(paste0(batchf,summaryf,"/scenarios_RTT_q05_dt.png"),
       unit="cm",
       dpi=300,
       width=30,
       height=15)

ggsave(paste0(batchf,summaryf,"/scenarios_RTT_q05_dt.svg"),
       unit="cm",
       dpi=300,
       width=30,
       height=15)


ggplot(data=log_attributes_l_runsum_l %>% filter(quant_n=="50%",scenarioid %in% c("O","A","B1"))) +
  geom_boxplot(aes(x=factor(round(step/365,2)) , y = quant/7,colour=scenarioname)) +
  facet_wrap(.~scenarioname_,nrow=1)+
  labs(x="Time (years)", y="RTT median waiting time (weeks)",
       footnote = "Distribution across replications")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        legend.position="none")+
  theme(text = element_text(size = 16))

ggsave(paste0(batchf,summaryf,"/scenariosOAB2_RTT_q05_dt.png"),
       unit="cm",
       dpi=300,
       width=30,
       height=15)

ggsave(paste0(batchf,summaryf,"/scenariosOAB2_RTT_q05_dt.svg"),
       unit="cm",
       dpi=300,
       width=30,
       height=15)


mon_audit$scenarioname <- factor(mon_audit$scenarioname,levels=scenariosorder)
mon_audit$scenarioname_ <- unlist(lapply(strwrap(mon_audit$scenarioname,
                                                                 width=20, simplify=FALSE), paste, 
                                                         collapse="\n"))
mon_audit$scenarioname <- factor(mon_audit$scenarioname,levels=scenariosorder)
mon_audit$scenarioname_ <- factor(mon_audit$scenarioname_,levels=scenariosorder_)


ggplot(data=mon_audit) +
  geom_boxplot(aes(x=factor(round(time/365,1)) , y = priority.3.patients.waiting,colour=scenario)) +
  facet_wrap(.~scenarioname_,nrow=2)+
  labs(x="Audit timepoint (years)", y="RTT waiting list size (no.)",
       footnote = "Distribution across replications")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        legend.position="none")+
  theme(text = element_text(size = 16))+
  scale_x_discrete(breaks=c(5,6,7,8,50))

ggsave(paste0(batchf,summaryf,"/scenarios_RTT_WLmean_dt.png"),
       unit="cm",
       dpi=300,
       width=30,
       height=15)

ggsave(paste0(batchf,summaryf,"/scenarios_RTT_WLmean_dt.svg"),
       unit="cm",
       dpi=300,
       width=30,
       height=15)


ggplot(data=mon_audit) +
  geom_boxplot(aes(x=factor(round(time/365,1)) , y = resources.occupied,colour=scenario)) +
  facet_wrap(.~scenarioname_,nrow=2)+
  labs(x="Audit timepoint (years)", y="Resources occupied (no. slots)",
       footnote = "Distribution across replications")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1),
        legend.position="none")+
  theme(text = element_text(size = 16))+
  scale_x_discrete(breaks=c(5,6,7,8,50))

ggsave(paste0(batchf,summaryf,"/scenarios_resourcesoccupied_dt.png"),
       unit="cm",
       dpi=300,
       width=30,
       height=15)

ggsave(paste0(batchf,summaryf,"/scenarios_resourcesoccupied_dt.svg"),
       unit="cm",
       dpi=300,
       width=30,
       height=15)


#### Summary KPIs tabulated ####


mon_kpi <- data.frame() ; i=1
for (subf in scenarios){
  scenarioid <- scenariosid[i]
  scenarioname <- scenariosshorthand[i]
  mon_kpi_i <-  read.csv(paste0("./",batchf,subf,"/batch_kpi.csv"))
  mon_kpi_i <- mon_kpi_i %>% mutate(scenario=subf, scenarioid=scenarioid,scenarioname=scenarioname)
  mon_kpi <- mon_kpi %>% bind_rows(mon_kpi_i)
  i = i+1
}

mon_kpi$scenarioname <- factor(mon_kpi$scenarioname,levels=scenariosorder)


write.csv(mon_kpi,paste0("./",batchf,summaryf,"kpi_bundle.csv"))

kpinow = "RTT_WL_end"
titlenow = "Mean number of RTT patients waiting\n (at end of observation period)"
ggplot(data = mon_kpi %>% filter(KPI==kpinow), aes(x=scenarioname, y = KPI_mean, group=KPI)) +
  geom_point(colour="#1e81b0",size=3) +
  geom_errorbar(aes(ymin=KPI_LCI,ymax=KPI_UCI),width=.3,
                position=position_dodge(0.05),colour="#1e81b0")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        legend.position="bottom")+
  labs(y="KPI, across batches",x="What-if scenario",title=titlenow)


ggsave(paste0("./",batchf,summaryf,"scenarios_",kpinow,".png"),
       unit="cm",
       dpi=300,
       width=15,
       height=15)
ggsave(paste0("./",batchf,summaryf,"scenarios_",kpinow,".svg"),
       unit="cm",
       dpi=300,
       width=15,
       height=15)


kpinow = "RTT_mean"
titlenow = "Waiting time for first outpatient\n (mean weeks over last year)"
ggplot(data = mon_kpi %>% filter(KPI==kpinow), aes(x=scenarioname, y = KPI_mean/7, group=KPI)) +
  geom_point(colour="#1e81b0",size=3) +
  geom_errorbar(aes(ymin=KPI_LCI/7,ymax=KPI_UCI/7),width=.3,
                position=position_dodge(0.05),colour="#1e81b0")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        legend.position="bottom")+
  labs(y="KPI, across batches",x="What-if scenario",title=titlenow)


ggsave(paste0("./",batchf,summaryf,"scenarios_",kpinow,".png"),
       unit="cm",
       dpi=300,
       width=15,
       height=15)
ggsave(paste0("./",batchf,summaryf,"scenarios_",kpinow,".svg"),
       unit="cm",
       dpi=300,
       width=15,
       height=15)


kpinow = "RTT_q0.5"
titlenow = "Waiting time for first outpatient\n (median weeks over last year)"
ggplot(data = mon_kpi %>% filter(KPI==kpinow), aes(x=scenarioname, y = KPI_mean/7, group=KPI)) +
  geom_point(colour="#1e81b0",size=3) +
  geom_errorbar(aes(ymin=KPI_LCI/7,ymax=KPI_UCI/7),width=.3,
                position=position_dodge(0.05),colour="#1e81b0")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        legend.position="bottom")+
  labs(y="KPI, across batches",x="What-if scenario",title=titlenow)


ggsave(paste0("./",batchf,summaryf,"scenarios_",kpinow,".png"),
       unit="cm",
       dpi=300,
       width=15,
       height=15)
ggsave(paste0("./",batchf,summaryf,"scenarios_",kpinow,".svg"),
       unit="cm",
       dpi=300,
       width=15,
       height=15)

kpinow = "resources occupied"
titlenow = "Resources (slots) occupied\n (at end of observation period)"
ggplot(data = mon_kpi %>% filter(KPI==kpinow), aes(x=scenarioname, y = KPI_mean, group=KPI)) +
  geom_point(colour="#1e81b0",size=3) +
  geom_hline(yintercept=43,colour="red",linetype="dashed")+
  geom_errorbar(aes(ymin=KPI_LCI,ymax=KPI_UCI),width=.3,
                position=position_dodge(0.05),colour="#1e81b0")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
        legend.position="bottom")+
  labs(y="KPI, across batches",x="What-if scenario",title=titlenow) +
  expand_limits(x = 0, y = 0)


ggsave(paste0("./",batchf,summaryf,"scenarios_",kpinow,".png"),
       unit="cm",
       dpi=300,
       width=15,
       height=15)
ggsave(paste0("./",batchf,summaryf,"scenarios_",kpinow,".svg"),
       unit="cm",
       dpi=300,
       width=15,
       height=15)



#### Convergence analysis ####

setwd("C:/Users/marti/OneDrive/Documents/GitHub/au-HSMA-DOP/DES-rheum")

#subf = 'out_sand'
#subf = 'out_central_15_pifu50'
subf = 'bundleV3/out_central_rep50'
mon_app <- read.csv(paste0("./",subf,"/batch_mon_appointments.csv"))
mon_app <- mon_app %>% mutate(end_q=start_q+q_time)

mon_audit <- read.csv(paste0("./",subf,"/batch_mon_audit.csv"))
mon_audit_red <- mon_audit %>% select(time,priority.3.patients.waiting,resources.occupied,rep)
mon_audit_red$scenario = "baseline"
mon_audit_red_l <- mon_audit_red %>% pivot_longer(cols=c("priority.3.patients.waiting","resources.occupied"),
                                                  names_to = "quant_n",
                                                  values_to  ="quant")



### per-rep kpis

t_max = mon_app$start_q %>% max()
t_window = 365
mystep = 28*4

log_attributes_l_runsum_l <- mon_app %>%
  filter(end_q > t_max - t_window, priority==3) %>%
  group_by(rep) %>%
  summarise(
    quant_n = c(names(quantile(q_time,c(0.25,0.5,0.75,0.9),na.rm=T)),"mean"),
    quant=c(quantile(q_time,c(0.25,0.5,0.75,0.9),na.rm=T),mean(q_time,na.rm=T)))





### inter-rep (batch) kpis

conv_df <- data.frame()
conv_reps <- c(3,5,10,15,20,25,30,35,40,45,50)

alpha <- 0.05
for (i in conv_reps){
  
  
  mon_kpis_i <- log_attributes_l_runsum_l %>%
    filter(rep < i) %>% # starts at 0 so <
    group_by(quant_n) %>%
    summarise(mean = mean(quant),
              lower = mean(quant) - qt(1- alpha/2, (n() - 1))*sd(quant)/sqrt(n()),
              upper = mean(quant) + qt(1- alpha/2, (n() - 1))*sd(quant)/sqrt(n())) %>%
    mutate(nreps = i)
  
  mon_kpis_wl_i <- mon_audit_red_l %>%
    filter(time==max(time)) %>%
    group_by(quant_n) %>%
    summarise(mean = mean(quant),
              lower = mean(quant) - qt(1- alpha/2, (n() - 1))*sd(quant)/sqrt(n()),
              upper = mean(quant) + qt(1- alpha/2, (n() - 1))*sd(quant)/sqrt(n())) %>%
    mutate(nreps = i)
  
  conv_df <- conv_df %>% bind_rows(mon_kpis_i) %>% bind_rows(mon_kpis_wl_i)
  
}


conv_df <- conv_df %>% mutate(quant_n = ifelse(quant_n == "priority.3.patients.waiting","Waiting-list-size",quant_n))

sel_vars <- c("Waiting-list-size","resources.occupied","WT-50%","WT-mean")

ggplot(data = conv_df %>% filter(quant_n %in% sel_vars) , aes(x=nreps,y=mean,color=quant_n)) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=lower,ymax=upper),width=.3,
                position=position_dodge(0.05))+
  facet_wrap(~quant_n,
             ncol=1,
             scales="free_y")+
  labs(x="Number of replications",y="KPI summary")+
  expand_limits(x = 0, y = 0)+
  theme(legend.position = "none")


ggsave(paste0(subf,"/Convergence_analysis_Abs.png"),
       unit="cm",
       dpi=300,
       width=15,
       height=15)

ggsave(paste0(subf,"/Convergence_analysis_Abs.svg"),
       unit="cm",
       dpi=300,
       width=15,
       height=15)

ggplot(data = conv_df %>% filter(quant_n %in% sel_vars) , aes(x=nreps,y=0,color=quant_n)) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=lower-mean,ymax=upper-mean),width=.3,
                position=position_dodge(0.05))+
  facet_wrap(~quant_n,
             ncol=1,
             scales="free_y")
  #expand_limits(x = 0, y = 0)

ggplot(data = conv_df %>% filter(quant_n %in% sel_vars) , aes(x=nreps,y=0,color=quant_n)) +
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=(lower-mean)/mean*100,ymax=(upper-mean)/mean*100),width=.3,
                position=position_dodge(0.05))+
  labs(x="Number of replications",y="KPI summary - relative error")+
  facet_wrap(~quant_n,
             ncol=1,
             scales="free_y")+
  labs(y="Relative error (%)")+
  theme(legend.position = "none")
#expand_limits(x = 0, y = 0)


ggsave(paste0(subf,"/Convergence_analysis_RE.png"),
       unit="cm",
       dpi=300,
       width=15,
       height=15)

ggsave(paste0(subf,"/Convergence_analysis_RE.svg"),
       unit="cm",
       dpi=300,
       width=15,
       height=15)
