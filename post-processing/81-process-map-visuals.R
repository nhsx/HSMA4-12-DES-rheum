# Author: Martina Fonseca, NHSE
# Date: 27/09/2022
# Contact: martina.fonseca@nhs.net
#
# R post-processing script to create illustrative/stylistic process map and animations from simpy model logs
# This should not be used for a very large simulation (e.g. 1000+ patients), as it will incur memory issues.
#
# # NOTE: Raw simpy model outputs (appointment and patient logs) not provided in repo due to file size
# Script can be adapted to show outputs from new user runs
#
# More commented version TBA
#
###############################################################################

library(here)
library(bupaR)
library(processmapR)
library(tidyverse)
library(processanimateR)
library(webshot)
library(htmlwidgets)
mydate <- "20220921"

bundle="bundleVanim_"
case="out_central_rep1_5w3s_intarr1_pifu20_inter100"
case="out_central_rep1_5w3s_intarr1"

attrib_log <- read.csv(here::here(bundle,case, "batch_mon_appointments.csv"))

attrib_logres <- read.csv(here::here(bundle,case, "appt_result.csv"))

myid = 5
myrep = 0

forevent <- attrib_logres %>%
  filter(rep==myrep) %>%
  mutate(
    end_q = start_q + q_time,
    activity_instance=1:nrow(.),
    resource_id = "Outpatient wait",
    status="start")

my_origin = "2017-04-01"
forevent$time <- as.POSIXct(forevent$start_q*60*60*24,origin=my_origin)

# Discriminate PIFU after decision
forevent <- forevent %>% group_by(P_ID) %>% mutate(aux=end_q-min(end_q),aux2=ifelse(aux<365,99999,aux))
forevent <- forevent %>% mutate(type1=ifelse(type!="Traditional",type,ifelse(aux>365,"Traditional FU","Traditional FU - Year1")))#,
                               # type1 = ifelse(pathway=="PIFU"&type1=="TFU","TFU-Y1",type1))

forevent <- forevent %>% group_by(P_ID) %>%
  mutate(type1 = ifelse(aux2==min(aux2) & type=="Traditional","Traditional FU - Year1",type1))

forevent <- forevent %>% group_by(P_ID) %>%
  mutate(type1 = ifelse(type=="Traditional" & pathway=="PIFU","Traditional FU - Year1",type1))

forevent_comp <- forevent %>% mutate(status="complete")
forevent_comp$time <- as.POSIXct(forevent$end_q*60*60*24,origin=my_origin)

forevent_appt <- forevent_comp %>% filter(type %in% c('First','First-only')) %>% mutate(status="start") %>%
  mutate(resource_id = "Consultant",
         time = time, # s
         activity_instance = activity_instance+max(forevent$activity_instance),
         type1=paste0(type1," appointment"))

forevent_appt_end <- forevent_appt %>% mutate(status="start") %>%
  mutate(time = time + 15*60)


forevent_pifudec <- forevent_comp %>% filter(type=="Traditional",pathway=="PIFU") %>% group_by(P_ID) %>% filter(aux==max(aux)) %>% ungroup()

forevent_pifudec <- forevent_pifudec %>%
  mutate(status="start") %>%
  mutate(resource_id = "Consultant",
         time = time+0.01, # s
         activity_instance = activity_instance+2*max(forevent$activity_instance),
         type1=paste0(pathway, " pathway"),
         type="PIFU")
  

forevent_preferral <- forevent %>% filter(type %in% c('First','First-only')) %>% mutate(status="complete") %>%
  mutate(resource_id = "Consultant",
         time = time - 15*60, # s
         activity_instance = activity_instance+3*max(forevent$activity_instance),
         type1="Pre-referral")

forevent <- bind_rows(forevent,forevent_comp,forevent_appt,forevent_appt_end,forevent_pifudec,forevent_preferral)

forevent <- forevent %>% mutate(type1 = ifelse(type1 %in% c('First','First-only'),"   Referral to Treatment wait   ",type1))

forevent <- forevent %>% arrange(time)


forevent_ori <- forevent

# forevent <- forevent_ori %>%
#   #filter(start_q>365*5) %>%
#   group_by(P_ID) %>%
#   filter(min(start_q)>5*365) %>%
#   #filter(min(end_q)<365*8-365/2) %>%  # exclude those entering system in last 6M
#   filter(min(end_q)<365*6) %>%  # exclude those entering system after year 6
#   ungroup()

 # forevent <- forevent_ori %>%
 #   filter(start_q>365*4)
 
 
 forevent <- forevent_ori %>%
   filter(start_q>365*4) %>%
   group_by(P_ID) %>%
   filter(min(end_q)<365*8-365*0.5) %>%  # exclude those entering system in last 6M
   ungroup()

event_DES <- eventlog(forevent,
                      case_id="P_ID",
                      activity_id="type1",
                      activity_instance_id = "activity_instance",
                      timestamp="time",
                      lifecycle_id="status",
                      resource_id="resource_id"
                      )


event_DES %>% summary()

event_DES %>% process_map()

mapping(event_DES)

n_activities(event_DES)

activity_labels(event_DES)

activities(event_DES)# %>% View()

event_DES$typef <- factor(event_DES$type,levels=c("First","First-only","PIFU","Traditional"))

animate_process(event_DES,
                mode = "absolute",
                epsilon_time=0.01,
                jitter=0,
                duration=15,
                legend="color",
                initial_time=0,
                mapping = token_aes(color = token_scale("typef", 
                                                        scale = "ordinal", 
                                                        range = RColorBrewer::brewer.pal(4, "Set1")[c(1,2,4)])))

saveviewer <- function(aux,mymetric){
  saveWidget(aux,"temp.html",selfcontained=FALSE)
  webshot("temp.html",file=here::here(bundle,case,paste0("map-",mymetric,"viz.png")))
}



### Process maps - static


det=5
aa <- event_DES #%>%
  #filter_trim(start_activities = "New", end_activities =  c("New",paste0("FUp ",1:det)))

mymetric = "median"
aux <- aa %>% process_map(performance(median, "days")) ; aux
#htmlwidgets::saveWidget(aux, file = "./Process mapping/outputs/apr17_1000_perf_median_.html")
saveviewer(aux,mymetric)


mymetric="absolute"
aux <- aa %>% process_map();aux
#saveviewer(aux,mycohort,mymetric,det)

mymetric="relative"
aux<- aa %>% process_map(type = frequency(mymetric));aux
#saveviewer(aux,mycohort,mymetric,det)

mymetric="relative_case"
aux <- aa %>% process_map(type = frequency(mymetric));aux
#saveviewer(aux,mycohort,mymetric,det)


### Animated

saveviewerhtml <- function(aux,mymetric="absolute"){
  
  saveWidget(aux,file=here::here(bundle,case,paste0("map-",mymetric,"rep",myrep,"id",myid,"_",Sys.Date(),".html")),selfcontained=FALSE)
  
}

saveviewerhtml_self <- function(aux,mymetric="absolute"){
  
  saveWidget(aux,file=here::here(bundle,case,paste0("map-",mymetric,"rep",myrep,"id",myid,"_",Sys.Date(),"_self.html")),selfcontained=TRUE)
  
}


ap_aa <- animate_process(event_DES,
                         mode = "absolute",
                         epsilon_time=0.01,
                         jitter=0,
                         duration=15,
                         legend="color",
                         initial_time=0,
                         mapping = token_aes(color = token_scale("type", 
                                                                 scale = "ordinal", 
                                                                 range = RColorBrewer::brewer.pal(4, "Set1")[c(1,2,4)])))
ap_aa

saveviewerhtml(ap_aa)
