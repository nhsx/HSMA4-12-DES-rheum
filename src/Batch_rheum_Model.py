#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# v0.01
# Last changed: 25/02/2022
# Pathway: rheumatology
# refactor with day steps and slots per day - try something more like https://qualitysafety.bmj.com/content/qhc/23/5/373.full.pdf

# Time unit: day

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import simpy
import random
import pandas as pd
import numpy as np
import csv
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import scipy.stats as st

from src.helpers import mean_confidence_interval
from src.initialisers import g
#from src.patient import FOPA_Patient
from src.rheum_Model import rheum_Model

       

# Class for Batch runs / replications of the model
class Batch_rheum_model:
    
    # Initialise Class for Batch run model. Instantiate g.
    def __init__(self,in_res=5 , in_inter_arrival=1, in_prob_pifu=0.6, in_path_horizon_y=3,audit_interval=7,in_savepath="temp/",in_FOavoidable=0,in_interfu_perc=0.6):
        
        self.batch_mon_appointments = pd.DataFrame()
        self.batch_mon_audit = pd.DataFrame()
        self.batch_mon_app_kpit = pd.DataFrame()
        self.batch_kpi = pd.DataFrame()
        self.savepath=in_savepath
        self.g = g(in_res,in_inter_arrival,in_prob_pifu, in_path_horizon_y, audit_interval,in_FOavoidable=in_FOavoidable,in_interfu_perc=in_interfu_perc) # instance of global variables
        
    
    # Method to read csv results to constructor instance dataframe
    def read_logs_to_self(self):        
        # Read in results back into self dataframe
        self.trial_results_df = pd.read_csv(self.savepath +"patient_result2.csv")
        self.batch_mon_appointments = pd.read_csv(self.savepath +"appt_result.csv")
        self.batch_mon_audit = pd.read_csv(self.savepath + "batch_mon_audit_ls.csv")
    
    
    # Plotting an overview of behaviour at audit timepoints (across reps)    
    def plot_audit_reps(self):
        
        t_warm = self.g.warm_duration
        
        fig_q = plt.figure(figsize=(12,12))
        sns.violinplot(x='type',y='q_time',data=self.batch_mon_appointments[self.batch_mon_appointments['start_q']>t_warm],hue='type')


        # Other plot
        #fig_q2 = plt.figure(figsize=(12,12))
        #sns.boxplot(x='time',y='priority 1 patients waiting',data=self.batch_mon_audit)  

        # Audit plot
        audit_kpis = ['priority 3 patients waiting','priority 1 patients waiting','priority 2 patients waiting','resources occupied']
        pricolors = ['orange','skyblue','green','gray']
        audit_kpis_name = ['RTT appointment waits','Traditional appointment waits','PIFU appointment waits','Resources (slots) occupied']
        batch_mon_audit_melted = pd.melt(self.batch_mon_audit,id_vars=['time','rep'],value_vars=audit_kpis,var_name='audit_KPI')
              
        fig, axes = plt.subplots(len(audit_kpis),1, figsize=(20, 20), sharey=False)
        #fig.suptitle('Audit point KPIs')        
        i=0
        for kpi in audit_kpis:
            ax=sns.boxplot(ax=axes[i],x='time',y='value',data=batch_mon_audit_melted[batch_mon_audit_melted['audit_KPI']==kpi],color=pricolors[i])
            ax.set_xticklabels(ax.get_xticklabels(),rotation = 30)
            if i==len(audit_kpis)-1:
                ax.set_xlabel("Audit timepoint (days)", fontsize = 20)
            ax.set_ylabel(audit_kpis_name[i], fontsize = 20)
            i+=1
            
        return fig, fig_q
    
    # Computing headline/core KPIs on queuing time, resources and waiting list size (batch / inter-replication)
    def headline_KPI(self,window_tail=365):
        

        max_time = self.g.warm_duration +  self.g.obs_duration
        batch_mon_appointments = self.batch_mon_appointments.copy()
        batch_mon_appointments['end_q']=batch_mon_appointments['start_q']+batch_mon_appointments['q_time']
        batch_mon_appointments = batch_mon_appointments[batch_mon_appointments['start_q'] > max_time - window_tail]
        batch_KPIs = batch_mon_appointments[batch_mon_appointments['priority']==3].copy() # only first referrals
        
        batch_KPIs_q = batch_KPIs[['rep','q_time']].groupby(['rep']).quantile([0.50,0.75,0.92,0.95]).reset_index()
        batch_KPIs_q=batch_KPIs_q.rename(columns={'level_1':'KPI','q_time':'value'})
        batch_KPIs_q['KPI'] = 'RTT_q'+ batch_KPIs_q['KPI'].astype(str)
        
        batch_KPIs_mu=  batch_KPIs[['rep','q_time']].groupby(['rep']).mean().reset_index()
        batch_KPIs_mu=batch_KPIs_mu.rename(columns={'q_time':'value'})
        batch_KPIs_mu['KPI']='RTT_mean'
        
        batch_KPIs_n=  batch_KPIs[['rep','q_time']].groupby(['rep']).size().reset_index()
        batch_KPIs_n=batch_KPIs_n.rename(columns={0:'value'})
        batch_KPIs_n['KPI']='RTT_n'
        
        
        batch_mon_audit = self.batch_mon_audit.copy()
        batch_KPIs_wl = batch_mon_audit[batch_mon_audit['time']==max(batch_mon_audit['time'])].rename(columns={'priority 3 patients waiting':'RTT_WL_end'})
        batch_KPIs_wl= batch_KPIs_wl[['rep','RTT_WL_end','resources occupied']]
        batch_KPIs_wl = pd.melt(batch_KPIs_wl,id_vars=['rep'],var_name='KPI')
        
        
        batch_kpi_rep = pd.concat([batch_KPIs_q,batch_KPIs_mu,batch_KPIs_n,batch_KPIs_wl])
        
        a = batch_kpi_rep.groupby(['KPI']).agg({'value': lambda x: mean_confidence_interval(x)[0]}).rename(columns={'value':'KPI_mean'})
        b = batch_kpi_rep.groupby(['KPI']).agg({'value': lambda x: mean_confidence_interval(x)[1]}).rename(columns={'value':'KPI_LCI'})
        c = batch_kpi_rep.groupby(['KPI']).agg({'value': lambda x: mean_confidence_interval(x)[2]}).rename(columns={'value':'KPI_UCI'})
        
        batch_kpi = pd.concat([a,b,c],axis=1)
        
        self.batch_kpi = batch_kpi
        
        return batch_kpi

                
    
    # Plotting an overview of queueing time appointment KPI behaviour (across reps and by time intervals)
    def plot_monappKPI_reps(self,step=365/4):
        
        batch_mon_appointments = self.batch_mon_appointments.copy()

        batch_mon_appointments['end_q'] = batch_mon_appointments['start_q']+batch_mon_appointments['q_time']
        
        batch_mon_appointments['interval'] = (batch_mon_appointments['end_q']//step)*step
        
        batch_mon_appointments = batch_mon_appointments[['rep','priority','interval','q_time']] # mfadd 11/12/2023

        batch_mon_app_kpit = batch_mon_appointments.groupby(['rep','priority','interval']).quantile([0.50],numeric_only=True).reset_index()
        
        batch_mon_app_kpitmu =batch_mon_appointments
        batch_mon_app_kpitmu['level_3']='mean'
        batch_mon_app_kpitmu = batch_mon_appointments.groupby(['rep','priority','interval','level_3']).mean().reset_index()
                
        batch_mon_app_kpit = pd.concat([batch_mon_app_kpit,batch_mon_app_kpitmu])
        
        batch_mon_app_kpit = batch_mon_app_kpit.rename(columns={'level_3':'KPI'})
        
        batch_mon_app_kpiseen = batch_mon_appointments.groupby(['rep','priority','interval']).size().reset_index(name='q_time')
        batch_mon_app_kpiseen['KPI']='Seen'
        batch_mon_app_kpit = pd.concat([batch_mon_app_kpit,batch_mon_app_kpiseen])
        
        #batch_mon_app_kpit.to_csv('batch_mon_app_kpit.csv')
        
        self.batch_mon_app_kpit = batch_mon_app_kpit
        
        ## Create boxplot
        #fig_qtime = plt.figure(figsize=(12,12))
        KPI_types = batch_mon_app_kpit['KPI'].unique()[:-1]
        i=0
        #axes = fig_qtime.add_subplot(1,len(KPI_types),1)
        fig, axes = plt.subplots(len(KPI_types),1, figsize=(20, 20), sharey=True)
        fig.suptitle('RTT queueing time (days)',fontsize=20)
       
        
        for KPI_t in KPI_types:
            
            data_now = batch_mon_app_kpit[batch_mon_app_kpit['priority']==3]
            #ax=sns.boxplot(ax=axes[i],x='interval',y='q_time',data=batch_mon_app_kpit[batch_mon_app_kpit['KPI']==KPI_t],hue='priority')
           # ax=sns.boxplot(ax=axes[i],x='interval',y='q_time',data=batch_mon_app_kpit[(batch_mon_app_kpit['KPI']==KPI_t) & (batch_mon_app_kpit['priority']==3)])
            ax=sns.boxplot(ax=axes[i],x='interval',y='q_time',data=data_now[data_now['KPI']==KPI_t],hue='priority') 
            ax.set_xticklabels(ax.get_xticklabels(),rotation = 80)
            ax.set_xlabel("Simulation time interval (days)", fontsize = 20)
            ax.set_ylabel("KPI" + str(KPI_t), fontsize = 20)   
            #ax.set_title("Queueing time")
            i+=1
        
        i=0
        prilist = [3,1,2]
        prinames = ['RTT seen','Traditional F/U seen','PIFU seen']
        pricolors = ['orange','skyblue','green']
        fig2, axes2 = plt.subplots(3,1, figsize=(20, 20), sharey=True)
        fig2.suptitle('Number seen per quarter',fontsize=20)
        for pri in prilist:
            data_now = batch_mon_app_kpit[batch_mon_app_kpit['KPI']=='Seen'].copy()
            ax2=sns.lineplot(ax=axes2[i],x='interval',y='q_time',data=data_now[data_now['priority']==pri],color=pricolors[i],
                             markers=True,marker='o',markersize=16)
            #ax2.set_xticklabels(ax2.get_xticklabels(),rotation = 80)
            if i==2:
                ax2.set_xlabel("Simulation time interval (days)", fontsize = 20)
            ax2.set_ylabel(prinames[i], fontsize = 20) 
            i=i+1
     
        return fig, fig2
        
    
    # Method to run replications. Calls run method
    def run_reps(self,reps):
        
        for run in range(reps):
            
            if self.g.debug  and self.g.debuglevel>=0:
                print (f"Run {run+1} of {reps}")
            
            # Instance of rheumatology model
            my_ed_model = rheum_Model(run,
                                      in_res=self.g.number_of_slots,
                                      in_inter_arrival=self.g.wl_inter,
                                      in_prob_pifu=self.g.prob_pifu,
                                      in_path_horizon_y=self.g.max_fuopa_tenor_y,
                                      audit_interval=self.g.audit_interval,
                                      repid = run,
                                      savepath = self.savepath,
                                      in_FOavoidable = self.g.in_FOavoidable,
                                      in_interfu_perc = self.g.interfu_perc) # create instance of rheumatology model (constructor init)
            
            
            start=datetime.now()
            if run<reps-1:
                my_ed_model.run() # apply custom method run to instance
                
            else:
                chart_output_lastrep, text_output_lastrep, quant_output_lastrep = my_ed_model.run()
            if self.g.debug  and self.g.debuglevel>=0:
                    scenario1run = datetime.now()-start
                    print(f"Run-time of Run {run+1}: {scenario1run}")
                
            print ()
            
            # Load up audit results of replication
            e_results = my_ed_model.g.results # audit counts . patients waiting per time
            self.batch_mon_audit = pd.concat([self.batch_mon_audit, e_results]) # Replication audit result added to Batch audit result df

            # Load up appointment log results of replication
            if self.g.loglinesave:
                self.read_logs_to_self() # read from file
            else:
                e_appt_queuing_result= my_ed_model.g.appt_queuing_results # replication appointments monitor in memory
                self.batch_mon_appointments = pd.concat([self.batch_mon_appointments,e_appt_queuing_result]) # append for batch
             
            # keep only post warm-up queue starts (rather q finish????)
            self.batch_mon_appointments = self.batch_mon_appointments[self.batch_mon_appointments['start_q']>self.g.warm_duration]
            
        ### Batch summaries (plots, KPIs...)    
        fig_audit_reps, fig_q_audit_reps = self.plot_audit_reps() # generate audit plots
        
        fig_monappKPI_reps, fig_monappKPIn_reps = self.plot_monappKPI_reps(step=365/4) # generate KPI over time plots
        
        self.headline_KPI(365) # generate core/headline KPIs
        
        return fig_audit_reps, chart_output_lastrep, text_output_lastrep, quant_output_lastrep, fig_q_audit_reps, fig_monappKPI_reps, fig_monappKPIn_reps
        
    # Save aggregate logs (cross-replication)    
    def save_logs(self):
        
       self.batch_mon_appointments.to_csv(self.savepath + 'batch_mon_appointments.csv')
       self.batch_mon_audit.to_csv(self.savepath + 'batch_mon_audit.csv')
       self.batch_mon_app_kpit.to_csv(self.savepath + 'batch_mon_app_kpit.csv')
       self.batch_kpi.to_csv(self.savepath + 'batch_kpi.csv')
        #my_batch_model.plot_monappKPI_reps(step=28*2)


            
