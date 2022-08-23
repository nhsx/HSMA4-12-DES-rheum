# -*- coding: utf-8 -*-
"""
Created on Mon May  9 16:50:31 2022

@author: martina.fonseca_nhsx
"""

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import simpy_rheum_v0031 as rheum
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import datetime
import os
import csv
import random

random.seed(9001)

scriptrun_flag = True
batch_flag = True
reps=30
#savepath = 'out_sand/'
#savepath = 'bundleV2/out_central_15_pifu20_interpifu20/'
#savepath = 'bundleV2/out_central_3_base/'
#savepath = 'bundleV2/out_10years_rep3_base_capplus1_distexp/'
#savepath = 'bundleV2/out_10years_rep1_base_capplus5_12mshock100/'

#V1 vs V3 . V3 with 5+3 instead of 3+4 warm+obs. PIFU starts being offered at end of wamr-up to 'eligible at stock fup' patients not new 'stock new' patients.
#savepath = 'bundleV3/out_central_rep30_pifu20_interpifu20/'
savepath = 'bundleV3_/out_central_rep30_pifu10/'
start=datetime.now()
file1 = savepath + 'patient_result2.csv'
file2 = savepath + 'appt_result.csv'
file3 = savepath + 'batch_mon_audit_ls.csv'
#file4 = savepath + 'batch_kpis.csv'


# Check whether the specified path exists or not
isExist = os.path.exists(savepath)
if not isExist:
  # Create a new directory because it does not exist 
  os.makedirs(savepath)
  print("The new directory is created!")

rheum.Trial_Results_initiate(file1,file2,file3)

# Create a file to store trial results, and write the column headers
with open(savepath + "trial_results.csv", "w") as f:
    writer = csv.writer(f, delimiter=",")
    column_headers = ["Run", "Mean_Q_Time_FOPA",
                      "Mean_Q_Time_FUOPA",
                      "Mean_Q_Time_Total"]
    writer.writerow(column_headers)


if scriptrun_flag:
    if batch_flag:
        
        
        intarr = 1/6 # 1/6
        in_prob_pifu = 0.1 # 0
        in_path_horizon_y=3 #3
        #in_FOavoidable = 0.43
        in_FOavoidable = 0
        in_interfu_perc = 0.6 # base should be 0.6!
        #in_interfu_perc = 0.2 # base should be 0.6!
        g_defaults = rheum.g()
        audit_interval = 28
        #cap = 1/intarr * ((1+ 3*12/4.5)*0.6 + 1*0.4)
        #cap  = 1/intarr * ((1 + in_path_horizon_y / g_defaults.mean_interOPA *365) * (1-g_defaults.prob_firstonly) + 1 * g_defaults.prob_firstonly)# * 7/5
        cap  = 1/intarr * ((2 + in_path_horizon_y / g_defaults.mean_interOPA *365) * (1-g_defaults.prob_firstonly) + 2 * g_defaults.prob_firstonly)# * 7/5
        cap_diff = 0
        #cap = 22
        #cap_diff=0
        
        my_batch_model = rheum.Batch_rheum_model(in_res=np.round(cap,0) + cap_diff,
                                                 in_inter_arrival=intarr,
                                                 in_prob_pifu=in_prob_pifu,
                                                 in_path_horizon_y=in_path_horizon_y,
                                                 audit_interval=audit_interval,
                                                 in_savepath = savepath,
                                                 in_FOavoidable = in_FOavoidable,
                                                 in_interfu_perc=in_interfu_perc)
        
        fig_audit_reps, chart_output_lastrep, text_output_lastrep, quant_output_lastrep, fig_q_audit_reps,fig_monappKPI_reps, fig_monappKPIn_reps = my_batch_model.run_reps(reps=reps)
        
        
        my_batch_model.headline_KPI(365)
        
        mydf = my_batch_model.headline_KPI(365)
        
        my_batch_model.save_logs()
        
       # my_batch_model.batch_kpi.to_csv(savepath + 'batch_kpi.csv')
        
        print(my_batch_model.batch_kpi)
        
        #my_batch_model.batch_mon_appointments.to_csv(savepath + 'batch_mon_appointments.csv')
        #my_batch_model.batch_mon_audit.to_csv(savepath + 'batch_mon_audit.csv')
        #my_batch_model.plot_monappKPI_reps(step=28*2)
        
        
# =============================================================================
#         ## Create mean, quantile KPIs per replication, chunked by intervals.
#         batch_mon_appointments = my_batch_model.batch_mon_appointments
#         
#         step = 28
#         batch_mon_appointments['interval'] = (batch_mon_appointments['start_q']//step)*step
#         
#         batch_mon_app_kpit = batch_mon_appointments.groupby(['rep','priority','interval']).quantile([0.25,0.50,0.75]).reset_index()
#         
#         batch_mon_app_kpitmu =batch_mon_appointments.copy()
#         batch_mon_app_kpitmu['level_3']='mean'
#         batch_mon_app_kpitmu = batch_mon_appointments.groupby(['rep','priority','interval','level_3']).mean().reset_index()
#                 
#         batch_mon_app_kpit = pd.concat([batch_mon_app_kpit,batch_mon_app_kpitmu])
#         
#         
#         batch_mon_app_kpit = batch_mon_app_kpit.rename(columns={'level_3':'KPI'})
#         
#         batch_mon_app_kpiseen = batch_mon_appointments.groupby(['rep','priority','interval']).size().reset_index(name='q_time')
#         batch_mon_app_kpiseen['KPI']='Seen'
#         
#         batch_mon_app_kpit = pd.concat([batch_mon_app_kpit,batch_mon_app_kpiseen])
#         
#         ## Create boxplot
#         #fig_qtime = plt.figure(figsize=(12,12))
#         KPI_types = batch_mon_app_kpit['KPI'].unique()[:-1]
#         i=0
#         #axes = fig_qtime.add_subplot(1,len(KPI_types),1)
#         fig, axes = plt.subplots(len(KPI_types),1, figsize=(20, 20), sharey=True)
#         fig.suptitle('Queueing time KPIs per replication/interval')
#         for KPI_t in KPI_types:
#             
#             ax=sns.boxplot(ax=axes[i],x='interval',y='q_time',data=batch_mon_app_kpit[batch_mon_app_kpit['KPI']==KPI_t],hue='priority')
#             ax.set_xticklabels(ax.get_xticklabels(),rotation = 80)
#             ax.set_xlabel("Simulation time interval (days)", fontsize = 20)
#             ax.set_ylabel("KPI: " + str(KPI_t), fontsize = 20)   
#             #ax.set_title("Queueing time")
#             i+=1
#         
#         
#         fig2, axes2 = plt.subplots(figsize=(20, 20))
#         ax2=sns.boxplot(ax=axes2,x='interval',y='q_time',data=batch_mon_app_kpit[batch_mon_app_kpit['KPI']=='Seen'],hue='priority')
#         ax2.set_xticklabels(ax2.get_xticklabels(),rotation = 80)
#         ax2.set_xlabel("Simulation time interval (days)", fontsize = 14)
#         ax2.set_ylabel("Number seen", fontsize = 20)   
# =============================================================================
        
        
# =============================================================================
#         batch_mon_appointments = pd.DataFrame()
#         batch_mon_audit = pd.DataFrame()
#         
#         for run in range(reps):
#             print (f"Run {run+1} of {reps}")
#             
#             intarr = 1
#             cap = 1/intarr * (1+ 3*12/4.5)
#             
#             my_ed_model = rheum.rheum_Model(run,
#                                             in_res=np.round(cap,0),
#                                             in_inter_arrival=intarr,
#                                             in_prob_pifu=0.5,
#                                             in_path_horizon_y=3,
#                                             audit_interval=28
#                                         ) # create instance of rheumatology model (constructor init)
#             my_ed_model.run() # apply custom method run to instance
#             print ()
#             
#             p = rheum.rheum_Model.chart(my_ed_model)
# 
#             e_results = my_ed_model.g.results # audit counts . patients waiting per time
#             e_results["rep"] = run
# 
#             e_resultsdf = my_ed_model.results_df # initial df with time each patient waited for first OPA
# 
#             e_appt_queuing_result= my_ed_model.g.appt_queuing_results # appointments monitor
#             e_appt_queuing_result["rep"] = run
# 
#             batch_mon_appointments = pd.concat([batch_mon_appointments,e_appt_queuing_result])
#             batch_mon_audit = pd.concat([batch_mon_audit, e_results])
#             
#             t_warm = my_ed_model.g.warm_duration
#             
#             sns.boxplot(x='type',y='q_time',data=batch_mon_appointments[batch_mon_appointments['start_q']>t_warm],hue='type')
# 
#             audit_kpis = ['priority 1 patients waiting','priority 2 patients waiting','priority 3 patients waiting','resources occupied']
#             batch_mon_audit_melted = pd.melt(batch_mon_audit,id_vars=['time','rep'],value_vars=audit_kpis,var_name='audit_KPI')
# 
#             sns.boxplot(x='time',y='priority 1 patients waiting',data=batch_mon_audit)
#             
#             fig, axes = plt.subplots(len(audit_kpis),1, figsize=(20, 20), sharey=True)
#             fig.suptitle('Audit point KPIs')
#             
#             i=0
#             for kpi in audit_kpis:
#                 ax=sns.boxplot(ax=axes[i],x='time',y='value',data=batch_mon_audit_melted[batch_mon_audit_melted['audit_KPI']==kpi])
#                 ax.set_xticklabels(ax.get_xticklabels(),rotation = 30)
#                 if i==len(audit_kpis)-1:
#                     ax.set_xlabel("Audit timepoint (days)", fontsize = 20)
#                 ax.set_ylabel(kpi, fontsize = 20)
#                 i+=1
#         
# =============================================================================
        
        # Once the trial is complete, we'll create an instance of the
        # Trial_Result_Calculator class and run the print_trial_results method
        my_trial_results_calculator = rheum.Trial_Results_Calculator(savepath)
        my_trial_results_calculator.print_trial_results(savepath)
    
    else:
        #my_ed_model = rheum.rheum_Model(0,32,(365/4590),0.5)
        
        intarr = 1
        cap = 1/intarr * (1+ 3*12/4.5)
        
        my_ed_model = rheum.rheum_Model(0,
                                        in_res=np.round(cap,0),
                                        in_inter_arrival=intarr,
                                        in_prob_pifu=0.5,
                                        in_path_horizon_y=3,
                                        audit_interval=7
                                        )
        my_ed_model.run()
        
scenario2run = datetime.now()-start
print(f"Overall run-time of {scenario2run}")        
     

    
a = pd.DataFrame({"P_ID":[99],
                                    "App_ID":[99],
                                    "priority":[99],
                                    "type":[99],
                                    "q_time":[99],
                                    "start_q":[99],
                                    "DNA":[99]})
