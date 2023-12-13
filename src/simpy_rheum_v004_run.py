""" Python script (.py) to run the experiment. Similar can be achieved with the .ipynb file. """
import os
import csv
import random
import warnings
from datetime import datetime
import numpy as np
warnings.simplefilter(action='ignore', category=FutureWarning)


os.chdir('../') ## go up one dir
import src.Batch_rheum_Model as rheum ##
from src.helpers import Trial_Results_initiate ##
from src.initialisers import g ##

random.seed(9001)

scriptrun_flag = True # True to save each log line by line (more efficient)
reps=30 # Number of model replications | Baseline: 30 replications
outputdir = 'outputs/'
savepath = 'out_sand/'
savepath = outputdir + savepath
start=datetime.now()
file1 = savepath + 'patient_result2.csv'
file2 = savepath + 'appt_result.csv'
file3 = savepath + 'batch_mon_audit_ls.csv'
#file4 = savepath + 'batch_kpis.csv'

# Check whether the output directory exists or not ##
isExist = os.path.exists(outputdir)
if not isExist:
    # Create a new directory because it does not exist
    os.makedirs(outputdir)
    print("The new directory is created!")

# Check whether the specified path exists or not
isExist = os.path.exists(savepath)
if not isExist:
    # Create a new directory because it does not exist
    os.makedirs(savepath)
    print("The new directory is created!")

Trial_Results_initiate(file1,file2,file3)

# Create a file to store trial results, and write the column headers
with open(savepath + "trial_results.csv", "w",encoding="cp1252") as f:
    writer = csv.writer(f, delimiter=",")
    column_headers = ["Run", "Mean_Q_Time_FOPA",
                      "Mean_Q_Time_FUOPA",
                      "Mean_Q_Time_Total"]
    writer.writerow(column_headers)


if scriptrun_flag:

    intarr = 1/6 # Inter-arrival rate (in days) | Baseline: 1/6 days, i.e. 6 per day
    in_path_horizon_y=3 # Patient follow-up horizon, simplification on how long each non first-only pathway lasts (years) | Baseline: 3

    #in_prob_pifu = 0
    in_prob_pifu = 0.2 # PIFU proportion - probability of PIFU pathway for non first-only pathways | Baseline: 0 | Scenarios: 0.1,0.2,0.3,0.5

    in_FOavoidable = 0 # A&G proportion - proportion of first-only pathways avoidable via A&G | Baseline: 0 | Scenario: 0.43
    #in_FOavoidable = 0.43

    in_interfu_perc = 0.6 # Percentage increase in inter-appointment interval with PIFU (vs traditional), i.e. 0.6 means 60% longer interval | Baseline: 0.6 | Scenarios: 0.6, 0.2
    #in_interfu_perc = 0.2 # More conservative. Used in scenario C

    g_defaults = g()
    audit_interval = 28 # audit timepoint (in simulation days)

    cap  = 1/intarr * ((2 + in_path_horizon_y / g_defaults.mean_interOPA *365) * (1-g_defaults.prob_firstonly) + 2 * g_defaults.prob_firstonly)# * Heuristic of daily slots (365 days) needed to deal with steady-state model demand
    cap_diff = 0 # Increment or decrement in daily slots to be applied to the heuristic | Baseline: 0

    # Set up constructor instance
    my_batch_model = rheum.Batch_rheum_model(in_res=np.round(cap,0) + cap_diff,
                                             in_inter_arrival=intarr,
                                             in_prob_pifu=in_prob_pifu,
                                             in_path_horizon_y=in_path_horizon_y,
                                             audit_interval=audit_interval,
                                             in_savepath = savepath,
                                             in_FOavoidable = in_FOavoidable,
                                             in_interfu_perc=in_interfu_perc)

    # Run model
    fig_audit_reps, chart_output_lastrep, text_output_lastrep, quant_output_lastrep, fig_q_audit_reps,fig_monappKPI_reps, fig_monappKPIn_reps = my_batch_model.run_reps(reps=reps)

#
# =============================================================================

    # Some output KPIs generated and printed
    #my_batch_model.headline_KPI(365)
    #mydf = my_batch_model.headline_KPI(365)
    #my_batch_model.save_logs()
    print(my_batch_model.batch_kpi)

# Model performance - time
scenario2run = datetime.now()-start
print(f"Overall run-time of {scenario2run}")
     