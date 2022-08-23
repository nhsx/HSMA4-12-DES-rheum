# from numpy.lib.arraysetops import ediff1d
import streamlit as st
import numpy as np
import simpy_rheum_v0031 as rheum
from datetime import datetime
import os
import csv
import random

st.write('| Toy tool of backlog rheumatology outpatient Discrete Event Simulation Model. The effect of Patient Initiated Follow-up (PIFU) and Advice & Guidance (A&G) can be simulated. The runs may take 5-10 minutes, 3 simulation replications are used so caution is needed - more are used in report examples.')


st.title('Rheumatology PIFU Queueing Simulation - main scenario')

savepath = 'out_streamlit_S1/'

# Check whether the specified path exists or not
isExist = os.path.exists(savepath)
if not isExist:
  # Create a new directory because it does not exist 
  os.makedirs(savepath)
  print("The new directory is created!")
file1 = savepath + 'patient_result2.csv'
file2 = savepath + 'appt_result.csv'
file3 = savepath + 'batch_mon_audit_ls.csv'

rheum.Trial_Results_initiate(file1,file2,file3)

# Create a file to store trial results, and write the column headers
with open(savepath + "trial_results.csv", "w") as f:
    writer = csv.writer(f, delimiter=",")
    column_headers = ["Run", "Mean_Q_Time_FOPA",
                      "Mean_Q_Time_FUOPA",
                      "Mean_Q_Time_Total"]
    writer.writerow(column_headers)

#col1, col2, col3, col4, col5 = st.columns(3)
col1, col2, col3,col4 = st.columns(4)
in_daily_arrivals = col1.slider(
    'Average daily referrals',1.0, 6.0, 6.0,step=0.1)
in_prob_pifu = col2.slider('PIFU - proportion that goes to PIFU (%) *', 0,100,0,step=1)/100
in_interfu_perc = col3.slider('PIFU - increase in inter-appointment interval (%)', 0.0,100.0,60.0,step=1.0)/100
in_FOavoidable = col4.slider('A&G - first-only pathways avoidable (%)', 0,100,0,step=1)/100

g_defaults = rheum.g()

in_path_horizon_y = 3
audit_interval = 28*2
#in_FOavoidable = 0
#in_interfu_perc = 0.6
cap  = np.round(in_daily_arrivals * ((2 + in_path_horizon_y / g_defaults.mean_interOPA *365) * (1-g_defaults.prob_firstonly) + 2 * g_defaults.prob_firstonly))# * Heuristic of daily slots (365 days) needed to deal with steady-state model demand

st.write(f"Capacity: without interventions, {np.round(cap,0)} is the approx heuristic number of daily slots that ensures demand-supply balance, given {in_daily_arrivals} daily arrivals, a {g_defaults.prob_firstonly*100}% probability of being a first-only pathway (internal assumption), a {in_path_horizon_y} year follow-up horizon for non first-only pathways (internal assumption) and a mean {g_defaults.mean_interOPA} day traditional inter-appointment interval (internal assumption).")

col5, col6, col7,col8 = st.columns(4)
cap_diff = col5.slider('Adjustment in daily slots', -2,2,0,step=1)
in_res = cap + cap_diff

rheum_model = rheum.Batch_rheum_model(in_res= in_res,
                                #in_inter_arrival=in_inter_arrival / (24*60),
                                in_inter_arrival = 1/in_daily_arrivals,
                                in_prob_pifu=in_prob_pifu,
                                in_path_horizon_y=in_path_horizon_y,
                                audit_interval = audit_interval,
                                in_savepath = savepath,
                                in_FOavoidable = in_FOavoidable,
                                in_interfu_perc=in_interfu_perc
                                )



st.title('Rheumatology PIFU Queueing Simulation - scenario 2')

savepath_S2 = 'out_streamlit_S2/'

# Check whether the specified path exists or not
isExist = os.path.exists(savepath_S2)
if not isExist:
  # Create a new directory because it does not exist 
  os.makedirs(savepath_S2)
  print("The new directory is created!")
file1 = savepath_S2 + 'patient_result2.csv'
file2 = savepath_S2 + 'appt_result.csv'
file3 = savepath_S2 + 'batch_mon_audit_ls.csv'

rheum.Trial_Results_initiate(file1,file2,file3)


#colS2_1, colS2_2, colS2_3, colS2_4, colS2_5 = st.columns(3)
colS2_1, colS2_2, colS2_3, colS2_4 = st.columns(4)
in_daily_arrivals_S2 = colS2_1.slider(
    'Average daily referrals.',1.0, 6.0, 6.0,step=0.1)
in_prob_pifu_S2 = colS2_2.slider('PIFU - proportion that goes to PIFU (%) *.', 0,100,50,step=1)/100
in_interfu_perc_S2 = colS2_3.slider('PIFU - increase in inter-appointment interval (%).', 0,100,60,step=1)/100
in_FOavoidable_S2 = colS2_4.slider('A&G - first-only pathways avoidable (%).', 0,100,0,step=1)/100

in_path_horizon_y_S2 = 3
audit_interval_S2 = 28*2

cap_S2  = np.round(in_daily_arrivals_S2 * ((2 + in_path_horizon_y_S2 / g_defaults.mean_interOPA *365) * (1-g_defaults.prob_firstonly) + 2 * g_defaults.prob_firstonly))# * Heuristic of daily slots (365 days) needed to deal with steady-state model demand

st.write(f"Capacity: without interventions, {np.round(cap_S2,0)} is the approx heuristic number of daily slots that ensures demand-supply balance, given {in_daily_arrivals_S2} daily arrivals, a {g_defaults.prob_firstonly*100}% probability of being a first-only pathway (internal assumption), a {in_path_horizon_y_S2} year follow-up horizon for non first-only pathways (internal assumption) and a mean {g_defaults.mean_interOPA} day traditional inter-appointment interval (internal assumption).")

colS2_5, colS2_6, colS2_7,colS2_8 = st.columns(4)
cap_diff_S2 = colS2_5.slider('Adjustment in daily slots.',-2,2,0,step=1)
in_res_S2 = cap_S2 + cap_diff_S2

st.write(' >[*] among those non first-only pathways, i.e. with follow-up')

rheum_model_S2 = rheum.Batch_rheum_model(in_res = in_res_S2,
                                   #in_inter_arrival =in_inter_arrival_S2/ (24*60),
                                   in_inter_arrival = 1/in_daily_arrivals_S2,
                                   in_prob_pifu=in_prob_pifu_S2,
                                   in_path_horizon_y=in_path_horizon_y_S2,
                                   audit_interval = audit_interval_S2,
                                   in_savepath = savepath_S2,
                                   in_FOavoidable = in_FOavoidable_S2,
                                   in_interfu_perc=in_interfu_perc_S2
                                   )


st.title('Runs')


with st.sidebar:
    main_button = st.button('Run rheum model - run main scenario',key="1")
    tworuns_button = st.button('Run rheum model - run both',key="2")


if main_button:
    # Get results
    start=datetime.now()
    random.seed(9001)
    fig_audit_reps, chart, text , quant, fig_q_audit_reps,fig_monappKPI_reps,fig_monappKPIn_reps = rheum_model.run_reps(2)
    scenario1run = datetime.now()-start
    print(f"Run-time of {scenario1run}")
    # Show chart
    st.write('---')
    st.write('Core KPIs')
    #st.write(rheum_model.batch_kpi.iloc[[0,1,2,3,7],])
    st.write(rheum_model.batch_kpi)
    st.write('---')
    st.write('Waiting time KPIs')
    st.pyplot(fig_monappKPI_reps)
    st.write('---')
    st.write('Audit KPIs')
    st.pyplot(fig_audit_reps)
    st.write('---')
    st.write('Patients seen KPIs')
    st.pyplot(fig_monappKPIn_reps)

    #st.write('---')
    #st.write('Single replication example:')
    # Show chart
    #st.pyplot(chart)
    # Display results (text)
    # for t in text:
    #     st.write(t)
    #st.write(quant)
    st.write(scenario1run)


if tworuns_button:
    
    st.subheader("Main scenario")
    # Get results
    start=datetime.now()
    random.seed(9001)
    fig_audit_reps, chart, text , quant, fig_q_audit_reps,fig_monappKPI_reps, fig_monappKPIn_reps = rheum_model.run_reps(2)
    scenario1run = datetime.now()-start
    print(f"Run-time of {scenario1run}")
    # Show chart
    st.write('---')
    st.write('Core KPIs')
    #st.write(rheum_model.batch_kpi.iloc[[0,1,2,3,7],])
    st.write(rheum_model.batch_kpi)
    st.write('---')
    st.write('Waiting time KPIs')
    st.pyplot(fig_monappKPI_reps)
    st.write('---')
    st.write('Audit KPIs')
    st.pyplot(fig_audit_reps)
    st.write('---')
    st.write('Patients seen KPIs')
    st.pyplot(fig_monappKPIn_reps)
    #st.write('---')
    #st.write('Single replication example:')
    # Show chart
    #st.pyplot(chart)
    # Display results (text)
    # for t in text:
    #     st.write(t)
    #st.write(quant)
    st.write(scenario1run)
    
    st.subheader("Scenario 2")
    # Get results
    start=datetime.now()
    fig_audit_reps_S2, chart_S2, text_S2 , quant_S2, fig_q_audit_reps_S2,fig_monappKPI_reps_S2, fig_monappKPIn_reps_S2 = rheum_model_S2.run_reps(2)
    scenario2run = datetime.now()-start
    print(f"Run-time of {scenario2run}")
    st.write('---')
    st.write('Core KPIs')
    #st.write(rheum_model_S2.batch_kpi.iloc[[0,1,2,3,7],])
    st.write(rheum_model_S2.batch_kpi)
    # Show chart
    st.write('---')
    st.write('Waiting time KPIs')
    st.pyplot(fig_monappKPI_reps_S2)
    st.write('---')
    st.write('Audit KPIs')
    st.pyplot(fig_audit_reps_S2)
    st.write('---')
    st.write('Patients seen KPIs')
    st.pyplot(fig_monappKPIn_reps_S2)

    #st.write('---')
    #st.write('Single rpelication example:')
    # Show chart
    #st.pyplot(chart_S2)
    # Display results (text)
    # for t in text_S2:
    #     st.write(t)
    #st.write(quant_S2)
    st.write(scenario2run)


