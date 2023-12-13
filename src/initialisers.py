""" includes initialising functions class (globals)"""

import pandas as pd
import numpy as np


class g:
    """ Class to store global parameter values.
    
    Not used optimally - should be revisited. Currently code creates instances per replication, should likely be global
    """
    mean_interOPA = np.round(4.5 * 365/12,0) # mean days inbetween appointments (traditional pathway) - exponential distribution
    interOPA_tri = np.round(np.array([3,6,4.5])*365/12) # days inbetween appointments (traditional pathway) - low, high , mode for triangular distribution
    t_decision = 1 * 365 # days, time mark for pathway PIFU decision / stratification (from first OPA appointment)
    obs_duration=365*3 # observation period or window (days) for simulation, in addition to warm-up period

    debug=True # whether to print to console
    debuglevel = 1 # level of debug prints - 1 as lowest ; 4 for most detailed


    def __init__(self,in_res=5,in_inter_arrival=1,in_prob_pifu=0,in_path_horizon_y=3,audit_interval=7,in_reps=1,repid=1,savepath='temp',in_FOavoidable=0,in_interfu_perc=0.6):
        """ Initialise global parameter values."""

        self.prob_firstonly = 0.35 # % of rheumatology RTT patients have no follow-ups | Baseline: ~35% with no follow-ups
        self.number_of_runs=in_reps # no replications
        self.prob_pifu=in_prob_pifu # PIFU proportion - probability of PIFU pathway for non first-only pathways[%]
        self.wl_inter = in_inter_arrival # inter-arrival time [days]
        self.number_of_slots = in_res # number of daily slots [slots]
        self.in_FOavoidable = in_FOavoidable # A&G proportion - proportion of first-only pathways avoidable via A&G [%]

        self.warm_duration=365*5 # Warm-up period or window [days] for simulation
        self.max_fuopa_tenor_y = in_path_horizon_y # Patient follow-up horizon [years], simplification on how long each non first-only pathway lasts (years)
        self.max_fuopa_tenor = in_path_horizon_y * 365 # Patient follow-up horizon [days]
        self.appt_counter = 0 # Counter for number of appointments [appointments], initialised
        self.DNA_pifu_pro=0.07 # Did not attend probability (PIFU pathways) [%]
        self.DNA_tra_pro=0.077 # Did not attend probability (traditional pathways) [%]
        self.interfu_perc = in_interfu_perc # Percentage increase in inter-appointment interval with PIFU (vs traditional), i.e. 0.6 means 60% longer interval
        self.mean_interPIFU = np.round(np.sum(self.interOPA_tri)/3 * (1+self.interfu_perc),0) # resultant mean days inbetween appointments for PIFU pathways (from traditional triangular mean) [days]. Used in exponential distribution
        self.PIFUbigbang = False # [boolean] Whether, when PIFU starts being used, it is offered to all eligible patients when they visit (e.g. those already followed up for years) - big-bang - or only new eligible patients

        self.unavail_on = False # [boolean] Whther to use resource unavailability functionality
        self.unavail_byshock = False # [boolean] for now model only for unavailability by single shock period OR by periodic (e.g. weekends). Can be improved in future.

        # Parameters for a shock to available resource (one-off) - if unavail_byshock = True and unavail_on = True
        self.unavail_shock_tmin = 365 * 3 # [days] Time in simulation days from which shock starts
        self.unavail_shock_period = np.floor(365/12 * 12) # [days] Period in simulation days for which the shock lasts e.g. from half-March (emergency response) to half-August (letter on Aug2020 for NHS response)
        self.unavail_shock_nrslots = int(np.floor(in_res*1)) # [slots] No of slots unavailable during shock

        # Parameters for periodic resource unavailability (e.g. weekends, leave) - if unavail_byshock = False and unavail_on = True
        self.unavail_slot=2 # [days] Time window of unavailability
        self.unavail_freq_slot=5 # [days] Time window of availability
        self.unavail_nrslots = int(np.floor(in_res*1)) # [slots] No of slots unavailable during each unavail period. Baseline: 100% of slots unavailable

        self.savepath = savepath # [string] Save path for outputs
        self.repid = repid # [integer] Id of current replication (within batch)
        self.loglinesave = True # if true saves each line to file, if false creates dataframe that stays in memory (former found to be more efficient)
        self.audit_time = []
        self.audit_interval = audit_interval # time step for audit metrics [simulation days]
        self.audit_patients_waiting = [] # vector of patients waiting at audit timepoints. populated in perform_audit
        self.audit_patients_waiting_p1 = [] # vector of priority 1 patients waiting at audit timepoints. populated in perform_audit
        self.audit_patients_waiting_p2 = [] # vector of priority 2 patients waiting at audit timepoints. populated in perform_audit
        self.audit_patients_waiting_p3 = [] # vector of priority 3 patients waiting at audit timepoints. populated in perform_audit
        self.audit_patients_in_system = [] # vector of patients in system. populated in perform_audit (lifetime to discharge...)
        self.audit_resources_used = [] # vector of resources in use at audit timepoints. populated in perform_audit
        self.patients_waiting = 0 # auxiliary counter for patients waiting
        self.patients_waiting_by_priority = [0,0,0] # auxiliary counters for patients waiting by priority
        self.results = pd.DataFrame() # populated in build_audit_results. audit results in one DF
        self.appt_queuing_results = (pd.DataFrame(
            columns=['P_ID','App_ID','priority','type',"pathway",'q_time','start_q','DNA'])) # populated in pathway (FOPA_Patiet)
        self.appt_queuing_results.set_index("App_ID", inplace=True) # reset index

    def change_reps(self,reps):
        """ Change number of replications for batch run """
        self.number_of_runs = reps

    def give_reps(self):
        """ Give number of replications for batch run """
        return self.number_of_runs
