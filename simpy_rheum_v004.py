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

# Code to compute (small sample) confidence interval, taken from web.
def mean_confidence_interval(data, confidence=0.95):
            a = 1.0 * np.array(data)
            n = len(a)
            m, se = np.mean(a), st.sem(a)
            h = se * st.t.ppf((1 + confidence) / 2., n-1)
            return m, m-h, m+h

# Class to store global parameter values.  Not used optimally - should be revisited. Currently code creates instances per replication,
# should likely be global
class g:
    mean_interOPA = np.round(4.5 * 365/12,0) # mean days inbetween appointments (traditional pathway) - exponential distribution
    interOPA_tri = np.round(np.array([3,6,4.5])*365/12) # days inbetween appointments (traditional pathway) - low, high , mode for triangular distribution
    t_decision = 1 * 365 # days, time mark for pathway PIFU decision / stratification (from first OPA appointment)
    obs_duration=365*3 # observation period or window (days) for simulation, in addition to warm-up period
    
    debug=True # whether to print to console
    debuglevel = 1 # level of debug prints - 1 as lowest ; 4 for most detailed

        
        def __init__(self,in_res,in_inter_arrival,in_prob_pifu,in_path_horizon_y,audit_interval,in_FOavoidable,in_interfu_perc):
        self.prob_firstonly = 0.35 # % of rheumatology RTT patients have no follow-ups | Baseline: ~35% with no follow-ups
       # self.number_of_runs=in_reps # no replications
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
        
       # self.savepath = savepath # [string] Save path for outputs
       # self.repid = repid # [integer] Id of current replication (within batch)
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
    
    # Change number of replications
    def change_reps(self,reps):
        self.number_of_runs = reps


class patient_blocker:

    def __init__(self,blockerid):        
        self.blockerid = blockerid
    
# Class representing our RTT patients entering the secondary care rheumatology pathway.
class FOPA_Patient:

    # The following dictionaries store patients
    all_patients = {}

    # Initialised parameters
    def __init__(self, p_id,prob_pifu,in_path_horizon,DNA_pifu_pro,DNA_tra_pro):
        self.apptype = []
        self.id = p_id # patient id
        self.prob_pifu = prob_pifu # PIFU probability
        self.q_time_fopa = float("nan") # queueing time first outpatient (deprecated)
        self.q_time_fuopa = float("nan") # queueing time follow-up outpatient (deprecated)
        self.q_times_fuopa = [] # queueing times follow-up appointments (deprecated)
        self.topifu = False # Whether patient is on PIFU pathway or not. Initialised to False
        self.priority = 3 # Patient priority. Initialise at lower priority
        self.max_fuopa = 999 # Maximum follow up outpatients per pathway (not used, see below)
        self.max_fuopa_tenor = in_path_horizon # Maximum horizon for follow-up per pathway (as days since first OPA) [days]
        self.used_fuopa = 0 # Counter for no f/up appointments used
        self.ls_appt = [] # Vector to hold ids of appointments for pathway
        self.DNA_pifu_pro=DNA_pifu_pro # DNA probability (traditional appt)
        self.DNA_tra_pro=DNA_tra_pro # DNA probability (PIFU appt)
        self.tradition_dna=False # instantiate DNA for traditional type as False
        self.pifu_dna=False # instantiated DNA for PIFU type as False
        self.FOavoided = False # instantiate whether first outpatient avoided fully as False
        self.df_appt_to_add = (pd.DataFrame(
            columns=['P_ID','App_ID','priority','type',"pathway",'q_time','start_q'])) # Initialise dataframe to hold appointment log
        self.ls_appt_to_add = [] # Initialise list to hold appointment info 
        self.ls_patient_to_add = [] # Initialise list to hold pathway info
        self.RTT_sub = 0 # Initialise a variable that can 'scramble' further priority of first outpatient - priority increment (to increase variance) [Commented out]
    
    # Method to decide and assign at random PIFU faith
    def triage_decision(self):
        if random.random() < self.prob_pifu:
            self.topifu = True
            self.type = "PIFU"
        else:
            self.type = "TFU"
            #self.priority = 2 # assign 2nd highest priority in the sense that those on traditional already are scheduled in/planned so not really 'moveable'
    
    # Method to change priority to '2', i.e. PIFU
    def give_pifu_priority(self):
        self.priority =  2 # assign 2nd highest priority in the sense that those on traditional already are scheduled in/planned so not really 'moveable'
        
    # Method to change priority to '1', i.e. traditional / booked
    def give_tfu_priority(self):
        self.priority = 1
    
    # Method to assign 'first-only' pathway and appointment status or 'first' (of more)
    def assign_firstonly(self,prob_firstonly):
        if random.random() < prob_firstonly:
            self.type = "First-only" # single appointment
            self.apptype = "First-only"
            self.max_fuopa_tenor = 0 # no follow-ups
        else:
           self.type = "First" # first leading to long-term follow-up
           self.apptype = "First"
            
    # Method to decide and assign at random DNA faith of appointment (PIFU)
    def decision_DNA_pifu(self):
        if  random.random() < self.DNA_pifu_pro:
            self.pifu_dna = True
             
    # Method to decide and assign at random DNA faith of appointment (first ; traditional)
    def decision_DNA_tradtion(self):
        if  random.random() < self.DNA_tra_pro:
            self.tradition_dna = True 
    
    # Method to further scramble priority of first appointments (rationale: add variability to RTT waiting list times).
    # Currently not in use, needs development
    def sub_RTT_priority(self):
        # to add some spread to the RTT list
        
        # Bernoulli p=0.5
        #self.RTT_sub = random.randint(0, 1)
        
        # Bernoulli p chosen
        #if  random.random() < 0.3:
        #    self.RTT_sub=1
        
        # Sampling from options, equiprobable
        #self.RTT_sub = np.random.choice([0,1,2,3,4,5,6,7,8,9], 1, replace=False)[0]
        
        #print()
        pass
    
    # Method to decide and assign at random first-only pathway or not faith    
    def avoidable_firstonly(self,in_FOavoidable):
        # A&G, whether a first-only appointment can be avoided
        if random.random() < in_FOavoidable:            
            self.FOavoided = True
        
# Class representing our overall model of the rheumatology outpatient clinic
class rheum_Model:
    # Here, the constructor sets up the SimPy environment, sets a patient
    # counter to 0 (which we'll use for assigning patient IDs), and sets up
    # our resources (here appointment slot units (symbolycally 15 min), with capacity given by
    # the number stored in the g class)

    def __init__(self, run_number, in_res , in_inter_arrival, in_prob_pifu, in_path_horizon_y,audit_interval,repid, savepath,in_FOavoidable,in_interfu_perc):
        self.env = simpy.Environment() # instance of environment

        self.g = g(in_res,in_inter_arrival,in_prob_pifu, in_path_horizon_y, audit_interval, in_FOavoidable = in_FOavoidable,in_interfu_perc=in_interfu_perc) # instance of global variables for this replication
        self.repid=repid
        self.patient_counter = 0 # patient counter instantiated to 0
        self.block_counter = 0 # block counter instantiated to 0 (to control that right no of unavailable slots are enforced)
        
        self.consultant = simpy.PriorityResource(self.env, capacity=self.g.number_of_slots) # set up resources, i.e. appointment slot units (assume 1 unit - 15 min slot)

        self.run_number = run_number # [integer] run number id
        self.savepath = savepath # [string] savepath
        
        self.mean_q_time_total = pd.DataFrame() # [running but deprecated]
        self.results_df = pd.DataFrame() # [running but deprecated]
        self.results_df["P_ID"] = [] # [running but deprecated]
        self.results_df["Q_Time_fopa"] = [] # [running but deprecated]
        self.results_df["Q_Time_fuopa"] = [] # [running but deprecated]
        self.results_df.set_index("P_ID", inplace=True) # [running but deprecated]
        
    # A method that generates patients arriving for the RTT outpatient 'clinic'
    def generate_wl_arrivals(self):
        # Keep generating indefinitely (until the simulation ends)
        while True:
            # Increment the patient counter by 1
            self.patient_counter += 1
            
            # Create a new patient - an instance of the FOPA_Patient 
            # class, and give the patient an ID determined by the patient
            # counter, its PIFU prob, its follow-up tenor, its DNA probabilities
            wp = FOPA_Patient(self.patient_counter,self.g.prob_pifu,self.g.max_fuopa_tenor, self.g.DNA_pifu_pro,self.g.DNA_tra_pro)
            
            # Add patient to dictionary of patients
            FOPA_Patient.all_patients[wp.id] = wp
            
            # Get the SimPy environment to run the attend_OPA method
            # with this patient
            self.env.process(self.attend_OPA(wp))
            
            # Randomly sample the time to the next patient arriving for the
            # RTT outpatient 'clinic'.  The details of patient and pathway are stored in the g replication instance.
            #sampled_interarrival = int(np.round(random.expovariate(1.0 / self.g.wl_inter),0))
            sampled_interarrival = random.expovariate(1.0 / self.g.wl_inter)
            
            # Freeze this function until that time has elapsed
            yield self.env.timeout(sampled_interarrival)


    
    # A method to obstruct a single slot
    def obstruct_slot(self,patient_blocker,unavail_timeperiod):
            
        # Once this time has elapsed, request a slot wiyh priority
        # of -1 (so that we know this will get the top priority, as none
        # of our pathway appointment requests will have a negative priority), and hold them
        # for the specified unavailability amount of time
        with self.consultant.request(priority=-1) as req:
            # Freeze the function until the request can be met (this
            # ensures that the slot will finish before becoming unavailabe)
            yield req
            
            yield self.env.timeout(unavail_timeperiod)
                
    # A method to obstruct multiple slots            
    def obstruct_slots(self):
        
        # If unavailability is single shock period
        if self.g.unavail_byshock:
            
            # Let the period-to-unavailability pass
            yield self.env.timeout(self.g.unavail_shock_tmin)
            
            # Iterate over number of slots that needs blocking
            for i in range(self.g.unavail_shock_nrslots):
                
                self.block_counter += 1
                slot_block = patient_blocker(self.block_counter)
                
                # Get the SimPy environment to run the obstruct_slot method with this slot block
                self.env.process(self.obstruct_slot(slot_block,self.g.unavail_shock_period))
                
                if self.g.debug and self.g.debuglevel>=3:
                    print ("Appointment will not be able to book at",
                           f"{self.env.now + self.g.unavail_freq_slot:.1f}")
        
        # If unavailability is periodic                     
        else:
            
            # Run indefinitely till end of simulation
            while True:
                # Iterate over number of slots that needs blocking 
                for i in range(self.g.unavail_nrslots):
                    
                    self.block_counter += 1
                    slot_block = patient_blocker(self.block_counter)
                    
                    # Get the SimPy environment to run the obstruct_slot method with this slot block
                    self.env.process(self.obstruct_slot(slot_block,self.g.unavail_slot))
                
                    if self.g.debug and self.g.debuglevel>=3:
                        print ("Appointment will not be able to book at",
                               f"{self.env.now + self.g.unavail_freq_slot:.1f}")
                
                # Freeze the function for the time period during which no unavailability
                yield self.env.timeout(self.g.unavail_freq_slot)
                
            
    # A method that models the processes / RTT patient pathway for attending the outpatient rheumatology clinic.
    # The method needs to be passed an RTT patient who will go through these processes
    def attend_OPA(self, patient):
                        
        patient.assign_firstonly(self.g.prob_firstonly) # Assign whether first-only pathway       
        patient.sub_RTT_priority() # add some variability to priority within RTT queue (increment to its '3' priority)        
        patient.avoidable_firstonly(self.g.in_FOavoidable) # Assign, if 'first-only', whether pathway is avoided or not (e.g. A&G)
        
        # If first-only AND avoidance from A&G AND past warm-up period
        if patient.type == "First-only" and patient.FOavoided and self.env.now > self.g.warm_duration:
            # if first-only pathway and avoidable through A&G and current day within intervention/study period, skip anything further for patient
            if self.g.debug and self.g.debuglevel>=2:
                print(f"Patient {patient.id} had first outpatient avoided. Not added to log.")
          
        # ELSE
        else:            

            ############################################
            ### First outpatient appointment ####
            ############################################        
            
            # Record the time the patient started queuing for the first outpatient
            start_q_fopa = self.env.now
            self.g.appt_counter +=1 # increment
            patient.ls_appt.append(self.g.appt_counter) # append
            self.g.patients_waiting += 1 # increment
            self.g.patients_waiting_by_priority[patient.priority-1] += 1 # increment
            
            # Request a slot
            with self.consultant.request(priority = patient.priority + patient.RTT_sub) as req:
                # Freeze the function until the request for a slot can be met
                yield req
                
                # if non-first-only pathway
                if patient.type != "First-only":                
                # determine already their PIFU faith                
                    if self.g.PIFUbigbang: # if 'big-bang' ('stock'), i.e. PIFU applied to all FY cohorts / pathways, draw PIFU for all                    
                        patient.triage_decision()
                        
                    else:
                        if self.env.now > (self.g.warm_duration - self.g.t_decision): # else, draw PIFU only if they are a new pathway entering PIFU eligibility (1 year follow-up) from after warm-up period                    
                            patient.triage_decision()
                        else:
                            patient.type = "TFU" # if pathway started pre warm-up, keep them on traditional
        
                # reduce patients waiting counts
                self.g.patients_waiting_by_priority[patient.priority-1] -= 1 # decrement
                self.g.patients_waiting -= 1 # decrement
                
                # Record the time the patient finished queuing for a consultant
                end_q_fopa = self.env.now
    
                 # Calculate the time this patient spent queuing for the consultant (FU) and
                # store in the patient's attribute
                patient.q_time_fopa = end_q_fopa - start_q_fopa
                
                if self.g.debug and self.g.debuglevel>=2:
                    print(f"Req {patient.ls_appt[-1]}: Patient {patient.id} queued {np.round(patient.q_time_fopa,2)} days for 1st app. Priority {patient.priority}")
                
                # Freeze this function until the day time unit has elapsed
                #yield self.env.timeout(1) # freeze for one time-unit (a day) - that same slot will only be available the next day
                yield self.env.timeout(2) # freeze for two time-units (two dayz) - that same slot will only be available in two days (simplification/ discretisation to deal with first outpatient being ~30 min, so 2 of our slot units)
                
                patient.decision_DNA_tradtion() # Decide whether this is a DNA or not
                
                if patient.tradition_dna==True:
                    if self.g.debug and self.g.debuglevel>=2:
                        print(f" Patient {patient.id} queued {np.round(patient.q_time_fopa,2)} days for 1st app. Priority {patient.priority} but didn't attend the appointment")
    
                else:
                    if self.g.debug and self.g.debuglevel>=2:
                        print(f" Patient {patient.id} queued {np.round(patient.q_time_fopa,2)} days for 1st app. Priority {patient.priority}")
                    
    
                # Line/list to add to appointment log held in memory (df) or saved (csv)         
                patient.ls_appt_to_add = [patient.id, patient.ls_appt[-1],patient.priority,patient.apptype,patient.type,patient.q_time_fopa,start_q_fopa,patient.tradition_dna,self.repid]
                patient.ls_patient_to_add = [patient.id , patient.q_time_fopa,999,self.repid] # deprecated
                    
                # Whether to save each log line or hold in memory by appending
                if self.g.loglinesave:
                    
                    with open(self.savepath +"appt_result.csv", "a") as f:
                         writer = csv.writer(f, delimiter=",")
                         writer.writerow(patient.ls_appt_to_add)
                    
                    if start_q_fopa > self.g.warm_duration: # don't save things in warm-up period
                        with open(self.savepath +"patient_result2.csv", "a") as f:
                             writer = csv.writer(f, delimiter=",")
                             writer.writerow(patient.ls_patient_to_add)
                                              
                else:
                    patient.df_appt_to_add = pd.DataFrame( columns = ["P_ID","App_ID","priority","type","pathway","q_time","start_q","DNA","rep"] ,
                                                          data=[patient.ls_appt_to_add]) # row list to row dataframe
                    patient.df_appt_to_add.set_index("App_ID", inplace=True)                    
                    self.g.appt_queuing_results=self.g.appt_queuing_results.append(patient.df_appt_to_add)
                    
                    df_to_add = pd.DataFrame( columns = ["P_ID","Q_time_fopa","Q_time_fuopa","rep"], data =[patient.ls_patient_to_add])
                    df_to_add.set_index("P_ID", inplace=True)                                    
                    if start_q_fopa > self.g.warm_duration: # don't save things in warm-up period           
                        self.results_df = self.results_df.append(df_to_add)
                                   
            patient.give_tfu_priority() # Assign traditional follow-up priority to subsequent requests
            
            ############################################
            ### Traditional Follow-up appointments ####
            ############################################
            
            while self.env.now - end_q_fopa < patient.max_fuopa_tenor: # While within pathway horizon / tenor
            #while patient.used_fuopa < patient.max_fuopa:

                # Determine time till next needing F/U
                #sampled_interfu_duration = int(random.expovariate(1.0 / g.mean_interOPA)) # integer only (days)
                sampled_interfu_duration = int(random.triangular(g.interOPA_tri[0],g.interOPA_tri[1],g.interOPA_tri[2]))
                # Freeze this function until time has elapsed
                yield self.env.timeout(sampled_interfu_duration)
                
                self.g.patients_waiting += 1 # increment
                self.g.patients_waiting_by_priority[patient.priority-1] += 1 # increment
                patient.used_fuopa+=1 # count the follow-up outpatient
                
                self.g.appt_counter +=1 # increment
                patient.ls_appt.append(self.g.appt_counter) # append
                start_q_fuopa = self.env.now # current time
                # Request a slot for follow-up
                with self.consultant.request(priority = patient.priority) as req:
                    # Freeze the function until the request for a slot can be met
                    yield req
                                   
                    # reduce patients waiting counts
                    self.g.patients_waiting_by_priority[patient.priority-1] -= 1 # decrement
                    self.g.patients_waiting -= 1
    
                    # Record the time the patient finished queuing for a follow-up slot
                    end_q_fuopa = self.env.now
                                   
                    # Calculate the time this patient spent queuing for a slot and
                    # store in the patient's attribute
                    patient.q_times_fuopa.append(end_q_fuopa - start_q_fuopa) # add latest followup
                    #print(patient.q_time_fuopa)
                    
                    if patient.used_fuopa ==1 : # deprecated, not relevant
                        patient.q_time_fuopa = end_q_fuopa - start_q_fuopa
                    
                    patient.decision_DNA_tradtion() # Determine DNA faith
                    if patient.tradition_dna==True:
                        if self.g.debug  and self.g.debuglevel>=2:
                            print(f"Req {patient.ls_appt[-1]}: Patient {patient.id} queued {np.round(end_q_fuopa - start_q_fuopa,2)} for app {patient.used_fuopa}. Priority {patient.priority}")
                    else:
                        if self.g.debug  and self.g.debuglevel>=2:
                            print(f"Req {patient.ls_appt[-1]}: Patient {patient.id} queued {np.round(end_q_fuopa - start_q_fuopa,2)} for app {patient.used_fuopa}. Priority {patient.priority} but DNAd")
                    
                            
                    # Freeze this function until the day time unit has elapsed
                    yield self.env.timeout(1) # freeze for one time-unit (a day) - that same slot will only be available the next day
                         
                    
                    # Add to appointment log or save            
                    patient.ls_appt_to_add = [patient.id, patient.ls_appt[-1],patient.priority,"Traditional",patient.type,patient.q_time_fuopa,start_q_fuopa,patient.tradition_dna,self.repid]                            
                    if self.g.loglinesave:                    
                        with open(self.savepath +"appt_result.csv", "a") as f:
                             writer = csv.writer(f, delimiter=",")
                             writer.writerow(patient.ls_appt_to_add)
    
                    else:                         
                        patient.df_appt_to_add = pd.DataFrame( columns = ["P_ID","App_ID","priority","type","pathway","q_time","start_q","DNA","rep"] ,
                                                              data=[patient.ls_appt_to_add])
                        patient.df_appt_to_add.set_index("App_ID", inplace=True)
                        self.g.appt_queuing_results=self.g.appt_queuing_results.append(patient.df_appt_to_add)
    
                
                # If current simulation time is beyond warm-up , and if current time exceeds timing for PIFU eligilibity to be adequate for this patient / pathway
                if self.env.now - end_q_fopa > self.g.t_decision  and self.env.now > self.g.warm_duration:
                    # if patient is PIFU pathway assigned, 'break' from traditional appointments to enable PIFU appointment cycle below
                    if patient.topifu:
                        break
                
    
            # If patient is PIFU pathway assigned (will only get to this portion of code if 'break' from traditional appointment cycle)    
            if patient.topifu:
                
                if self.g.debug  and self.g.debuglevel>=2:
                    print(f"Patient {patient.id} PIFU. Follows {patient.used_fuopa} traditional apps.")
                
                patient.give_pifu_priority() # Assign PIFU priority to all further slot requests
                #print(f"Patient {patient.id} entered PIFU. Has {patient.used_fuopa} traditional apps. Priority {patient.priority}")
                
                while patient.topifu: # while true (indefinitely while simulation running, will break if not within follow-up horizon)
                    
                    patient.used_fuopa+=1 # count the follow-up outpatient
                    # Determine time till next needing PIF/U
                    sampled_interpifu_duration = int(np.round(random.expovariate(1.0 / self.g.mean_interPIFU),0)) # integer only (days)
                    # Freeze this function until time has elapsed (inter-pifu)
                    yield self.env.timeout(sampled_interpifu_duration)
                    
                    self.g.appt_counter +=1 # increment
                    patient.ls_appt.append(self.g.appt_counter) # append
                    self.g.patients_waiting += 1 # increment
                    self.g.patients_waiting_by_priority[patient.priority-1] += 1 # increment
    
    
                    start_q_pifuopa = self.env.now            
                    # Request slit
                    with self.consultant.request(priority = patient.priority) as req:
                        # Freeze the function until the request for a slot can be met
                        yield req
    
                        # reduce patients waiting counts
                        self.g.patients_waiting_by_priority[patient.priority-1] -= 1 # decrement
                        self.g.patients_waiting -= 1 # decrement
                        
                        # Record the time the patient finished queuing for a consultant
                        end_q_pifuopa = self.env.now
                        
                        # Calculate the time this patient spent queuing for a consultant and
                        # store in the patient's attribute
                        patient.q_time_pifuopa = end_q_pifuopa - start_q_pifuopa
                        
                        # Freeze this function until the day time unit has elapsed
                        yield self.env.timeout(1) # freeze for one time-unit (a day) - that same slot will only be available the next day
                    
                        patient.decision_DNA_pifu() # Determine DNA status of appointment
                        if patient.pifu_dna==True:
                            if self.g.debug  and self.g.debuglevel>=2:
                                print(f"Req {patient.ls_appt[-1]}: Patient {patient.id} queued {np.round(patient.q_time_pifuopa,2)} days for app {patient.used_fuopa} - PIFU. Priority {patient.priority} but did not attend")
                        else:
                             if self.g.debug  and self.g.debuglevel>=2:
                                print(f"Req {patient.ls_appt[-1]}: Patient {patient.id} queued {np.round(patient.q_time_pifuopa,2)} days for app {patient.used_fuopa} - PIFU. Priority {patient.priority}")
                        
                            
                        # Add to appointment log or save            
                        patient.ls_appt_to_add = [patient.id, patient.ls_appt[-1],patient.priority,"PIFU",patient.type,end_q_pifuopa - start_q_pifuopa,start_q_pifuopa,patient.pifu_dna,self.repid]                            
                        if self.g.loglinesave:                    
                            with open(self.savepath +"appt_result.csv", "a") as f:
                                 writer = csv.writer(f, delimiter=",")
                                 writer.writerow(patient.ls_appt_to_add) 
                        
                        else:
                            patient.df_appt_to_add = pd.DataFrame( columns = ["P_ID","App_ID","priority","type","pathway","q_time","start_q","DNA","rep"] ,
                                                                  data=[patient.ls_appt_to_add])
                            patient.df_appt_to_add.set_index("App_ID", inplace=True)
                            self.g.appt_queuing_results=self.g.appt_queuing_results.append(patient.df_appt_to_add)
                    
                    # break if time elapsed since first appointment exceeds follow-up horizon    
                    if self.env.now - end_q_fopa > patient.max_fuopa_tenor:
                        break
             
            else:
                
                if self.g.debug  and self.g.debuglevel>=2:
                    print(f"Patient {patient.id} not PIFU. Follows {patient.used_fuopa} traditional apps.")
            
     
            # Delete patient (removal from patient dictionary removes only
                # reference to patient and Python then automatically cleans up)
            del FOPA_Patient.all_patients[patient.id]

    

    def build_audit_results(self):
        """Compiles single run audit results into single dataframe held in g.results """
        
        self.g.results['time'] = self.g.audit_time
        
        self.g.results['patients in system'] = (
            self.g.audit_patients_in_system)
        
        self.g.results['all patients waiting'] = (
            self.g.audit_patients_waiting)
        
        self.g.results['priority 1 patients waiting'] = (
            self.g.audit_patients_waiting_p1)
        
        self.g.results['priority 2 patients waiting'] = (
            self.g.audit_patients_waiting_p2)
        
        self.g.results['priority 3 patients waiting'] = (
            self.g.audit_patients_waiting_p3)
        
        self.g.results['resources occupied'] = (
            self.g.audit_resources_used)
        
        self.g.results['rep']=self.repid


    def calculate_mean_q_time(self):
        """Deprecated. A method that calculates the average quueing time for the nurse.
        We can call this at the end of each run """
        
        self.results_df['Q_Time_total']=self.results_df.sum(axis=1)
        self.mean_q_time_total = self.results_df.mean(axis=0)
        
    """ Deprecated. A method to write run results to file.  Here, we write the run number
    # against the the calculated mean queuing time for the nurse across
    # patients in the run.  Again, we can call this at the end of each run"""
    def write_run_results(self):
        with open(self.savepath +"trial_results.csv", "a") as f:
            writer = csv.writer(f, delimiter=",")
            results_to_write = [self.run_number,
                                self.mean_q_time_total['Q_time_fopa'],
                                self.mean_q_time_total['Q_time_fuopa'],
                                self.mean_q_time_total['Q_Time_total']]
            print(results_to_write)
            writer.writerow(results_to_write)
    

    """ Plot results relevant to one run """
    def chart(self):
        # plot results at end of run #
        self.results_df = self.results_df.sort_index()
        #self.g.appt_queuing_results = self.g.appt_queuing_results.sort_values(by=['start_q'])
        self.g.appt_queuing_results = self.g.appt_queuing_results.sort_index()

        # Define fig size
        fig = plt.figure(figsize=(20,12))
        
        # Figure 1: first opa patient level results
        ax1 = fig.add_subplot(2,4,1) # 1 row, 4 cols, pos 1
        x = self.results_df.index
        y = self.results_df['Q_time_fopa']

        ax1.plot(x,y,marker='^',color='k')
        ax1.set_xlabel('Patient')
        ax1.set_ylabel('Queuing time')
        ax1.legend()
        ax1.set_title('Queuing time for first outpatient')
        ax1.grid(True,which='both',lw=1,ls='',c='.75')

        # Figure 2: appointment level results - queueing time vs appointment id
        ax22 = fig.add_subplot(242)  # 1 row, 4 cols, chart position 2
        x = self.g.appt_queuing_results.index
        # Chart loops through 3 priorites
        markers = ['o', 'x', '^']
        for priority in range(1, 4):
            x = (self.g.appt_queuing_results[self.g.appt_queuing_results['priority'] == priority].index)
            
            y = (self.g.appt_queuing_results
                 [self.g.appt_queuing_results['priority'] == priority]['q_time'])
            
            ax22.plot(x, y, marker=markers[priority - 1], label='Priority ' + str(priority))
        ax22.set_xlabel('Appointment')
        ax22.set_ylabel('Queuing time')
        ax22.legend()
        ax22.set_title('Queuing time by type')
        ax22.grid(True, which='both', lw=1, ls='--', c='.75')

        
        di = {1:"Follow-up", 2:"Follow-up",3:"RTT"}
        step = 365
        mon_appointments = self.g.appt_queuing_results
        mon_appointments['interval'] = np.round(((mon_appointments['start_q']+mon_appointments['q_time'])//step)*step,0)
        mon_appointments['ttype']=mon_appointments['priority']
        mon_appointments.replace({"ttype": di},inplace=True)
        ax32 = fig.add_subplot(2,4,3)
        sns.violinplot(ax=ax32,x='interval',y='q_time',data=mon_appointments,hue='ttype',palette="Set2",split=False)
        ax32.set_xticklabels(ax32.get_xticklabels(),rotation = 30)


        # Figure 4: Staff usage
        ax3 = fig.add_subplot(244)  # 1 row, 4 cols, chart position 4
        x = self.g.results['time']
        y = self.g.results['resources occupied']
        ax3.plot(x, y, label='Resources occupied')
        ax3.set_xlabel('Time')
        ax3.set_ylabel('Resources occupied')
        ax3.set_title('Resources occupied')
        ax3.grid(True, which='both', lw=1, ls='--', c='.75')

        # Figure 5-7: appointment level results - queueing time vs app id - facetwrap
        mylabels = ['Traditional F/U appt waiting','PIFU appt waiting','First outpatient appt waiting','All']
        markers = ['o', 'x', '^']
        colors = ['r','g','k']
        # for priority in range(1, 4):
        #     ax22 = fig.add_subplot(3,4,4 + priority)  # 1 row, 4 cols, chart position 2
            
        #     x = self.g.appt_queuing_results.index
        #     # Chart loops through 3 priorites
            
        #     x = (self.g.appt_queuing_results[self.g.appt_queuing_results['priority'] == priority].index) 
        #     y = (self.g.appt_queuing_results
        #          [self.g.appt_queuing_results['priority'] == priority]['q_time'])
            
        #     ax22.plot(x, y, marker=markers[priority - 1], linestyle="none",label=mylabels[priority-1],color = colors[priority - 1])
        #     ax22.set_xlabel('Appointment nr')
        #     ax22.set_ylabel('Queuing time')
        #     ax22.legend()
        #     ax22.set_title('Queuing time by type')
        #     ax22.grid(True, which='both', lw=1, ls='--', c='.75')
        
        
        # Figure 5: System level queuing results - number waiting vs audit time point # change to bar chart?
        ax2 = fig.add_subplot(245)
        x = self.g.results['time']
        y1 = self.g.results['priority 1 patients waiting']
        y2 = self.g.results['priority 2 patients waiting']
        y3 = self.g.results['priority 3 patients waiting']
        y4 = self.g.results['all patients waiting']
        ax2.plot(x, y1, marker='o', linestyle="none",label='Traditional F/U appt waiting')
        ax2.plot(x, y2, marker='x',linestyle="none", label='PIFU appt waiting')
        ax2.plot(x, y3, marker='^',linestyle="none", label='First outpatient appt waiting')
        ax2.plot(x, y4, marker='s',linestyle="none", label='All')
        ax2.set_xlabel('Time')
        ax2.set_ylabel('Appointment waits')
        ax2.legend()
        ax2.set_title('No. Appointment waits by type')
        ax2.grid(True, which='both', lw=1, ls='--', c='.75')
        
        myvars = ['priority 1 patients waiting','priority 2 patients waiting','priority 3 patients waiting']
        for priority in range(1, 4):
            ax32 = fig.add_subplot(2,4,5 + priority)  # 1 row, 4 cols, chart position 2
            
            x = self.g.results['time']
            y1 = self.g.results[myvars[priority-1]]
            y4 = self.g.results['all patients waiting']
            ax32.plot(x, y1, marker=markers[priority-1], linestyle="none",label=mylabels[priority-1],color=colors[priority-1])
            ax32.plot(x, y4, marker='',linestyle="-", label='All')
            ax32.set_xlabel('Time')
            ax32.set_ylabel('Appointment waits')
            ax32.legend()
            ax32.set_title('No. Appointment waits')
            ax32.grid(True, which='both', lw=1, ls='--', c='.75')

        # Adjust figure spacing
        fig.tight_layout(pad=2)
        
        

        # return fig
        return fig

    """ single run summaries for streamlit (quartiles) - NEED CHECKING, NOT RUNNING AS SHOULD """
    def summarise(self):
        """Produces displayed text summary of model run"""
        
        self.g.appt_queuing_results['priority']=pd.to_numeric(self.g.appt_queuing_results['priority'])
        self.g.appt_queuing_results['q_time']=pd.to_numeric(self.g.appt_queuing_results['q_time'])
        
        # Put text in string to return
        text = []
        text.append('APPOINTMENT-CENTERED METRICS:')
        text.append ('Lower quartile time in system by priority:')
        quant=self.g.appt_queuing_results[['priority','q_time']].groupby('priority').quantile(0.25)
        
        text.append (self.g.appt_queuing_results[['priority','q_time']].convert_dtypes().groupby('priority').quantile(0.25))
        text.append ('Median time in system by priority:')
        text.append (self.g.appt_queuing_results[['priority','q_time']].groupby('priority').quantile(0.50))
        text.append ('Upper quartile time in system by priority:')
        text.append (self.g.appt_queuing_results[['priority','q_time']].groupby('priority').quantile(0.75))
        text.append ('Maximum time in system by priority:')
        text.append (self.g.appt_queuing_results[['priority','q_time']].groupby('priority').quantile(1))
        text.append('---')
        text.append ('SYSTEM-CENTRED METRICS - audit points:')
        text.append (self.g.results.describe().drop('time', axis=1))

        return text, quant

    def perform_audit(self):
        """Monitors modelled system at regular intervals (as defined by audit 
        interval in self.g)"""

        # Delay before first aurdit if length of warm-up
        yield self.env.timeout(self.g.warm_duration)
        
        # The trigger repeated audits
        while True:
            # Record time
            self.g.audit_time.append(self.env.now)
            if self.g.debug  and self.g.debuglevel>=1:
                print("")
                print(f"-- Audit Day {self.g.audit_time[-1]}")
                print(f"--- Patients waiting: First: {self.g.patients_waiting_by_priority[2]}, PIFU: {self.g.patients_waiting_by_priority[1]}, Traditional: {self.g.patients_waiting_by_priority[0]}")
                print(f"--- Slots used: {self.consultant.count}")
                print("--")
      
            ## alternative with save to file                        
            if self.g.loglinesave:
                ls_audit_to_add = [self.env.now, len(FOPA_Patient.all_patients), self.g.patients_waiting, self.g.patients_waiting_by_priority[0], self.g.patients_waiting_by_priority[1], self.g.patients_waiting_by_priority[2],self.consultant.count,self.repid]                                               
                with open(self.savepath +"batch_mon_audit_ls.csv", "a") as f:
                     writer = csv.writer(f, delimiter=",")
                     writer.writerow(ls_audit_to_add)
                     
            else:
                #Record patients waiting by referencing global variables
                self.g.audit_patients_waiting.append(self.g.patients_waiting)
                
                (self.g.audit_patients_waiting_p1.append
                  (self.g.patients_waiting_by_priority[0]))
                
                (self.g.audit_patients_waiting_p2.append
                  (self.g.patients_waiting_by_priority[1]))
                
                (self.g.audit_patients_waiting_p3.append
                  (self.g.patients_waiting_by_priority[2]))
                
                # Record patients waiting by asking length of dictionary of all patients 
                # (another way of doing things)
                self.g.audit_patients_in_system.append(len(FOPA_Patient.all_patients))
                # Record resources occupied (consultant)
                self.g.audit_resources_used.append(self.consultant.count)                
            
            # Trigger next audit after interval
            yield self.env.timeout(self.g.audit_interval)


    # The run method starts up the entity generators, and tells SimPy to start
    # running the environment for the duration specified in the g class. After
    # the simulation has run, it calls the methods that calculate run
    # results, and the method that writes these results to file
    def run(self,repid=1):
        
        
        # Start processes: entity generators and audit
        self.env.process(self.generate_wl_arrivals())
        
        if self.g.unavail_on:
            self.env.process(self.obstruct_slots())
        
        self.env.process(self.perform_audit())
        
        # Run simulation
        self.env.run(until=self.g.obs_duration + self.g.warm_duration)
        
        # End of simulation run. Build and save results.

        # Load Results log - audit
        if self.g.loglinesave:
            self.g.results = pd.read_csv(self.savepath + "batch_mon_audit_ls.csv") # read from csv
            self.g.results = (self.g.results[self.g.results['rep'] == self.repid])
        else:
            self.build_audit_results() # assemple from lists in memory
        
        # Load Results log - patient 
        if self.g.loglinesave:    
            self.results_df = pd.read_csv(self.savepath +"patient_result2.csv")
            
        # Load Results log - patient 
        if self.g.loglinesave:    
            self.g.appt_queuing_results = pd.read_csv(self.savepath +"appt_result.csv")
        
        # Calculate run results (aggregate). Run but deprecated in favour of batch methods
        self.calculate_mean_q_time()
        
        # Write run results to file. Run but deprecated in favour of batch methods
        self.write_run_results()

        # Get a chart of results
        chart_output = self.chart()
        
        # Get text summary of results
        [text_output,quant_output] = self.summarise()

        return chart_output, text_output, quant_output


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
        
        batch_mon_app_kpit = batch_mon_appointments.groupby(['rep','priority','interval']).quantile([0.50]).reset_index()
        
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


# Method to create files that will hold logs - RTT patient (1), appointment (2), audit (3)
def Trial_Results_initiate(file1,file2,file3):
    
    # Create a file to store trial results, and write the column headers
        with open(file1, "w") as f:
            writer = csv.writer(f, delimiter=",")
            column_headers = ["P_ID", "Q_time_fopa" , "Q_time_fuopa" , "rep"      
                           ]
            writer.writerow(column_headers)
        with open(file2, "w") as f:
            writer = csv.writer(f, delimiter=",")
            column_headers = ["P_ID", "Appt_ID","priority","type","pathway","q_time","start_q","DNA" , "rep"                   
                           ]
            writer.writerow(column_headers)
            
        with open(file3,"w") as f:
            writer = csv.writer(f, delimiter=",")
            column_headers = ['time','patients in system','all patients waiting','priority 1 patients waiting','priority 2 patients waiting','priority 3 patients waiting','resources occupied','rep']
            writer.writerow(column_headers)
            
        # with open(file4,"w") as f:
        #     writer = csv.writer(f, delimiter=",")
        #     column_headers = ['KPI','KPI_mean','KPI_LCI','KPI_UCI']
        #     writer.writerow(column_headers)
            
