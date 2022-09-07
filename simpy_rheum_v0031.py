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

# Class to store global parameter values.  We don't create an instance of this
# class - we just refer to the class blueprint itself to access the numbers
# inside
class g:
    #wl_inter = 1.48 * 214 * 0.5 # 214 trusts  .inter-arrival time (minutes)
    #mean_FOPA = 15 # literature: 15 minutes per appointment (f/u?)
    #mean_FUOPA = 15 # literature: 15 minutes per appointment
    #mean_PIFUOPA = 15 # may be longer....
    mean_interOPA = np.round(4.5 * 365/12,0) # days inbetween appointments
    interOPA_tri = np.round(np.array([3,6,4.5])*365/12) # triangular - low, high , mode
    #mean_interOPA = np.round(365/12,0) # days inbetween appointments
    t_decision = 1 * 365 # days . Time mark for PIFU decision / stratification (from first OPA)
    
    obs_duration=365*3 # observation (days)
    
    debug=True
    debuglevel = 1

        
    def __init__(self,in_res=5,in_inter_arrival=1,in_prob_pifu=0,in_path_horizon_y=3,audit_interval=7,in_reps=1,repid=1,savepath='temp',in_FOavoidable=0,in_interfu_perc=0.6):
        
        self.prob_firstonly = 0.35 # ~35% of rheumatology patients have no follow-ups
        self.number_of_runs=in_reps
        self.prob_pifu=in_prob_pifu # probability of being PIFU patient
        self.wl_inter = in_inter_arrival
        self.number_of_slots = in_res # slots
        self.in_FOavoidable = in_FOavoidable
        
        self.warm_duration=365*5 # warm-up (days)
        self.max_fuopa_tenor_y = in_path_horizon_y # years
        self.max_fuopa_tenor = in_path_horizon_y * 365 # days
        self.appt_counter = 0
        self.DNA_pifu_pro=0.07
        self.DNA_tra_pro=0.077
        self.interfu_perc = in_interfu_perc
        self.mean_interPIFU = np.round(np.sum(self.interOPA_tri)/3 * (1+self.interfu_perc),0) # 
        self.PIFUbigbang = False # whether, when PIFU starts being used, it is offered to all eligible patients when they visit (e.g. those already followed up for years) - big-bang - or only new eligible patients

        self.unavail_on = False
        self.unavail_byshock = False # for now model only for unavailability by shock OR by periodic. Can be improved in future.        
        # Parameters for a shock to available resource (one-off)
        self.unavail_shock_tmin = 365 * 3
        self.unavail_shock_period = np.floor(365/12 * 12) # e.g. from half-March (emergency response) to half-August (letter on Aug2020 for NHS response)
        self.unavail_shock_nrslots = int(np.floor(in_res*1))
        
        # Parameters for periodic resource unavailability (e.g. weekends, leave) - if unavail_byshock = False
        self.unavail_slot=2
        self.unavail_freq_slot=5
        self.unavail_nrslots = int(np.floor(in_res*1)) # 100% of slots unavailable
        
        self.savepath = savepath
        self.repid = repid
        self.loglinesave = True # if true saves each line to file, if false creates dataframe that stays in memory
        self.audit_time = []
        self.audit_interval = audit_interval # time step for audit (days)
        self.audit_patients_waiting = [] # populated in perform_audit
        self.audit_patients_waiting_p1 = [] # populated in perform_audit
        self.audit_patients_waiting_p2 = [] # populated in perform_audit
        self.audit_patients_waiting_p3 = [] # populated in perform_audit
        self.audit_patients_in_system = [] # populated in perform_audit (lifetime to discharge...)
        self.audit_resources_used = [] # populated in perform_audit
        self.patients_waiting = 0 # counter
        self.patients_waiting_by_priority = [0,0,0] # counter
        self.results = pd.DataFrame() # populated in build_audit_results. audit results in one DF
        self.appt_queuing_results = (pd.DataFrame(
            columns=['P_ID','App_ID','priority','type',"pathway",'q_time','start_q','DNA'])) # populated in pathway
        self.appt_queuing_results.set_index("App_ID", inplace=True)
        
    def change_reps(self,reps):
        self.number_of_runs = reps


class patient_blocker:

    def __init__(self,blockerid):        
        self.blockerid = blockerid
    
# Class representing our patients coming in for the weight loss clinic.
# This time we've added another attribute, that will store the calculated
# queuing time for the nurse for each patient (each instance of this class)
class FOPA_Patient:

    # The following dictionaries store patients
    all_patients = {}

    def __init__(self, p_id,prob_pifu,in_path_horizon,DNA_pifu_pro,DNA_tra_pro):
        self.apptype = []
        self.id = p_id
        self.prob_pifu = prob_pifu
        self.q_time_fopa = float("nan")
        self.q_time_fuopa = float("nan")
        self.q_times_fuopa = []
        self.q_time_nurse=float("nan")
        self.topifu = False
        self.priority = 3 # start with lower priority
        self.max_fuopa = 999 # placeholder maximum follow up outpatients
        self.max_fuopa_tenor = in_path_horizon # max horizon for follow-up (as days since first OPA)
        self.used_fuopa = 0
        self.ls_appt = []
        self.DNA_pifu_pro=DNA_pifu_pro
        self.DNA_tra_pro=DNA_tra_pro
        self.tradition_dna=False # instantiate DNA for traditional type
        self.pifu_dna=False # instantiated DNA for PIFU type
        self.FOavoided = False
        self.df_appt_to_add = (pd.DataFrame(
            columns=['P_ID','App_ID','priority','type',"pathway",'q_time','start_q']))
        self.ls_appt_to_add = []
        self.ls_patient_to_add = []
        self.RTT_sub = 0
        
    def triage_decision(self):
        if random.random() < self.prob_pifu:
            self.topifu = True
            self.type = "PIFU"
        else:
            self.type = "TFU"
            #self.priority = 2 # assign 2nd highest priority in the sense that those on traditional already are scheduled in/planned so not really 'moveable'
            
    def give_pifu_priority(self):
        self.priority =  2 # assign 2nd highest priority in the sense that those on traditional already are scheduled in/planned so not really 'moveable'
        
    def give_tfu_priority(self):
        self.priority = 1
        
    def assign_firstonly(self,prob_firstonly):
        if random.random() < prob_firstonly:
            self.type = "First-only" # single appointment
            self.apptype = "First-only"
            #self.type = "First" # single appointment
            self.max_fuopa_tenor = 0 # no follow-ups
        else:
           self.type = "First" # first leading to long-term follow-up
           self.apptype = "First"
            
            
    def decision_DNA_pifu(self):
        if  random.random() < self.DNA_pifu_pro:
            self.pifu_dna = True
             #self.priority = 2 # assign 2nd highest priority in the sense that those on traditional already are scheduled in/planned so not really 'moveable'
    def decision_DNA_tradtion(self):
        if  random.random() < self.DNA_tra_pro:
            self.tradition_dna = True 
            
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
        
    def avoidable_firstonly(self,in_FOavoidable):
        # A&G, whether a first-only appointment can be avoided
        if random.random() < in_FOavoidable:            
            self.FOavoided = True
        
# Class representing our model of the rheumatology outpatient clinic
class rheum_Model:
    # Here, the constructor sets up the SimPy environment, sets a patient
    # counter to 0 (which we'll use for assigning patient IDs), and sets up
    # our resources (here just a nurse resource, with capacity given by
    # the number stored in the g class)
    # We also now set up a Pandas DataFrame to store run results, an
    # attribute storing the run number, which gets passed in to the instance
    # of the class when it's instantiated, and an attribute to store the mean
    # queuing time for the nurse across patients in this run of the model.
    def __init__(self, run_number, in_res=2 , in_inter_arrival=(365/4590), in_prob_pifu=0.6, in_path_horizon_y=3,audit_interval=1,repid=1, savepath='temp',in_FOavoidable=0,in_interfu_perc=0.6):
        self.env = simpy.Environment() # instance of environment

        self.g = g(in_res,in_inter_arrival,in_prob_pifu, in_path_horizon_y, audit_interval, repid = repid, in_FOavoidable = in_FOavoidable,in_interfu_perc=in_interfu_perc) # instance of global variables
        
        self.patient_counter = 0
        self.block_counter = 0
        
        self.consultant = simpy.PriorityResource(self.env, capacity=self.g.number_of_slots)

        self.run_number = run_number
        self.savepath = savepath
        
        self.mean_q_time_total = pd.DataFrame()
        self.results_df = pd.DataFrame()
        self.results_df["P_ID"] = []
        self.results_df["Q_Time_fopa"] = []
        self.results_df["Q_Time_fuopa"] = []
        self.results_df.set_index("P_ID", inplace=True)
        
    # A method that generates patients arriving for the weight loss clinic
    def generate_wl_arrivals(self):
        # Keep generating indefinitely (until the simulation ends)
        while True:
            # Increment the patient counter by 1
            self.patient_counter += 1
            
            # Create a new patient - an instance of the Weight_Loss_Patient 
            # class, and give the patient an ID determined by the patient
            # counter
            wp = FOPA_Patient(self.patient_counter,self.g.prob_pifu,self.g.max_fuopa_tenor, self.g.DNA_pifu_pro,self.g.DNA_tra_pro)
            
            # Add patient to dictionary of patients
            FOPA_Patient.all_patients[wp.id] = wp
            
            # Get the SimPy environment to run the attend_wl_clinic method
            # with this patient
            self.env.process(self.attend_OPA(wp))
            
            # Randomly sample the time to the next patient arriving for the
            # weight loss clinic.  The mean is stored in the g class.
            #sampled_interarrival = int(np.round(random.expovariate(1.0 / self.g.wl_inter),0))
            sampled_interarrival = random.expovariate(1.0 / self.g.wl_inter)
            
            # Freeze this function until that time has elapsed
            yield self.env.timeout(sampled_interarrival)


    

    def obstruct_slot(self,patient_blocker,unavail_timeperiod):
            
        # Once this time has elapsed, request an ed doctor with a priority
        # of -1 (so that we know this will get the top priority, as none
        # of our patients will have a negative priority), and hold them
        # for the specified unavailability amount of time
        with self.consultant.request(priority=-1) as req:
            # Freeze the function until the request can be met (this
            # ensures that the doctor will finish what they're doing first)
            yield req
            
            yield self.env.timeout(unavail_timeperiod)
                
                
    def obstruct_slots(self):
        
        if self.g.unavail_byshock:
            
            yield self.env.timeout(self.g.unavail_shock_tmin)
            
            for i in range(self.g.unavail_shock_nrslots):
                
                self.block_counter += 1
                
                slot_block = patient_blocker(self.block_counter)
                
                # Get the SimPy environment to run the attend_obstruct_slot method with this slot block
                self.env.process(self.obstruct_slot(slot_block,self.g.unavail_shock_period))
                
                #self.obstruct_slot(self.g.unavail_shock_period)
                
                if self.g.debug and self.g.debuglevel>=3:
                    print ("Appointment will not be able to book at",
                           f"{self.env.now + self.g.unavail_freq_slot:.1f}")
                             
        else:
            
            while True:        
                for i in range(self.g.unavail_nrslots):
                    
                    self.block_counter += 1
                    slot_block = patient_blocker(self.block_counter)
                    
                    # Get the SimPy environment to run the attend_obstruct_slot method with this slot block
                    self.env.process(self.obstruct_slot(slot_block,self.g.unavail_slot))
                
                    if self.g.debug and self.g.debuglevel>=3:
                        print ("Appointment will not be able to book at",
                               f"{self.env.now + self.g.unavail_freq_slot:.1f}")
                
                # Freeze the function for the unavailability frequency period.
                # This could represent time on shift (the time the doctor is
                # available).  Here, we use a fixed value to reflect a set shift,
                # but you could also sample from a probability distribution.
                yield self.env.timeout(self.g.unavail_freq_slot)
                
            
    # A method that models the processes for attending the outpatient rheumatology clinic.
    # The method needs to be passed a patient who will go through these processes
    def attend_OPA(self, patient):
                        
        patient.assign_firstonly(self.g.prob_firstonly)        
        patient.sub_RTT_priority() # add some variability to priority within RTT queue (0 or 1 increment to its '3' priority)        
        patient.avoidable_firstonly(self.g.in_FOavoidable) # if 'first-only', whether it is avoided or not (e.g. A&G)
        
        if patient.type == "First-only" and patient.FOavoided and self.env.now > self.g.warm_duration:
            # if first-only pathway and avoidable through A&G and current day within intervention/study period, skip anything further for patient
            if self.g.debug and self.g.debuglevel>=2:
                print(f"Patient {patient.id} had first outpatient avoided. Not added to log.")
                
        else:            

                        
            #print(f"RTT increment:{patient.RTT_sub}")
            ############################################
            ### First outpatient appointment ####
            ############################################        
            
            # Record the time the patient started queuing for the first outpatient
            start_q_fopa = self.env.now
            self.g.appt_counter +=1
            patient.ls_appt.append(self.g.appt_counter)
            self.g.patients_waiting += 1
            self.g.patients_waiting_by_priority[patient.priority-1] += 1
            
            # Request a consultant
            with self.consultant.request(priority = patient.priority + patient.RTT_sub) as req:
                # Freeze the function until the request for a consultant can be met
                yield req
                
                
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
                self.g.patients_waiting -= 1
                
                # Record the time the patient finished queuing for a consultant
                end_q_fopa = self.env.now
    
                 # Calculate the time this patient spent queuing for the consultant (FU) and
                # store in the patient's attribute
                patient.q_time_fopa = end_q_fopa - start_q_fopa
                
                if self.g.debug and self.g.debuglevel>=2:
                    print(f"Req {patient.ls_appt[-1]}: Patient {patient.id} queued {np.round(patient.q_time_fopa,2)} days for 1st app. Priority {patient.priority}")
                
                # Freeze this function until the day time unit has elapsed
                #yield self.env.timeout(1) # freeze for one time-unit (a day) - that same slot will only be available the next day
                yield self.env.timeout(2) # freeze for two time-units (two dayz) - that same slot will only be available in two days
                
                if patient.tradition_dna==True:
                    if self.g.debug and self.g.debuglevel>=2:
                        print(f" Patient {patient.id} queued {np.round(patient.q_time_fopa,2)} days for 1st app. Priority {patient.priority} but didn't attend the appointment")
    
                else:
                    if self.g.debug and self.g.debuglevel>=2:
                        print(f" Patient {patient.id} queued {np.round(patient.q_time_fopa,2)} days for 1st app. Priority {patient.priority}")
                
                
                patient.decision_DNA_tradtion()
    
    
                # Add to appointment log or save            
                patient.ls_appt_to_add = [patient.id, patient.ls_appt[-1],patient.priority,patient.apptype,patient.type,patient.q_time_fopa,start_q_fopa,patient.tradition_dna,self.g.repid]
                patient.ls_patient_to_add = [patient.id , patient.q_time_fopa,999,self.g.repid]
                            
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
                                                          data=[patient.ls_appt_to_add])
                    patient.df_appt_to_add.set_index("App_ID", inplace=True)
                    
                    self.g.appt_queuing_results=self.g.appt_queuing_results.append(patient.df_appt_to_add)
                    
                    df_to_add = pd.DataFrame( columns = ["P_ID","Q_time_fopa","Q_time_fuopa","rep"], data =[patient.ls_patient_to_add])

                    df_to_add.set_index("P_ID", inplace=True)
                                    
                    if start_q_fopa > self.g.warm_duration: # don't save things in warm-up period           
                        self.results_df = self.results_df.append(df_to_add)
                        
            # time between first and first follow-up needed    
            #patient.priority = 1 # priority should now be max (in the sense that appointments are planned / scheduled)
            patient.give_tfu_priority()
            
            ############################################
            ### Traditional Follow-up appointments ####
            ############################################
            
            while self.env.now - end_q_fopa < patient.max_fuopa_tenor:
            #while patient.used_fuopa < patient.max_fuopa:
                    
                
                
                # Determine time till next needing F/U
                #sampled_interfu_duration = int(random.expovariate(1.0 / g.mean_interOPA)) # integer only (days)
                sampled_interfu_duration = int(random.triangular(g.interOPA_tri[0],g.interOPA_tri[1],g.interOPA_tri[2]))
                # Freeze this function until time has elapsed
                yield self.env.timeout(sampled_interfu_duration)
                
                self.g.patients_waiting += 1
                self.g.patients_waiting_by_priority[patient.priority-1] += 1
                patient.used_fuopa+=1 # count the follow-up outpatient
                
                self.g.appt_counter +=1
                patient.ls_appt.append(self.g.appt_counter)
                start_q_fuopa = self.env.now
                # Request a consultant for follow-up
                with self.consultant.request(priority = patient.priority) as req:
                    # Freeze the function until the request for a consultant can be met
                    yield req
                                   
                    # reduce patients waiting counts
                    self.g.patients_waiting_by_priority[patient.priority-1] -= 1 # decrement
                    self.g.patients_waiting -= 1
    
                    # Record the time the patient finished queuing for a nurse
                    end_q_fuopa = self.env.now
                                   
                    # Calculate the time this patient spent queuing for a consultant and
                    # store in the patient's attribute
                    patient.q_times_fuopa.append(end_q_fuopa - start_q_fuopa) # add latest followup
                    #print(patient.q_time_fuopa)
                    
                    if patient.used_fuopa ==1 :
                        patient.q_time_fuopa = end_q_fuopa - start_q_fuopa
                    
                    patient.decision_DNA_tradtion()
                    if patient.tradition_dna==True:
                        if self.g.debug  and self.g.debuglevel>=2:
                            print(f"Req {patient.ls_appt[-1]}: Patient {patient.id} queued {np.round(end_q_fuopa - start_q_fuopa,2)} for app {patient.used_fuopa}. Priority {patient.priority}")
                    else:
                        if self.g.debug  and self.g.debuglevel>=2:
                            print(f"Req {patient.ls_appt[-1]}: Patient {patient.id} queued {np.round(end_q_fuopa - start_q_fuopa,2)} for app {patient.used_fuopa}. Priority {patient.priority} but DNAd")
                    
                            
                    # Freeze this function until the day time unit has elapsed
                    yield self.env.timeout(1) # freeze for one time-unit (a day) - that same slot will only be available the next day
                    
                    
                    
                    # Add to appointment log or save            
                    patient.ls_appt_to_add = [patient.id, patient.ls_appt[-1],patient.priority,"Traditional",patient.type,patient.q_time_fuopa,start_q_fuopa,patient.tradition_dna,self.g.repid]                            
                    if self.g.loglinesave:                    
                        with open(self.savepath +"appt_result.csv", "a") as f:
                             writer = csv.writer(f, delimiter=",")
                             writer.writerow(patient.ls_appt_to_add)
    
                    else:                         
                        patient.df_appt_to_add = pd.DataFrame( columns = ["P_ID","App_ID","priority","type","pathway","q_time","start_q","DNA","rep"] ,
                                                              data=[patient.ls_appt_to_add])
                        patient.df_appt_to_add.set_index("App_ID", inplace=True)
                        self.g.appt_queuing_results=self.g.appt_queuing_results.append(patient.df_appt_to_add)
    
                
                
                if self.env.now - end_q_fopa > self.g.t_decision  and self.env.now > self.g.warm_duration:
                    if patient.topifu:
                        break
                
    
                
            if patient.topifu:
                
                if self.g.debug  and self.g.debuglevel>=2:
                    print(f"Patient {patient.id} PIFU. Follows {patient.used_fuopa} traditional apps.")
                
                patient.give_pifu_priority() # intermediate
                #print(f"Patient {patient.id} entered PIFU. Has {patient.used_fuopa} traditional apps. Priority {patient.priority}")
                
                while patient.topifu: # while true
                    
                    patient.used_fuopa+=1 # count the follow-up outpatient
                    # Determine time till next needing PIF/U
                    sampled_interpifu_duration = int(np.round(random.expovariate(1.0 / self.g.mean_interPIFU),0)) # integer only (days)
                    # Freeze this function until time has elapsed
                    yield self.env.timeout(sampled_interpifu_duration)
                    
                    self.g.appt_counter +=1
                    patient.ls_appt.append(self.g.appt_counter)
                    self.g.patients_waiting += 1
                    self.g.patients_waiting_by_priority[patient.priority-1] += 1
    
    
                    start_q_pifuopa = self.env.now            
                    # Request consultant
                    with self.consultant.request(priority = patient.priority) as req:
                        # Freeze the function until the request for a consultant can be met
                        yield req
    
                        # reduce patients waiting counts
                        self.g.patients_waiting_by_priority[patient.priority-1] -= 1 # decrement
                        self.g.patients_waiting -= 1
                        
                        # Record the time the patient finished queuing for a consultant
                        end_q_pifuopa = self.env.now
                        
                        # Calculate the time this patient spent queuing for a consultant and
                        # store in the patient's attribute
                        patient.q_time_pifuopa = end_q_pifuopa - start_q_pifuopa
                        
                        # Freeze this function until the day time unit has elapsed
                        yield self.env.timeout(1) # freeze for one time-unit (a day) - that same slot will only be available the next day
                    
                        patient.decision_DNA_pifu()
                        if patient.pifu_dna==True:
                            if self.g.debug  and self.g.debuglevel>=2:
                                print(f"Req {patient.ls_appt[-1]}: Patient {patient.id} queued {np.round(patient.q_time_pifuopa,2)} days for app {patient.used_fuopa} - PIFU. Priority {patient.priority} but did not attend")
                        else:
                             if self.g.debug  and self.g.debuglevel>=2:
                                print(f"Req {patient.ls_appt[-1]}: Patient {patient.id} queued {np.round(patient.q_time_pifuopa,2)} days for app {patient.used_fuopa} - PIFU. Priority {patient.priority}")
                        
                            
                        # Add to appointment log or save            
                        patient.ls_appt_to_add = [patient.id, patient.ls_appt[-1],patient.priority,"PIFU",patient.type,end_q_pifuopa - start_q_pifuopa,start_q_pifuopa,patient.pifu_dna,self.g.repid]                            
                        if self.g.loglinesave:                    
                            with open(self.savepath +"appt_result.csv", "a") as f:
                                 writer = csv.writer(f, delimiter=",")
                                 writer.writerow(patient.ls_appt_to_add) 
                        
                        else:
                            patient.df_appt_to_add = pd.DataFrame( columns = ["P_ID","App_ID","priority","type","pathway","q_time","start_q","DNA","rep"] ,
                                                                  data=[patient.ls_appt_to_add])
                            patient.df_appt_to_add.set_index("App_ID", inplace=True)
                            self.g.appt_queuing_results=self.g.appt_queuing_results.append(patient.df_appt_to_add)
                    
                        
                    if self.env.now - end_q_fopa > patient.max_fuopa_tenor:
                        break
             
            else:
                
                if self.g.debug  and self.g.debuglevel>=2:
                    print(f"Patient {patient.id} not PIFU. Follows {patient.used_fuopa} traditional apps.")
            
            # Store the start and end queue times alongside the patient ID in
            # the Pandas DataFrame of the GP_Surgery_Model class
            
            # Add appts to model log, only if greater than warm up
            #self.g.appt_queuing_results = self.g.appt_queuing_results.append(patient.df_appt_to_add[patient.df_appt_to_add['start_q']>self.g.warm_duration])
            
    
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
        
        self.g.results['rep']=self.g.repid

    # A method that calculates the average quueing time for the nurse.  We can
    # call this at the end of each run
    def calculate_mean_q_time(self):
        
        self.results_df['Q_Time_total']=self.results_df.sum(axis=1)
        self.mean_q_time_total = self.results_df.mean(axis=0)
        
    # A method to write run results to file.  Here, we write the run number
    # against the the calculated mean queuing time for the nurse across
    # patients in the run.  Again, we can call this at the end of each run
    def write_run_results(self):
        with open(self.savepath +"trial_results.csv", "a") as f:
            writer = csv.writer(f, delimiter=",")
            results_to_write = [self.run_number,
                                self.mean_q_time_total['Q_time_fopa'],
                                self.mean_q_time_total['Q_time_fuopa'],
                                #self.mean_q_time_total['Q_Time_eddoc'],
                                #self.mean_q_time_total['Q_Time_acudoc'],
                                self.mean_q_time_total['Q_Time_total']]
            print(results_to_write)
            writer.writerow(results_to_write)
    

    # plot results - one run
    def chart(self):
        # plot results at end of run #
        self.results_df = self.results_df.sort_index()
        #self.g.appt_queuing_results = self.g.appt_queuing_results.sort_values(by=['start_q'])
        self.g.appt_queuing_results = self.g.appt_queuing_results.sort_index()

        # Define fig size
        fig = plt.figure(figsize=(12,12))
        
        # Figure 1: first opa patient level results
        ax1 = fig.add_subplot(3,4,1) # 1 row, 4 cols, pos 1
        x = self.results_df.index
        y = self.results_df['Q_time_fopa']

        ax1.plot(x,y,marker='^',color='k')
        ax1.set_xlabel('Patient')
        ax1.set_ylabel('Queuing time')
        ax1.legend()
        ax1.set_title('Queuing time for first outpatient')
        ax1.grid(True,which='both',lw=1,ls='',c='.75')

        # Figure 2: appointment level results
        ax22 = fig.add_subplot(342)  # 1 row, 4 cols, chart position 2
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

        # Figure 3: System level queuing results # change to bar chart?
        ax2 = fig.add_subplot(343)
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


        # Figure 4: Staff usage
        ax3 = fig.add_subplot(344)  # 1 row, 4 cols, chart position 4
        x = self.g.results['time']
        y = self.g.results['resources occupied']
        ax3.plot(x, y, label='Resources occupied')
        ax3.set_xlabel('Time')
        ax3.set_ylabel('Resources occupied')
        ax3.set_title('Resources occupied')
        ax3.grid(True, which='both', lw=1, ls='--', c='.75')

        # Figure 5-7: appointment level results
        mylabels = ['Traditional F/U appt waiting','PIFU appt waiting','First outpatient appt waiting','All']
        markers = ['o', 'x', '^']
        colors = ['r','g','k']
        for priority in range(1, 4):
            ax22 = fig.add_subplot(3,4,4 + priority)  # 1 row, 4 cols, chart position 2
            
            x = self.g.appt_queuing_results.index
            # Chart loops through 3 priorites
            
            x = (self.g.appt_queuing_results[self.g.appt_queuing_results['priority'] == priority].index) 
            y = (self.g.appt_queuing_results
                 [self.g.appt_queuing_results['priority'] == priority]['q_time'])
            
            ax22.plot(x, y, marker=markers[priority - 1], linestyle="none",label=mylabels[priority-1],color = colors[priority - 1])
            ax22.set_xlabel('Appointment nr')
            ax22.set_ylabel('Queuing time')
            ax22.legend()
            ax22.set_title('Queuing time by type')
            ax22.grid(True, which='both', lw=1, ls='--', c='.75')

        myvars = ['priority 1 patients waiting','priority 2 patients waiting','priority 3 patients waiting']
        for priority in range(1, 4):
            ax32 = fig.add_subplot(3,4,8 + priority)  # 1 row, 4 cols, chart position 2
            
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

    # summaries for streamlit (quartiles)
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
                ls_audit_to_add = [self.env.now, len(FOPA_Patient.all_patients), self.g.patients_waiting, self.g.patients_waiting_by_priority[0], self.g.patients_waiting_by_priority[1], self.g.patients_waiting_by_priority[2],self.consultant.count,self.g.repid]                                               
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

        if self.g.loglinesave:
            self.g.results = pd.read_csv(self.savepath + "batch_mon_audit_ls.csv")
            self.g.results = (self.g.results[self.g.results['rep'] == self.g.repid])
        else:
            self.build_audit_results()
        
        if self.g.loglinesave:    
            self.results_df = pd.read_csv(self.savepath +"patient_result2.csv")
        
        # Calculate run results (aggregate)
        self.calculate_mean_q_time()
        
        # Write run results to file
        self.write_run_results()

        # Get a chart of results
        chart_output = self.chart()
        
        # Get text summary of results
        [text_output,quant_output] = self.summarise()

        return chart_output, text_output, quant_output


# Class for batch runs
class Batch_rheum_model:
    
    def __init__(self,in_res=5 , in_inter_arrival=1, in_prob_pifu=0.6, in_path_horizon_y=3,audit_interval=7,in_savepath="temp/",in_FOavoidable=0,in_interfu_perc=0.6):
        
        self.batch_mon_appointments = pd.DataFrame()
        self.batch_mon_audit = pd.DataFrame()
        self.batch_mon_app_kpit = pd.DataFrame()
        self.batch_kpi = pd.DataFrame()
        self.savepath=in_savepath
        self.g = g(in_res,in_inter_arrival,in_prob_pifu, in_path_horizon_y, audit_interval,in_FOavoidable=in_FOavoidable,in_interfu_perc=in_interfu_perc) # instance of global variables
        
    
    
    def read_logs_to_self(self):        
        # Read in results back into self dataframe
        self.trial_results_df = pd.read_csv(self.savepath +"patient_result2.csv")
        self.batch_mon_appointments = pd.read_csv(self.savepath +"appt_result.csv")
        self.batch_mon_audit = pd.read_csv(self.savepath + "batch_mon_audit_ls.csv")
    
    
    # Plotting an overview of behaviour at audit timepoints (across reps)    
    def plot_audit_reps(self):
        
        t_warm = self.g.warm_duration
        
        fig_q = plt.figure(figsize=(12,12))
        sns.boxplot(x='type',y='q_time',data=self.batch_mon_appointments[self.batch_mon_appointments['start_q']>t_warm],hue='type')

        # Audit plot
        audit_kpis = ['priority 3 patients waiting','priority 1 patients waiting','priority 2 patients waiting','resources occupied']
        pricolors = ['orange','skyblue','green','gray']
        audit_kpis_name = ['RTT appointment waits','Traditional appointment waits','PIFU appointment waits','Resources (slots) occupied']
        batch_mon_audit_melted = pd.melt(self.batch_mon_audit,id_vars=['time','rep'],value_vars=audit_kpis,var_name='audit_KPI')
        sns.boxplot(x='time',y='priority 1 patients waiting',data=self.batch_mon_audit)        
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
        
        
        # fig2, axes2 = plt.subplots(figsize=(20, 20))
        # ax2=sns.boxplot(ax=axes2,x='interval',y='q_time',data=batch_mon_app_kpit[batch_mon_app_kpit['KPI']=='Seen'],hue='priority')
        # ax2.set_xticklabels(ax2.get_xticklabels(),rotation = 80)
        # ax2.set_xlabel("Simulation time interval (days)", fontsize = 14)
        # ax2.set_ylabel("Number seen", fontsize = 20)
        
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
        
    
    def run_reps(self,reps):
        
        for run in range(reps):
            
            if self.g.debug  and self.g.debuglevel>=0:
                print (f"Run {run+1} of {reps}")
            
            
            
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
            
            #p = rheum_Model.chart(my_ed_model)
            e_results = my_ed_model.g.results # audit counts . patients waiting per time
            self.batch_mon_audit = pd.concat([self.batch_mon_audit, e_results])
            
            #e_resultsdf = my_ed_model.results_df # initial df with time each patient waited for first OPA

            if self.g.loglinesave:
                self.read_logs_to_self()
            else:
                e_appt_queuing_result= my_ed_model.g.appt_queuing_results # appointments monitor
                                
                #e_appt_queuing_result["rep"] = run
                self.batch_mon_appointments = pd.concat([self.batch_mon_appointments,e_appt_queuing_result])
             
            # keep only post warm-up queue starts (rather q finish????)
            self.batch_mon_appointments = self.batch_mon_appointments[self.batch_mon_appointments['start_q']>self.g.warm_duration]
            
            
        fig_audit_reps, fig_q_audit_reps = self.plot_audit_reps()
        
        fig_monappKPI_reps, fig_monappKPIn_reps = self.plot_monappKPI_reps(step=365/4)
        
        self.headline_KPI(365)
        
        return fig_audit_reps, chart_output_lastrep, text_output_lastrep, quant_output_lastrep, fig_q_audit_reps, fig_monappKPI_reps, fig_monappKPIn_reps
        
        
    def save_logs(self):
        
       self.batch_mon_appointments.to_csv(self.savepath + 'batch_mon_appointments.csv')
       self.batch_mon_audit.to_csv(self.savepath + 'batch_mon_audit.csv')
       self.batch_mon_app_kpit.to_csv(self.savepath + 'batch_mon_app_kpit.csv')
       self.batch_kpi.to_csv(self.savepath + 'batch_kpi.csv')
        #my_batch_model.plot_monappKPI_reps(step=28*2)



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
            
            

# Class to store, calculate and manipulate trial results in a Pandas DataFrame
class Trial_Results_Calculator:
    # The constructor creates a new Pandas DataFrame, and stores this as an
    # attribute of the class instance
    def __init__(self,savepath):
        self.trial_results_df = pd.DataFrame()
        
        
        
    # A method to read in the trial results (that we wrote out elsewhere in the
    # code) and print them for the user
    def print_trial_results(self,savepath):
        print ("TRIAL RESULTS")
        print ("-------------")
        
        # Read in results from each run into our DataFrame
        self.trial_results_df = pd.read_csv(savepath + "trial_results.csv")
        
        # Take average over runs
        trial_mean_q_time_total = (
            self.trial_results_df["Mean_Q_Time_Total"].mean())
        trial_mean_q_time_fopa = (
            self.trial_results_df["Mean_Q_Time_FOPA"].mean())
        
        print ("Mean Time Queueing over Trial :",
               f"{trial_mean_q_time_total:.2f}")
        print ("Mean Time Queueing for First Outpatient :",
               f"{trial_mean_q_time_fopa:.2f}")



