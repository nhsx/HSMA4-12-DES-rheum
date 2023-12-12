import random
import pandas as pd
import numpy as np

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