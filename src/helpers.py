import scipy.stats as st
import numpy as np
import csv

# Code to compute (small sample) confidence interval, taken from web.
def mean_confidence_interval(data, confidence=0.95):
            a = 1.0 * np.array(data)
            n = len(a)
            m, se = np.mean(a), st.sem(a)
            h = se * st.t.ppf((1 + confidence) / 2., n-1)
            return m, m-h, m+h


class patient_blocker:

    def __init__(self,blockerid):        
        self.blockerid = blockerid


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