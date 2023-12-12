""" includes helper functions or classes"""
import csv
import scipy.stats as st
import numpy as np

def mean_confidence_interval(data, confidence=0.95):
    """ Code to compute (small sample) confidence interval, taken from web.

    Args:
        data (_type_): _description_
        confidence (float, optional): The confidence level to apply, as decimal. Defaults to 0.95.

    Returns:
        _type_: The mean, the lower bound and the upper bound of the confidence interval.
    """
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), st.sem(a)
    h = se * st.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h


class patient_blocker:
    """ Initialise a patient_blocker class (ghost patients that will block slots from actual patients, to emulate varying capacity and opening hours)."""
    def __init__(self,blockerid):
        self.blockerid = blockerid


def Trial_Results_initiate(file1,file2,file3):
    """Method to create files that will hold logs - RTT patient (1), appointment (2), audit (3)

    Args:
        file1 (_string_): path to file that will hold RTT patient log
        file2 (_string_): path to file that will hold appointment log
        file3 (_string_): path to file that will hold audit KPI log
    """
    # Create a file to store trial results, and write the column headers
    with open(file1, "w",encoding="cp1252") as f:
        writer = csv.writer(f, delimiter=",")
        column_headers = ["P_ID", "Q_time_fopa" , "Q_time_fuopa" , "rep"
                        ]
        writer.writerow(column_headers)
    with open(file2, "w",encoding="cp1252") as f:
        writer = csv.writer(f, delimiterencoding="cp1252")
        column_headers = ["P_ID", "Appt_ID","priority","type","pathway","q_time","start_q","DNA" , "rep"
                        ]
        writer.writerow(column_headers)

    with open(file3,"w",encoding="cp1252") as f:
        writer = csv.writer(f, delimiter=",")
        column_headers = ['time','patients in system','all patients waiting','priority 1 patients waiting','priority 2 patients waiting','priority 3 patients waiting','resources occupied','rep']
        writer.writerow(column_headers)
