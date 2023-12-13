""" Module includes model (single replication) and utilities for plotting and saving"""
import csv
import random
import simpy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from src.patient import FOPA_Patient
from src.helpers import patient_blocker
from src.initialisers import g


class rheum_Model:
    """Class representing our overall model of the rheumatology outpatient clinic.

    # Here, the constructor sets up the SimPy environment, sets a patient
    # counter to 0 (which we'll use for assigning patient IDs), and sets up
    # our resources (here appointment slot units (symbolycally 15 min), with capacity given by
    # the number stored in the g class)
    """

    def __init__(self, run_number, in_res=2 , in_inter_arrival=(365/4590), in_prob_pifu=0.6, in_path_horizon_y=3,audit_interval=1,repid=1, savepath='temp',in_FOavoidable=0,in_interfu_perc=0.6):
        """Initialise rhematology outpatient clinic model.

        Args:
            run_number (_type_): The present run number.
            in_res (int, optional): The number of daily slots [slots]. Defaults to 2.
            in_inter_arrival (tuple, optional): The inter-arrival time [days]. Defaults to (365/4590).
            in_prob_pifu (float, optional): PIFU proportion - probability of PIFU pathway for non first-only pathways[%]. Defaults to 0.6.
            in_path_horizon_y (int, optional): Patient follow-up horizon [years], simplification on how long each non first-only pathway lasts (years). Defaults to 3.
            audit_interval (int, optional): Time step for audit metrics [simulation days]. Defaults to 1.
            repid (int, optional): id of current replication (within batch). Defaults to 1.
            savepath (str, optional): Save path for outputs. Defaults to 'temp'.
            in_FOavoidable (float, optional): A&G proportion - proportion of first-only pathways avoidable via A&G [%]. Defaults to 0.
            in_interfu_perc (float, optional): Percentage increase in inter-appointment interval with PIFU (vs traditional), i.e. 0.6 means 60% longer interval. Defaults to 0.6.
        """
        self.env = simpy.Environment() # instance of environment

        self.g = g(in_res,in_inter_arrival,in_prob_pifu, in_path_horizon_y, audit_interval, repid = repid, in_FOavoidable = in_FOavoidable,in_interfu_perc=in_interfu_perc) # instance of global variables for this replication

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

    def generate_wl_arrivals(self):
        """A method that generates patients arriving for the RTT outpatient 'clinic'"""

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


    def obstruct_slot(self,unavail_timeperiod):
        """  A method to obstruct a single slot (emulate unavailability)

        Args:
            unavail_timeperiod (_type_): Time that the slot is obstructed.

        Yields:
            _type_: A given timed out period where the slot is unavailable
        """

        # Once this time has elapsed, request a slot wiyh priority
        # of -1 (so that we know this will get the top priority, as none
        # of our pathway appointment requests will have a negative priority), and hold them
        # for the specified unavailability amount of time
        with self.consultant.request(priority=-1) as req:
            # Freeze the function until the request can be met (this
            # ensures that the slot will finish before becoming unavailabe)
            yield req

            yield self.env.timeout(unavail_timeperiod)

    def obstruct_slots(self):
        """ A method to obstruct multiple slots (emulate unavailability)"""

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

    def attend_OPA(self, patient):
        """    A method that models the processes / RTT patient pathway for attending the outpatient rheumatology clinic.

        The method needs to be passed an RTT patient who will go through these processes

        Args:
            patient (_FOPA_Patient class_): An instantiated object of class FOPA_Patient
        """

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

                if patient.tradition_dna:
                    if self.g.debug and self.g.debuglevel>=2:
                        print(f" Patient {patient.id} queued {np.round(patient.q_time_fopa,2)} days for 1st app. Priority {patient.priority} but didn't attend the appointment")

                else:
                    if self.g.debug and self.g.debuglevel>=2:
                        print(f" Patient {patient.id} queued {np.round(patient.q_time_fopa,2)} days for 1st app. Priority {patient.priority}")

                # Line/list to add to appointment log held in memory (df) or saved (csv)
                patient.ls_appt_to_add = [patient.id, patient.ls_appt[-1],patient.priority,patient.apptype,patient.type,patient.q_time_fopa,start_q_fopa,patient.tradition_dna,self.g.repid]
                patient.ls_patient_to_add = [patient.id , patient.q_time_fopa,999,self.g.repid] # deprecated

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
                    if patient.tradition_dna:
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
                        if patient.pifu_dna:
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

        self.g.results['rep']=self.g.repid


    def calculate_mean_q_time(self):
        """Deprecated. A method that calculates the average quueing time for the nurse.
        We can call this at the end of each run """

        self.results_df['Q_Time_total']=self.results_df.sum(axis=1)
        self.mean_q_time_total = self.results_df.mean(axis=0)

    """ Deprecated. A method to write run results to file.  Here, we write the run number
    against the the calculated mean queuing time for the nurse across
    patients in the run.  Again, we can call this at the end of each run"""
    def write_run_results(self):
        with open(self.savepath +"trial_results.csv", "a") as f:
            writer = csv.writer(f, delimiter=",")
            results_to_write = [self.run_number,
                                self.mean_q_time_total['Q_time_fopa'],
                                self.mean_q_time_total['Q_time_fuopa'],
                                self.mean_q_time_total['Q_Time_total']]
            print(results_to_write)
            writer.writerow(results_to_write)


    def chart(self):
        """ Plot results relevant to one run """
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
            x = self.g.appt_queuing_results[self.g.appt_queuing_results['priority'] == priority].index

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

    def summarise(self):
        """ single run summaries for streamlit (quartiles) - NEED CHECKING, NOT RUNNING AS SHOULD """
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
        """Monitors modelled system at regular intervals (as defined by audit interval in self.g)"""

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


    def run(self,repid=1):
        """  Run method to do a single run of the model.

        The run method starts up the entity generators, and tells SimPy to start
        running the environment for the duration specified in the g class. After
        the simulation has run, it calls the methods that calculate run
        results, and the method that writes these results to file

        Args:
            repid (int, optional): The replication id. Defaults to 1.

        Returns:
            chart_output: Chart output for streamlit
            text_output: Text output for streamlit
            quant_output: KPI output for streamlit
            Other outputs are stored within object (self) rather than returned.
        """

        # Start processes: entity generators and audit
        self.env.process(self.generate_wl_arrivals())

        # Check for unavailable feature use or not. If so, create slot obstructor generator.
        if self.g.unavail_on:
            self.env.process(self.obstruct_slots())

        # Generator to perform audit at fixed intervals
        self.env.process(self.perform_audit())

        # Run simulation
        self.env.run(until=self.g.obs_duration + self.g.warm_duration)

        # End of simulation run. Build and save results.

        # Load Results log - audit
        if self.g.loglinesave:
            self.g.results = pd.read_csv(self.savepath + "batch_mon_audit_ls.csv") # read from csv
            self.g.results = self.g.results[self.g.results['rep'] == self.g.repid]
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
