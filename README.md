![image](https://user-images.githubusercontent.com/69211911/186386851-61f28c91-c9de-4018-be58-49d7fdc80385.png)


## Context 
PIFU is where, based on clinical steer, a patient goes on a pathway where they can arrange their own follow-up appointments as-and-when-needed rather than attending regular scheduled appointments. It means for instance that time spent on low clinical value appointments can be repurposed for managing complex patients and upstream care. The interplay between the downstream efficiencies and the upstream referral-to-treatment waiting list, as well as changes to planned levels of resourcing and allocation, can however be hard to grasp.

As such, in this project we propose to use Discrete Event Simulation - where individual patients’ pathway movement through a system with constrained resources is simulated and monitored - to better understand the role that PIFU and other interventions can play in addressing the elective (referral-to-treatment) backlog and supporting patient care.


## Project objectives
*	Use Discrete Event Simulation to better understand the role patient initiated follow-up (PIFU) can play in addressing the (outpatient and elective) referral-to-treatment (RTT) backlog, with particular emphasis on Rheumatology
*	Advocate and demonstrate the value of pathway and behavioural modelling in central decision-making for digital pathway redesign, elective recovery and the LTP personalisation agenda to support better patient care, experience and system resilience


## Main model

Can be run from `simpy_rheum_v004_run.py` in Python (spyder, VS Code). The main `simpy` model implementation is in `simpy_rheum_v004.py`. 

## Toy tool

A toy interactive app leveraging streamlit.

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://nhsx-hsma4-12-des-rheum-streamlit-model-mf-v0031-batch-e91w2z.streamlitapp.com/)
