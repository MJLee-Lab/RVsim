# RVsim
This function simulates population-level responses to DNA damage, both with and without activation of cell death. The simulated population sizes are parameterized around experimentally observed growth and death rates, and these are used to calculate relative viability (RV) and relative sensitivity to cell death. These metrics can be used to determine how much (or how little) the observed cell death rate can contribute to traditional evaluations of drug sensitivity using population size (Honeywell et al. 2023, bioRxiv).

# Data collection
* Evaluation of the effect of cell death on relative viability requires a minimum of three experimentally observed parameters: 1) doubling time of untreated cells, 2) doubling time of treated cells, and 3) the death rate of treated cells. These features can be derived simultaneously using the GRADE framework (MJLee-Lab/GRADE), or through separate measurement of live cells (counted over time) and measurement of death rate in FLICK (Richards et al., STAR Protocols, 2021).

# RVsim structure
* The RVsim function has 6 inputs, ordered as shown below:

      RVsim(untreated_tau, treated_tau, treated_dr, finaltime, dose_curve_tp, ec50)

    **untreated_tau** = untreated cell growth rate, in doublings per hour
    
    **treated_tau** = treated cell growth rate, in doublings per hour
    
    **treated_dr** = death rate of treated cells, in % death per hour
    
    **finaltime** = assay length, in hours
    
    **dose_curve_tp** = timepoint for plotted dose curve, in hours
    
    **ec50** = simulated EC50 for dose curves, value from 1-8 (default = 4)
    

* The first 5 inputs for the function are required, with the 6th (ec50) available as an optional input to change the aesthetics of the simulated dose curves.

# Running RVsim
* To calculate the drug-induced death rate, call the function RVsim:

      RV_timecourse = RVsim(untreated_tau, treated_tau, treated_dr, finaltime, dose_curve_tp, ec50)

* Example parameters are shown below for a condition where cells undergo complete growth arrest and activate low levels of cell death in response to DNA damage

    **untreated_tau** = 0.0417
    
    **treated_tau** = 0
    
    **treated_dr** = 0.3
    
    **finaltime** = 168
    
    **dose_curve_tp** = 72
    
    **ec50** = 4
    

* To run example data:

      RV_timecourse = RVsim(0.0417, 0, 0.3, 168, 72, 4)

* This yields a table, ‘RV_timecourse’, which contains the simulated population size over time, and the fitted RV dose curves for each timepoint. 

* The function will also produce two plots: 1) a plot of the relative sensitivity to cell death, with the input death rate highlighted, and 2) the simulated RV at the selected timepoint, showing the RV in the presence and absence of cell death.
