# IDED
**Intra-Dimensional/Extra-Dimensional set-shifting task**

This code belongs to the paper "Age effects on electrocortical processing during intradimensional/extradimensional set shifting" by Darna et al. (2023), Preprint will be available soon

## Requirements
This code was developed and run using the following software:
- Windows 10 Professional 64-bit
- [MATLAB R2021b](https://de.mathworks.com/help/matlab/release-notes.html) (Version 9.11)
- [MATLAB Statistics and ML Toolbox](https://de.mathworks.com/products/statistics.html) (Version 12.2)

To run the paradigm:
- [Psychtoolbox](http://psychtoolbox.org/) (Version 3.0.18)
- [Mex-File Plug-in IO64](http://apps.usd.edu/coglab/psyc770/IO64.html) (for triggers)

For EEG analysis:
- [Brain Vision Analyzer](https://www.brainproducts.com/downloads/analyzer/) (Version 2.2)
- [FieldTrip](https://www.fieldtriptoolbox.org/download/) (Version 20211122)

For statistics:
- [RStudio](https://posit.co/download/rstudio-desktop/) (R version 4.2.2)

## The IDED Paradigm:
In order to be able to run the IDED paradigm on your computer copy [ExpCode](https://github.com/margdarna/IDED/tree/main/ExpCode) on your computer.

### Instructions:
- Open [IDED.m](https://github.com/margdarna/IDED/blob/main/ExpCode/IDED.m)
- Modify 'project_dir' to your project directory. Important: [ExpCode](https://github.com/margdarna/IDED/tree/main/ExpCode) should be located in project_dir
- Modify 'raw_data_dir' to your directory where the raw behavioural data should be saved.
- Run the script.
- In the pop-up window: add subject number from 001 and 500.

### Additional information:
- The [IDED.m](https://github.com/margdarna/IDED/blob/main/ExpCode/IDED.m) script creates a subject folder called "SXXX", where XXX corresponds to the subject number. The folder includes three subfolders.
  - "Behavior": Here, behavior files are automatically saved during and after the experiment. More details about these files are provided in [IDED.m](https://github.com/margdarna/IDED/blob/main/ExpCode/IDED.m).
  - "EEG": Here, the EEG files can be manually saved.
  - "Exp_Files": Here, the preallocated subject_stimulus file is saved again as backup.
- The script does not allow double use of subject numbers. in other words, if a subject folder with a specific number already exists, you will be asked to input a different subject number. Exception is subject number 500.
- Subject number 500 is reserved for testing, a.k.a whenever you use this subject nnumber its folder will be overwritten. This is useful, when you want to test if a specific aspect of the code is working and you do not wish to save any actual data.
- The [IDED.m](https://github.com/margdarna/IDED/blob/main/ExpCode/IDED.m) requires the [Subject_stimuli](https://github.com/margdarna/IDED/tree/main/ExpCode/Subject_Stimuli) files. These are pregenerated files that have predetermined pseudorandomized order of trials.
- Here, I provide the subject stimuli files for 50 subjects. If you wish to create more or new ones use the function [allocate_stimuli_IDED.m](https://github.com/margdarna/IDED/blob/main/ExpCode/allocate_stimuli_IDED.m). To run it, please remember to modify the directories accordingly. Also adapt the variable "total_subj" to the total subjects number you wish to use. Be careful: this function overrides the currently available files in the [Subject_stimuli](https://github.com/margdarna/IDED/tree/main/ExpCode/Subject_Stimuli) folder.
- The [IDED.m](https://github.com/margdarna/IDED/blob/main/ExpCode/IDED.m) already includes a set of 20 practice trials that include all types of shifts at least one time and feedback. These always appear first before the actual paradigm.
- You can always stop the paradigm by pressing the Esc button while stimuli are being presented.
- Currently, trigger lines are commented out, because the code does not run if io64 is not available on your device. If you want to use triggers, e.g. for the purposes of EEG, you would need to remove the comments from those lines. These are the lines that concern the variable 'trig'.
- The instruction text is currently in German. Feel free to adapt.

## MATLAB Scripts
The MATLAB Scripts require the folder [Functions](https://github.com/margdarna/IDED/tree/main/Functions) to work correctly. If the folder is saved in the same folder as the other MATLAB scripts function
IDED Analysis, Behavioral Results:

### Behavioral results
The extraction of behavioral results as .csv table for further statistical analysis is performed with the script [summary_stat](https://github.com/margdarna/IDED/blob/main/summary_stat.m).
The script reads participant information from the table 'Protocol.xlsx'.

### EEG Analysis:
After pre-processing in [Brain Vision Analyzer](https://www.brainproducts.com/downloads/analyzer/) the raw files are opened in MATLAB using [FieldTrip](https://www.fieldtriptoolbox.org/download/).

Here is a brief description of the scripts:
- [IDED_ERP_Analysis.m](https://github.com/margdarna/IDED/blob/main/IDED_ERP_Analysis.m): Creates epochs, calculates subject ERP averages per conditions, extracts P300 average as .csv, plots ERP
- [IDED_ERP_Analysis_Latency.m](https://github.com/margdarna/IDED/blob/main/IDED_ERP_Analysis_Latency.m): Uses maximum peak method for latency calculation, extracts P300 latency as .csv
- [IDED_TFR_Analysis.m](https://github.com/margdarna/IDED/blob/main/IDED_TFR_Analysis.m): creates epochs, calculates subject specific TFR, extracts theta power as .csv, plots results, creates movies
- [IDED_single_trial_shepherd.m](https://github.com/margdarna/IDED/blob/main/IDED_single_trial_shepherd.m): perform shepher's Pi correlation, extracts z-values as .csv

## Statistics:
The statistic scripts are performed in sequential order from 1 to 8. 

Here is a brief description of the scripts:
- [1_check_summary_stats.R](https://github.com/margdarna/IDED/blob/main/Statistics/1_check_summary_stats.R): check the differences between summary statistics characteristics between the two groups (age, education years, MWTB performance, gender distribution).
- [2_ASST_survival.R](https://github.com/margdarna/IDED/blob/main/Statistics/2_ASST_survival.R): statistical analysis of the ASST by estimating a survival curve using the Kaplan-Meier method.
- [3_two_way_ANOVA_EMCO_RT.R](https://github.com/margdarna/IDED/blob/main/Statistics/3_two_way_ANOVA_EMCO_RT.R): mixed 2x3 ANOVA of reaction times.
- [4_two_way_ANOVA_EMCO_Error.R](https://github.com/margdarna/IDED/blob/main/Statistics/4_two_way_ANOVA_EMCO_Error.R): mixed 2x3 ANOVA of error rates.
- [5_two_way_ANOVA_EMCO_P300latency.R](https://github.com/margdarna/IDED/blob/main/Statistics/5_two_way_ANOVA_EMCO_P300latency.R): mixed 2x3 ANOVA of P300 latency.
- [6_two_way_ANOVA_EMCO_P300amplitude.R](https://github.com/margdarna/IDED/blob/main/Statistics/6_two_way_ANOVA_EMCO_P300amplitude.R): mixed 2x3 ANOVA of P300 amplitude.
- [7_two_way_ANOVA_EMCO_theta250power.R](https://github.com/margdarna/IDED/blob/main/Statistics/7_two_way_ANOVA_EMCO_theta250power.R): mixed 2x3 ANOVA of frontocentral theta power.
- [8_IDED_singletrl_theta_RT_Shepherd.R](https://github.com/margdarna/IDED/blob/main/Statistics/8_IDED_singletrl_theta_RT_Shepherd.R): mixed 2x3 ANOVA of z-values (a.k.a transformed shepherd's pi from single trial correlations between reaction time and theta power).

## Contact information
For any questions, please contact:

margarita.darna@lin-magdeburg.de

Margarita Darna

PhD Candidate

Department Behavioral Neurology

Working Group Synapse-Brain-Cognition

Leibniz Institut for Neurobiology

Brenneckestrasse 6

39118 Magdeburg, Germany
