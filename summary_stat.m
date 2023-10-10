% This script collects summary statistics such as reaction time and error
% rate of participants that solve the IDED task. These then get saved for
% further processing for statistics purposes
%
% Created on:    01/02/2022
% Last modified: 01/09/2023
% Created by: Margarita Darna
% margarita.darna@lin-magdeburg.de

%%
%----------------------------------------------------------------------
%                  Prepare workspace and directories
%----------------------------------------------------------------------
% Cleaning workspace and command window
clear;clc;
% Setting up needed directories
dirs = {};
% change project_dir accordingly
dirs.proj_dir = 'C:/your_project_directory/';  
dirs.dt_dir         = strcat (dirs.proj_dir, 'Data/');
dirs.exp_dir        = strcat (dirs.proj_dir, 'IDED_v1_Analysis/');
dirs.raw_dt_dir     = strcat(dirs.dt_dir, 'Raw_data/');
dirs.derived_dt_dir = strcat(dirs.exp_dir, 'Data/Derived_data/');
dirs.analysis_dir   = strcat(dirs.exp_dir, 'Functions/');
dirs.output_dir     = strcat(dirs.exp_dir, 'Output/');
dirs.prepr_dir      = strcat(dirs.output_dir, '1_Preprocessing/');

% adding analysis path and subfolders
addpath(genpath(dirs.exp_dir));

%%
%----------------------------------------------------------------------
%            Load excel table with Subject Information
%----------------------------------------------------------------------
subj_info = readtable(strcat(dirs.exp_dir, 'Protocol.xlsx'));
subj_info = subj_info(subj_info.Excluded== 0,:);
subs      = subj_info.Pseudonym;
age_group = subj_info.age_cohort;
is_young  = categorical(age_group) == 'young';
is_old  = categorical(age_group) == 'old';
subj_young = subs(is_young);
subj_old = subs(~is_young);
age= subj_info.age;
gender = subj_info.Gender;

% display some summary information
fprintf('Young:\nAge:   %.1f +- %.1f\n', mean(age(is_young)), std(age(is_young)))
fprintf('N:     %i (Female: %i)\n', numel(subs(is_young)), sum(gender(is_young) == 2))
fprintf('Education years:  %.1f +- %.1f\n', mean(subj_info.Education(is_young), 'omitnan'), std(subj_info.Education(is_young), 'omitnan'))
fprintf('PSS10: %.1f +- %.1f\n', mean(subj_info.PSS10(is_young), 'omitnan'), std(subj_info.PSS10(is_young), 'omitnan'))
fprintf('MWTB:  %.1f +- %.1f\n\n', mean(subj_info.MWTB(is_young), 'omitnan'), std(subj_info.MWTB(is_young), 'omitnan'))
fprintf('Old:\nAge:   %.1f +- %.1f\n', mean(age(is_old)), std(age(is_old)))
fprintf('N:     %i (Female: %i)\n', numel(subs(is_old)), sum(gender(is_old) == 2))
fprintf('Education years:  %.1f +- %.1f\n', mean(subj_info.Education(is_old), 'omitnan'), std(subj_info.Education(is_old), 'omitnan'))
fprintf('PSS10: %.1f +- %.1f\n', mean(subj_info.PSS10(is_old), 'omitnan'), std(subj_info.PSS10(is_old), 'omitnan'))
fprintf('MWTB:  %.1f +- %.1f\n', mean(subj_info.MWTB(is_old), 'omitnan'), std(subj_info.MWTB(is_old), 'omitnan'))
fprintf('MMSE:  %.1f +- %.1f\n\n', mean(subj_info.MMSE(is_old), 'omitnan'), std(subj_info.MMSE(is_old), 'omitnan'))

%% 
%----------------------------------------------------------------------
%      Create individual performance struct and load directories
%----------------------------------------------------------------------
task = "IDED";

perf_all = {};
perf_all.subs = subs;

for n = 1: numel(task)
    fprintf('\nLoading files for %s:\n-----------------------------------\n', task(n));

    if task(n) == "IDED"
        perf_all.IDED.perf_dirs = get_perf_dirs(dirs.raw_dt_dir, dirs.derived_dt_dir, subs, task(n));
    else
        fprintf("Task %s not found!!!\n-----------------------------------\n", task(n))
    end
end

%% 
%------------------------------------------------------------------------
%                   Retrieve summary statistics
%------------------------------------------------------------------------
perf_all = get_perf_all(dirs, task, perf_all);

% print some overall results
fprintf('IDED Young:\nReaction Time: %.0f +- %.0f ms\n', mean(perf_all.IDED.RT.median_all(is_young)) * 1000, std(perf_all.IDED.RT.median_all(is_young)) * 1000)
fprintf('Error rate:    %.2f +- %.2f %% \n\n', mean(perf_all.IDED.error.error_all(is_young)) * 100, std(perf_all.IDED.error.error_all(is_young)) * 100)
fprintf('IDED Old:\nReaction Time: %.0f +- %.0f ms\n', mean(perf_all.IDED.RT.median_all(is_old)) * 1000, std(perf_all.IDED.RT.median_all(is_old)) * 1000)
fprintf('Error rate:    %.2f +- %.2f %% \n\n', mean(perf_all.IDED.error.error_all(is_old)) * 100, std(perf_all.IDED.error.error_all(is_old)) * 100)

fprintf('IDED All:\nReaction Time: %.0f +- %.0f ms\n', mean(perf_all.IDED.RT.median_all) * 1000, std(perf_all.IDED.RT.median_all) * 1000)
fprintf('Skewness: %.2f +- %.2f \n', mean(perf_all.IDED.RT.skewness_all), std(perf_all.IDED.RT.skewness_all))
fprintf('Error rate:    %.2f +- %.2f %% \n', mean(perf_all.IDED.error.error_all) * 100, std(perf_all.IDED.error.error_all) * 100)
fprintf('Error rate:    %.2f +- %.2f %% \n\n', mean(perf_all.IDED.error.error_all) * 100, std(perf_all.IDED.error.error_all) * 100)


%%
%------------------------------------------------------------------------
%                       Create table
%------------------------------------------------------------------------
% Reaction time
% use second repeat as control condition
% separate by condition
dv = [perf_all.IDED.RT.median_repeat2(is_young,:); perf_all.IDED.RT.median_repeat2(is_old,:); ...
    perf_all.IDED.RT.median_ID(is_young,:); perf_all.IDED.RT.median_ID(is_old,:); ...
    perf_all.IDED.RT.median_ED(is_young,:); perf_all.IDED.RT.median_ED(is_old,:)];
subj_num = [1:numel(subs) 1:numel(subs) 1:numel(subs)]';

X = [dv subj_num];

% Create table to be saved
T = array2table(X);
% Default heading for the columns will be A1, A2 and so on. 
% You can assign the specific headings to your table in the following manner
T.Properties.VariableNames = {'dv','subj_num'};
T.between = [repmat('young', sum(is_young), 1); repmat(' old ', sum(~is_young), 1) ;...
    repmat('young', sum(is_young), 1); repmat(' old ', sum(~is_young), 1);
    repmat('young', sum(is_young), 1); repmat(' old ', sum(~is_young), 1)];
T.within = [repmat('repeat', sum(is_young), 1); repmat('repeat', sum(~is_young), 1); ...
    repmat('  ID  ', sum(is_young), 1); repmat('  ID  ', sum(~is_young), 1);
    repmat('  ED  ', sum(is_young), 1); repmat('  ED  ', sum(~is_young), 1)];
T.gender =  [gender(is_young)' - 1; gender(is_old)' - 1; gender(is_young)' - 1; gender(is_old)' - 1; ...
    gender(is_young)' - 1; gender(is_old)' - 1;];
writetable(T, strcat(dirs.output_dir, '4_Stats\stat_RT_repeat_ID_ED.csv'));


% Error rate
% use second repeat as control condition
dv = [perf_all.IDED.error.error_repeat2(is_young,:); perf_all.IDED.error.error_repeat2(is_old,:); ...
    perf_all.IDED.error.error_ID(is_young,:); perf_all.IDED.error.error_ID(is_old,:); ...
    perf_all.IDED.error.error_ED(is_young,:); perf_all.IDED.error.error_ED(is_old,:)];
subj_num = [1:numel(subs) 1:numel(subs) 1:numel(subs)]';

X = [dv subj_num];

% Create table to be saved
T = array2table(X);
T.Properties.VariableNames = {'dv','subj_num'};
T.between = [repmat('young', sum(is_young), 1); repmat(' old ', sum(~is_young), 1) ;...
    repmat('young', sum(is_young), 1); repmat(' old ', sum(~is_young), 1);
    repmat('young', sum(is_young), 1); repmat(' old ', sum(~is_young), 1)];
T.within = [repmat('repeat', sum(is_young), 1); repmat('repeat', sum(~is_young), 1); ...
    repmat('  ID  ', sum(is_young), 1); repmat('  ID  ', sum(~is_young), 1);
    repmat('  ED  ', sum(is_young), 1); repmat('  ED  ', sum(~is_young), 1)];
T.gender =  [gender(is_young)' - 1; gender(is_old)' - 1; gender(is_young)' - 1; gender(is_old)' - 1; ...
    gender(is_young)' - 1; gender(is_old)' - 1;];
writetable(T, strcat(dirs.output_dir, '4_Stats\stat_error_repeat2_ID_ED.csv'));
