% This script performs single trial analysis for data obtained from the
% IDED task of the A05 project of CRC1436
%
% Created on:    04/11/2022
% Last modified: 01/09/2023
% Created by: Margarita Darna
% margarita.darna@lin-magdeburg.de

%----------------------------------------------------------------------
%               0 - Prepare workspace and directories
%----------------------------------------------------------------------

% Cleaning workspace and command window
clear; clc;

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
addpath(genpath(dirs.analysis_dir));

% making sure that fieldtrip is called correctly
ft_defaults

% define the subjects
subj_info = readtable(strcat(dirs.analysis_dir, 'Protocol.xlsx'));
subj_info = subj_info(subj_info.Excluded== 0,:);
subs      = subj_info.Pseudonym;
age_group = subj_info.age_cohort;
is_young  = categorical(age_group) == 'young';
subj_young = subs(is_young);
subj_old = subs(~is_young);

% Calculate average age for each age group
tblstats = grpstats(subj_info(:,2:3),"age_cohort", ["mean" "std"]);
% count gender in each group
[GC_young,GR_young] = groupcounts(subj_info(is_young,:).Gender);
[GC_old,GR_old] = groupcounts(subj_info(~is_young,:).Gender);

%%
%----------------------------------------------------------------------
%   1 - Load reaction times and frequency power for each participant
%----------------------------------------------------------------------

for i = 1:height(subs)
    % Define subject
    subj = subs{i};

    % print information in command window
    fprintf("************************\nStarting %s\n************************\n", subj)

    % Load perf.mat of that subject to retrieve RTs
    load(sprintf('%s%s/%s_task-IDED_perf.mat', dirs.derived_dt_dir, subj, subj))

    % load frequency power of specific subject to retrieve theta power
    load(sprintf('%s2_Subj_Avg/%s_subj_avg_stimpres_TFR_trial.mat', dirs.output_dir, subj))

    %%theta%%
    % averaging over frontal theta, but keeping trials
    cfg = [];
    cfg.channel     = {'F1', 'F2', 'Fz', 'FC1', 'FC2', 'FCz'};
    cfg.avgoverchan = 'yes';
    cfg.avgovertime = 'yes';
    cfg.frequency   =  [4 8];
    cfg.avgoverfreq =  'yes';
    cfg.latency     =  [0.25 0.50];
    TFR             = ft_selectdata(cfg, alltrl_TFR_trial);
    theta_alltrl    = TFR.powspctrm;
    TFR             = ft_selectdata(cfg, repeat1_TFR_trial);
    theta_repeat1    = TFR.powspctrm;
    TFR             = ft_selectdata(cfg, repeat2_TFR_trial);
    theta_repeat2    = TFR.powspctrm;
    TFR             = ft_selectdata(cfg, ID_TFR_trial);
    theta_ID        = TFR.powspctrm;
    TFR             = ft_selectdata(cfg, ED_TFR_trial);
    theta_ED        = TFR.powspctrm;
    TFR             = ft_selectdata(cfg, last_TFR_trial);
    theta_last      = TFR.powspctrm;

    % calculate baseline and calculate difference
    cfg.latency     =  [-0.4 -0.20];
    TFR_bl          = ft_selectdata(cfg, alltrl_TFR_trial);
    theta_alltrl    = theta_alltrl - TFR_bl.powspctrm;
    TFR_bl          = ft_selectdata(cfg, repeat1_TFR_trial);
    theta_repeat1    = theta_repeat1 - TFR_bl.powspctrm;
    TFR_bl          = ft_selectdata(cfg, repeat2_TFR_trial);
    theta_repeat2    = theta_repeat2 - TFR_bl.powspctrm;
    TFR_bl          = ft_selectdata(cfg, ID_TFR_trial);
    theta_ID        = theta_ID - TFR_bl.powspctrm;
    TFR_bl          = ft_selectdata(cfg, ED_TFR_trial);
    theta_ED        = theta_ED - TFR_bl.powspctrm;
    TFR_bl          = ft_selectdata(cfg, last_TFR_trial);
    theta_last      = theta_last- TFR_bl.powspctrm;
    
   % Choose trials from perf that have not been rejected in EEG analysis
    % only include epochs with correct answers and epochs without artefacts
    % epochs with correct answers have been already selected in trialinfo
    % from previous script
    % also exclude outliers
    ind            = ismember(perf.resp_mat.trial, alltrl_TFR_trial.trialinfo(:,3));
    RT_alltrl      = perf.resp_mat.RT(ind);

    ind            = ismember(perf.resp_mat.trial, repeat1_TFR_trial.trialinfo(:,3));
    RT_repeat1      = perf.resp_mat.RT(ind);

    ind            = ismember(perf.resp_mat.trial, repeat2_TFR_trial.trialinfo(:,3));
    RT_repeat2      = perf.resp_mat.RT(ind);

    ind            = ismember(perf.resp_mat.trial, ID_TFR_trial.trialinfo(:,3));
    RT_ID          = perf.resp_mat.RT(ind);
 
    ind            = ismember(perf.resp_mat.trial, ED_TFR_trial.trialinfo(:,3));
    RT_ED          = perf.resp_mat.RT(ind);

    ind            = ismember(perf.resp_mat.trial, last_TFR_trial.trialinfo(:,3));
    RT_last        = perf.resp_mat.RT(ind);

    % perform shepherd's Pi correlation
    % alltrl
    [Pi, p, ~] = Shepherd(theta_alltrl, RT_alltrl);
    % save information in array
    Pi_alltrl(i) = Pi;
    p_alltrl(i)  = p;
    %outlier_alltrl = o1;
    % calculate Fisher's z
    z_alltrl(i) = 0.5 * log((1 + Pi) / (1 - Pi));

    % repeat1
    [Pi, p, ~] = Shepherd(theta_repeat1, RT_repeat1);
    % save information in array
    Pi_repeat1(i) = Pi;
    p_repeat1(i)  = p;
    %outlier_repeat1 = o1;
    % calculate Fisher's z
    z_repeat1(i) = 0.5 * log((1 + Pi) / (1 - Pi));

    % repeat2
    [Pi, p, ~] = Shepherd(theta_repeat2, RT_repeat2);
    % save information in array
    Pi_repeat2(i) = Pi;
    p_repeat2(i)  = p;
    %outlier_repeat2 = o1;
    % calculate Fisher's z
    z_repeat2(i) = 0.5 * log((1 + Pi) / (1 - Pi));
    
    % ID
    [Pi, p, ~] = Shepherd(theta_ID, RT_ID);
    % save information in array
    Pi_ID(i) = Pi;
    p_ID(i)  = p;
    %outlier_ID = o1;
    % calculate Fisher's z
    z_ID(i) = 0.5 * log((1 + Pi) / (1 - Pi));

    % ED
    [Pi, p, ~] = Shepherd(theta_ED, RT_ED);
    % save information in array
    Pi_ED(i) = Pi;
    p_ED(i)  = p;
    %outlier_ED = o1;
    % calculate Fisher's z
    z_ED(i) = 0.5 * log((1 + Pi) / (1 - Pi));

    % last
    [Pi, p, ~] = Shepherd(theta_last, RT_last);
    % save information in array
    Pi_last(i) = Pi;
    p_last(i)  = p;
    %outlier_last = o1;
    % calculate Fisher's z
    z_last(i) = 0.5 * log((1 + Pi) / (1 - Pi));
   
    clear ind TFR RT_repeat1 RT_repeat2 RT_ID RT_ED RT_last theta_repeat1 theta_repeat2 ...
        theta_ID theta_ED theta_last TFR_bl Pi p
end

% save data
save(strcat(dirs.output_dir, '6_Single_Trial_Analysis\IDED_theta250_single_trial_shepherd.mat'),"Pi_alltrl", ...
    "Pi_repeat1", "Pi_repeat2", "Pi_ID", "Pi_ED" , "Pi_last", "p_alltrl", "p_repeat1", "p_repeat2", "p_ID", ...
    "p_ED", "p_last", "z_alltrl", "z_repeat1", "z_repeat2", "z_ID", "z_ED", "z_last");

%%
%----------------------------------------------------------------------
%                       2 - Generate matrix
%----------------------------------------------------------------------
load(strcat(dirs.output_dir, '6_Single_Trial_Analysis\IDED_theta250_single_trial_shepherd.mat'))

% get age group from subject info
subs      = subj_info.Pseudonym;
age_group = subj_info.age_cohort;
is_young  = categorical(age_group) == 'young';

% generate matrix for ANOVA
% depending on condition
% depending variable
dv = [z_repeat2'; z_ID'; z_ED'];
% subject number
subj_num = [1:numel(subs) 1:numel(subs) 1:numel(subs)]';
% generate table
X = [dv subj_num];
T = array2table(X);
% Assign headings to the table
T.Properties.VariableNames = {'dv', 'subj_num'};
% add between and within variables
T.is_young = [is_young; is_young; is_young];
T.within = [repmat('repeat', numel(subs), 1); repmat('  ID  ', numel(subs), 1); ...
    repmat('  ED  ', numel(subs), 1)];
% save table as csv
writetable(T, strcat(dirs.output_dir, '6_Single_Trial_Analysis\IDED_theta250_single_trial_shepherd_repeat2_ID_ED.csv'));
