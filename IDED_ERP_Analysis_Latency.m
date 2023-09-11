% EEG Analysis script for Fieldtrip for single subjects and grand averages
%
% This script opens preprocessed eeg data from brain vision analyzer,
% and gives out basic analytics
% Based on tutorial from fieldtriptoolbox.org
%
% Written by: Margarita Darna
% Created on: 06. December 2022
% Last modified on: 06. Dezember 2022

%%
%----------------------------------------------------------------------
%                  Prepare workspace and directories
%----------------------------------------------------------------------
clear;clc;
% Setting up needed directories
dirs = {};
% change project_dir accordingly
dirs.proj_dir = 'C:/your_project_directory/';   
dirs.dt_dir         = strcat (dirs.proj_dir, 'Data/');
dirs.exp_dir        = strcat (dirs.proj_dir, 'IDED_v1_Analysis/');
dirs.raw_dt_dir     = strcat(dirs.dt_dir, 'Raw_data/');
dirs.derived_dt_dir = strcat(dirs.dt_dir, 'Derived_data/IDED/');
dirs.analysis_dir   = strcat(dirs.exp_dir, 'Analysis/');
dirs.output_dir     = strcat(dirs.exp_dir, 'Output/');
dirs.prepr_dir      = strcat(dirs.output_dir, '1_Preprocessing/');

% adding analysis path and subfolders
addpath(genpath(dirs.analysis_dir));

% define the subjects
% big epochs for the purpose of time frequency analysis
prestim = 1.25;
poststim = 3;

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
%%
%----------------------------------------------------------------------
%                         Retrieve ERPs
%----------------------------------------------------------------------
% load grandaverages
% Load single subject data to extract the channels and time range of
% interest
load(strcat(dirs.output_dir, '2_Subj_Avg\young_subj_avg_stimpres_ERP'));
load(strcat(dirs.output_dir, '2_Subj_Avg\old_subj_avg_stimpres_ERP'));

% P300 - P3b - parietal
% prepare the cfg
cfg = [];
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'no';
cfg.channel     = {'Pz', 'P1', 'P2', 'POz', 'PO3', 'PO4'};%{'Pz'}; 
cfg.latency     =  [0.3 0.6];

for i = 1:numel(young_repeat1_stimpres_ERP)
fprintf('************************\n%s\n************************\n', subj_young{i});
young_repeat1_stimpres{i} = ft_selectdata(cfg, young_repeat1_stimpres_ERP{i});
young_repeat2_stimpres{i} = ft_selectdata(cfg, young_repeat2_stimpres_ERP{i});
young_shift_stimpres{i}  = ft_selectdata(cfg, young_shift_stimpres_ERP{i});
young_ID_stimpres{i}     = ft_selectdata(cfg, young_ID_stimpres_ERP{i});
young_ED_stimpres{i}     = ft_selectdata(cfg, young_ED_stimpres_ERP{i});
young_last_stimpres{i}   = ft_selectdata(cfg, young_last_stimpres_ERP{i});
young_alltrl_stimpres{i} = ft_selectdata(cfg, young_alltrl_stimpres_ERP{i});
end

for i = 1:numel(old_repeat1_stimpres_ERP)
fprintf('************************\n%s\n************************\n', subj_old{i});
old_repeat1_stimpres{i}  = ft_selectdata(cfg, old_repeat1_stimpres_ERP{i});
old_repeat2_stimpres{i}  = ft_selectdata(cfg, old_repeat2_stimpres_ERP{i});
old_shift_stimpres{i}   = ft_selectdata(cfg, old_shift_stimpres_ERP{i});
old_ID_stimpres{i}      = ft_selectdata(cfg, old_ID_stimpres_ERP{i});
old_ED_stimpres{i}      = ft_selectdata(cfg, old_ED_stimpres_ERP{i});
old_last_stimpres{i}    = ft_selectdata(cfg, old_last_stimpres_ERP{i});
old_alltrl_stimpres{i}  = ft_selectdata(cfg, old_alltrl_stimpres_ERP{i});
end
%%
%----------------------------------------------------------------------
%               Using maximum peak method for latency calculation
%----------------------------------------------------------------------

for i = 1:numel(young_repeat1_stimpres)
    young_alltrl_latency(i) = get_latency_maxpeak(young_alltrl_stimpres{i}.avg, 0.3, 500);
    young_repeat1_latency(i) = get_latency_maxpeak(young_repeat1_stimpres{i}.avg, 0.3, 500);
    young_repeat2_latency(i) = get_latency_maxpeak(young_repeat2_stimpres{i}.avg, 0.3, 500);
    young_ID_latency(i) = get_latency_maxpeak(young_ID_stimpres{i}.avg, 0.3, 500);
    young_ED_latency(i) = get_latency_maxpeak(young_ED_stimpres{i}.avg, 0.3, 500);
end

for i = 1:numel(old_repeat1_stimpres)
    old_alltrl_latency(i) = get_latency_maxpeak(old_alltrl_stimpres{i}.avg, 0.3, 500);
    old_repeat1_latency(i) = get_latency_maxpeak(old_repeat1_stimpres{i}.avg, 0.3, 500);
    old_repeat2_latency(i) = get_latency_maxpeak(old_repeat2_stimpres{i}.avg, 0.3, 500);
    old_ID_latency(i) = get_latency_maxpeak(old_ID_stimpres{i}.avg, 0.3, 500);
    old_ED_latency(i) = get_latency_maxpeak(old_ED_stimpres{i}.avg, 0.3, 500);
end

% Get mean all tral latency and standard deviation
% Calculate mean latency of the P300 component to include in amplitude
% analysis
ind_P300p_young = mean(young_alltrl_latency);
ind_P300p_old = mean(old_alltrl_latency);

% Save variables to be loaded in IDED_ERP_Analysis
save("ind_P300.mat", 'ind_P300p_young', 'ind_P300p_old')
% Display result
fprintf('IDED young:\n P300 latency: %.0f +- %.0f ms\n', ind_P300p_young  * 1000, std(young_alltrl_latency)*1000)
fprintf('IDED old:\n P300 latency: %.0f +- %.0f ms\n', ind_P300p_old  * 1000, std(old_alltrl_latency)*1000)

% generate matrix for ANOVA
% depending variable
% dv = [young_repeat1_latency; old_repeat1_latency; young_ID_latency; old_ID_latency; ...
%       young_ED_latency; old_ED_latency];
dv = [young_repeat2_latency old_repeat2_latency young_ID_latency old_ID_latency ...
      young_ED_latency old_ED_latency];
% subject number
subj_num = [1:numel(subs) 1:numel(subs) 1:numel(subs)];
% generate table
X = [dv' subj_num'];
T = array2table(X);
% Assign headings to the table
T.Properties.VariableNames = {'dv', 'subj_num'};
% add between and within variables
T.between = [repmat('young', sum(is_young), 1); repmat(' old ', sum(~is_young), 1); repmat('young', sum(is_young), 1); ...
    repmat(' old ', sum(~is_young), 1); repmat('young', sum(is_young), 1); repmat(' old ', sum(~is_young), 1)];
T.within = [repmat('repeat', sum(is_young), 1); repmat('repeat', sum(~is_young), 1); repmat('  ID  ', sum(is_young), 1); ...
    repmat('  ID  ', sum(~is_young), 1); repmat('  ED  ', sum(is_young), 1); repmat('  ED  ', sum(~is_young), 1)];
% save table as csv
writetable(T, strcat(dirs.output_dir, '4_Stats\stat_ERP_P300p_latency.csv'));
