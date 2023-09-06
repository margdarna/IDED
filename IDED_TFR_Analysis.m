% EEG Analysis script for Fieldtrip
%
% This script opens preprocessed eeg data from brain vision analyzer,
% and gives out basic analytics
% Based on tutorial from fieldtriptoolbox.org
%
% Written by: Margarita Darna
% Created on: 03. August 2022
% Last modified on: 01. September 2023

%%
%----------------------------------------------------------------------
%               0 - Prepare workspace and directories
%----------------------------------------------------------------------
clear;clc;
% Setting up needed directories
dirs = {};
%dirs.proj_dir = 'C:\Users\mdarna\Documents\PhD\A05-SFB1436\';    
dirs.proj_dir       = '//linstore01/home/mdarna/PhD/A05-SFB1436/';    
dirs.dt_dir         = strcat (dirs.proj_dir, 'Data/');
dirs.exp_dir        = strcat (dirs.proj_dir, 'IDED_v1_Analysis/');
dirs.raw_dt_dir     = strcat(dirs.dt_dir, 'Raw_data/');
dirs.derived_dt_dir = strcat(dirs.dt_dir, 'Derived_data/IDED/');
dirs.analysis_dir   = strcat(dirs.exp_dir, 'Analysis/');
dirs.output_dir     = strcat(dirs.exp_dir, 'Output/');
dirs.prepr_dir      = strcat(dirs.output_dir, '1_Preprocessing/');

% adding analysis path and subfolders
addpath(genpath(dirs.analysis_dir));

% making sure that fieldtrip is called correctly
ft_defaults

% define the subjects
% big epochs for the purpose of time frequency analysis
prestim = 1.25;
poststim = 3;

%%
%----------------------------------------------------------------------
%            Load excel table with Subject Information
%----------------------------------------------------------------------
subj_info = readtable(strcat(dirs.analysis_dir, 'Protocol.xlsx'));
subj_info = subj_info(subj_info.Excluded== 0,:);
subs      = subj_info.Pseudonym;
age_group = subj_info.age_cohort;
is_young  = categorical(age_group) == 'young';
subj_young = subs(is_young);
subj_old = subs(~is_young);
age= subj_info.age;
gender = subj_info.Gender;

% Calculate average age for each age group
tblstats = grpstats(subj_info(:,2:3),"age_cohort", ["mean" "std"]);
% count gender in each group
[GC_young,GR_young] = groupcounts(subj_info(is_young,:).Gender);
[GC_old,GR_old] = groupcounts(subj_info(~is_young,:).Gender);

%%
%----------------------------------------------------------------------
%            1 - Load single subject EEG and prepare epochs
%----------------------------------------------------------------------

for i = 1:height(subs)
% ensure that we don't mix up subjects
clear subjectdata
subj = subs{i};

% print information in command windo
fprintf("************************\nStarting %s\n************************\n", subj)

% define the filenames, parameters and other information that is subject
% specific
subjectdata.subj = subj;
% save all subject relevant information in subjectdata

% read dataset
cfg = [];
cfg.dataset = [dirs.derived_dt_dir, subj, '_task-IDED_eeg_prepr.dat'];
data_eeg = ft_preprocessing(cfg);

% select different trial types
% Port signalling:
% Fixation Cross:            10
% Stimulus Presentation:
%             practice:      210 + Target position
%                  pre:      220 + Target position
%               repeat:      230 + Target position
%                   ID:      240 + Target position
%                   ED:      250 + Target position
%                 last:      290 + Target position %% actual trigger: 35 or 36   
% Button press:              30  + correctness
% Feedback(practice trials): 40  + Feedback
% Break start:               50  + number of break
% Break end:                 60  + number of break

% extract epochs, stimulus locked
% here we are using a custom function that also rejects trials with
% artifacts and extracts trial number
cfg = [];
cfg.dataset               = [dirs.derived_dt_dir, subj, '_task-IDED_eeg_prepr.dat'];
cfg.trialfun              = 'IDED_stimpres_trialfun'; % custom function
cfg.trialdef.prestim      = prestim; % in seconds
cfg.trialdef.poststim     = poststim; % in seconds
cfg.trialdef.eventtype    = 'Stimulus';
cfg.trialdef.eventvalue   = {'S231', 'S232', 'S 35', 'S 36', 'S241', 'S242', 'S251', 'S252'};
cfg.trialdef.trigger_type = 'nexttrigger';
cfg_stimpres              = ft_definetrial(cfg);

% Change trial number in subject S023 because recording started after first
% break, Trial number taken out of perf struct of S023
% Recording started sometime after first break, I calculated the trials
% based on the second break. Number has been validated, do not change!
if strcmp(subj, 'S023')
    cfg_stimpres.trl(:,6) =  cfg_stimpres.trl(:,6)+ 163;
end

% create datasets
% stimulus locked
data_stimpres   = ft_redefinetrial(cfg_stimpres, data_eeg);

% remove baseline from stimulus locked trials
cfg = [];
cfg.demean          = 'yes';
cfg.baselinewindow  = [-0.25 0];
data_stimpres       = ft_preprocessing(cfg, data_stimpres);

% use original dataset (with no additional filtering) for TFR
data_TFR = data_stimpres;

% load previous subject data and save changes
load(strcat(dirs.prepr_dir, subj, '_IDED_ERP_ft_prepr_1.mat'), "subjectdata");

% save subject specific information
subjectdata.data_eeg = data_eeg;
subjectdata.data_TFR = data_TFR;

save(strcat(dirs.prepr_dir, subj, '_IDED_TFR_ft_prepr_1.mat'), "subjectdata", '-v7.3');
fprintf(strcat(subj, ' succesfully saved!\n\n'))
end
clear data_eeg data_stimpres data_TFR cfg_stimpres

%% 
%----------------------------------------------------------------------
%            2 - Calculate subject averages for each condition
%----------------------------------------------------------------------
for i = 1:height(subs)
% ensure that we don't mix up subjects
clear subjectdata
subj = subs{i};

% print information in command windo
fprintf("************************\nStarting %s\n************************\n", subj)

% load mat file that was prepared in the previous section
load(strcat(dirs.prepr_dir, subj, '_IDED_TFR_ft_prepr_1.mat'), "subjectdata");

data_TFR = subjectdata.data_TFR;

 % define cfg for TFR
 cfg              = [];                                             
 cfg.output       = 'pow';
 cfg.method       = 'wavelet';
 cfg.foi          = 1:1:60;            % start and ending frequency with steps
 cfg.toi          = -prestim:0.05:poststim;    % time window "slides" from prestim to poststim in steps of 0.05 sec (50 ms)
 cfg.width        = 7;
 cfg.gwidth       = 3;

% perform freqanalysis (TFR) for each trial type
% here we only choose trials with correct responses: data_TFR.trialinfo(:,2) == 31

% shift trials
shift_ind  = (data_TFR.trialinfo(:,1)==241 | data_TFR.trialinfo(:,1)==242 | ...
             data_TFR.trialinfo(:,1)==251 | data_TFR.trialinfo(:,1)==252) ...
             & data_TFR.trialinfo(:,2) == 31;
cfg.trials = find(shift_ind);
cfg.keeptrials  = 'no';
shift_TFR  = ft_freqanalysis(cfg, data_TFR);
% keep trials for single trial analysis
cfg.keeptrials  = 'yes';
shift_TFR_trial  = ft_freqanalysis(cfg, data_TFR);

% repeat1 trials - only choose the first trial after a shift
repeat1_ind  = [false; shift_ind(1:end - 1,:)];
cfg.trials  = find(repeat1_ind & (data_TFR.trialinfo(:,1)==231 |...
              data_TFR.trialinfo(:,1)==232) & data_TFR.trialinfo(:,2) == 31); % choose only correct responses
cfg.keeptrials  = 'no';
repeat1_TFR  = ft_freqanalysis(cfg, data_TFR);
% keep trials for single trial analysis
cfg.keeptrials  = 'yes';
repeat1_TFR_trial  = ft_freqanalysis(cfg, data_TFR);

% repeat2 trials - only choose the second trial after a shift
repeat2_ind  = [false; false; shift_ind(1:end - 2,:)];
cfg.trials  = find(repeat2_ind & (data_TFR.trialinfo(:,1)==231 |...
              data_TFR.trialinfo(:,1)==232) & data_TFR.trialinfo(:,2) == 31); % choose only correct responses
cfg.keeptrials  = 'no';
repeat2_TFR  = ft_freqanalysis(cfg, data_TFR);
% keep trials for single trial analysis
cfg.keeptrials  = 'yes';
repeat2_TFR_trial  = ft_freqanalysis(cfg, data_TFR);

% last trials
cfg.trials = find((data_TFR.trialinfo(:,1)==35 | data_TFR.trialinfo(:,1)==36)& ...
             data_TFR.trialinfo(:,2) == 31);
cfg.keeptrials  = 'no';
last_TFR   = ft_freqanalysis(cfg, data_TFR);
% keep trials for single trial analysis
cfg.keeptrials  = 'yes';
last_TFR_trial  = ft_freqanalysis(cfg, data_TFR);

% ID trials
cfg.trials = find((data_TFR.trialinfo(:,1)==241 | data_TFR.trialinfo(:,1)==242)& ...
             data_TFR.trialinfo(:,2) == 31);
cfg.keeptrials  = 'no';
ID_TFR     = ft_freqanalysis(cfg, data_TFR);
% keep trials for single trial analysis
cfg.keeptrials  = 'yes';
ID_TFR_trial  = ft_freqanalysis(cfg, data_TFR);

% ED trials
cfg.trials = find((data_TFR.trialinfo(:,1)==251 | data_TFR.trialinfo(:,1)==252)& ...
             data_TFR.trialinfo(:,2) == 31);
cfg.keeptrials  = 'no';                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    R     = ft_freqanalysis(cfg, data_TFR);
ED_TFR     = ft_freqanalysis(cfg, data_TFR);
% keep trials for single trial analysis
cfg.keeptrials  = 'yes';
ED_TFR_trial  = ft_freqanalysis(cfg, data_TFR);

% all trials with correct responses, exclude practice trials and trials
% after break
cfg.trials = find(data_TFR.trialinfo(:,2) == 31);
cfg.keeptrials  = 'no';
alltrl_TFR = ft_freqanalysis(cfg, data_TFR);
% keep trials for single trial analysis
cfg.keeptrials  = 'yes';
alltrl_TFR_trial  = ft_freqanalysis(cfg, data_TFR);

% for TFR: normalize each condition using the average of all trials
cfg               = [];
cfg.baseline      = [-0.4 -0.2];
cfg.baselinetype  =  'db';
cfg.parameter     = 'powspctrm';
shift_TFR_norm    = ft_freqbaseline_from_avg(cfg, shift_TFR , alltrl_TFR);
repeat1_TFR_norm  = ft_freqbaseline_from_avg(cfg, repeat1_TFR , alltrl_TFR);
repeat2_TFR_norm  = ft_freqbaseline_from_avg(cfg, repeat2_TFR , alltrl_TFR);
last_TFR_norm    = ft_freqbaseline_from_avg(cfg, last_TFR , alltrl_TFR);
ID_TFR_norm      = ft_freqbaseline_from_avg(cfg, ID_TFR , alltrl_TFR);
ED_TFR_norm      = ft_freqbaseline_from_avg(cfg, ED_TFR , alltrl_TFR);
alltrl_TFR_norm  = ft_freqbaseline_from_avg(cfg, alltrl_TFR , alltrl_TFR);

% save subject data
save(strcat(dirs.output_dir, '2_Subj_Avg\', subj, '_subj_avg_stimpres_TFR.mat'), "repeat1_TFR", "repeat2_TFR", "shift_TFR", ...
    "last_TFR", "ID_TFR", "ED_TFR", "alltrl_TFR", "repeat1_TFR_norm", "repeat2_TFR_norm", "shift_TFR_norm", "last_TFR_norm", ...
    "ID_TFR_norm", "ED_TFR_norm", "alltrl_TFR_norm");
save(strcat(dirs.output_dir, '2_Subj_Avg\', subj, '_subj_avg_stimpres_TFR_trial.mat'), "repeat1_TFR_trial", "repeat2_TFR_trial",...
    "shift_TFR_trial", "last_TFR_trial", "ID_TFR_trial", "ED_TFR_trial", "alltrl_TFR_trial");

% clear unnecessary variables
clear repeat1_TFR repeat2_TFR shift_TFR last_TFR ID_TFR ED_TFR  alltrl_TFR repeat1_TFR_norm repeat2_TFR_norm shift_TFR_norm last_TFR_norm ...
    ID_TFR_norm ED_TFR_norm  alltrl_TFR_norm data_TFR

fprintf(strcat('Subject', subj, ' succesful!\n\n'));
end
%% 
%----------------------------------------------------------------------
%            3 - Generate struct array to save all subject averages
%----------------------------------------------------------------------
for i = 1:height(subs)
    subj = subs{i};
    load(strcat(dirs.output_dir, '2_Subj_Avg\', subj, '_subj_avg_stimpres_TFR.mat'))
    % save the struct arrays in cell array for all subjects
    % using the normalized data
    all_repeat1_stimpres_TFR{i} = repeat1_TFR_norm;
    all_repeat2_stimpres_TFR{i} = repeat2_TFR_norm;
    all_shift_stimpres_TFR{i}  = shift_TFR_norm;
    all_last_stimpres_TFR{i}   = last_TFR_norm;
    all_ID_stimpres_TFR{i}     = ID_TFR_norm;
    all_ED_stimpres_TFR{i}     = ED_TFR_norm;
    all_alltrl_stimpres_TFR{i} = alltrl_TFR_norm;
    % for the purposes of baseline comparison, we take the pure alltrl
    % before normalization
    all_alltrl_stimpres_TFR_no_norm{i} = alltrl_TFR;
end

% save the struct array
save(strcat(dirs.output_dir, '2_Subj_Avg\all_subj_avg_stimpres_TFR'),"all_repeat1_stimpres_TFR","all_repeat2_stimpres_TFR", ...
    "all_shift_stimpres_TFR", "all_last_stimpres_TFR", "all_ID_stimpres_TFR", "all_ED_stimpres_TFR", ...
    "all_alltrl_stimpres_TFR", "all_alltrl_stimpres_TFR_no_norm");

% separate arrays by age group
young_repeat1_stimpres_TFR        = all_repeat1_stimpres_TFR(is_young);
old_repeat1_stimpres_TFR          = all_repeat1_stimpres_TFR(~is_young);
young_repeat2_stimpres_TFR        = all_repeat2_stimpres_TFR(is_young);
old_repeat2_stimpres_TFR          = all_repeat2_stimpres_TFR(~is_young);
young_shift_stimpres_TFR          = all_shift_stimpres_TFR(is_young);
old_shift_stimpres_TFR            = all_shift_stimpres_TFR(~is_young);
young_ID_stimpres_TFR             = all_ID_stimpres_TFR(is_young);
old_ID_stimpres_TFR               = all_ID_stimpres_TFR(~is_young);
young_ED_stimpres_TFR             = all_ED_stimpres_TFR(is_young);
old_ED_stimpres_TFR               = all_ED_stimpres_TFR(~is_young);
young_alltrl_stimpres_TFR         = all_alltrl_stimpres_TFR(is_young);
old_alltrl_stimpres_TFR           = all_alltrl_stimpres_TFR(~is_young);
young_last_stimpres_TFR           = all_last_stimpres_TFR(is_young);
old_last_stimpres_TFR             = all_last_stimpres_TFR(~is_young);
young_alltrl_stimpres_TFR_no_norm = all_alltrl_stimpres_TFR_no_norm(is_young);
old_alltrl_stimpres_TFR_no_norm   = all_alltrl_stimpres_TFR_no_norm(~is_young);

% save the separated struct arrays
save(strcat(dirs.output_dir, '2_Subj_Avg\young_subj_avg_stimpres_TFR'), "young_repeat1_stimpres_TFR", "young_repeat2_stimpres_TFR", ...
    "young_shift_stimpres_TFR", "young_last_stimpres_TFR", "young_ID_stimpres_TFR", "young_ED_stimpres_TFR", ...
    "young_alltrl_stimpres_TFR", "young_alltrl_stimpres_TFR_no_norm");

save(strcat(dirs.output_dir, '2_Subj_Avg\old_subj_avg_stimpres_TFR'), "old_repeat1_stimpres_TFR", "old_repeat2_stimpres_TFR", ...
    "old_shift_stimpres_TFR", "old_last_stimpres_TFR", "old_ID_stimpres_TFR", "old_ED_stimpres_TFR", ...
    "old_alltrl_stimpres_TFR", "old_alltrl_stimpres_TFR_no_norm");

% clear unnecasarry variables
clear all_repeat1_stimpres_TFR all_repeat2_stimpres_TFR all_shift_stimpres_TFR all_ID_stimpres_TFR all_ED_stimpres_TFR ...
    all_alltrl_stimpres_TFR all_last_stimpres_TFR all_alltrl_stimpres_TFR_no_norm 

%%
%----------------------------------------------------------------------
%            4 - Calculate grand average TFR for each condition
%---------------------------------------------------------------------- 
% TFR
% load subject averages
load(strcat(dirs.output_dir, '2_Subj_Avg\young_subj_avg_stimpres_TFR'));
load(strcat(dirs.output_dir, '2_Subj_Avg\old_subj_avg_stimpres_TFR'));
load(strcat(dirs.output_dir, '2_Subj_Avg\all_subj_avg_stimpres_TFR'));

% calculate grand average for each condition
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'powspctrm';
cfg.keepindividual = 'yes';

TFR.grandavg                      = ft_freqgrandaverage(cfg, all_alltrl_stimpres_TFR{:});
TFR.grandavg_young                = ft_freqgrandaverage(cfg, young_alltrl_stimpres_TFR{:});
TFR.grandavg_old                  = ft_freqgrandaverage(cfg, old_alltrl_stimpres_TFR{:});
TFR.grandavg_repeat1              = ft_freqgrandaverage(cfg, all_repeat1_stimpres_TFR{:});
TFR.grandavg_repeat2              = ft_freqgrandaverage(cfg, all_repeat2_stimpres_TFR{:});
TFR.grandavg_shift                = ft_freqgrandaverage(cfg, all_shift_stimpres_TFR{:});
TFR.grandavg_ID                   = ft_freqgrandaverage(cfg, all_ID_stimpres_TFR{:});
TFR.grandavg_ED                   = ft_freqgrandaverage(cfg, all_ED_stimpres_TFR{:});
TFR.grandavg_last                 = ft_freqgrandaverage(cfg, all_last_stimpres_TFR{:});
TFR.grandavg_repeat1_young        = ft_freqgrandaverage(cfg, young_repeat1_stimpres_TFR{:});
TFR.grandavg_repeat1_old          = ft_freqgrandaverage(cfg, old_repeat1_stimpres_TFR{:});
TFR.grandavg_repeat2_young        = ft_freqgrandaverage(cfg, young_repeat2_stimpres_TFR{:});
TFR.grandavg_repeat2_old          = ft_freqgrandaverage(cfg, old_repeat2_stimpres_TFR{:});
TFR.grandavg_shift_young          = ft_freqgrandaverage(cfg, young_shift_stimpres_TFR{:});
TFR.grandavg_shift_old            = ft_freqgrandaverage(cfg, old_shift_stimpres_TFR{:});
TFR.grandavg_ID_young             = ft_freqgrandaverage(cfg, young_ID_stimpres_TFR{:});
TFR.grandavg_ID_old               = ft_freqgrandaverage(cfg, old_ID_stimpres_TFR{:});
TFR.grandavg_ED_young             = ft_freqgrandaverage(cfg, young_ED_stimpres_TFR{:});                 
TFR.grandavg_ED_old               = ft_freqgrandaverage(cfg, old_ED_stimpres_TFR{:});
TFR.grandavg_last_young           = ft_freqgrandaverage(cfg, young_last_stimpres_TFR{:});
TFR.grandavg_last_old             = ft_freqgrandaverage(cfg, old_last_stimpres_TFR{:});
TFR.grandavg_alltrl_no_norm_young = ft_freqgrandaverage(cfg, young_alltrl_stimpres_TFR_no_norm{:});
TFR.grandavg_alltrl_no_norm_old   = ft_freqgrandaverage(cfg, old_alltrl_stimpres_TFR_no_norm{:});
%
save(strcat(dirs.output_dir, '3_Grandavg\grandavg_TFR_stimpres.mat'),"TFR", '-v7.3');

%%
%----------------------------------------------------------------------
%            5 - Collect average power for specific frequencies
%---------------------------------------------------------------------- 

% load grandaverages
load(strcat(dirs.output_dir, '3_Grandavg\grandavg_TFR_stimpres.mat'));

% define cfg
cfg = [];
%cfg.baseline     = [-0.30 -0.25];
%cfg.baselinetype = 'absolute';
cfg.xlim         = [-0.25 1.5];
cfg.showlabels   = 'yes';
cfg.layout       = 'acticap-64ch-standard2-EMCO.mat';
cfg.colorbar     = 'yes';

% theta
% check if have differences in latency
cfg.ylim         = [4 8];
cfg.zlim         = [-2.5 2.5];
figure; ft_multiplotTFR(cfg, TFR.grandavg);
figure; ft_multiplotTFR(cfg, TFR.grandavg_young);
figure; ft_multiplotTFR(cfg, TFR.grandavg_old);
% topoplot
cfg.xlim = [-0.1:0.1:0.7];
cfg.colorbar     = 'yes';
figure; ft_topoplotTFR(cfg,TFR.grandavg);

% calculate average theta power for young and old on the FCz channel to see
% if we have a differnt timeline
theta_power_fr       = mean(squeeze(mean(TFR.grandavg.powspctrm(:, 11, 4:8, :))));
theta_power_young_fr = mean(squeeze(mean(TFR.grandavg_young.powspctrm(:, 11, 4:8, :))));
theta_power_old_fr   = mean(squeeze(mean(TFR.grandavg_old.powspctrm(:, 11, 4:8, :))));
% calculate average theta power for young and old on the FCz channel to see
% if we have a differnt timeline
theta_power_p       = mean(squeeze(mean(TFR.grandavg.powspctrm(:, 63, 4:8, :))));
theta_power_young_p = mean(squeeze(mean(TFR.grandavg_young.powspctrm(:, 63, 4:8, :))));
theta_power_old_p   = mean(squeeze(mean(TFR.grandavg_old.powspctrm(:, 63, 4:8, :))));

figure()
plot(TFR.grandavg.time, theta_power_fr)
xlim([-0.25 2]);
ylim([0 4]);
hold on
plot(TFR.grandavg.time, theta_power_p)
figure()
plot(TFR.grandavg.time, theta_power_young_fr)
hold on
plot(TFR.grandavg.time, theta_power_young_p)
xlim([-0.25 2]);
ylim([0 4]);
figure()
plot(TFR.grandavg.time, theta_power_old_fr)
hold on
plot(TFR.grandavg.time, theta_power_old_p)
xlim([-0.25 2]);
ylim([0 4]);

% prepare cfg
cfg = [];
cfg.channel     = {'F1', 'F2', 'Fz', 'FC1', 'FC2', 'FCz'}; % {'Fz'};
cfg.avgoverchan = 'yes';
cfg.latency     =  [0.25 0.50];
cfg.avgovertime = 'yes';
cfg.frequency   =  [4 8];
cfg.avgoverfreq =  'yes';
% young
for i = 1:numel(young_repeat1_stimpres_TFR)
fprintf('************************\n%s\n************************\n', subj_young{i});
young_repeat1_stimpres{i} = ft_selectdata(cfg, young_repeat1_stimpres_TFR{i});
young_repeat2_stimpres{i} = ft_selectdata(cfg, young_repeat2_stimpres_TFR{i});
young_shift_stimpres{i} = ft_selectdata(cfg, young_shift_stimpres_TFR{i});
young_ID_stimpres{i} = ft_selectdata(cfg, young_ID_stimpres_TFR{i});
young_ED_stimpres{i} = ft_selectdata(cfg, young_ED_stimpres_TFR{i});
young_last_stimpres{i} = ft_selectdata(cfg, young_last_stimpres_TFR{i});
young_alltrl_stimpres{i} = ft_selectdata(cfg, young_alltrl_stimpres_TFR{i});
end
% old
for i = 1:numel(old_repeat1_stimpres_TFR)
fprintf('************************\n%s\n************************\n', subj_old{i});
old_repeat1_stimpres{i} = ft_selectdata(cfg, old_repeat1_stimpres_TFR{i});
old_repeat2_stimpres{i} = ft_selectdata(cfg, old_repeat2_stimpres_TFR{i});
old_shift_stimpres{i} = ft_selectdata(cfg, old_shift_stimpres_TFR{i});
old_ID_stimpres{i} = ft_selectdata(cfg, old_ID_stimpres_TFR{i});
old_ED_stimpres{i} = ft_selectdata(cfg, old_ED_stimpres_TFR{i});
old_last_stimpres{i} = ft_selectdata(cfg, old_last_stimpres_TFR{i});
old_alltrl_stimpres{i} = ft_selectdata(cfg, old_alltrl_stimpres_TFR{i});
end

% pool all subject averages together using the grandaverage function
% for now am just focusing on repeat and ID and ED trials
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'powspctrm';
cfg.keepindividual = 'yes';

theta250.young_repeat1 = ft_freqgrandaverage(cfg, young_repeat1_stimpres{:});
theta250.young_repeat2 = ft_freqgrandaverage(cfg, young_repeat2_stimpres{:});
theta250.young_last    = ft_freqgrandaverage(cfg, young_last_stimpres{:});
theta250.young_ID      = ft_freqgrandaverage(cfg, young_ID_stimpres{:});
theta250.young_ED      = ft_freqgrandaverage(cfg, young_ED_stimpres{:});
theta250.young_alltrl  = ft_freqgrandaverage(cfg, young_alltrl_stimpres{:});
theta250.old_repeat1   = ft_freqgrandaverage(cfg, old_repeat1_stimpres{:});
theta250.old_repeat2   = ft_freqgrandaverage(cfg, old_repeat2_stimpres{:});
theta250.old_last      = ft_freqgrandaverage(cfg, old_last_stimpres{:});
theta250.old_ID        = ft_freqgrandaverage(cfg, old_ID_stimpres{:});
theta250.old_ED        = ft_freqgrandaverage(cfg, old_ED_stimpres{:});
theta250.old_alltrl    = ft_freqgrandaverage(cfg, old_alltrl_stimpres{:});

% Showcase average theta for all trials
fprintf('IDED young:\nMean Theta (250 ms - 500 ms):\n %.2f +- %.2f dB\n', mean(theta250.young_alltrl.powspctrm), std(theta250.young_alltrl.powspctrm))
fprintf('IDED old:\nMean Theta (250 ms - 500 ms):\n %.2f +- %.2f dB\n', mean(theta250.old_alltrl.powspctrm), std(theta250.old_alltrl.powspctrm))

% save variables
save(strcat(dirs.output_dir, '4_Stats\stats_TFR_theta250.mat'), "young_repeat1_stimpres", "young_repeat2_stimpres", "young_shift_stimpres", ...
    "young_ID_stimpres", "young_ED_stimpres", "young_last_stimpres", "young_alltrl_stimpres", "old_repeat1_stimpres", "old_repeat2_stimpres", ...
    "old_shift_stimpres", "old_ID_stimpres", "old_ED_stimpres", "old_last_stimpres", "old_alltrl_stimpres","theta250");

%%
%----------------------------------------------------------------------
%            6 - Generate matrix for ANOVA in RStudio
%---------------------------------------------------------------------- 
% load variables
load(strcat(dirs.output_dir, '4_Stats\stats_TFR_theta250.mat'), "theta250")

% Theta 250 to 500 ms
% generate matrix for ANOVA
% depending variable
% dv = [theta250.young_repeat1.powspctrm; theta250.old_repeat1.powspctrm; theta250.young_ID.powspctrm; ...
%     theta250.old_ID.powspctrm; theta250.young_ED.powspctrm; theta250.old_ED.powspctrm];

dv = [theta250.young_repeat2.powspctrm; theta250.old_repeat2.powspctrm; theta250.young_ID.powspctrm; ...
    theta250.old_ID.powspctrm; theta250.young_ED.powspctrm; theta250.old_ED.powspctrm];
% subject number
subj_num = [1:numel(subs) 1:numel(subs) 1:numel(subs)]';
% generate table
X = [dv subj_num];
T = array2table(X);
% Assign headings to the table
T.Properties.VariableNames = {'dv', 'subj_num'};
% add between and within variables
T.between = [repmat('young', sum(is_young), 1); repmat(' old ', sum(~is_young), 1); repmat('young', sum(is_young), 1); ...
    repmat(' old ', sum(~is_young), 1); repmat('young', sum(is_young), 1); repmat(' old ', sum(~is_young), 1)];
T.within = [repmat('repeat', sum(is_young), 1); repmat('repeat', sum(~is_young), 1); repmat('  ID  ', sum(is_young), 1); ...
    repmat('  ID  ', sum(~is_young), 1); repmat('  ED  ', sum(is_young), 1); repmat('  ED  ', sum(~is_young), 1)];
T.gender =  [gender(is_young)' - 1; gender(~is_young)' - 1; gender(is_young)' - 1; gender(~is_young)' - 1; ...
    gender(is_young)' - 1; gender(~is_young)' - 1;];
% save table as csv
writetable(T, strcat(dirs.output_dir, '4_Stats\stat_TFR_theta250.csv'));

% Theta 250 to 500 ms - all trials
% generate matrix for ANOVA
% depending variable
dv = [theta250.young_alltrl.powspctrm; theta250.old_alltrl.powspctrm];
% subject number
subj_num = [1:numel(subs)]';
% generate table
X = [dv subj_num];
T = array2table(X);
% Assign headings to the table
T.Properties.VariableNames = {'dv', 'subj_num'};
% add between and within variables
T.between = [repmat('young', sum(is_young), 1); repmat(' old ', sum(~is_young), 1)];
T.age    = [subj_info.age(is_young); subj_info.age(~is_young)];
T.gender =  [gender(is_young)' - 1; gender(~is_young)' - 1];
% save table as csv
writetable(T, strcat(dirs.output_dir, '4_Stats\stat_TFR_theta250_alltrl.csv'));

%% 6 - Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load grandaverages
load(strcat(dirs.output_dir, '3_Grandavg\grandavg_TFR_stimpres.mat'));

% theta250
cfg = [];
cfg.channel     = {'F1', 'F2', 'Fz', 'FC1', 'FC2', 'FCz'};
cfg.latency     =  [-0.25 1];
cfg.frequency   =  [4 8];
cfg.avgoverfreq = 'yes';
cfg.avgoverchan = 'yes';

theta250.grandavg              = ft_selectdata(cfg, TFR.grandavg);
theta250.grandavg_repeat1_young = ft_selectdata(cfg, TFR.grandavg_repeat1_young);
theta250.grandavg_repeat2_young = ft_selectdata(cfg, TFR.grandavg_repeat2_young);
theta250.grandavg_ID_young     = ft_selectdata(cfg, TFR.grandavg_ID_young);
theta250.grandavg_ED_young     = ft_selectdata(cfg, TFR.grandavg_ED_young);
theta250.grandavg_repeat1_old   = ft_selectdata(cfg, TFR.grandavg_repeat1_old);
theta250.grandavg_repeat2_old   = ft_selectdata(cfg, TFR.grandavg_repeat2_old);
theta250.grandavg_ID_old       = ft_selectdata(cfg, TFR.grandavg_ID_old);
theta250.grandavg_ED_old       = ft_selectdata(cfg, TFR.grandavg_ED_old);

% collect everything together to make plotting easier
theta250_young = squeeze([theta250.grandavg_repeat2_young.powspctrm; theta250.grandavg_ID_young.powspctrm; ...
    theta250.grandavg_ED_young.powspctrm]);
theta250_old = squeeze([theta250.grandavg_repeat2_old.powspctrm; theta250.grandavg_ID_old.powspctrm; ...
    theta250.grandavg_ED_old.powspctrm]);

% generating the group variable
within_val={'repeat' 'ID' 'ED'};

% young
ind = [repmat(1, sum(is_young), 1); repmat(2, sum(is_young), 1); repmat(3, sum(is_young), 1)];
within = within_val(ind);

figure()
g = gramm('x', theta250.grandavg.time*1000, 'y', theta250_young, 'linestyle', within, 'lightness', within);
g.geom_polygon('x',{[250 500 500 250]}, 'y', {[-0.5 -0.5 2.5 2.5]}, 'alpha', 0.1)
g.stat_summary('type', 'ci'); % plot with 95% confidence interval
g.set_title(sprintf('Theta (4-8 Hz) young\n(F1, F2, Fz, FC1, FC2, FCz)'));
g.set_names('color','Condition','x','Time (ms)','y','Power (dB)', 'linestyle', 'Condition', 'lightness','Condition');
g.axe_property('XLim',[-250 1000],'YLim',[-0.5 2.5]);
g.set_color_options('lightness_range',[80 30], 'legend','merge');
g.set_order_options('linestyle',0, 'lightness', 0);
g.draw();

% old
ind = [repmat(1, sum(~is_young), 1); repmat(2, sum(~is_young), 1); repmat(3, sum(~is_young), 1)];
within = within_val(ind);
figure()
g = gramm('x', theta250.grandavg.time*1000, 'y', theta250_old, 'linestyle', within, 'lightness', within);
g.geom_polygon('x',{[250 500 500 250]}, 'y', {[-0.5 -0.5 2.5 2.5]}, 'alpha', 0.1)
g.stat_summary('type', 'ci'); % plot with 95% confidence interval
g.set_title(sprintf('Theta (4-8 Hz) old\n(F1, F2, Fz, FC1, FC2, FCz)'));
g.set_names('color','Condition','x','Time (ms)','y','Power (dB)', 'linestyle', 'Condition', 'lightness','Condition');
g.axe_property('XLim',[-250 1000],'YLim',[-0.5 2.5]);
g.set_color_options('lightness_range',[80 30],'hue_range',[175 175], 'legend','merge');
g.set_order_options('linestyle',0, 'lightness', 0);
g.draw();

% plot theta 250 topography

% define cfg
cfg = [];
cfg.baseline     = [-0.30 -0.25];
cfg.layout       = 'acticap-64ch-standard2-EMCO.mat';
cfg.colorbar     = 'yes';
cfg.ylim         = [4 8];
cfg.zlim         = [-3.5 3.5];
% topoplot
cfg.xlim = [0.25, 0.5];
cfg.colorbar     = 'yes';
figure; ft_topoplotTFR(cfg,TFR.grandavg_repeat2_young);
figure; ft_topoplotTFR(cfg,TFR.grandavg_ID_young);
figure; ft_topoplotTFR(cfg,TFR.grandavg_ED_young);
figure; ft_topoplotTFR(cfg,TFR.grandavg_repeat2_old);
figure; ft_topoplotTFR(cfg,TFR.grandavg_ID_old);
figure; ft_topoplotTFR(cfg,TFR.grandavg_ED_old);

%% Create movie for presentations

% load grandaverages
load(strcat(dirs.output_dir, '3_Grandavg\grandavg_TFR_stimpres.mat'));

% theta
% Create png for every frame and save it
cfg = [];
cfg.layout       = 'acticap-64ch-standard2-EMCO.mat';
cfg.ylim         = [4 8];
cfg.zlim         = [-2.5 2.5];
cfg.xlim         = [time time];
cfg.colorbar     = 'yes';
cfg.comment      = 'no';
cfg.interactive  = 'no';

k = 1;
for time = -0.400:0.05:1
cfg.xlim         = [time time];

% overall grandaverages
% ft_topoplotTFR(cfg,TFR.grandavg_young);
% saveas(gcf, sprintf('%s8_movies/Theta_topo_young_%i.png', dirs.output_dir, k));
% ft_topoplotTFR(cfg,TFR.grandavg_old);
% saveas(gcf, sprintf('%s8_movies/Theta_topo_old_%i.png', dirs.output_dir, k));

% young separated by condition
ft_topoplotTFR(cfg,TFR.grandavg_repeat2_young);
saveas(gcf, sprintf('%s8_movies/Theta_topo_young_repeat_%i.png', dirs.output_dir, k));
ft_topoplotTFR(cfg,TFR.grandavg_ID_young);
saveas(gcf, sprintf('%s8_movies/Theta_topo_young_ID_%i.png', dirs.output_dir, k));
ft_topoplotTFR(cfg,TFR.grandavg_ED_young);
saveas(gcf, sprintf('%s8_movies/Theta_topo_young_ED_%i.png', dirs.output_dir, k));

% old separated by condition
ft_topoplotTFR(cfg,TFR.grandavg_repeat2_old);
saveas(gcf, sprintf('%s8_movies/Theta_topo_old_repeat_%i.png', dirs.output_dir, k));
ft_topoplotTFR(cfg,TFR.grandavg_ID_old);
saveas(gcf, sprintf('%s8_movies/Theta_topo_old_ID_%i.png', dirs.output_dir, k));
ft_topoplotTFR(cfg,TFR.grandavg_ED_old);
saveas(gcf, sprintf('%s8_movies/Theta_topo_old_ED_%i.png', dirs.output_dir, k));

k = k + 1;
end
