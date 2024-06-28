% EEG Analysis script for Fieldtrip for single subjects and grand averages
%
% This script opens preprocessed eeg data from brain vision analyzer,
% and gives out basic analytics
% Based on tutorial from fieldtriptoolbox.org
%
% Written by: Margarita Darna
% Created on: 03. August 2022
% Last modified on: 30. November 2023

%% --------------------------------------------------------------------
%               0 - Prepare workspace and directories
%----------------------------------------------------------------------
clear;clc;
% Setting up needed directories
dirs = {};
dirs.proj_dir       = 'your_directory';
dirs.dt_dir         = strcat (dirs.proj_dir, 'Data/');
dirs.exp_dir        = strcat (dirs.proj_dir, 'IDED_Analysis/');
dirs.raw_dt_dir     = strcat(dirs.dt_dir, 'Raw_data/');
dirs.derived_dt_dir = strcat(dirs.dt_dir, 'Derived_data/IDED/');
dirs.analysis_dir   = strcat(dirs.exp_dir, 'Functions/');
dirs.output_dir     = strcat(dirs.exp_dir, 'Output/');
dirs.prepr_dir      = strcat(dirs.output_dir, '1_Preprocessing/');

% adding analysis path and subfolders
addpath(genpath(dirs.analysis_dir));

% making sure that fieldtrip is called correctly
ft_defaults

% load the subjects
subj_info = readtable(strcat(dirs.exp_dir, 'Protocol.xlsx'));
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

%% --------------------------------------------------------------------
%               1 - Load single subject EEG and prepare epochs
%----------------------------------------------------------------------
for i = 1 : height(subs)
   
    % define subject
    subj = subs{i};

    % print information in command window
    fprintf("************************\nStarting %s\n************************\n", subj)
    
    % load subjectdata as they were created in script IDED_ERP_Analysis
    load(strcat(dirs.prepr_dir, subj, '_IDED_ERP_ft_prepr_1.mat'), "subjectdata")

    % redefine trials to be response locked - use negative offset
    cfg = [];
    cfg.offset = -subjectdata.data_ERP.trialinfo(:,4); % negative offset, so that the axis moves to the front
    data_ERP_response_wholetimewindow = ft_redefinetrial(cfg, subjectdata.data_ERP);

    % choose the same time window for all trials
    cfg = [];
    cfg.toilim = [-1 1.5]; % prestim = 1 s; poststim = 1.5 s
    data_ERP_response = ft_redefinetrial(cfg, data_ERP_response_wholetimewindow );

    % save subject specific information
    subjectdata.data_ERP_response = data_ERP_response;

    fprintf('Saving %s ... ', subj)
    save(strcat(dirs.prepr_dir, subj, '_IDED_ERP_ft_prepr_1.mat'), "subjectdata", '-v7.3');
    fprintf('succesful!\n\n')

    % remove unnecesary variables
    clear subjectdata data_ERP_response_wholetimewindow  data_ERP_response
end

%% --------------------------------------------------------------------
%          2a - Calculate subject average for each condition
%----------------------------------------------------------------------
for i = 1:height(subs)
    % ensure that we don't mix up subjects by removing subject data
    clear subjectdata

    % define subject
    subj = subs{i};

    % print information in command window
    fprintf("************************\nStarting %s\n************************\n", subj)

    % check if the process has already been done
    %subj_avg_done = exist(strcat(dirs.output_dir, '2_Subj_Avg\', subj, '_subj_avg_response_ERP.mat'), 'file');
    subj_avg_done = 1;
    if subj_avg_done == 2
        fprintf('%s_IDED_ERP_ft_prepr_1.mat File already exists\n\n', subj)
    else
        % load mat file that was prepared in the previous section
        load(strcat(dirs.prepr_dir, subj, '_IDED_ERP_ft_prepr_1.mat'), "subjectdata");
        data_ERP = subjectdata.data_ERP_response;

        % perform timelockanalysis (ERP) for each trial type
        % here we only choose trials with correct responses: data_ERP.trialinfo(:,2) == 31

        % shift trials
        shift_ind  = (data_ERP.trialinfo(:,1)==241 | data_ERP.trialinfo(:,1)==242 | ...
            data_ERP.trialinfo(:,1)==251 | data_ERP.trialinfo(:,1)==252) ...
            & data_ERP.trialinfo(:,2) == 31;
        cfg        = [];
        cfg.trials = find(shift_ind);
        shift_ERP  = ft_timelockanalysis(cfg, data_ERP);

        % repeat trials - only choose the first trial after a shift - this will be
        % our first repeat: repeat1
        repeat1_ind  = [false; shift_ind(1:end - 1,:)];
        cfg         = [];
        cfg.trials  = find(repeat1_ind & (data_ERP.trialinfo(:,1)==231 |...
            data_ERP.trialinfo(:,1)==232) & data_ERP.trialinfo(:,2) == 31); % choose only correct responses
        repeat1_ERP  = ft_timelockanalysis(cfg, data_ERP);

        % repeat trials - only choose the second trial after a shift - this will be
        % our second repeat: repeat2
        repeat2_ind  = [false; false; shift_ind(1:end - 2,:)];
        cfg         = [];
        cfg.trials  = find(repeat2_ind & (data_ERP.trialinfo(:,1)==231 |...
            data_ERP.trialinfo(:,1)==232) & data_ERP.trialinfo(:,2) == 31); % choose only correct responses
        repeat2_ERP  = ft_timelockanalysis(cfg, data_ERP);

        % last trials
        cfg        = [];
        cfg.trials = find((data_ERP.trialinfo(:,1)==35 | data_ERP.trialinfo(:,1)==36)& ...
            data_ERP.trialinfo(:,2) == 31);
        last_ERP   = ft_timelockanalysis(cfg, data_ERP);

        % ID trials
        cfg        = [];
        cfg.trials = find((data_ERP.trialinfo(:,1)==241 | data_ERP.trialinfo(:,1)==242)& ...
            data_ERP.trialinfo(:,2) == 31);
        ID_ERP     = ft_timelockanalysis(cfg, data_ERP);

        % ED trials
        cfg        = [];
        cfg.trials = find((data_ERP.trialinfo(:,1)==251 | data_ERP.trialinfo(:,1)==252)& ...
            data_ERP.trialinfo(:,2) == 31);
        ED_ERP     = ft_timelockanalysis(cfg, data_ERP);

        % all trials with correct responses, exclude practice trials and trials
        % after break
        cfg        = [];
        cfg.trials = find(data_ERP.trialinfo(:,2) == 31);
        alltrl_ERP = ft_timelockanalysis(cfg, data_ERP);

        % these are some interactive plots to plot individual subject ERPS
        % it can be modified accordingly
        %plot ERP
        % cfg = [];
        % cfg.layout = 'acticap-64ch-standard2-EMCO';
        % cfg.interactive = 'yes';
        % cfg.showoutline = 'yes';
        % ft_multiplotER(cfg, repeat_ERP, shift_ERP)

        % save subject data
        save(strcat(dirs.output_dir, '2_Subj_Avg\', subj, '_subj_avg_response_ERP.mat'), "repeat1_ERP", "repeat2_ERP", ...
            "shift_ERP", "last_ERP", "ID_ERP", "ED_ERP", "alltrl_ERP", '-v7.3');

        % clear unnecessary variables
        clear repeat1_ERP repeat2_ERP shift_ERP last_ERP ID_ERP ED_ERP alltrl_ERP data_ERP

        fprintf(strcat('Subject', subj, ' succesful!\n\n'));
    end
end

%% --------------------------------------------------------------------
%        2b - Generate struct array tp save all subject averages
%----------------------------------------------------------------------
% add all subject information in one struct array
for i = 1:height(subs)

    % define subjecta
    subj = subs{i};

    % load subject average
    load(strcat(dirs.output_dir, '2_Subj_Avg\', subj, '_subj_avg_response_ERP.mat'))

    % save the struct arrays in cell array for all subjects
    % ERP
    all_repeat1_response_ERP{i} = repeat1_ERP;
    all_repeat2_response_ERP{i} = repeat2_ERP;
    all_shift_response_ERP{i}  = shift_ERP;
    all_last_response_ERP{i}   = last_ERP;
    all_ID_response_ERP{i}     = ID_ERP;
    all_ED_response_ERP{i}     = ED_ERP;
    all_alltrl_response_ERP{i} = alltrl_ERP;

    clear repeat1_ERP repeat2_ERP shift_ERP last_ERP ID_ERP ED_ERP alltrl_ERP data_ERP
end

% save the struct array
save(strcat(dirs.output_dir, '2_Subj_Avg\all_subj_avg_response_ERP'),...
    "all_repeat1_response_ERP",  "all_repeat2_response_ERP", "all_shift_response_ERP", "all_last_response_ERP", "all_ID_response_ERP", ...
    "all_ED_response_ERP", "all_alltrl_response_ERP", "subj_info", "tblstats");

% separate arrays by age group
young_repeat1_response_ERP = all_repeat1_response_ERP(is_young);
old_repeat1_response_ERP   = all_repeat1_response_ERP(~is_young);
young_repeat2_response_ERP = all_repeat2_response_ERP(is_young);
old_repeat2_response_ERP   = all_repeat2_response_ERP(~is_young);
young_shift_response_ERP  = all_shift_response_ERP(is_young);
old_shift_response_ERP    = all_shift_response_ERP(~is_young);
young_ID_response_ERP     = all_ID_response_ERP(is_young);
old_ID_response_ERP       = all_ID_response_ERP(~is_young);
young_ED_response_ERP     = all_ED_response_ERP(is_young);
old_ED_response_ERP       = all_ED_response_ERP(~is_young);
young_alltrl_response_ERP = all_alltrl_response_ERP(is_young);
old_alltrl_response_ERP   = all_alltrl_response_ERP(~is_young);
young_last_response_ERP   = all_last_response_ERP(is_young);
old_last_response_ERP     = all_last_response_ERP(~is_young);

% save the separated struct arrays
save(strcat(dirs.output_dir, '2_Subj_Avg\young_subj_avg_response_ERP'),"young_repeat1_response_ERP", ...
    "young_repeat2_response_ERP", "young_shift_response_ERP", "young_last_response_ERP", "young_ID_response_ERP", ...
    "young_ED_response_ERP", "young_alltrl_response_ERP", "subj_info", "tblstats");

save(strcat(dirs.output_dir, '2_Subj_Avg\old_subj_avg_response_ERP'), "old_repeat1_response_ERP", ...
    "old_repeat2_response_ERP", "old_shift_response_ERP", "old_last_response_ERP", "old_ID_response_ERP", ...
    "old_ED_response_ERP", "old_alltrl_response_ERP", "subj_info", "tblstats");

% clear unnecasarry variables
clear all_repeat1_response_ERP all_repeat2_response_ERP all_shift_response_ERP all_ID_response_ERP all_ED_response_ERP all_alltrl_response_ERP ...
    all_last_response_ERP

%% --------------------------------------------------------------------
%              3 - Calculate grand average ERP for each condition
%----------------------------------------------------------------------
% ERP
% load subject averages
load(strcat(dirs.output_dir, '2_Subj_Avg\young_subj_avg_response_ERP'));
load(strcat(dirs.output_dir, '2_Subj_Avg\old_subj_avg_response_ERP'));
load(strcat(dirs.output_dir, '2_Subj_Avg\all_subj_avg_response_ERP'));

% calculate grand average for each condition for the purposes of figures
cfg = [];
cfg.channel        = 1:65;
cfg.latency        = 'all';
cfg.parameter      = 'avg';
cfg.keepindividual = 'yes'; % important to keep to be able to calculate the confidence interval

ERP.grandavg               = ft_timelockgrandaverage(cfg, all_alltrl_response_ERP{:});
ERP.grandavg_young         = ft_timelockgrandaverage(cfg, young_alltrl_response_ERP{:});
ERP.grandavg_old           = ft_timelockgrandaverage(cfg, old_alltrl_response_ERP{:});
ERP.grandavg_repeat1       = ft_timelockgrandaverage(cfg, all_repeat1_response_ERP{:});
ERP.grandavg_repeat2       = ft_timelockgrandaverage(cfg, all_repeat2_response_ERP{:});
ERP.grandavg_shift         = ft_timelockgrandaverage(cfg, all_shift_response_ERP{:});
ERP.grandavg_repeat1_young = ft_timelockgrandaverage(cfg, young_repeat1_response_ERP{:});
ERP.grandavg_repeat1_old   = ft_timelockgrandaverage(cfg, old_repeat1_response_ERP{:});
ERP.grandavg_repeat2_young = ft_timelockgrandaverage(cfg, young_repeat2_response_ERP{:});
ERP.grandavg_repeat2_old   = ft_timelockgrandaverage(cfg, old_repeat2_response_ERP{:});
ERP.grandavg_shift_young   = ft_timelockgrandaverage(cfg, young_shift_response_ERP{:});
ERP.grandavg_shift_old     = ft_timelockgrandaverage(cfg, old_shift_response_ERP{:});
ERP.grandavg_ID_young      = ft_timelockgrandaverage(cfg, young_ID_response_ERP{:});
ERP.grandavg_ID_old        = ft_timelockgrandaverage(cfg, old_ID_response_ERP{:});
ERP.grandavg_ED_young      = ft_timelockgrandaverage(cfg, young_ED_response_ERP{:});
ERP.grandavg_ED_old        = ft_timelockgrandaverage(cfg, old_ED_response_ERP{:});
ERP.grandavg_last_young    = ft_timelockgrandaverage(cfg, young_last_response_ERP{:});
ERP.grandavg_last_old      = ft_timelockgrandaverage(cfg, old_last_response_ERP{:});

% save in matfile
save(strcat(dirs.output_dir, '3_Grandavg\grandavg_ERP_response.mat'), "ERP", '-v7.3');

%% --------------------------------------------------------------------
%          4- Find ERP components that differ from zero
%----------------------------------------------------------------------

% load grandaverages
load(strcat(dirs.output_dir, '3_Grandavg\grandavg_ERP_response.mat'));

% Visualize average ERP to find the peaks of interest
cfg = [];
cfg.showlabels  = 'yes';
cfg.layout      = 'acticap-64ch-standard2-EMCO.mat';
cfg.xlim        = [-0.5 0.5];
figure; ft_multiplotER(cfg, ERP.grandavg)

figure; ft_multiplotER(cfg, ERP.grandavg_young, ERP.grandavg_old)

%% --------------------------------------------------------------------
%          5- Collect average amplitudes in the chosen latencies
%----------------------------------------------------------------------
% load subject averages
load(strcat(dirs.output_dir, '2_Subj_Avg\young_subj_avg_response_ERP'));
load(strcat(dirs.output_dir, '2_Subj_Avg\old_subj_avg_response_ERP'));
load(strcat(dirs.output_dir, '2_Subj_Avg\all_subj_avg_response_ERP'));

% a) Motor preparation
% prepare the cfg
cfg = [];
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';
cfg.channel     = {'Pz', 'P1', 'P2', 'POz', 'PO3', 'PO4'};%{'Pz'};
cfg.latency_young     =  [-0.3 -0.1];
cfg.latency_old       =  [-0.3 -0.1];


for i = 1:numel(young_repeat1_response_ERP)
cfg.latency = cfg.latency_young;
young_repeat1_response{i} = ft_selectdata(cfg, young_repeat1_response_ERP{i});
young_repeat2_response{i} = ft_selectdata(cfg, young_repeat2_response_ERP{i});
young_shift_response{i} = ft_selectdata(cfg, young_shift_response_ERP{i});
young_ID_response{i} = ft_selectdata(cfg, young_ID_response_ERP{i});
young_ED_response{i} = ft_selectdata(cfg, young_ED_response_ERP{i});
young_last_response{i} = ft_selectdata(cfg, young_last_response_ERP{i});
young_alltrl_response{i} = ft_selectdata(cfg, young_alltrl_response_ERP{i});
end

for i = 1:numel(old_repeat1_response_ERP)
cfg.latency = cfg.latency_old;
old_repeat1_response{i} = ft_selectdata(cfg, old_repeat1_response_ERP{i});
old_repeat2_response{i} = ft_selectdata(cfg, old_repeat2_response_ERP{i});
old_shift_response{i} = ft_selectdata(cfg, old_shift_response_ERP{i});
old_ID_response{i} = ft_selectdata(cfg, old_ID_response_ERP{i});
old_ED_response{i} = ft_selectdata(cfg, old_ED_response_ERP{i});
old_last_response{i} = ft_selectdata(cfg, old_last_response_ERP{i});
old_alltrl_response{i} = ft_selectdata(cfg, old_alltrl_response_ERP{i});
end

% pool all subject averages together using the grandaverage function
% for now am just focusing on repeat and ID and ED trials
cfg = [];
cfg.channel        = 'all';
cfg.latency        = 'all';
cfg.parameter      = 'avg';
cfg.keepindividual = 'yes';

component_avg.young_repeat1 = ft_timelockgrandaverage(cfg, young_repeat1_response{:});
component_avg.young_repeat2 = ft_timelockgrandaverage(cfg, young_repeat2_response{:});
component_avg.young_ID     = ft_timelockgrandaverage(cfg, young_ID_response{:});
component_avg.young_ED     = ft_timelockgrandaverage(cfg, young_ED_response{:});
component_avg.young_last   = ft_timelockgrandaverage(cfg, young_last_response{:});
component_avg.old_repeat1   = ft_timelockgrandaverage(cfg, old_repeat1_response{:});
component_avg.old_repeat2   = ft_timelockgrandaverage(cfg, old_repeat2_response{:});
component_avg.old_ID       = ft_timelockgrandaverage(cfg, old_ID_response{:});
component_avg.old_ED       = ft_timelockgrandaverage(cfg, old_ED_response{:});
component_avg.old_last     = ft_timelockgrandaverage(cfg, old_last_response{:});

% generate matrix for ANOVA
% depending variable

% choose second repeat trial
dv = [component_avg.young_repeat2.individual; component_avg.old_repeat2.individual; component_avg.young_ID.individual; component_avg.old_ID.individual; ...
      component_avg.young_ED.individual; component_avg.old_ED.individual];
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
%T.gender =  [gender(is_young)' - 1; gender(~is_young)' - 1; gender(is_young)' - 1; gender(~is_young)' - 1; ...
%    gender(is_young)' - 1; gender(~is_young)' - 1;];
% save table as csv
writetable(T, strcat(dirs.output_dir, '4_Stats\stat_', 'ERP', '_', 'motor_preparation', '.csv'));



% b) Time after motor movement
% prepare the cfg
cfg = [];
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';
cfg.channel     = {'F1', 'F2', 'Fz', 'FC1', 'FC2', 'FCz'};%{'Pz'};
cfg.latency_young     =  [0 0.2];
cfg.latency_old       =  [0 0.2];

for i = 1:numel(young_repeat1_response_ERP)
cfg.latency = cfg.latency_young;
young_repeat1_response{i} = ft_selectdata(cfg, young_repeat1_response_ERP{i});
young_repeat2_response{i} = ft_selectdata(cfg, young_repeat2_response_ERP{i});
young_shift_response{i} = ft_selectdata(cfg, young_shift_response_ERP{i});
young_ID_response{i} = ft_selectdata(cfg, young_ID_response_ERP{i});
young_ED_response{i} = ft_selectdata(cfg, young_ED_response_ERP{i});
young_last_response{i} = ft_selectdata(cfg, young_last_response_ERP{i});
young_alltrl_response{i} = ft_selectdata(cfg, young_alltrl_response_ERP{i});
end

for i = 1:numel(old_repeat1_response_ERP)
cfg.latency = cfg.latency_old;
old_repeat1_response{i} = ft_selectdata(cfg, old_repeat1_response_ERP{i});
old_repeat2_response{i} = ft_selectdata(cfg, old_repeat2_response_ERP{i});
old_shift_response{i} = ft_selectdata(cfg, old_shift_response_ERP{i});
old_ID_response{i} = ft_selectdata(cfg, old_ID_response_ERP{i});
old_ED_response{i} = ft_selectdata(cfg, old_ED_response_ERP{i});
old_last_response{i} = ft_selectdata(cfg, old_last_response_ERP{i});
old_alltrl_response{i} = ft_selectdata(cfg, old_alltrl_response_ERP{i});
end

% pool all subject averages together using the grandaverage function
% for now am just focusing on repeat and ID and ED trials
cfg = [];
cfg.channel        = 'all';
cfg.latency        = 'all';
cfg.parameter      = 'avg';
cfg.keepindividual = 'yes';

component_avg.young_repeat1 = ft_timelockgrandaverage(cfg, young_repeat1_response{:});
component_avg.young_repeat2 = ft_timelockgrandaverage(cfg, young_repeat2_response{:});
component_avg.young_ID     = ft_timelockgrandaverage(cfg, young_ID_response{:});
component_avg.young_ED     = ft_timelockgrandaverage(cfg, young_ED_response{:});
component_avg.young_last   = ft_timelockgrandaverage(cfg, young_last_response{:});
component_avg.old_repeat1   = ft_timelockgrandaverage(cfg, old_repeat1_response{:});
component_avg.old_repeat2   = ft_timelockgrandaverage(cfg, old_repeat2_response{:});
component_avg.old_ID       = ft_timelockgrandaverage(cfg, old_ID_response{:});
component_avg.old_ED       = ft_timelockgrandaverage(cfg, old_ED_response{:});
component_avg.old_last     = ft_timelockgrandaverage(cfg, old_last_response{:});

% generate matrix for ANOVA
% depending variable

% choose second repeat trial
dv = [component_avg.young_repeat2.individual; component_avg.old_repeat2.individual; component_avg.young_ID.individual; component_avg.old_ID.individual; ...
      component_avg.young_ED.individual; component_avg.old_ED.individual];
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
%T.gender =  [gender(is_young)' - 1; gender(~is_young)' - 1; gender(is_young)' - 1; gender(~is_young)' - 1; ...
%    gender(is_young)' - 1; gender(~is_young)' - 1;];
% save table as csv
writetable(T, strcat(dirs.output_dir, '4_Stats\stat_', 'ERP', '_', 'response_0to200ms', '.csv'));

%% --------------------------------------------------------------------
%                       6 - Plot results
%----------------------------------------------------------------------

% Calculate average of posterior channels
% prepare the cfg
cfg = [];
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'no';
cfg.channel     = {'Pz', 'P1', 'P2', 'POz', 'PO3', 'PO4'};%{'Pz'};


young_repeat_post = ft_selectdata(cfg, ERP.grandavg_repeat2_young);
young_ID_post     = ft_selectdata(cfg, ERP.grandavg_ID_young);
young_ED_post = ft_selectdata(cfg, ERP.grandavg_ED_young);

old_repeat_post = ft_selectdata(cfg, ERP.grandavg_repeat2_old);
old_ID_post     = ft_selectdata(cfg, ERP.grandavg_ID_old);
old_ED_post = ft_selectdata(cfg, ERP.grandavg_ED_old);

% looking at frontal channels
cfg.channel     ={'F1', 'F2', 'Fz', 'FC1', 'FC2', 'FCz'};
young_repeat_fr = ft_selectdata(cfg, ERP.grandavg_repeat2_young);
young_ID_fr     = ft_selectdata(cfg, ERP.grandavg_ID_young);
young_ED_fr = ft_selectdata(cfg, ERP.grandavg_ED_young);

old_repeat_fr = ft_selectdata(cfg, ERP.grandavg_repeat2_old);
old_ID_fr     = ft_selectdata(cfg, ERP.grandavg_ID_old);
old_ED_fr = ft_selectdata(cfg, ERP.grandavg_ED_old);

% young
% generating the group variable
within_val={'repeat', '  ID  ' '  ED  '};
ind = [repmat(1, sum(is_young), 1); repmat(2, sum(is_young), 1); repmat(3, sum(is_young), 1)];
within = within_val(ind);

% posterior channels
y = squeeze([young_repeat_post.individual; young_ID_post.individual; young_ED_post.individual]);

figure()
g = gramm('x', young_repeat_post.time*1000, 'y', y, 'linestyle', within, 'lightness', within);
g.geom_polygon('x',{[-300 -100 -100 -300]}, 'y', {[-9 -9 9 9]}, 'alpha', 0.1)
g.stat_summary('type', 'ci'); % plot with 95% confidence interval
g.set_title(sprintf('Response Locked young\n(Pz, P1, P2, POz, PO3, PO4)'));
g.set_names('color','Condition','x','Time (ms)','y','Potential (uV)', 'linestyle', 'Condition', 'lightness','Condition');
g.axe_property('XLim',[-500 500], 'YLim',[-4 8]);
g.set_color_options('lightness_range',[80 30], 'legend','merge');
g.set_order_options('linestyle',0, 'lightness', 0);
g.draw();

% frontal channels
y = squeeze([young_repeat_fr.individual; young_ID_fr.individual; young_ED_fr.individual]);

figure()
g = gramm('x', young_repeat_post.time*1000, 'y', y, 'linestyle', within, 'lightness', within);
g.geom_polygon('x',{[0 200 200 0]}, 'y', {[-9 -9 9 9]}, 'alpha', 0.1)
g.stat_summary('type', 'ci'); % plot with 95% confidence interval
g.set_title(sprintf('Response Locked young\n(F1, F2, Fz, FC1, FC2, FCz)'));
g.set_names('color','Condition','x','Time (ms)','y','Potential (uV)', 'linestyle', 'Condition', 'lightness','Condition');
g.axe_property('XLim',[-500 500], 'YLim',[-4 8]);
g.set_color_options('lightness_range',[80 30], 'legend','merge');
g.set_order_options('linestyle',0, 'lightness', 0);
g.draw();

% old
% generating the group variable
within_val={'repeat', '  ID  ' '  ED  '};
ind = [repmat(1, sum(~is_young), 1); repmat(2, sum(~is_young), 1); repmat(3, sum(~is_young), 1)];
within = within_val(ind);

y = squeeze([old_repeat_post.individual; old_ID_post.individual; old_ED_post.individual]);

% posterior channels
figure()
g = gramm('x', old_repeat_post.time*1000, 'y', y, 'linestyle', within, 'lightness', within);
g.geom_polygon('x',{[-300 -100 -100 -300]}, 'y', {[-9 -9 9 9]}, 'alpha', 0.1)
g.stat_summary('type', 'ci'); % plot with 95% confidence interval
g.set_title(sprintf('Response Locked old\n(Pz, P1, P2, POz, PO3, PO4)'));
g.set_names('color','Condition','x','Time (ms)','y','Potential (uV)', 'linestyle', 'Condition', 'lightness','Condition');
g.axe_property('XLim',[-500 500], 'YLim',[-4 8]);
g.set_color_options('lightness_range',[80 30],'hue_range',[175 175], 'legend','merge');
g.set_order_options('linestyle',0, 'lightness', 0);
g.draw();

% frontal channels
y = squeeze([old_repeat_fr.individual; old_ID_fr.individual; old_ED_fr.individual]);

figure()
g = gramm('x', old_repeat_post.time*1000, 'y', y, 'linestyle', within, 'lightness', within);
g.geom_polygon('x',{[0 200 200 0]}, 'y', {[-9 -9 9 9]}, 'alpha', 0.1)
g.stat_summary('type', 'ci'); % plot with 95% confidence interval
g.set_title(sprintf('Response Locked old\n(F1, F2, Fz, FC1, FC2, FCz)'));
g.set_names('color','Condition','x','Time (ms)','y','Potential (uV)', 'linestyle', 'Condition', 'lightness','Condition');
g.axe_property('XLim',[-500 500], 'YLim',[-4 8]);
g.set_color_options('lightness_range',[80 30],'hue_range',[175 175], 'legend','merge');
g.set_order_options('linestyle',0, 'lightness', 0);
g.draw();


% %% performing t-test with moving time window
% load(strcat(dirs.output_dir, '2_Subj_Avg\young_subj_avg_response_ERP'));
% load(strcat(dirs.output_dir, '2_Subj_Avg\old_subj_avg_response_ERP'));
% load(strcat(dirs.output_dir, '2_Subj_Avg\all_subj_avg_response_ERP'));
% 
% % load grandaverages
% load(strcat(dirs.output_dir, '3_Grandavg\grandavg_ERP_response.mat'));
% 
% % define the parameters for the statistical comparison
% cfg = [];
% cfg.channel     = 'all';
% cfg.avgovertime = 'yes';
% cfg.parameter   = 'avg';
% cfg.method      = 'analytic';
% cfg.statistic   = 'ft_statfun_depsamplesT';
% cfg.alpha       = 0.05;
% cfg.correctm    = 'no';
% 
% Nsub = 39;
% cfg.design(1,:)  = [ones(1,Nsub) 2*ones(1,Nsub)];
% cfg.design(2,:)  = [1:Nsub 1:Nsub];
% cfg.ivar         = 1; % the 1st row in cfg.design contains the independent variable
% cfg.uvar         = 2; % the 2nd row in cfg.design contains the subject number
% 
% % perform moving window analysis
% % define window in which we perform analysis
% time_window = -0.5:0.01:0;
% % we will overlap the windows by 50 ms
% % prepare struct array
% stat = struct([]);
% for i = 1:numel(time_window)
% cfg.latency     = [time_window(i)-0.05 time_window(i)+0.05];
% stat{i} = ft_timelockstatistics(cfg, all_repeat2_response_ERP{:}, all_ED_response_ERP{:});
% end
% 
% for i = 1:numel(time_window)
% % make the plot
% cfg = [];
% cfg.style     = 'blank';
% cfg.layout    = 'acticap-64ch-standard2-EMCO.mat';
% cfg.highlight = 'on';
% cfg.highlightchannel = find(stat{i}.mask);
% cfg.comment   = 'no';
% figure; ft_topoplotER(cfg, ERP.grandavg_repeat2)
% title(sprintf('Time: %d:%d ms', time_window(i)*1000-50, time_window(i)*1000+50))
% end
% 
% % only young participants
% % define the parameters for the statistical comparison
% cfg = [];
% cfg.channel     = 'all';
% cfg.avgovertime = 'yes';
% cfg.parameter   = 'avg';
% cfg.method      = 'analytic';
% cfg.statistic   = 'ft_statfun_depsamplesT';
% cfg.alpha       = 0.05;
% cfg.correctm    = 'no';
% 
% Nsub = numel(subj_young);
% cfg.design(1,:)  = [ones(1,Nsub) 2*ones(1,Nsub)];
% cfg.design(2,:)  = [1:Nsub 1:Nsub];
% cfg.ivar         = 1; % the 1st row in cfg.design contains the independent variable
% cfg.uvar         = 2; % the 2nd row in cfg.design contains the subject number
% 
% % perform moving window analysis
% % define window in which we perform analysis
% time_window = -0.5:0.01:0;
% % we will overlap the windows by 50 ms
% % prepare struct array
% stat = struct([]);
% for i = 1:numel(time_window)
% cfg.latency     = [time_window(i)-0.05 time_window(i)+0.05];
% stat{i} = ft_timelockstatistics(cfg, young_repeat2_response_ERP{:}, young_ED_response_ERP{:});
% end
% 
% for i = 1:numel(time_window)
% % make the plot
% cfg = [];
% cfg.style     = 'blank';
% cfg.layout    = 'acticap-64ch-standard2-EMCO.mat';
% cfg.highlight = 'on';
% cfg.highlightchannel = find(stat{i}.mask);
% cfg.comment   = 'no';
% figure; ft_topoplotER(cfg, ERP.grandavg_repeat2_young)
% title(sprintf('Time: %d:%d ms', time_window(i)*1000-50, time_window(i)*1000+50))
% end
% 
% %%
% cfg = [];
% cfg.channel     = 'all';
% cfg.avgovertime = 'no';
% cfg.parameter   = 'avg';
% cfg.method      = 'analytic';
% cfg.statistic   = 'ft_statfun_depsamplesT';
% cfg.alpha       = 0.05;
% cfg.correctm    = 'no';
% 
% Nsub = numel(subj_young);
% cfg.design(1,:)  = [ones(1,Nsub) 2*ones(1,Nsub)];
% cfg.design(2,:)  = [1:Nsub 1:Nsub];
% cfg.ivar         = 1; % the 1st row in cfg.design contains the independent variable
% cfg.uvar         = 2; % the 2nd row in cfg.design contains the subject number
% cfg.latency     = [-0.3 0.0];
% stat = ft_timelockstatistics(cfg, young_repeat2_response_ERP{:}, young_ED_response_ERP{:});
% 
% 
% for i = 1:length(stat.mask)
% % make the plot
% cfg = [];
% cfg.style     = 'blank';
% cfg.layout    = 'acticap-64ch-standard2-EMCO.mat';
% cfg.highlight = 'on';
% cfg.highlightchannel = find(stat.mask(:,i));
% cfg.comment   = 'no';
% figure; ft_topoplotER(cfg, ERP.grandavg_repeat2_young)
% title(sprintf('Time: %d ms', round(stat.time(i)*1000)))
% end
% 
% %%
% cfg = [];
% cfg.channel     = 'all';
% cfg.avgovertime = 'no';
% cfg.parameter   = 'avg';
% cfg.method      = 'analytic';
% cfg.statistic   = 'ft_statfun_depsamplesT';
% cfg.alpha       = 0.05;
% cfg.correctm    = 'no';
% 
% Nsub = numel(subj_old);
% cfg.design(1,:)  = [ones(1,Nsub) 2*ones(1,Nsub)];
% cfg.design(2,:)  = [1:Nsub 1:Nsub];
% cfg.ivar         = 1; % the 1st row in cfg.design contains the independent variable
% cfg.uvar         = 2; % the 2nd row in cfg.design contains the subject number
% cfg.latency     = [-0.3 0.0];
% stat = ft_timelockstatistics(cfg, old_repeat2_response_ERP{:}, old_ED_response_ERP{:});
% 
% 
% for i = 1:length(stat.mask)
% % make the plot
% cfg = [];
% cfg.style     = 'blank';
% cfg.layout    = 'acticap-64ch-standard2-EMCO.mat';
% cfg.highlight = 'on';
% cfg.highlightchannel = find(stat.mask(:,i));
% cfg.comment   = 'no';
% figure; ft_topoplotER(cfg, ERP.grandavg_repeat2_old)
% title(sprintf('Time: %d ms', round(stat.time(i)*1000)))
% end
% 
