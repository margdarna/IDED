% EEG Analysis script for Fieldtrip for single subjects and grand averages
%
% This script opens preprocessed eeg data from brain vision analyzer,
% and gives out basic analytics
% Based on tutorial from fieldtriptoolbox.org
%
% Written by: Margarita Darna
% Created on: 03. August 2022
% Last modified on: 12. March 2024

%% --------------------------------------------------------------------
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
dirs.analysis_dir   = strcat(dirs.exp_dir, 'Functions/');
dirs.output_dir     = strcat(dirs.exp_dir, 'Output/');
dirs.prepr_dir      = strcat(dirs.output_dir, '1_Preprocessing/');

% adding analysis path and subfolders
addpath(genpath(dirs.analysis_dir));

% making sure that fieldtrip is called correctly
ft_defaults

% define epochs - big epochs for the purpose of time frequency analysis
prestim = 1.25; % prestimulus interval in seconds
poststim = 3;   % poststimulus interval in seconds

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
    % ensure that we don't mix up subjects - delete previously created
    % variable "subjectdata"
    clear subjectdata

    % define subject
    subj = subs{i};

    % print information in command window
    fprintf("************************\nStarting %s\n************************\n", subj)

    % save all subject relevant information in subjectdata
    subjectdata.subj = subj;    

    % Check if preprocessing has already been performed by checking if the matfile exists
    prepr_done = exist(strcat(dirs.prepr_dir, subj, '_IDED_ERP_ft_prepr_1.mat'), 'file');
    %prepr_done = 1;
    if prepr_done == 2
        fprintf('%s_IDED_ERP_ft_prepr_1.mat File already exists\n\n', subj)
    else
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

        % extract stimulus-locked epochs
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

        % create datasets by redifining the trials
        data_stimpres   = ft_redefinetrial(cfg_stimpres, data_eeg);

        % remove baseline from trials
        cfg = [];
        cfg.demean          = 'yes';
        cfg.baselinewindow  = repmat([-0.25 0], 434, 1);
        data_stimpres       = ft_preprocessing(cfg, data_stimpres);

        % lowpass filter for ERP analysis
        cfg = [];
        cfg.lpfilter  = 'yes';
        cfg.lpfreq    = 30;
        data_ERP = ft_preprocessing(cfg, data_stimpres);

        % save subject specific information
        subjectdata.data_eeg = data_eeg;
        subjectdata.data_ERP = data_ERP;
        
        % save subjectdata in mat file
        fprintf('Saving %s ... ', subj)
        save(strcat(dirs.prepr_dir, subj, '_IDED_ERP_ft_prepr_1.mat'), "subjectdata", '-v7.3');
        fprintf('succesful!\n\n')
    end
end
clear subjectdata data_eeg data_stimpres data_ERP cfg_stimpres

%% --------------------------------------------------------------------
%               1b - Check for saccades
%----------------------------------------------------------------------
1 - Load single subject EEG and prepare epochs
% 
% for i = 1 : height(subs)
%     % ensure that we don't mix up subjects by deleting subjectdata
%     clear subjectdata
% 
%     % define subject
%     subj = subs{i};
% 
%     % load subjectdata
%     load(strcat(dirs.prepr_dir, subj, '_IDED_ERP_ft_prepr_1.mat'), "subjectdata");
% 
%     % % detect movemebt with detece movement function
%     cfg = [];
%     cfg.method = 'velocity2D';
%     cfg.velocity2D.mindur  = 10;% a minimum saccade duration is around 20 ms
%     cfg.velocity2D.velther =100;
%     cfg.channel = {'HEOG'};
%     cfg.trials  = 'all';
%     [~, movement] = ft_detect_movement(cfg, subjectdata.data_ERP);
% 
%     % we need to figure out a way to exclude blinks from being detected
%     % as saccades
% 
%     % tranform to the correct format of events
%     saccades = [];
%     for k = 1:height(movement)
%         saccades(k).type     = 'Saccade';
%         saccades(k).value    = 'Saccade';
%         saccades(k).sample   = movement(k,1);
%         saccades(k).duration = movement(k,2) - movement(k,1);
%         saccades(k).offset   = [];
%     end
% 
%     % include saccades in event file
%     event = ft_read_event([dirs.derived_dt_dir, subj, '_task-IDED_eeg_prepr.dat']);
%     all_event = [event saccades];
% 
% %     cfg              = [];
% %     cfg.event        = all_event;
% %     cfg.channel      = {'VEOG', 'HEOG'};
% %     cfg.continuous   = 'no';
% %     cfg.allowoverlap = 'yes';
% %     ft_databrowser(cfg, subjectdata.data_ERP)
% 
%     % add saccade information in variable
%     data_ERP.saccadeinfo.saccades = saccades;
%     data_ERP.saccadeinfon_saccades = length(saccades);
% end

%% --------------------------------------------------------------------
%          2a - Calculate subject average for each condition
%----------------------------------------------------------------------
for i = 1:height(subs)

    % ensure that we don't mix up subjects by deleting subjectdata
    clear subjectdata

    % define subject
    subj = subs{i};

    % print information in command window
    fprintf("************************\nStarting %s\n************************\n", subj)

    % check if this step has already been performed
    subj_avg_done = exist(strcat(dirs.output_dir, '2_Subj_Avg\', subj, '_subj_avg_stimpres_ERP.mat'), 'file');
    %subj_avg_done = 1;
    if subj_avg_done == 2
        fprintf('%s_IDED_ERP_ft_prepr_1.mat File already exists\n\n', subj)
    else
        % load mat file that was prepared in the previous section
        load(strcat(dirs.prepr_dir, subj, '_IDED_ERP_ft_prepr_1.mat'), "subjectdata");
        data_ERP = subjectdata.data_ERP;

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
        % plot ERP
        % cfg = [];
        % cfg.layout = 'acticap-64ch-standard2-EMCO';
        % cfg.interactive = 'yes';
        % cfg.showoutline = 'yes';
        % ft_multiplotER(cfg, repeat_ERP, shift_ERP)

        % save subject data
        save(strcat(dirs.output_dir, '2_Subj_Avg\', subj, '_subj_avg_stimpres_ERP.mat'), "repeat1_ERP", "repeat2_ERP", ...
            "shift_ERP", "last_ERP", "ID_ERP", "ED_ERP", "alltrl_ERP");

        % clear unnecessary variables
        clear repeat1_ERP repeat2_ERP shift_ERP last_ERP ID_ERP ED_ERP alltrl_ERP data_ERP

        % print information in the command window
        fprintf(strcat('Subject', subj, ' succesful!\n\n'));
    end
end

%% --------------------------------------------------------------------
%        2b - Generate struct array tp save all subject averages
%----------------------------------------------------------------------
% add all subject information in one struct array
for i = 1:height(subs)
    % define subject
    subj = subs{i};

    % load subject mat file
    load(strcat(dirs.output_dir, '2_Subj_Avg\', subj, '_subj_avg_stimpres_ERP.mat'))

    % save the struct arrays in cell array for all subjects
    % ERP
    all_repeat1_stimpres_ERP{i} = repeat1_ERP;
    all_repeat2_stimpres_ERP{i} = repeat2_ERP;
    all_shift_stimpres_ERP{i}  = shift_ERP;
    all_last_stimpres_ERP{i}   = last_ERP;
    all_ID_stimpres_ERP{i}     = ID_ERP;
    all_ED_stimpres_ERP{i}     = ED_ERP;
    all_alltrl_stimpres_ERP{i} = alltrl_ERP;

end

% save the struct array
save(strcat(dirs.output_dir, '2_Subj_Avg\all_subj_avg_stimpres_ERP'),...
    "all_repeat1_stimpres_ERP",  "all_repeat2_stimpres_ERP", "all_shift_stimpres_ERP", "all_last_stimpres_ERP", "all_ID_stimpres_ERP", ...
    "all_ED_stimpres_ERP", "all_alltrl_stimpres_ERP", "subj_info", "tblstats");

% separate arrays by age group
young_repeat1_stimpres_ERP = all_repeat1_stimpres_ERP(is_young);
old_repeat1_stimpres_ERP   = all_repeat1_stimpres_ERP(~is_young);
young_repeat2_stimpres_ERP = all_repeat2_stimpres_ERP(is_young);
old_repeat2_stimpres_ERP   = all_repeat2_stimpres_ERP(~is_young);
young_shift_stimpres_ERP   = all_shift_stimpres_ERP(is_young);
old_shift_stimpres_ERP     = all_shift_stimpres_ERP(~is_young);
young_ID_stimpres_ERP      = all_ID_stimpres_ERP(is_young);
old_ID_stimpres_ERP        = all_ID_stimpres_ERP(~is_young);
young_ED_stimpres_ERP      = all_ED_stimpres_ERP(is_young);
old_ED_stimpres_ERP        = all_ED_stimpres_ERP(~is_young);
young_alltrl_stimpres_ERP  = all_alltrl_stimpres_ERP(is_young);
old_alltrl_stimpres_ERP    = all_alltrl_stimpres_ERP(~is_young);
young_last_stimpres_ERP    = all_last_stimpres_ERP(is_young);
old_last_stimpres_ERP      = all_last_stimpres_ERP(~is_young);

% save the separated struct arrays
save(strcat(dirs.output_dir, '2_Subj_Avg\young_subj_avg_stimpres_ERP'),"young_repeat1_stimpres_ERP", ...
    "young_repeat2_stimpres_ERP", "young_shift_stimpres_ERP", "young_last_stimpres_ERP", "young_ID_stimpres_ERP", ...
    "young_ED_stimpres_ERP", "young_alltrl_stimpres_ERP", "subj_info", "tblstats");

save(strcat(dirs.output_dir, '2_Subj_Avg\old_subj_avg_stimpres_ERP'), "old_repeat1_stimpres_ERP", ...
    "old_repeat2_stimpres_ERP", "old_shift_stimpres_ERP", "old_last_stimpres_ERP", "old_ID_stimpres_ERP", ...
    "old_ED_stimpres_ERP", "old_alltrl_stimpres_ERP", "subj_info", "tblstats");

% clear unnecasarry variables
clear all_repeat1_stimpres_ERP all_repeat2_stimpres_ERP all_shift_stimpres_ERP all_ID_stimpres_ERP all_ED_stimpres_ERP all_alltrl_stimpres_ERP ...
    all_last_stimpres_ERP

%% --------------------------------------------------------------------
%              3 - Calculate grand average ERP for each condition
%----------------------------------------------------------------------
% load subject averages
load(strcat(dirs.output_dir, '2_Subj_Avg\young_subj_avg_stimpres_ERP'));
load(strcat(dirs.output_dir, '2_Subj_Avg\old_subj_avg_stimpres_ERP'));
load(strcat(dirs.output_dir, '2_Subj_Avg\all_subj_avg_stimpres_ERP'));

% calculate grand average for each condition for the purposes of figures
cfg = [];
cfg.channel        = 1:65;
cfg.latency        = 'all';
cfg.parameter      = 'avg';
cfg.keepindividual = 'yes'; % important to keep to be able to calculate the confidence interval

ERP.grandavg               = ft_timelockgrandaverage(cfg, all_alltrl_stimpres_ERP{:});
ERP.grandavg_young         = ft_timelockgrandaverage(cfg, young_alltrl_stimpres_ERP{:});
ERP.grandavg_old           = ft_timelockgrandaverage(cfg, old_alltrl_stimpres_ERP{:});
ERP.grandavg_repeat1       = ft_timelockgrandaverage(cfg, all_repeat1_stimpres_ERP{:});
ERP.grandavg_repeat2       = ft_timelockgrandaverage(cfg, all_repeat2_stimpres_ERP{:});
ERP.grandavg_shift         = ft_timelockgrandaverage(cfg, all_shift_stimpres_ERP{:});
ERP.grandavg_ID            = ft_timelockgrandaverage(cfg, all_ID_stimpres_ERP{:});
ERP.grandavg_ED            = ft_timelockgrandaverage(cfg, all_ED_stimpres_ERP{:});
ERP.grandavg_repeat1_young = ft_timelockgrandaverage(cfg, young_repeat1_stimpres_ERP{:});
ERP.grandavg_repeat1_old   = ft_timelockgrandaverage(cfg, old_repeat1_stimpres_ERP{:});
ERP.grandavg_repeat2_young = ft_timelockgrandaverage(cfg, young_repeat2_stimpres_ERP{:});
ERP.grandavg_repeat2_old   = ft_timelockgrandaverage(cfg, old_repeat2_stimpres_ERP{:});
ERP.grandavg_shift_young   = ft_timelockgrandaverage(cfg, young_shift_stimpres_ERP{:});
ERP.grandavg_shift_old     = ft_timelockgrandaverage(cfg, old_shift_stimpres_ERP{:});
ERP.grandavg_ID_young      = ft_timelockgrandaverage(cfg, young_ID_stimpres_ERP{:});
ERP.grandavg_ID_old        = ft_timelockgrandaverage(cfg, old_ID_stimpres_ERP{:});
ERP.grandavg_ED_young      = ft_timelockgrandaverage(cfg, young_ED_stimpres_ERP{:});
ERP.grandavg_ED_old        = ft_timelockgrandaverage(cfg, old_ED_stimpres_ERP{:});
ERP.grandavg_last_young    = ft_timelockgrandaverage(cfg, young_last_stimpres_ERP{:});
ERP.grandavg_last_old      = ft_timelockgrandaverage(cfg, old_last_stimpres_ERP{:});

% save in matfile
save(strcat(dirs.output_dir, '3_Grandavg\grandavg_ERP_stimpres.mat'), "ERP", '-v7.3');

%% --------------------------------------------------------------------
%          4- Collect average amplitudes for specific EEG components
%----------------------------------------------------------------------
% load grandaverages
load(strcat(dirs.output_dir, '3_Grandavg\grandavg_ERP_stimpres.mat'));

% Visualize average ERP to find the peaks of interest
cfg = [];
cfg.showlabels  = 'yes';
cfg.layout      = 'acticap-64ch-standard2-EMCO.mat';
cfg.xlim        = [-0.25 1.5];
figure; ft_multiplotER(cfg, ERP.grandavg)

figure; ft_multiplotER(cfg, ERP.grandavg_young, ERP.grandavg_old)

% rertieving latencies from Latency calculation script - run this script
% now if it has not yet been performed
% IDED_ERP_Analysis_Lateny.m
load("ind_P300.mat")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% P300 - parietal electrodes based on distribution of grand average
% prepare the cfg
cfg = [];
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';
cfg.channel     = {'Pz', 'P1', 'P2', 'POz', 'PO3', 'PO4'};
cfg.latency_young     =  [ind_P300p_young-0.1 ind_P300p_young+0.1];
cfg.latency_old       =  [ind_P300p_old-0.1 ind_P300p_old+0.1];

% get component average
P300p = get_component_average(cfg, 'P300p', dirs);

%Generate matrix for ANOVA in RStudio

% P300p
create_R_matrix('P300p', 'ERP', subj_info, dirs)

%% --------------------------------------------------------------------
%                       6 - Plot results
%----------------------------------------------------------------------

% load grandaverages
load(strcat(dirs.output_dir, '3_Grandavg\grandavg_ERP_stimpres.mat'));

% load P300 index
load("ind_P300.mat")

% select only relevant channels and average across them for each individual
% P300
cfg = [];
cfg.channel     = {'Pz', 'P1', 'P2', 'POz', 'PO3', 'PO4'};
cfg.latency     =  [-0.25 1];
cfg.avgoverchan = 'yes';

P300p.grandavg               = ft_selectdata(cfg, ERP.grandavg);
P300p.grandavg_repeat1_young = ft_selectdata(cfg, ERP.grandavg_repeat1_young);
P300p.grandavg_repeat2_young = ft_selectdata(cfg, ERP.grandavg_repeat2_young);
P300p.grandavg_ID_young      = ft_selectdata(cfg, ERP.grandavg_ID_young);
P300p.grandavg_ED_young      = ft_selectdata(cfg, ERP.grandavg_ED_young);
P300p.grandavg_repeat1_old    = ft_selectdata(cfg, ERP.grandavg_repeat1_old);
P300p.grandavg_repeat2_old    = ft_selectdata(cfg, ERP.grandavg_repeat2_old);
P300p.grandavg_ID_old        = ft_selectdata(cfg, ERP.grandavg_ID_old);
P300p.grandavg_ED_old        = ft_selectdata(cfg, ERP.grandavg_ED_old);
P300p.grandavg_last_old      = ft_selectdata(cfg, ERP.grandavg_last_old);
P300p.grandavg_last_young    = ft_selectdata(cfg, ERP.grandavg_last_young);


P300_young = squeeze([P300p.grandavg_repeat2_young.individual; P300p.grandavg_ID_young.individual; ...
    P300p.grandavg_ED_young.individual]);

P300_old = squeeze([P300p.grandavg_repeat2_old.individual; P300p.grandavg_ID_old.individual; ...
    P300p.grandavg_ED_old.individual]);

% young
% generating the group variable
within_val={'repeat' 'ID' 'ED'};
ind = [repmat(1, sum(is_young), 1); repmat(2, sum(is_young), 1); repmat(3, sum(is_young), 1)];
within = within_val(ind);

figure()
g = gramm('x', P300p.grandavg.time*1000, 'y', P300_young, 'linestyle', within, 'lightness', within);
g.geom_polygon('x',{[ind_P300p_young*1000-100 ind_P300p_young*1000+100 ind_P300p_young*1000+100 ind_P300p_young*1000-100]},...
    'y', {[-9 -9 9 9]}, 'alpha', 0.1)
g.stat_summary('type', 'ci'); % plot with 95% confidence interval
g.set_title(sprintf('Average ERP young\n(Pz, P1, P2, POz, PO3, PO4)'));
g.set_names('color','Condition','x','Time (ms)','y','Potential (uV)', 'linestyle', 'Condition', 'lightness','Condition');
g.axe_property('XLim',[-250 1000], 'YLim',[-9 9]);
g.set_color_options('lightness_range',[80 30], 'legend','merge');
g.set_order_options('linestyle',0, 'lightness', 0);
g.draw();

% old
ind = [repmat(1, sum(~is_young), 1); repmat(2, sum(~is_young), 1); repmat(3, sum(~is_young), 1)];
within = within_val(ind);

figure()
g = gramm('x', P300p.grandavg.time*1000, 'y', P300_old, 'linestyle', within, 'lightness', within);
g.geom_polygon('x',{[ind_P300p_old*1000-100 ind_P300p_old*1000+100 ind_P300p_old*1000+100 ind_P300p_old*1000-100]},...
    'y', {[-9 -9 9 9]}, 'alpha', 0.1)
g.stat_summary('type', 'ci'); % plot with 95% confidence interval
g.set_title(sprintf('Average ERP old\n(Pz, P1, P2, POz, PO3, PO4)'));
g.set_names('color','Condition','x','Time (ms)','y','Potential (uV)', 'linestyle', 'Condition', 'lightness','Condition');
g.axe_property('XLim',[-250 1000], 'YLim',[-9 9]);
g.set_color_options('lightness_range',[80 30],'hue_range',[175 175], 'legend','merge');
g.set_order_options('linestyle',0, 'lightness', 0);
g.draw();
