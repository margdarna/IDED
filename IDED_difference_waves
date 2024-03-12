% EEG Analysis script for Fieldtrip for single subjects and grand averages
%
% This script open preprocessed ERPs and calculates difference waves
%
% Written by: Margarita Darna
% Created on: 01. March 2024
% Last modified on: 07. March 2024

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

% Get subject information
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
%               1 - Calculate difference waves
%----------------------------------------------------------------------
% load stimulus locked epoches
load(strcat(dirs.output_dir, '3_Grandavg\grandavg_ERP_stimpres.mat'));

% Calculate the difference wave
cfg           = [];
cfg.operation = 'x2-x1';
cfg.parameter = 'individual';

ID_repeat_difference_young = ft_math(cfg, ERP.grandavg_repeat2_young, ERP.grandavg_ID_young);
ED_repeat_difference_young = ft_math(cfg, ERP.grandavg_repeat2_young, ERP.grandavg_ED_young);

ID_repeat_difference_old = ft_math(cfg, ERP.grandavg_repeat2_old, ERP.grandavg_ID_old);
ED_repeat_difference_old = ft_math(cfg, ERP.grandavg_repeat2_old, ERP.grandavg_ED_old);

% Calculate difference wave for all
ID_repeat_difference_all = ft_math(cfg, ERP.grandavg_repeat2, ERP.grandavg_ID);
ED_repeat_difference_all = ft_math(cfg, ERP.grandavg_repeat2, ERP.grandavg_ED);

% calculate average
cfg.operation = '(x1 + x2)/2';
shift_repeat_difference_all = ft_math(cfg, ID_repeat_difference_all, ED_repeat_difference_all);

% plot ERP
cfg = [];
cfg.layout = 'acticap-64ch-standard2-EMCO';
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
cfg.xlim       = [-0.25 1];
cfg.ylim       = [-2.5 3];

figure(); ft_multiplotER(cfg, shift_repeat_difference_all)
% there are two components present. First a negative one followed by a
% positive one.
figure(); ft_multiplotER(cfg, ID_repeat_difference_young, ED_repeat_difference_young)
figure();ft_multiplotER(cfg, ID_repeat_difference_old, ED_repeat_difference_old)

%% --------------------------------------------------------------------
%          2- Collect average amplitudes in the chosen latencies
%----------------------------------------------------------------------
% Get latency of first component (negative)
% prepare the cfg
cfg = [];
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'no';
cfg.channel     = {'Pz', 'P1', 'P2', 'POz', 'PO3', 'PO4'};%{'Pz'};
cfg.latency     =  [0.3 0.6];

young_ID_repeat_dif_component1 = ft_selectdata(cfg, ID_repeat_difference_young);
young_ED_repeat_dif_component1 = ft_selectdata(cfg, ED_repeat_difference_young);
old_ID_repeat_dif_component1 = ft_selectdata(cfg, ID_repeat_difference_old);
old_ED_repeat_dif_component1 = ft_selectdata(cfg, ED_repeat_difference_old);

young_ID_repeat_dif_component1_latency = get_latency_maxpeak(-squeeze(young_ID_repeat_dif_component1.individual)', 0.3, 500);
young_ED_repeat_dif_component1_latency = get_latency_maxpeak(-squeeze(young_ED_repeat_dif_component1.individual)', 0.3, 500);

old_ID_repeat_dif_component1_latency = get_latency_maxpeak(-squeeze(old_ID_repeat_dif_component1.individual)', 0.3, 500);
old_ED_repeat_dif_component1_latency = get_latency_maxpeak(-squeeze(old_ED_repeat_dif_component1.individual)', 0.3, 500);

component1_latency_young = mean([young_ID_repeat_dif_component1_latency young_ED_repeat_dif_component1_latency]);
component1_latency_old = mean([old_ID_repeat_dif_component1_latency old_ED_repeat_dif_component1_latency]);

% Get latency of second component (positive)
% prepare the cfg
cfg.latency     =  [0.5 0.9];

young_ID_repeat_dif_component2 = ft_selectdata(cfg, ID_repeat_difference_young);
young_ED_repeat_dif_component2 = ft_selectdata(cfg, ED_repeat_difference_young);
old_ID_repeat_dif_component2 = ft_selectdata(cfg, ID_repeat_difference_old);
old_ED_repeat_dif_component2 = ft_selectdata(cfg, ED_repeat_difference_old);

young_ID_repeat_dif_component2_latency = get_latency_maxpeak(squeeze(young_ID_repeat_dif_component2.individual)', 0.5, 500);
young_ED_repeat_dif_component2_latency = get_latency_maxpeak(squeeze(young_ED_repeat_dif_component2.individual)', 0.5, 500);

old_ID_repeat_dif_component2_latency = get_latency_maxpeak(squeeze(old_ID_repeat_dif_component2.individual)', 0.5, 500);
old_ED_repeat_dif_component2_latency = get_latency_maxpeak(squeeze(old_ED_repeat_dif_component2.individual)', 0.5, 500);

component2_latency_young = mean([young_ID_repeat_dif_component2_latency young_ED_repeat_dif_component2_latency]);
component2_latency_old   = mean([old_ID_repeat_dif_component2_latency old_ED_repeat_dif_component2_latency]);

% extract mean amplitude for each component of the difference waves:
cfg = [];
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';
cfg.channel     = {'Pz', 'P1', 'P2', 'POz', 'PO3', 'PO4'};%{'Pz'};
cfg.keepchandim = 'no';
cfg.keeptimedim = 'no';

% component 1
% young
cfg.latency     =  [component1_latency_young-0.1 component1_latency_young+0.1];
young_ID_repeat_dif_component1_ampl = ft_selectdata(cfg, ID_repeat_difference_young);
young_ED_repeat_dif_component1_ampl = ft_selectdata(cfg, ED_repeat_difference_young);
% old
cfg.latency     =  [component1_latency_old-0.1 component1_latency_old+0.1];
old_ID_repeat_dif_component1_ampl = ft_selectdata(cfg, ID_repeat_difference_old);
old_ED_repeat_dif_component1_ampl = ft_selectdata(cfg, ED_repeat_difference_old);

% component 2
% young
cfg.latency     =  [component2_latency_young-0.1 component2_latency_young+0.1];
young_ID_repeat_dif_component2_ampl = ft_selectdata(cfg, ID_repeat_difference_young);
young_ED_repeat_dif_component2_ampl = ft_selectdata(cfg, ED_repeat_difference_young);
% old
cfg.latency     =  [component2_latency_old-0.1 component2_latency_old+0.1];
old_ID_repeat_dif_component2_ampl = ft_selectdata(cfg, ID_repeat_difference_old);
old_ED_repeat_dif_component2_ampl = ft_selectdata(cfg, ED_repeat_difference_old);

%% --------------------------------------------------------------------
%                  3 - Generate matrix for ANOVA
%----------------------------------------------------------------------

% component1
% depending variable
dv = [young_ID_repeat_dif_component1_ampl.individual; old_ID_repeat_dif_component1_ampl.individual;...
    young_ED_repeat_dif_component1_ampl.individual; old_ED_repeat_dif_component1_ampl.individual];
% subject number
subj_num = [1:numel(subs) 1:numel(subs)]';
% generate table
X = [dv subj_num];
T = array2table(X);
% Assign headings to the table
T.Properties.VariableNames = {'dv', 'subj_num'};
% add between and within variables
T.between = [repmat('young', sum(is_young), 1); repmat(' old ', sum(~is_young), 1); repmat('young', sum(is_young), 1); ...
    repmat(' old ', sum(~is_young), 1)];
T.within = [repmat('ID-repeat', sum(is_young), 1); repmat('ID-repeat', sum(~is_young), 1);
    repmat('ED-repeat', sum(is_young), 1); repmat('ED-repeat', sum(~is_young), 1)];
% save table as csv
writetable(T, strcat(dirs.output_dir, '4_Stats\stat_ERP_component1_shift-repeat.csv'));

% component2
% depending variable
dv = [young_ID_repeat_dif_component2_ampl.individual; old_ID_repeat_dif_component2_ampl.individual;...
    young_ED_repeat_dif_component2_ampl.individual; old_ED_repeat_dif_component2_ampl.individual];
% subject number
subj_num = [1:numel(subs) 1:numel(subs)]';
% generate table
X = [dv subj_num];
T = array2table(X);
% Assign headings to the table
T.Properties.VariableNames = {'dv', 'subj_num'};
% add between and within variables
T.between = [repmat('young', sum(is_young), 1); repmat(' old ', sum(~is_young), 1); repmat('young', sum(is_young), 1); ...
    repmat(' old ', sum(~is_young), 1)];
T.within = [repmat('ID-repeat', sum(is_young), 1); repmat('ID-repeat', sum(~is_young), 1);
    repmat('ED-repeat', sum(is_young), 1); repmat('ED-repeat', sum(~is_young), 1)];
% save table as csv
writetable(T, strcat(dirs.output_dir, '4_Stats\stat_ERP_component2_shift-repeat.csv'));

%% --------------------------------------------------------------------
%                       4 - Plot results
%----------------------------------------------------------------------
% Calculate average of posterior channels
% prepare the cfg
cfg = [];
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'no';
cfg.channel     = {'Pz', 'P1', 'P2', 'POz', 'PO3', 'PO4'};%{'Pz'};

young_ID_repeat_dif_post_avg = ft_selectdata(cfg, ID_repeat_difference_young);
young_ED_repeat_dif_post_avg = ft_selectdata(cfg, ED_repeat_difference_young);
old_ID_repeat_dif_post_avg   = ft_selectdata(cfg, ID_repeat_difference_old);
old_ED_repeat_dif_post_avg   = ft_selectdata(cfg, ED_repeat_difference_old);

% young
% generating the group variable
within_val={'ID-repeat' 'ED-repeat'};
ind = [repmat(1, sum(is_young), 1); repmat(2, sum(is_young), 1)];
within = within_val(ind);

y = squeeze([young_ID_repeat_dif_post_avg.individual; young_ED_repeat_dif_post_avg.individual]);

figure()
g = gramm('x', young_ID_repeat_dif_post_avg.time*1000, 'y', y, 'linestyle', within, 'lightness', within);
g.geom_polygon('x',{[component1_latency_young*1000-100 component1_latency_young*1000+100 ...
                    component1_latency_young*1000+100 component1_latency_young*1000-100]},...
    'y', {[-9 -9 9 9]}, 'alpha', 0.1)
g.geom_polygon('x',{[component2_latency_young*1000-100 component2_latency_young*1000+100 ...
                    component2_latency_young*1000+100 component2_latency_young*1000-100]},...
    'y', {[-9 -9 9 9]}, 'alpha', 0.1)
g.stat_summary('type', 'ci'); % plot with 95% confidence interval
g.set_title(sprintf('Difference wave young\n(Pz, P1, P2, POz, PO3, PO4)'));
g.set_names('color','Condition','x','Time (ms)','y','Potential (uV)', 'linestyle', 'Condition', 'lightness','Condition');
g.axe_property('XLim',[-250 1000], 'YLim',[-4 5]);
g.set_color_options('lightness_range',[80 30], 'legend','merge');
g.set_order_options('linestyle',0, 'lightness', 0);
g.draw();

% old
% generating the group variable
within_val={'ID-repeat' 'ED-repeat'};
ind = [repmat(1, sum(~is_young), 1); repmat(2, sum(~is_young), 1)];
within = within_val(ind);

y = squeeze([old_ID_repeat_dif_post_avg.individual; old_ED_repeat_dif_post_avg.individual]);

figure()
g = gramm('x', old_ID_repeat_dif_post_avg.time*1000, 'y', y, 'linestyle', within, 'lightness', within);
g.geom_polygon('x',{[component1_latency_old*1000-100 component1_latency_old*1000+100 ...
                    component1_latency_old*1000+100 component1_latency_old*1000-100]},...
    'y', {[-9 -9 9 9]}, 'alpha', 0.1)
g.geom_polygon('x',{[component2_latency_old*1000-100 component2_latency_old*1000+100 ...
                    component2_latency_old*1000+100 component2_latency_old*1000-100]},...
    'y', {[-9 -9 9 9]}, 'alpha', 0.1)
g.stat_summary('type', 'ci'); % plot with 95% confidence interval
g.set_title(sprintf('Difference wave old\n(Pz, P1, P2, POz, PO3, PO4)'));
g.set_names('color','Condition','x','Time (ms)','y','Potential (uV)', 'linestyle', 'Condition', 'lightness','Condition');
g.axe_property('XLim',[-250 1000], 'YLim',[-4 5]);
g.set_color_options('lightness_range',[80 30],'hue_range',[175 175], 'legend','merge');
g.set_order_options('linestyle',0, 'lightness', 0);
g.draw();

%% --------------------------------------------------------------------
%                   5 - Create topography movies
%----------------------------------------------------------------------

% a) Create png for every frame and save it
cfg = [];
cfg.layout       = 'acticap-64ch-standard2-EMCO.mat';
cfg.ylim         = [4 8];
cfg.zlim         = [-3 3];
cfg.colorbar     = 'yes';
cfg.comment      = 'no';
cfg.interactive  = 'no';

k = 1;
for time = -0.1:0.002:0.9
    cfg.xlim = [time time];

    % overall grandaverages
    % ft_topoplotTFR(cfg,TFR.grandavg_young);
    % saveas(gcf, sprintf('%s8_movies/Theta_topo_young_%i.png', dirs.output_dir, k));
    % ft_topoplotTFR(cfg,TFR.grandavg_old);
    % saveas(gcf, sprintf('%s8_movies/Theta_topo_old_%i.png', dirs.output_dir, k));

    % young difference waves separated by condition
    h = ft_topoplotTFR(cfg,ID_repeat_difference_young);
    saveas(gcf, sprintf('%s8_movies/ERP_topo_young_ID_repeat_dif_%i.png', dirs.output_dir, k));
    ft_topoplotTFR(cfg,ED_repeat_difference_young);
    saveas(gcf, sprintf('%s8_movies/ERP_topo_young_ED_repeat_dif_%i.png', dirs.output_dir, k));

    %old separated by condition
    ft_topoplotTFR(cfg,ID_repeat_difference_old);
    saveas(gcf, sprintf('%s8_movies/ERP_topo_old_ID_repeat_dif_%i.png', dirs.output_dir, k));
    ft_topoplotTFR(cfg,ED_repeat_difference_old);
    saveas(gcf, sprintf('%s8_movies/ERP_topo_old_ED_repeat_dif_%i.png', dirs.output_dir, k));

    k = k + 1;
    % close all figures
    close all
end

% b) prepare frames
frames_ID_repeat_dif_young = cell(601,1);
frames_ED_repeat_dif_young = cell(601,1);
frames_ID_repeat_dif_old   = cell(601,1);
frames_ED_repeat_dif_old   = cell(601,1);

for i = 1:601
    frames_ID_repeat_dif_young{i} = imread(sprintf('%s8_movies/ERP_topo_young_ID_repeat_dif_%i.png', dirs.output_dir, i));
    frames_ED_repeat_dif_young{i} = imread(sprintf('%s8_movies/ERP_topo_young_ED_repeat_dif_%i.png', dirs.output_dir, i));
    frames_ID_repeat_dif_old{i}   = imread(sprintf('%s8_movies/ERP_topo_old_ID_repeat_dif_%i.png', dirs.output_dir, i));
    frames_ED_repeat_dif_old{i}   = imread(sprintf('%s8_movies/ERP_topo_old_ED_repeat_dif_%i.png', dirs.output_dir, i));
end

% add text string with time info
time = -250:2:950;
position = [350 600];
for i = 1:601
    text_str = sprintf('%i ms', time(i));
    frames_ID_repeat_dif_young{i} = insertText(frames_ID_repeat_dif_young{i}, position, text_str, ...
        'FontSize', 24, 'BoxColor', 'white');

    frames_ED_repeat_dif_young{i} = insertText(frames_ED_repeat_dif_young{i}, position, text_str, ...
        'FontSize', 24, 'BoxColor', 'white');

    frames_ID_repeat_dif_old{i} = insertText(frames_ID_repeat_dif_old{i}, position, text_str, ...
        'FontSize', 24, 'BoxColor', 'white');

    frames_ED_repeat_dif_old{i} = insertText(frames_ED_repeat_dif_old{i}, position, text_str, ...
        'FontSize', 24, 'BoxColor', 'white');
end

% c) create the video writer
% young - ID-repeat
writerObj = VideoWriter(sprintf('%s8_movies/ERP_topo_young_ID_repeat_dif.avi', dirs.output_dir));
writerObj.FrameRate = 30;

% open the video writer
open(writerObj);
% write the frames to the video
for u = 1:length(frames_ID_repeat_dif_young)
    % convert the image to a frame
    frame = im2frame(frames_ID_repeat_dif_young{u});
    writeVideo(writerObj, frame);
end
%close the writer object
close(writerObj);


% young - ED-repeat
% create the video writer
writerObj = VideoWriter(sprintf('%s8_movies/ERP_topo_young_ED_repeat_dif.avi', dirs.output_dir));
writerObj.FrameRate = 30;

% open the video writer
open(writerObj);
% write the frames to the video
for u = 1:length(frames_ED_repeat_dif_young)
    % convert the image to a frame
    frame = im2frame(frames_ED_repeat_dif_young{u});
    writeVideo(writerObj, frame);
end
%close the writer object
close(writerObj);


% old - ID-repeat
% create the video writer
writerObj = VideoWriter(sprintf('%s8_movies/ERP_topo_old_ID_repeat_dif.avi', dirs.output_dir));
writerObj.FrameRate = 30;

% open the video writer
open(writerObj);
% write the frames to the video
for u = 1:length(frames_ID_repeat_dif_old)
    % convert the image to a frame
    frame = im2frame(frames_ID_repeat_dif_old{u});
    writeVideo(writerObj, frame);
end
%close the writer object
close(writerObj);


% old - ED-repeat
% create the video writer
writerObj = VideoWriter(sprintf('%s8_movies/ERP_topo_old_ED_repeat_dif.avi', dirs.output_dir));
writerObj.FrameRate = 30;

% open the video writer
open(writerObj);
% write the frames to the video
for u = 1:length(frames_ED_repeat_dif_old)
    % convert the image to a frame
    frame = im2frame(frames_ED_repeat_dif_old{u});
    writeVideo(writerObj, frame);
end
%close the writer object
close(writerObj);
