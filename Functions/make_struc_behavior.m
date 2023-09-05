function perf = make_struc_behavior (curr_dir, save_dir, subj, task)
% this funcions generates a struc variable out of all performance variables
% of one subject in one run of the chosen task
%
% INPUT:
%       curr_dir - string, directory where all_var.mat file is located
%       save_dir - string, directory of where perf.mat file will be saved
%       subj     - string, subject pseudonym (e.g. "S402")
%       taks     - strink, task name (e.g. "IDED") as it appears in the BIDS
%       conventional name of the data
%
% OUTPUT:
%       perf - struct, contains all information about subject performance in
%       one specific task
%
% Created on:    05/02/2022
% Last modified: 05/02/2022
%
% Created by: Margarita Darna
% margarita.darna@lin-magdeburg.de

% Loading all_var.mat
load(sprintf('%s%s_task-%s_all_var.mat', curr_dir, subj, task));

% Creating empty struct
perf = struct();

%% saving basic information:
perf.subj = subj;
perf.task = task;

%% saving trial information
if strcmp(task, "IDED")
    perf.all_trials = all_trials;
    perf.trial_decomposition = [practice_trials total_exp_trials]; % first practice then total_exp trials
    perf.trial_str = trial_str;
    % renaming first column
    perf.trial_str.Properties.VariableNames(1) = {'trial'};
    % adding variable is_break as last column
    try
    perf.trial_str.is_break = is_break;
    catch
        fprintf('%s: is_break not found in all_var.mat', subj);
    end
elseif strcmp(task, "ASST-D") || strcmp(task, "ASSTD")|| strcmp(task, "ASST") || strcmp(task, "ASST-G") || strcmp(task, "ASST-Gabor")|| strcmp(task, "ASSTG")
    perf.all_trials = trial - 1;
    perf.trial_decomposition = [8 perf.all_trials-8]; % predetermined: practice trials = 8
    %perf.trial_str = trial_str(1:perf.all_trials,:); 
    perf.stage_str = [stage_names stimulus_group target_dim target_num irr_dim distractor_num];
end

%% saving performance information
% creating column with trial numbers
perf.resp_mat = array2table(respMat(1:perf.all_trials,:));
 %%%% make sure no rows are skipped
perf.t_mat    = array2table(t_Mat(1: perf.all_trials,:));
if strcmp(task, "IDED")
    % resp_mat:
        % naming table columns
        perf.resp_mat.Properties.VariableNames = {'trial','targ_pos','response','RT'};
        % adding column that checks if response is correct
        perf.resp_mat.correct = perf.resp_mat.targ_pos == perf.resp_mat.response;
    % t_mat:
        % get the size
        t_size = size(perf.t_mat);
        % naming table columns
        if t_size(2) == 4
            perf.t_mat.Properties.VariableNames = {'cross', 'stim', 'press', 'black'};
        elseif t_size(2) == 3
            perf.t_mat.Properties.VariableNames = {'cross', 'stim', 'press'};
        end
elseif strcmp(task, "ASST-D") || strcmp(task, "ASSTD")|| strcmp(task, "ASST") || strcmp(task, "ASST-G") || strcmp(task, "ASST-Gabor")|| strcmp(task, "ASSTG")
    % resp_mat:       
        % naming table columns
        perf.resp_mat.Properties.VariableNames = {'trial', 'stage_num', 'targ_pos','response','RT', 'cons_corr'};
        % adding column that checks if response is correct
        perf.resp_mat.correct = perf.resp_mat.targ_pos == perf.resp_mat.response;
    % t_mat:
        % naming table columns
        perf.t_mat.Properties.VariableNames = {'trial', 'cross', 'stim', 'press', 'black', 'feedback'};
end

%% check if subject folder exists in save_dir
if ~isfolder(sprintf('%s%s', save_dir, subj))
   mkdir(sprintf('%s%s', save_dir, subj)) 
end
% saving struct in perf.mat file
save(sprintf('%s%s/%s_task-%s_perf.mat', save_dir, subj, subj, task), 'perf');
end