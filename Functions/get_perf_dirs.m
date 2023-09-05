function perf_dirs = get_perf_dirs(raw_dt_dir, derived_dt_dir, subs, task)
% This function goes through the raw data folder and chooses all subjects
% that are defined in subs. It locates the behavioral data for one specific
% task ("IDED" or "ASST" or "ASST_Gabor") and creates a string array called
% perf_dirs that includes the directories of the perf structures of each of the chosen subject.
% Every subject appears in one row
%
% INPUT:
%       raw_dt_dir  - string, directory of raw data
%       derived_dt_dir - string, directory of where derived data are saved
%       subs        - string array, pseodonyms of subs
%       task        - string, taks for which we want perf_all to be created
%       for
%
% OUTPUT:
%       perf_dirs	- string array, n x 1 cell array, where n: number of
%       subjects, in each cell the directory of perf.mat files for each
%       subject can be found
%
% Created on: 05/02/2022
% Last modified: 18/02/2022
%
% Created by: Margarita Darna
% margarita.darna@lin-magdeburg.de
perf_dirs = strings(numel(subs), 1);

for i = 1:numel(subs)
    subj = subs{i};
    % check if subject folder exists
    subj_dir = strcat(raw_dt_dir, subj);
    if isempty(dir(subj_dir))
        fprintf('%s: Subject folder not found\n', subj)
    else
        % check if Behavioral folder exists
        beh_dir = strcat(subj_dir, '/Behavior/');
        if ~exist(beh_dir, 'dir')
            fprintf('%s: Behavioral folder not found\n', subj)
        else
            % check if all_var matrix of chosen task is saved in behavioral folder
            if isempty(dir(sprintf('%s%s_task-%s_all_var.mat', beh_dir, subj, task)))
               fprintf('%s: Performance file for task %s not found\n', subj, task)               
            else
                % create perf.mat if does not already exist
                if isempty(dir(sprintf('%s%s/%s_task-%s_perf.mat', derived_dt_dir, subj, subj, task)))
                    try
                    % generate perf.mat using custom function
                        make_struc_behavior(beh_dir, derived_dt_dir, subj, task);
                    catch
                        fprintf('%s: Could not generate perf struct from all_var.mat\n', subj, task)
                        break
                    end
                end
                % try loading performance matrix
                try
                    load(sprintf('%s%s/%s_task-%s_perf.mat', derived_dt_dir, subj, subj, task))
                    if exist('perf', 'var')
                        % check if the perf matrix actually belongs to that
                        % participant
                        is_correct = subj == perf.subj;
                        if is_correct
                        % saving perf struc in perf_all
                        perf_dirs(i, 1) = sprintf('%s%s/%s_task-%s_perf.mat', derived_dt_dir, subj, subj, task);
                        fprintf('%s: Succesfull!!\n', subj)
                        clear('perf')
                        else
                            fprintf('%s: Performance matrix of %s is not assigned correctly. It was assigned to subject %s\n', ...
                                subj, task, perf.subj)
                        end
                    else
                        perf_dirs(i, 1) = '';
                        fprintf('%s: Could not find perf struct in perf_mat!\n', subj)
                    end
                catch
                    fprintf('%s: Could not load perf_mat for task %s\n', subj, task)
                end
            end
        end
    end
end
end