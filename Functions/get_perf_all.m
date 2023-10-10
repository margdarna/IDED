function perf_all = get_perf_all(dirs, task, perf_all)
% This function goes through all directories aff all chosen taks and
% calculates reaction times and error rates for each subject. These are
% then saved in the individual perf.mat files of each subject, but also all
% together in the perf_all file
%
% INPUT
%   dirs - struct, includes the names of all directories of this project
%   task - string array, includes the names of the tasks that should be
%   analyzed
%   perf_all - struct, intialized perf_all that already includes the
%   directories of the perf.mat files that should be loaded
%
%
% OUTPUT
%   perf_all - includes all performance parameters of all subjects
%

%% preparing empty arrays to store values that will be retrieved later
for n = 1: numel(task)
    if any(strcmp(task, "ASST-D")) || any(strcmp(task, "ASSTD"))|| any(strcmp(task, "ASST"))
        % ASST
        perf_all.ASSTD.RT.median_all             = nan(numel(perf_all.subs), 1);
        perf_all.ASSTD.RT.median_repeat          = nan(numel(perf_all.subs), 1);
        perf_all.ASSTD.RT.median_shift           = nan(numel(perf_all.subs), 1);
        perf_all.ASSTD.RT.median_first           = nan(numel(perf_all.subs), 1);
        perf_all.ASSTD.RT.median_second          = nan(numel(perf_all.subs), 1);
        perf_all.ASSTD.RT.median_last            = nan(numel(perf_all.subs), 1);
        perf_all.ASSTD.error.error_all           = nan(numel(perf_all.subs), 1);
        perf_all.ASSTD.error.error_in_stage      = nan(numel(perf_all.subs), 9);
        perf_all.ASSTD.error.trials_to_criterion = nan(numel(perf_all.subs), 9);
    elseif any(task(n) == "IDED")
        % IDED
        perf_all.IDED.RT.median_all    = nan(numel(perf_all.subs), 1);
        perf_all.IDED.RT.median_ID     = nan(numel(perf_all.subs), 1);
        perf_all.IDED.RT.median_ED     = nan(numel(perf_all.subs), 1);
        perf_all.IDED.RT.median_shift  = nan(numel(perf_all.subs), 1);
        perf_all.IDED.RT.median_last   = nan(numel(perf_all.subs), 1);
        perf_all.IDED.RT.median_repeat1 = nan(numel(perf_all.subs), 1);
        perf_all.IDED.RT.median_repeat2 = nan(numel(perf_all.subs), 1);
        perf_all.IDED.error.error_all  = nan(numel(perf_all.subs), 1);
        perf_all.IDED.error.error_ID         = nan(numel(perf_all.subs), 1);
        perf_all.IDED.error.error_ED         = nan(numel(perf_all.subs), 1);
        perf_all.IDED.error.error_shift      = nan(numel(perf_all.subs), 1);
        perf_all.IDED.error.error_last     = nan(numel(perf_all.subs), 1);
        perf_all.IDED.error.error_repeat1     = nan(numel(perf_all.subs), 1);
        perf_all.IDED.error.error_repeat2     = nan(numel(perf_all.subs), 1);
    end
end

%% Retrieving performance information
for n = 1:numel(task)
    fprintf('\nAnalyzing performance for %s:\n-----------------------------------\n', task(n));
    for i = 1: numel(perf_all.subs)
        if task(n) == "IDED"
            try
                % load perf struct for specific task
                load(perf_all.IDED.perf_dirs(i), 'perf')
            catch
                fprintf('%s: Perf.mat does not exist\n', perf_all.subs{i});
                continue
            end
            try
                % create indices, if not already there
                perf.trial_str.is_practice = perf.trial_str.trial_type == 'practice';
                perf.trial_str.is_ID       = perf.trial_str.trial_type == 'ID';
                perf.trial_str.is_ED       = perf.trial_str.trial_type == 'ED';
                % for repeat we only take the first trial after shift, so
                % that we have the same number of trials between shift and
                % repeat conditions
                perf.trial_str.is_repeat   =  perf.trial_str.is_ID | perf.trial_str.is_ED;
                perf.trial_str.is_repeat1   = [false; perf.trial_str.is_repeat(1:end - 1) ];
                % here we take the second trial after shift
                perf.trial_str.is_repeat2   = [false; false; perf.trial_str.is_repeat(1:end - 2) ];

                % all RTs (correct trials)
                x = perf.resp_mat.RT( ~perf.trial_str.is_practice  & perf.resp_mat.correct,:);
                % get skewness
                perf.RT.skewness_all = skewness(x);
                perf.RT.median_all = median(x, 'omitnan');
                %perf.RT.median_all = plot_his_RTs(task(n), perf_all.subs{i}, x, 'all RTs');
                
                % separated by shift type & exclude wrong trials
                x1 = perf.resp_mat.RT(perf.trial_str.is_ID & perf.resp_mat.correct,:);
                x2 = perf.resp_mat.RT(perf.trial_str.is_ED & perf.resp_mat.correct,:);
                x3 = perf.resp_mat.RT(perf.trial_str.is_repeat1 & perf.resp_mat.correct,:);
                x4 = perf.resp_mat.RT(perf.trial_str.is_repeat2 & perf.resp_mat.correct,:);

                perf.RT.median_ID = median(x1, 'omitnan');
                perf.RT.median_ED = median(x2, 'omitnan');
                perf.RT.median_repeat1 = median(x3, 'omitnan');
                perf.RT.median_repeat2 = median(x4, 'omitnan');

                %[perf.RT.median_ID, perf.RT.median_ED] = plot_his_RTs(task(n), perf_all.subs{i}, x1, 'ID', x2,  'ED');

                % separated by shift or last before shift
                x1 = perf.resp_mat.RT(perf.trial_str.is_last & perf.resp_mat.correct,:);
                x2 = perf.resp_mat.RT(perf.trial_str.is_shift  & perf.resp_mat.correct,:);

                perf.RT.median_last  = median(x1, 'omitnan');
                perf.RT.median_shift = median(x2, 'omitnan');

                %[perf.RT.median_last, perf.RT.median_shift] = plot_his_RTs(task(n), perf_all.subs{i}, x1, 'last', x2,  'shift');

                % reaction time for each shift condition (colour vs shape)
                % ID shifts
                x1 = perf.resp_mat.RT(perf.trial_str.shift_type == 'col_col' & perf.resp_mat.correct & perf.trial_str.trial_type ~= 'practice',:);
                x2 = perf.resp_mat.RT(perf.trial_str.shift_type == 'sh_sh' & perf.resp_mat.correct & perf.trial_str.trial_type ~= 'practice',:);
                % ED shifts
                x3 = perf.resp_mat.RT(perf.trial_str.shift_type == 'sh_col' & perf.resp_mat.correct & perf.trial_str.trial_type ~= 'practice',:);
                x4 = perf.resp_mat.RT(perf.trial_str.shift_type == 'col_sh' & perf.resp_mat.correct & perf.trial_str.trial_type ~= 'practice',:);
                % first repeats separated by colour as a control condition
                x5 = perf.resp_mat.RT(perf.trial_str.dimension == 'colour' & perf.resp_mat.correct & perf.trial_str.is_repeat1,:);
                x6 = perf.resp_mat.RT(perf.trial_str.dimension == 'shape' & perf.resp_mat.correct & perf.trial_str.is_repeat1,:);
                % second repeats separated by colour as a control condition
                x7 = perf.resp_mat.RT(perf.trial_str.dimension == 'colour' & perf.resp_mat.correct & perf.trial_str.is_repeat2,:);
                x8 = perf.resp_mat.RT(perf.trial_str.dimension == 'shape' & perf.resp_mat.correct & perf.trial_str.is_repeat2,:);

                perf.RT.median_IDcol     = median(x1, 'omitnan');
                perf.RT.median_IDsh      = median(x2, 'omitnan');
                perf.RT.median_EDcol     = median(x3, 'omitnan');
                perf.RT.median_EDsh      = median(x4, 'omitnan');
                perf.RT.median_repeat1col = median(x5, 'omitnan');
                perf.RT.median_repeat1sh  = median(x6, 'omitnan');
                perf.RT.median_repeat2col = median(x7, 'omitnan');
                perf.RT.median_repeat2sh  = median(x8, 'omitnan');

                % all errors
                perf.error.error_all =  1 - sum(perf.resp_mat{~perf.trial_str.is_practice, 5})/perf.trial_decomposition(2);

                % errors in different conditions
                perf.error.error_ID     =  1 - sum(perf.resp_mat{perf.trial_str.is_ID, 5})/sum(perf.trial_str.is_ID);
                perf.error.error_ED     =  1 - sum(perf.resp_mat{perf.trial_str.is_ED, 5})/sum(perf.trial_str.is_ED);
                perf.error.error_shift  =  1 - sum(perf.resp_mat{perf.trial_str.is_shift, 5})/sum(perf.trial_str.is_shift);
                perf.error.error_last   =  1 - sum(perf.resp_mat{perf.trial_str.is_last, 5})/sum(perf.trial_str.is_last);
                perf.error.error_repeat1 =  1 - sum(perf.resp_mat{perf.trial_str.is_repeat1, 5})/sum(perf.trial_str.is_repeat1);
                perf.error.error_repeat2 =  1 - sum(perf.resp_mat{perf.trial_str.is_repeat2, 5})/sum(perf.trial_str.is_repeat2);
                perf.error.error_IDcol  =  1 - sum(perf.resp_mat{perf.trial_str.shift_type == 'col_col' &...
                    perf.trial_str.trial_type ~= 'practice', 5})/sum(perf.trial_str.shift_type == 'col_col' & ...
                    perf.trial_str.trial_type ~= 'practice');
                perf.error.error_IDsh   =  1 - sum(perf.resp_mat{perf.trial_str.shift_type == 'sh_sh' &...
                    perf.trial_str.trial_type ~= 'practice', 5})/sum(perf.trial_str.shift_type == 'sh_sh' & ...
                    perf.trial_str.trial_type ~= 'practice');
                perf.error.error_EDcol  =  1 - sum(perf.resp_mat{perf.trial_str.shift_type == 'sh_col' &...
                    perf.trial_str.trial_type ~= 'practice', 5})/sum(perf.trial_str.shift_type == 'sh_col' & ...
                    perf.trial_str.trial_type ~= 'practice');
                perf.error.error_EDsh   =  1 - sum(perf.resp_mat{perf.trial_str.shift_type == 'col_sh' &...
                    perf.trial_str.trial_type ~= 'practice', 5})/sum(perf.trial_str.shift_type == 'col_sh' & ...
                    perf.trial_str.trial_type ~= 'practice');
                perf.error.error_repeat1col  =  1 - sum(perf.resp_mat{perf.trial_str.dimension == 'colour' &...
                    perf.trial_str.is_repeat1, 5})/sum(perf.trial_str.dimension == 'colour' & ...
                    perf.trial_str.is_repeat1);
                perf.error.error_repeat2col  =  1 - sum(perf.resp_mat{perf.trial_str.dimension == 'colour' &...
                    perf.trial_str.is_repeat2, 5})/sum(perf.trial_str.dimension == 'colour' & ...
                    perf.trial_str.is_repeat2);
                perf.error.error_repeat1sh  =  1 - sum(perf.resp_mat{perf.trial_str.dimension == 'shape' &...
                    perf.trial_str.is_repeat1, 5})/sum(perf.trial_str.dimension == 'shape' & ...
                    perf.trial_str.is_repeat1);
                perf.error.error_repeat2sh  =  1 - sum(perf.resp_mat{perf.trial_str.dimension == 'shape' &...
                    perf.trial_str.is_repeat2, 5})/sum(perf.trial_str.dimension == 'shape' & ...
                    perf.trial_str.is_repeat2);

                save(sprintf('%s%s/%s_task-%s_perf.mat', dirs.derived_dt_dir, perf_all.subs{i}, perf_all.subs{i}, task(n)), 'perf')

                % saving new parameters in perf_all
                perf_all.IDED.RT.median_all(i,1)         = perf.RT.median_all;
                perf_all.IDED.RT.skewness_all(i,1)       = perf.RT.skewness_all;
                perf_all.IDED.RT.median_ID(i,1)          = perf.RT.median_ID;
                perf_all.IDED.RT.median_ED(i,1)          = perf.RT.median_ED;
                perf_all.IDED.RT.median_shift(i,1)       = perf.RT.median_shift;
                perf_all.IDED.RT.median_last(i,1)        = perf.RT.median_last;
                perf_all.IDED.RT.median_repeat1(i,1)     = perf.RT.median_repeat1;
                perf_all.IDED.RT.median_repeat2(i,1)     = perf.RT.median_repeat2;
                perf_all.IDED.RT.median_IDcol(i,1)       = perf.RT.median_IDcol;
                perf_all.IDED.RT.median_IDsh(i,1)        = perf.RT.median_IDsh;
                perf_all.IDED.RT.median_EDcol(i,1)       = perf.RT.median_EDcol;
                perf_all.IDED.RT.median_EDsh(i,1)        = perf.RT.median_EDsh;
                perf_all.IDED.RT.median_repeat1col(i,1)  = perf.RT.median_repeat1col;
                perf_all.IDED.RT.median_repeat1sh(i,1)   = perf.RT.median_repeat1sh;
                 perf_all.IDED.RT.median_repeat2col(i,1) = perf.RT.median_repeat2col;
                perf_all.IDED.RT.median_repeat2sh(i,1)   = perf.RT.median_repeat2sh;
                perf_all.IDED.error.error_all(i,1)       = perf.error.error_all;
                perf_all.IDED.error.error_ID(i,1)        = perf.error.error_ID;
                perf_all.IDED.error.error_ED(i,1)        = perf.error.error_ED;
                perf_all.IDED.error.error_shift(i,1)     = perf.error.error_shift;
                perf_all.IDED.error.error_last(i,1)      = perf.error.error_last;
                perf_all.IDED.error.error_repeat1(i,1)   = perf.error.error_repeat1;
                perf_all.IDED.error.error_repeat2(i,1)   = perf.error.error_repeat2;
                perf_all.IDED.error.error_IDcol(i,1)     = perf.error.error_IDcol;
                perf_all.IDED.error.error_IDsh(i,1)      = perf.error.error_IDsh;
                perf_all.IDED.error.error_EDcol(i,1)     = perf.error.error_EDcol;
                perf_all.IDED.error.error_EDsh(i,1)      = perf.error.error_EDsh;
                perf_all.IDED.error.error_repeat1col(i,1)= perf.error.error_repeat1col;
                perf_all.IDED.error.error_repeat1sh(i,1) = perf.error.error_repeat1sh;
                perf_all.IDED.error.error_repeat2col(i,1)= perf.error.error_repeat2col;
                perf_all.IDED.error.error_repeat2sh(i,1) = perf.error.error_repeat2sh;

                clear('perf', 'x', 'x1', 'x2');
                fprintf('%s: Succesfull!\n', perf_all.subs{i});
            catch
                fprintf('%s: An error occured in graph generation\n', perf_all.subs{i});
                continue
            end
        elseif task(n) == "ASST-D" || task(n) == "ASSTD" || task(n) == "ASST"
            try
                % load perf struct for specific task
                load(perf_all.ASSTD.perf_dirs(i), 'perf')
            catch
                fprintf('%s: Perf.mat does not exist\n', perf_all.subs{i})
                continue
            end
            try
                % create indices if not already there
                %%%% temporary solution
                perf.resp_mat.is_last_cor    = perf.resp_mat.cons_corr == 8;
                % shifting it out by one to get shifts
                perf.resp_mat.is_shift = [false; perf.resp_mat.is_last_cor(1:end-1,:)];
                perf.resp_mat.is_practice = perf.resp_mat.stage_num == 1;
                perf.resp_mat.is_first_cor = perf.resp_mat.cons_corr == 1;
                perf.resp_mat.is_second_cor = perf.resp_mat.cons_corr == 2;

                % all RTs
                x = perf.resp_mat.RT( ~ perf.resp_mat.is_practice,:);

                perf.RT.median_all = median(x, 'omitnan');
                %perf.RT.median_all = plot_his_RTs(task(n), perf_all.subs{i}, x, 'all RTs');

                % separated by shift or no shift
                x1 = perf.resp_mat.RT(~perf.resp_mat.is_shift,:);
                x2 = perf.resp_mat.RT(perf.resp_mat.is_shift,:);

                perf.RT.median_repeat  = median(x1, 'omitnan');
                perf.RT.median_shift = median(x2, 'omitnan');

                %[perf.RT.median_repeat, perf.RT.median_shift] = plot_his_RTs(task(n), perf_all.subs{i}, x1, 'repeat', x2,  'shift');

                % first correct vs second correct
                x1 = perf.resp_mat.RT(perf.resp_mat.is_second_cor,:);
                x2 = perf.resp_mat.RT(perf.resp_mat.is_first_cor,:);

                perf.RT.median_second = median(x1, 'omitnan');
                perf.RT.median_first  = median(x2, 'omitnan');

                %[perf.RT.median_second, perf.RT.median_first] = plot_his_RTs(task(n), perf_all.subs{i}, x1, 'Second correct', x2,  'First Correct');

                % first correct vs last correct
                x1 = perf.resp_mat.RT(perf.resp_mat.is_last_cor,:);
                x2 = perf.resp_mat.RT(perf.resp_mat.is_first_cor,:);

                perf.RT.median_last = median(x1, 'omitnan');

                %perf.RT.median_last = plot_his_RTs(task(n), perf_all.subs{i}, x1, 'Last correct', x2,  'First Correct');

                save(sprintf('%s%s/%s_task-%s_perf.mat', dirs.derived_dt_dir, perf_all.subs{i}, perf_all.subs{i}, task(n)), 'perf')

            catch
                fprintf('%s: An error occured during graph generation of RTs\n', perf_all.subs{i})
                continue
            end

            try
                % saving new parameters in perf_all
                perf_all.ASSTD.RT.median_all(i,1)    = perf.RT.median_all;
                perf_all.ASSTD.RT.median_repeat(i,1) = perf.RT.median_repeat;
                perf_all.ASSTD.RT.median_shift(i,1)  = perf.RT.median_shift;
                perf_all.ASSTD.RT.median_first(i,1)  = perf.RT.median_first;
                perf_all.ASSTD.RT.median_second(i,1) = perf.RT.median_second;
                perf_all.ASSTD.RT.median_last(i,1)   = perf.RT.median_last;

                % Error rates
                perf.error.error_all = 1 - (sum(perf.resp_mat.correct)/perf.all_trials);
                perf.error.error_in_stage = nan(1, height(perf.stage_str));
                perf.error.trials_to_criterion  = nan(1, height(perf.stage_str));
                for k = 1: height(perf.stage_str)
                    ind = perf.resp_mat.stage_num == k;
                    x = perf.resp_mat(ind,:);
                    perf.error.error_in_stage(k) = 1 - (sum(x.correct)/height(x));
                    perf.error.trials_to_criterion(k)  = height(x);
                end

                % saving new parameters in perf_all
                perf_all.ASSTD.error.error_all(i,1)              = perf.error.error_all;
                perf_all.ASSTD.error.error_in_stage(i, :)      = perf.error.error_in_stage;
                perf_all.ASSTD.error.trials_to_criterion(i, :) = perf.error.trials_to_criterion;

                save(sprintf('%s%s/%s_task-%s_perf.mat', dirs.derived_dt_dir, perf_all.subs{i}, perf_all.subs{i}, task(n)), 'perf')

                clear('perf', 'x', 'ind');

                fprintf('%s: Succesfull!\n', perf_all.subs{i})
            catch
                fprintf('%s: An error occured during error calculation\n', perf_all.subs{i})
                continue
            end
        end
    end
end

%% Saving the results
clear('i', 'k', 'n', 'subs', 'task')
try
    save(sprintf('%sall/%s_perf_all.mat', dirs.derived_dt_dir, datetime('today', 'Format', 'yyyy-MM-dd')))
catch
    fprintf('Could not save perf_all.mat');
end
end
