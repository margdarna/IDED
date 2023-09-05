%----------------------------------------------------------------------
%                       Description
%----------------------------------------------------------------------

% This scripts allocates stimuli of IDED paradigm to the participants
% and saves those in mat files.
%
% Written by: Margarita Darna
% margarita.darna@lin-magdeburg.de
% Leibniz Institute for Neurobiology, Magdeburg
%
%
% Version 2:
% Fixed wrong allocation of is_shift and is_last
%
% Version 3:
% Fixed bug in stimulus presentation after break
% Created:     03.12.2021
% Last edited: 01.09.2023

%----------------------------------------------------------------------
%                       Preamble
%----------------------------------------------------------------------
clear; clc;

% Define directories
% from office PC
% project_dir = 'H:\PhD\A05-SFB1436\';
% from laptop
project_dir = '//linstore01/home/mdarna/PhD/A05-SFB1436/';
%project_dir = 'C:/Users/mdarna/Documents/PhD/A05-SFB1436/';
exp_dir = sprintf('%sCode/IDED/', project_dir);
subj_stim_dir = sprintf('%sSubject_Stimuli/', exp_dir);

 % load practice stimuli and add later on in arrays
 try
     load(sprintf('%spractice_stimuli.mat', subj_stim_dir));
 catch
     fprintf('Could not find practice stimuli mat file in chosen directory\n', subject);
     return;
 end

%----------------------------------------------------------------------
%                      Preparing Parameters
%----------------------------------------------------------------------
% define total number of trials per condition
total_ED     = 50;
total_ID     = total_ED;
trial_types  = categorical(["ID", "ED", "pre", "repeat", "practice"]);
% Indices 1 = intra-dimensional shift, 2 = extra-dimensional shift,  3 = pre-shift,
% 4 = repeat, 5 = practice
dimensions   = categorical(["colour", "shape"]);
% Indices: 1 = colour, 2 = shape
colours      = categorical(["red", "green", "blue", "yellow", "magenta", "cyan"]);
% Indices: 1 = red, 2 = green, 3 = blue, 4 = yellow, 5 = magenta, 6 = cyan
shapes       = categorical(["circle", "star", "diamond", "triangle", "pentagon", "cross"]);
% Indices: 1 = circle, 2 = star, 3 = diamond, 4 = triangle, 5 = pentagon, 6 = cross
shift_types  = categorical(["col_col", "sh_sh", "col_sh", "sh_col", "none"]);

% Indices: 
col_ind  = 1:6;
sh_ind   = 1:6;
shift_ind         = [ones(total_ID,1); repmat(2,total_ED,1)]; % 1 = ID, 2 = ED;

% % we want to equally divide number of trials in each set so we calculate
% that
trials_in_set_num = floor((numel(shift_ind)+6)/6);
% calculate number of shifts we get and see how much we are missing
missing_set_num = (numel(shift_ind)+6) - trials_in_set_num*6;
% add the remaining sets, starting with sets that have the lowest number of
% trials
trials_in_set_ind = [repmat(3,trials_in_set_num,1); repmat(4,trials_in_set_num,1); repmat(5,trials_in_set_num,1); ...
    repmat(6,trials_in_set_num,1); repmat(7,trials_in_set_num,1);  repmat(8,trials_in_set_num,1)];
for i = 1:missing_set_num
   trials_in_set_ind = [trials_in_set_ind; repmat(i+2,1,1)]; 
end
dimensions_ind    = [1, 2];
total_exp_trials  = sum(trials_in_set_ind);
practice_trials   = numel(trial_num_pr);
all_trials        = practice_trials + total_exp_trials;
trial_num         = 1:all_trials;
% take a break approximately after the calculation of trials_in_set_num
% (aka take a break after approximately an equal amount of sets)
take_break_at_set = [];
for i = 1:5
    take_break_at_set = [take_break_at_set i*trials_in_set_num];
end

%----------------------------------------------------------------------
%                       Stimulus allocation
%----------------------------------------------------------------------

for subject = 1:500
    run = 1;
    % Setting parameter for allocating
    allocate = true;
    while allocate
        % Setting parameters
        target_col    = [target_col_pr; categorical(nan(total_exp_trials,1))];
        target_sh     = [target_sh_pr; categorical(nan(total_exp_trials,1))];
        match_col     = [match_col_pr; categorical(nan(total_exp_trials,1))];
        match_sh      = [match_sh_pr; categorical(nan(total_exp_trials,1))];
        no_match_col  = [no_match_col_pr; categorical(nan(total_exp_trials,1))];
        no_match_sh   = [no_match_sh_pr; categorical(nan(total_exp_trials,1))];
        dimension     = [dimension_pr; categorical(nan(total_exp_trials,1))];
        shift_type    = [shift_type_pr; categorical(nan(total_exp_trials,1))];
        trial_type    = [trial_type_pr; categorical(nan(total_exp_trials,1))];
        set_trial     = [set_trial_pr; zeros(total_exp_trials,1)];
        is_shift      = [is_shift_pr; false(total_exp_trials,1)];
        is_last       = [is_last_pr; false(total_exp_trials,1)];
        is_break      = [false(all_trials, 1)];
        set = 1;
        trial = practice_trials + 1;

        % Shuffling arrays
        sh_shift_ind = Shuffle(shift_ind, 2);
        % add the break trials after break trials to make sure that we have
        % the correct number of elements in sh_shift_ind
        % 3 marks sets after breaks
        sh_shift_ind = [3; sh_shift_ind];
        for i = 1:5
        sh_shift_ind = [sh_shift_ind(1:take_break_at_set(i)-1); 3; sh_shift_ind(take_break_at_set(i):end)];
        end
        sh_trials_in_set_ind = Shuffle(trials_in_set_ind);
        while trial <= all_trials
            trial_in_set = 1;
            while trial_in_set <= sh_trials_in_set_ind(set)
                % Defining important variables
                shift = sh_shift_ind(set);
                if ismember(set, [1 take_break_at_set])
                    shift_type_ind = 5;
                    trial_type_ind = 3;
                elseif trial_in_set == 1
                    is_shift(trial) = true;
                    if shift == 1 % Intra-dimensional shift
                        trial_type_ind = 1;
                        % Checking previous dimension to save shift_type
                        if dimension_ind == 1
                            shift_type_ind = 1; % colour to colour
                        elseif dimension_ind == 2
                            shift_type_ind = 2; % shape to shape
                        end
                    elseif shift == 2 % Extra-dimensional shift
                        trial_type_ind = 2;
                        % Checking previous dimension to save shift_type
                        if dimension_ind == 1
                            shift_type_ind = 3; % colour to shape
                        elseif dimension_ind == 2
                            shift_type_ind = 4; % shape to color
                        end
                    end
                else
                    shift_type_ind = 5;
                    trial_type_ind = 4;
                end
                if trial_in_set == sh_trials_in_set_ind(set)
                    is_last(trial) = true;
                end

                % Defining targets
                if ismember(set, [1 take_break_at_set]) && trial_in_set == 1
                    % in first trial after break randomly assign target colour and shape
                    target_col_ind = randsample(col_ind, 1);
                    target_sh_ind = randsample(sh_ind, 1);
                elseif trial_in_set == 1
                    % change target completely (it should not have the colour or
                    % shape of the previous target or previous matching stimulus)
                    target_col_ind = randsample(col_ind(col_ind ~= target_col_ind & col_ind ~= match_col_ind ), 1);
                    target_sh_ind = randsample(sh_ind(sh_ind ~= target_sh_ind & sh_ind ~= match_sh_ind ), 1);
                end

                % Defining matching stimuls
                if ismember(set, [1 take_break_at_set]) && trial_in_set == 1
                    %randomly choose which dimension is the relevant one
                    dimension_ind = randsample(dimensions_ind,1);
                    if dimension_ind == 1 % colour
                        % match by colour
                        match_col_ind = target_col_ind;
                        % different shape
                        match_sh_ind = randsample(sh_ind(sh_ind ~= target_sh_ind), 1);
                    elseif dimension_ind == 2 % shape
                        % different colour
                        match_col_ind  = randsample(col_ind(col_ind ~= target_col_ind), 1);
                        % match by shape
                        match_sh_ind = target_sh_ind;
                    else
                        fprintf("ERROR: shift not found in first trial\n")
                        return;
                    end
                else
                    if trial_in_set == 1 && shift == 2 % Extra-dimensional shift
                        % dimension changes
                        dimension_ind = dimensions_ind(dimensions_ind ~= dimension_ind);
                    end
                    if dimension_ind == 1
                        % match by colour
                        match_col_ind = target_col_ind;
                        % different shape from target and previous matching  stimulus
                        match_sh_ind = randsample(sh_ind(sh_ind ~= match_sh_ind & sh_ind ~= target_sh_ind), 1);
                    elseif dimension_ind == 2
                        % different colour from target and previous matching  stimulus
                        match_col_ind = randsample(col_ind(col_ind ~= match_col_ind & col_ind ~= target_col_ind), 1);
                        % match by shape
                        match_sh_ind = target_sh_ind;
                    else
                        fprintf('ERROR: dimension_ind not found\n')
                        return;
                    end
                end

                % Defining non-matching stimulus
                if ismember(set, [1 take_break_at_set])
                    no_match_col_ind = randsample(col_ind(col_ind ~= match_col_ind & col_ind ~= target_col_ind), 1);
                    no_match_sh_ind  = randsample(sh_ind(sh_ind ~= match_sh_ind & sh_ind ~= target_sh_ind), 1);
                else
                    no_match_col_ind = randsample(col_ind(col_ind ~= match_col_ind & col_ind ~= target_col_ind ...
                        & col_ind ~= no_match_col_ind), 1);
                    no_match_sh_ind  = randsample(sh_ind(sh_ind ~= match_sh_ind & sh_ind ~= target_sh_ind ...
                        & sh_ind ~= no_match_sh_ind), 1);
                end

                % checking if we have a break
                if ismember(set, take_break_at_set) && trial_in_set == 1
                    is_break(trial - 1) = true;
                end

                % Completing arrays
                target_col(trial)   = colours(target_col_ind);
                target_sh(trial)    = shapes(target_sh_ind);
                match_col(trial)    = colours(match_col_ind);
                match_sh(trial)     = shapes(match_sh_ind);
                no_match_col(trial) = colours(no_match_col_ind);
                no_match_sh(trial)  = shapes(no_match_sh_ind);
                dimension(trial)    = dimensions(dimension_ind);
                shift_type(trial)   = shift_types(shift_type_ind);
                trial_type(trial)   = trial_types(trial_type_ind);
                set_trial(trial)    = trial_in_set;

                % increasing index of trials
                trial_in_set = trial_in_set + 1;
                trial = trial + 1;
            end
            set = set + 1;
        end

        % Check with count cats if we have a balanced set of shif types, if any of
        % the appearances for any shift type is smaller than 25 then we repeat
        % stimulus allocation for this participant
        if min(countcats(shift_type)) < total_ED/2 + 1 % adding one because of the practice trials
            allocate = true; % allocate stimuli again
            run = run + 1;
            % break after 1000 unsuccesfull runs
        else
            allocate = false;
            fprintf('S%d succesful!!\n Saving data... \n', subject)
        end
    end

    % Creating table with variables
    trial_str = table(trial_num', set_trial, trial_type, shift_type, dimension, target_col, target_sh, ...
        match_col, match_sh, no_match_col, no_match_sh, is_shift, is_last, is_break);

    % save mat file
    save(sprintf('%sS%03d_stimuli.mat', subj_stim_dir, subject), 'trial_num', 'set_trial', 'trial_type', 'shift_type', ...
        'dimension', 'target_col', 'target_sh', 'match_col', 'match_sh', 'no_match_col', 'no_match_sh', 'is_shift', ...
        'is_last', 'is_break', 'trial_str', 'practice_trials', 'total_exp_trials', 'all_trials');
end