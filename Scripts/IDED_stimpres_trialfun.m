function trl = IDED_stimpres_trialfun(cfg)
% this function requires the following fields to be specified
% cfg.dataset
% cfg.trialdef.eventtype
% cfg.trialdef.eventvalue
% cfg.trialdef.prestim
% cfg.trialdef.poststim
% cfg.trialdef.trigger_type can be 'nexttrigger', 'previoustrigger'
% we use 'nexttrigger' when we want trials to be stimulus locked
% we use 'previoustrigger' when we want trials to be response locked

% define header and events
hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

% define variables
trl = [];
art_ind = false (length(event), 1);
stim_ind = false (length(event), 1);
artifact = [];
rejected = 0;

% locate Stimulus and 'Bad Interval'events
for i=1:length(event)
    if strcmp(event(i).type, 'Bad Interval')
        art_ind(i) = true;
        artifact(i,:) = [event(i).sample, event(i).sample+event(i).duration]; % [artifact_start artifact_finish]
    elseif strcmp(event(i).type, 'Stimulus')
        stim_ind(i) = true;
    end
end

% count practice trials, to make check if we have more practice trials than
% expected
n_practice = 0;
n_trial = 0;
for i=1:length(event)
    if strcmp(event(i).value, 'S211') || strcmp(event(i).value, 'S212')
        n_practice = n_practice + 1;
    end
    if strcmp(event(i).value, 'S 10')
        n_trial = n_trial + 1;
    end
end

% find the trials add add trial number
% find presentation of fixation cross (S10)
if n_trial == 599
    n = 0;
else
    % if n_trial is not what we expected give warning
    fprintf('WARNING!! Number of trials does not much expected number!\n Checking if more practice trials were run\n')
    % if n_practice is higher than 20 it means that the participant required
    % more than one round of practice trials. In order to have compatibility
    % with perf.mat we change the starting n. This way we can ensure
    % trial-to-trial match between the eeg and the perf data. (Practice
    % trials are rejected anyway)
    if n_practice > 20
        % calculate how many additional trials we have
        extra_practice_trl = n_trial - 599;
        n = 0 - extra_practice_trl; % count the last twenty, because they are always counted
        fprintf('Extra practice trials removed succesfully!\n')
    else
        fprintf('WARNING!! Could not find reason of unexpected number of trials. Please check manually!\n Number of trials: %d\nStill defining n as normal:', ...
            n_trial)
        n = 0;
    end
end
for i=1:length(event)
    if strcmp(event(i).value, 'S 10')        
        n = n + 1;
    end
    event(i).trial_num = n;
end

% save index for stimuli and artifacts
artifact = artifact(art_ind,:);
stimulus = event(stim_ind);

% Calculate timing between events
for i = 1:length(stimulus)
    if i == 1
        stimulus(i).dt_sample = 0;
    else
        stimulus(i).dt_sample = stimulus(1,i).sample - stimulus(1,i-1).sample;
    end
end

% reject trials with artifacts (as defined in Brain Vision Analyzer)
% reject epochs with marked events 'Bad Interval'
% add more information in trl definition
for i=1:length(stimulus)
    if strcmp(stimulus(i).type, cfg.trialdef.eventtype)
        % it is a trigger, see whether it has the right value
        if ismember(stimulus(i).value, cfg.trialdef.eventvalue)
            % add this to the trl definition
            begsample    = stimulus(i).sample - cfg.trialdef.prestim*hdr.Fs;
            endsample    = stimulus(i).sample + cfg.trialdef.poststim*hdr.Fs - 1;
            offset       = -cfg.trialdef.prestim*hdr.Fs;
            trigger      = stimulus(i).value; % remember the trigger (=condition) for each trial
            trial_num    = stimulus(i).trial_num; % get trial number from table
                        % depending whether we want the nexttrigger (when stimulus
            % locked) or previous trigger (when response locked), define
            % adjtrigger
            if strcmp(cfg.trialdef.trigger_type, 'nexttrigger')
                if i == length(stimulus)
                    adjtrigger  = 'NaN';
                else
                    adjtrigger  = stimulus(i + 1).value; % the condition of the next trial
                end
                % Reaction time in sample points (not in ms!!)
                RT = stimulus(i + 1).dt_sample;
            elseif  strcmp(cfg.trialdef.trigger_type, 'previoustrigger')
                if i == 1
                    adjtrigger  = 'NaN';
                else
                    adjtrigger  = stimulus(i - 1).value; % the previous event
                end
                RT = stimulus(i).dt_sample; % this is for event locked data
            end
            % check if the trial overlaps with an artifact
            % check if any artifact starts before epoch begins but finishes after epoch begins 
            ind_beg = begsample >= artifact(:,1) & begsample <= artifact(:,2);
            % check if any artifact starts before epoch ends but finishes after epoch ends 
            ind_end = endsample >= artifact(:,1) & endsample <= artifact(:,2);
            % check if any artifact starts after epoch starts but finishes before epoch finishes            
            ind_in = begsample <= artifact(:,1) & endsample >= artifact(:,2);
            if ~any(ind_beg) && ~any(ind_end) && ~any(ind_in)
                % include trials without artifacts
                % columns are: start end offset trigger response trial_num
                trl(end+1, :) = [round([begsample endsample offset])  ...
                str2double(trigger(2:end)) str2double(adjtrigger(2:end)) trial_num RT];
            else
                % reject trials with artifact
                rejected = rejected + 1;
            end
            
        end
    end
end
% inform user about number of rejected trials
fprintf('rejected %i trials with artifacts \n', rejected)
end
