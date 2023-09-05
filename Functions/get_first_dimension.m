function dimensions = get_first_dimension(ASSTD_stim_dir, ASSTG_stim_dir)
% prepare table to save everything
% create emtpy vectors
Subject = strings(500,1);
ASSTD   = strings(500,1);
ASSTG   = strings(500,1);

for subject = 1:500

Subject(subject) = sprintf('S%03i', subject);
load(sprintf('%sS%03i_stimuli.mat', ASSTD_stim_dir, subject), 'target_dim');
ASSTD(subject) = target_dim(1);
load(sprintf('%sS%03i_stimuli.mat', ASSTG_stim_dir, subject), 'target_dim');
ASSTG(subject) = target_dim(1);
clear 'target_dim'
dimensions = [Subject, ASSTD, ASSTG];
end
end