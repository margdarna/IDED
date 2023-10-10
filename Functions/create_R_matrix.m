function create_R_matrix(component_name, analysis, subj_info, dirs)
% create R matrix for variable
% component_name: component name as character vector, e.g. 'P300'
% analysis: type of analysis to be performed as character vector, e.g. 'ERP' or 'TFR'
% subj_info: subj_info array, as created in original scripts.
% dirs: struc with experiment directories as generated in original scripts.

% load mat file of that component
load(strcat(dirs.output_dir, '4_Stats\stats_', analysis, '_', component_name, '.mat'), "component_avg")

% get age group from subject info
subs      = subj_info.Pseudonym;
age_group = subj_info.age_cohort;
is_young  = categorical(age_group) == 'young';

dem_info =  readtable(strcat(dirs.raw_dt_dir, 'Demographics.xlsx'));
% arrange the age in the order of subjects
for i = 1:numel(subs)
    subj = subs{i};
    %age(i) = dem_info.Age(strcmp(dem_info.Pseudonym,subj), :);
    gender(i) = dem_info.Gender(strcmp(dem_info.Pseudonym,subj), :);
end

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
T.gender =  [gender(is_young)' - 1; gender(~is_young)' - 1; gender(is_young)' - 1; gender(~is_young)' - 1; ...
    gender(is_young)' - 1; gender(~is_young)' - 1;];
% save table as csv
writetable(T, strcat(dirs.output_dir, '4_Stats\stat_', analysis, '_', component_name, '.csv'));
