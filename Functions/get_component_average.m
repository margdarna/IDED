function component_avg = get_component_average(cfg, component_name, dirs)

% Load single subject data to extract the channels and time range of
% interest
load(strcat(dirs.output_dir, '2_Subj_Avg\young_subj_avg_stimpres_ERP'));
load(strcat(dirs.output_dir, '2_Subj_Avg\old_subj_avg_stimpres_ERP'));

for i = 1:numel(young_repeat1_stimpres_ERP)
cfg.latency = cfg.latency_young;
young_repeat1_stimpres{i} = ft_selectdata(cfg, young_repeat1_stimpres_ERP{i});
young_repeat2_stimpres{i} = ft_selectdata(cfg, young_repeat2_stimpres_ERP{i});
young_shift_stimpres{i} = ft_selectdata(cfg, young_shift_stimpres_ERP{i});
young_ID_stimpres{i} = ft_selectdata(cfg, young_ID_stimpres_ERP{i});
young_ED_stimpres{i} = ft_selectdata(cfg, young_ED_stimpres_ERP{i});
young_last_stimpres{i} = ft_selectdata(cfg, young_last_stimpres_ERP{i});
young_alltrl_stimpres{i} = ft_selectdata(cfg, young_alltrl_stimpres_ERP{i});
end

for i = 1:numel(old_repeat1_stimpres_ERP)
cfg.latency = cfg.latency_old;
old_repeat1_stimpres{i} = ft_selectdata(cfg, old_repeat1_stimpres_ERP{i});
old_repeat2_stimpres{i} = ft_selectdata(cfg, old_repeat2_stimpres_ERP{i});
old_shift_stimpres{i} = ft_selectdata(cfg, old_shift_stimpres_ERP{i});
old_ID_stimpres{i} = ft_selectdata(cfg, old_ID_stimpres_ERP{i});
old_ED_stimpres{i} = ft_selectdata(cfg, old_ED_stimpres_ERP{i});
old_last_stimpres{i} = ft_selectdata(cfg, old_last_stimpres_ERP{i});
old_alltrl_stimpres{i} = ft_selectdata(cfg, old_alltrl_stimpres_ERP{i});
end

% pool all subject averages together using the grandaverage function
% for now am just focusing on repeat and ID and ED trials
cfg = [];
cfg.channel        = 'all';
cfg.latency        = 'all';
cfg.parameter      = 'avg';
cfg.keepindividual = 'yes';

component_avg.young_repeat1 = ft_timelockgrandaverage(cfg, young_repeat1_stimpres{:});
component_avg.young_repeat2 = ft_timelockgrandaverage(cfg, young_repeat2_stimpres{:});
component_avg.young_ID     = ft_timelockgrandaverage(cfg, young_ID_stimpres{:});
component_avg.young_ED     = ft_timelockgrandaverage(cfg, young_ED_stimpres{:});
component_avg.young_last   = ft_timelockgrandaverage(cfg, young_last_stimpres{:});
component_avg.old_repeat1   = ft_timelockgrandaverage(cfg, old_repeat1_stimpres{:});
component_avg.old_repeat2   = ft_timelockgrandaverage(cfg, old_repeat2_stimpres{:});
component_avg.old_ID       = ft_timelockgrandaverage(cfg, old_ID_stimpres{:});
component_avg.old_ED       = ft_timelockgrandaverage(cfg, old_ED_stimpres{:});
component_avg.old_last     = ft_timelockgrandaverage(cfg, old_last_stimpres{:});

% save variables
save(strcat(dirs.output_dir, '4_Stats\stats_ERP_', component_name, '.mat'), "young_repeat1_stimpres", ...
    "young_repeat2_stimpres","young_shift_stimpres", "young_ID_stimpres", "young_ED_stimpres", "young_last_stimpres", ...
    "young_alltrl_stimpres", "old_repeat1_stimpres", "old_repeat2_stimpres", "old_shift_stimpres", "old_ID_stimpres", ...
    "old_ED_stimpres", "old_last_stimpres", "old_alltrl_stimpres","component_avg");