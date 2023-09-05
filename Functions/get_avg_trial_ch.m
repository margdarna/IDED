function avg = get_avg_trial_ch(data, ch_ind)
for i = 1:numel(data)
    data_avg = data{i}.avg(ch_ind, :);
    % compute average of each column
    avg(i,:) = mean(data_avg, 1);
end