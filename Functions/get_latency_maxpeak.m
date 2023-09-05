function [latency, max_value] = get_latency_maxpeak(ERP, time_start, sample_frequency)
    ERP_avg = ERP;
    [max_value, ind_ERP] = max(ERP_avg); % retrieve the frame where maximum is found
    % translate frame into data point
    latency = time_start + ind_ERP/sample_frequency;
end