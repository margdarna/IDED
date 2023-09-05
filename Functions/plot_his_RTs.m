function [median_x1, median_x2] = plot_his_RTs(task, subj, x1, value1, x2,  value2)
% This function plots histograms of reaction time distribution of one
% subject.
%
% INPUT:
%   task    -   string, task name
%   subj    -   string, subject pseudonym
%   x1      -   array, values that should be plotted
%   value1  -   string, name given to array x1 (will be presented on graph)
%   x2      -   array, optional, second values that should be plotted on top of x1
%   value2  -   string, optional, name given to array x2 (will be presented on graph)
%
% OUTPUT
%  figure
%  median_x1 - median value for x1
%  median:x2 - median value for x2
% Created on    05/02/2022
% Last modified 05/02/2022
%
% Created by: Margarita Darna
% margarita.darna@lin-magdeburg.de

% Check if both x1 and x2 were input
if exist('x2') == 1

    % calculate median
    median_x1 = median(x1, 'omitnan');% calculate median
    median_x2 = median(x2, 'omitnan');

    figure()
    hold on
    histogram(x1, 0:0.05:4, 'Normalization','probability');
    histogram(x2, 0:0.05:4, 'Normalization','probability');
    plot([median_x1 median_x1], [0 1], 'b');
    plot([median_x2 median_x2], [0 1], 'r');
    txt1 = {sprintf('%s: Median = %.3f', value1, median_x1)};
    txt2 = {sprintf('%s: Median = %.3f', value2, median_x2)};
    text(1.2, 0.2, txt1,  'Color', 'b')
    text(1.2, 0.18, txt2,  'Color', 'r')
    xlim([0, 2]);
    ylim([0 0.3])
    xlabel('RT (s)', 'FontSize', 14);
    ylabel('Bin Count', 'FontSize', 14);
    title(sprintf('%s - %s', task, subj), 'FontSize', 14);
    hold off
else
    % calculate median
    median_x1 = median(x1, 'omitnan');
    % plot histogram of RTs
    figure()
    hold on
    h =histogram(x1, 0:0.05:4, 'Normalization','probability');
    plot([median_x1 median_x1], [0 1], 'b');
    txt = {sprintf('%s: Median = %.3f', value1, median_x1)};
    text(1.2, 0.2, txt, 'Color', 'b')
    xlim([0, 4]);
    ylim([0 0.3])
    xlabel('RT (s)', 'FontSize', 14);
    ylabel('Bin Count', 'FontSize', 14);
    title(sprintf('%s - %s', task, subj), 'FontSize', 14);
    hold off
end


