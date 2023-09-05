function plot_ERP_custom(dataset1, dataset2, time)
bounds1 = std(dataset1, [], 1);
bounds2 = std(dataset2, [], 1);
%linecolor = [0 0 0 0.08];
%plot(time, data_all.avg_young', 'Color', linecolor);
hold on
boundedline(time, mean(dataset1, 1), bounds1, 'alpha', 'b'); % alpha makes bounds transparent
boundedline(time, mean(dataset2, 1), bounds2, 'alpha', 'r'); % alpha makes bounds transparent
xlabel('Time (ms)');
ylabel('Potential (μV)');
ax = gca(); % this Gets the Current Axis so we can set properties
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
%ax.TickDir = 'out';
% Remove the box around the plot
box off;
% And move the x-axis label to underneath the axis:
ax.XLabel.Position(2) = 0;