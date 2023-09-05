function latency = get_latency(ERP)
    total_area = trapz(ERP.time,ERP.avg_cor);
if total_area == 0
    % some participants have only relative ERPs, they slowly go down.
    % to combat that I am here calculating the negative area under the x
    % curve and finding the middle. This should be an appropriate solution
    % to this problem
     ERP.avg_cor = ERP.avg;
     ERP.avg_cor(ERP.avg > 0)  = 0;
     total_area = trapz(ERP.time,ERP.avg);     
end
for k = 2:numel(ERP.time)-1
    area = trapz(ERP.time(1:k),ERP.avg_cor(1:k));
    difference(k-1) = total_area/2 - area;
end
% calculate difference from entire area
[~, ind]= min(abs(difference));
latency = ERP.time(ind + 1);
end
