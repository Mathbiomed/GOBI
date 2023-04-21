clc;
clear;
close all;

%% parameter
period = 2*pi;
time_interval = 1/30;
tspan = linspace(0,period, period/time_interval+1);

time_interval_pre = pi/4;
zt_pre = linspace(0,period, period/time_interval_pre+1);
zt = linspace(0,period, period/time_interval+1);
num_data = 100;

%% create time series
y_total = {};
for i = 1:num_data
    % create source
    z_pre = rand([length(zt_pre),1]) - 0.5;
    z_tmp = spline(zt_pre,z_pre,zt);
    y_0 = rand([1,3]) - 0.5;
    
    % create other timeseries
    y_tmp = zeros(length(tspan),4);
    [t0, y0] = ode15s(@(t,y) sample(t,y,zt,z_tmp), tspan, y_0);
    y_tmp(:,2:4) = y0;
    y_tmp(:,1) = z_tmp;
    
    %normalizing
    y1 = y_tmp;
    y2 = y_tmp;
    for j = 1:4
        y1(:,j) = y1(:,j) - min(y1(:,j));
        y2(:,j) = y1(:,j)/max(y1(:,j));
    end
    y_total{end+1} = y2;
end



%% plot time series
figure(1)
plot(t0, y2(:,1), 'k','linewidth',2)
hold on
plot(t0, y2(:,2), 'r','linewidth',2)
hold on
plot(t0, y2(:,3), 'g','linewidth',2)
hold on
plot(t0, y2(:,4), 'b','linewidth',2)

%% save data
t = tspan.';
filename = 'sample_timeseries';
save(filename, 'y_total', 't', 'time_interval','num_data')