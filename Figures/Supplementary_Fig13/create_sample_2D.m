clc;
clear;
close all;

%% parameter
period = 1;
time_interval = 1/100;
tspan = linspace(0,period, period/time_interval+1);

% create input signal
time_interval_pre = 1/4;
bt_pre = linspace(0,period, period/time_interval_pre+1);
bt = linspace(0,period, period/time_interval+1);
ct_pre = linspace(0,period, period/time_interval_pre+1);
ct = linspace(0,period, period/time_interval+1);
num_data = 100;

%% create time series
y_total = {};
for i = 1:num_data
    % create source
    b_pre = rand([length(bt_pre),1])*2-1;
    b_tmp = spline(bt_pre,b_pre,bt);
    c_pre = rand([length(ct_pre),1])*2-1;
    c_tmp = spline(ct_pre,c_pre,ct);
    y_0 = rand(1)*2-1;
    
    % create other timeseries
    y_tmp = zeros(length(tspan),3);
    [t0, y0] = ode45(@(t,y) example_2D(t,y,bt,b_tmp,ct,c_tmp), tspan, y_0);
    y_tmp(:,1) = y0;
    y_tmp(:,2) = b_tmp;
    y_tmp(:,3) = c_tmp;
    
    %normalizing
    y1 = y_tmp;
    y2 = y_tmp;
    for j = 1:3
        y1(:,j) = y1(:,j) - min(y1(:,j));
        y2(:,j) = y1(:,j)/max(y1(:,j));
    end
    y_total{end+1} = y2;
end



%% plot time series
figure(1)
plot(t0, y2(:,1), 'r','linewidth',2)
hold on
plot(t0, y2(:,2), 'g','linewidth',2)
hold on
plot(t0, y2(:,3), 'b','linewidth',2)

%% save data

t = tspan.';
filename = '2D_timeseries';
save(filename, 'y_total', 't', 'time_interval','num_data')