clc;
clear;
close all;

%% parameter
period = 2*pi;
time_interval = 1/30;
tspan = linspace(0,period, period/time_interval+1);

time_interval_pre = pi/2; % pi/2 for SFL & IFL, pi/4 for CFL
zt_pre = linspace(0,period, period/time_interval_pre+1);
zt = linspace(0,period, period/time_interval+1);
num_data = 100;
num_component = 3;

%% create time series
y_total_1 = {};
y_total_2 = {};
y_total_3 = {};

for i = 1:num_data
    % create source
    z_pre = rand([length(zt_pre),1]) - 0.5;
    z_tmp = spline(zt_pre,z_pre,zt);
    y_0 = rand([1,2]) - 0.5;
    
    % create other timeseries
    y_tmp = zeros(length(tspan),3);
    [t0, y0] = ode45(@(t,y) IFL(t,y,zt,z_tmp), tspan, y_0);
    y_tmp(:,2:3) = y0;
    y_tmp(:,1) = z_tmp;
    
    %normalizing
    y1 = y_tmp;
    y2 = y_tmp;
    for j = 1:3
        y1(:,j) = y1(:,j) - min(y1(:,j));
        y2(:,j) = y1(:,j)/max(y1(:,j));
    end
    y_total_1{end+1} = y2;
    
    % create other timeseries
    y_tmp = zeros(length(tspan),3);
    [t0, y0] = ode45(@(t,y) CFL(t,y,zt,z_tmp), tspan, y_0);
    y_tmp(:,2:3) = y0;
    y_tmp(:,1) = z_tmp;
    
    %normalizing
    y1 = y_tmp;
    y2 = y_tmp;
    for j = 1:3
        y1(:,j) = y1(:,j) - min(y1(:,j));
        y2(:,j) = y1(:,j)/max(y1(:,j));
    end
    y_total_2{end+1} = y2;
    
    % create other timeseries
    y_tmp = zeros(length(tspan),3);
    [t0, y0] = ode45(@(t,y) SFL(t,y,zt,z_tmp), tspan, y_0);
    y_tmp(:,2:3) = y0;
    y_tmp(:,1) = z_tmp;
    
    %normalizing
    y1 = y_tmp;
    y2 = y_tmp;
    for j = 1:3
        y1(:,j) = y1(:,j) - min(y1(:,j));
        y2(:,j) = y1(:,j)/max(y1(:,j));
    end
    y_total_3{end+1} = y2;
end

%% save data

t = tspan.';
y_total = y_total_1;
filename = 'IFL_timeseries';
%save(filename, 'y_total', 't', 'time_interval','num_data')

t = tspan.';
y_total = y_total_2;
filename = 'CFL_timeseries';
%save(filename, 'y_total', 't', 'time_interval','num_data')

y_total = y_total_3;
filename = 'SFL_timeseries';
save(filename, 'y_total', 't', 'time_interval','num_data')
