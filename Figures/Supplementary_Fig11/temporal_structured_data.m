clc;
clear;
close all;

addpath('../GOBI') 

%% parameter
range = 100; % range of time-series
time_interval = 1/20; % 1/time_interval represents sampling rate
range_interval = 25;
num_data = 100;

%% variables: time, parameters
% time
tspan = linspace(0,range, range/time_interval+1);
t_pre = tspan;

% parameter
range_length = range_interval / time_interval;
increase_parameter = [0:1/range_length:1];
decrease_parameter = 1 - increase_parameter;

alpha_pre = [ones(1,range_length), decrease_parameter, zeros(1,range_length*2)];
beta_pre  = [zeros(1,range_length), increase_parameter, ones(1,range_length), zeros(1,range_length)];
gamma_pre = [zeros(1,range_length*3), ones(1,range_length+1)];

% time for source
time_interval_pre = 1/2;
xt_pre = linspace(0,range, range/time_interval_pre+1);
zt_pre = xt_pre;

%% create time series
y_total = {};
for i = 1:num_data
    % create source (X, Z)
    x_pre = rand([length(xt_pre),1]);
    z_pre = rand([length(zt_pre),1]);
    x_source = spline(xt_pre,x_pre,tspan);
    z_source = spline(zt_pre,z_pre,tspan);
    x_source = x_source - min(x_source);
    z_source = z_source - min(z_source);
    x_source = x_source / max(x_source);
    z_source = z_source / max(z_source);
    
    %% initial condition
    y_initial = rand(1); % initial condition for Y
    
    %% simulate ODE
    y_ori = zeros(length(tspan),3);
    y_ori(:,1) = x_source;
    y_ori(:,3) = z_source;

    [t1, y1] = ode45(@(t,y) temporal_total(t,y,t_pre, alpha_pre, beta_pre, gamma_pre, x_source, z_source), tspan, y_initial);
    y_ori(:,2) = y1;
    if ~isempty(find(y1 == 0))
        i = i-1;
        continue
    end
    % save time series at the variable 'y_total'
    y_total{end+1} = y_ori;
end

%% save data
filename = 'data_temporal';
save(filename, 'y_total', 'tspan', 'time_interval','range', 'range_interval', 'num_data', 'alpha_pre', 'beta_pre', 'gamma_pre', 't_pre')
