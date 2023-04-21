clc;
clear;
close all;

%% parameter
range = 1;
time_interval = 1/100;
tspan = linspace(0,range, range/time_interval+1);
num_data = 100;

%% parameter for source (Z)
time_interval_pre = 1/4; % fix 5 points
zt_pre = linspace(0,range, range/time_interval_pre+1);
zt = linspace(0,range, range/time_interval+1);

%% create time series
y_total = {};
for i = 1:num_data
    % create source (Z)
    %randomly sampled 5 points and connect them to generate Z 
    z_pre = rand([length(zt_pre),1])-0.5;
    z_tmp = spline(zt_pre,z_pre,zt);
    y_0 = rand([1,2])-0.5; % initial condition for X and Y
    
    % create other timeseries
    y_tmp = zeros(length(tspan),3);
    [t0, y0] = ode45(@(t,y) sample(t,y,zt,z_tmp), tspan, y_0);
    y_tmp(:,2:3) = y0;
    y_tmp(:,1) = z_tmp;
    
    %normalizing
    y1 = y_tmp;
    y2 = y_tmp;
    for j = 1:3
        y1(:,j) = y1(:,j) - min(y1(:,j));
        y2(:,j) = y1(:,j)/max(y1(:,j));
    end
    y_total{end+1} = y2;
end

%% figure: Multiple time-series
figure(1)
idx = 1;
t = tspan.';
op = 1;
cx = [236,0,140]/255;
cy = [0,174,239]/255;
cz = [88,89,91]/255;
cx = [cx, op];
cy = [cy, op];
cz = [cz, op];
for i = idx
    y_tmp = cell2mat(y_total(i));
    plot(t,y_tmp(:,1),'color',cz)
    hold on
    plot(t,y_tmp(:,2),'color',cx)
    hold on
    plot(t,y_tmp(:,3),'color',cy)
    hold on
end
xlim([0,range])
ylim([0,1])
xticks([0,range])
yticks([0,1])

%% save time-series
filename = 'sample_timeseries';
save(filename, 'y_total', 't', 'time_interval','num_data')