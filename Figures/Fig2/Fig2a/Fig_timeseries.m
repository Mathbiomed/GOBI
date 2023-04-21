clc;
clear;
close all;

%% load time-series
load('sample_timeseries')

%% plot time-series
op = 1;
cx = [236,0,140]/255;
cy = [0,174,239]/255;
cz = [88,89,91]/255;
cx = [cx, op];
cy = [cy, op];
cz = [cz, op];

figure(1)
idx = 1;  % change idx to plot various time-series
for i = idx
    y_tmp = cell2mat(y_total(i));
    plot(t,y_tmp(:,1),'color',cz)
    hold on
    plot(t,y_tmp(:,2),'color',cx)
    hold on
    plot(t,y_tmp(:,3),'color',cy)
    hold on
end

xlim([0,1])
ylim([0,1])
xticks([0,1])
yticks([0,1])