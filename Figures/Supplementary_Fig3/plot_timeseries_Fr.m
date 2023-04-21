clc;
clear;
close all;

%% load data
load('timeseries_1')
y1 = y_total;

load('timeseries_5')
y2 = y_total;

load('timeseries_10')
y3 = y_total;

%% plot time series
figure(1)
for i = 1:1
    y_tmp = cell2mat(y1(i));
    
    plot(t, y_tmp(:,1), 'k')
    hold on
    plot(t, y_tmp(:,2), 'r')
    hold on
    plot(t, y_tmp(:,3), 'b')
    hold on
end

xlim([0,10])
xticks([0,10])
ylim([0,1])
yticks([0,1])

figure(2)
for i = 1:5
    y_tmp = cell2mat(y2(i));
    
    plot(t, y_tmp(:,1), 'k')
    hold on
    plot(t, y_tmp(:,2), 'r')
    hold on
    plot(t, y_tmp(:,3), 'b')
    hold on
end

xlim([0,10])
xticks([0,10])
ylim([0,1])
yticks([0,1])

figure(3)
for i = 1:10
    y_tmp = cell2mat(y3(i));
    
    plot(t, y_tmp(:,1), 'k')
    hold on
    plot(t, y_tmp(:,2), 'r')
    hold on
    plot(t, y_tmp(:,3), 'b')
    hold on
end

xlim([0,10])
xticks([0,10])
ylim([0,1])
yticks([0,1])