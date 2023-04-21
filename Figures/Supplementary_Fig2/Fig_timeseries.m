clc;
clear;
close all;

load('sample_timeseries')

op = 1;

cx = [236,0,140]/255;
cy = [0,174,239]/255;
cz = [88,89,91]/255;
cx = [cx, op];
cy = [cy, op];
cz = [cz, op];

figure(1)
for i = 1
    y_tmp = cell2mat(y_total(i));
    plot(t,y_tmp(:,1),'k')
    hold on
    plot(t,y_tmp(:,2),'b')
    hold on
    plot(t,y_tmp(:,3),'r')
    hold on
    plot(t,y_tmp(:,4),'g')
    hold on
end

xlim([0,pi*2])
ylim([0,1])
xticks([0,pi*2])
yticks([0,1])