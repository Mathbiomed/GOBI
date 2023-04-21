 clc;
clear;
close all;
%% law data
data = readtable('data.xlsx');
data = table2array(data);

t = data(:,1);
t1 = t;
y = data(:,2:end);

t = t - t(1);
data = data(1:end,2:end);

% for i = 1:length(data(1,:))
%     data(:,i) = data(:,i) - min(data(:,i));
%     data(:,i) = data(:,i) / max(data(:,i));
% end

figure(1)
plot(t,data(:,1),'k')
hold on
plot(t,data(:,2),'r')
hold on
plot(t,data(:,3),'g')
hold on
plot(t,data(:,4),'b')

xlim([0, t(end)])
ylim([0,1])
xticks([0,t(end)])
yticks([0,1])
%% fit data
t_fit = linspace(0, 2.45, 2.45*100+1).';

y_tmp = data;
y_fit = zeros(length(t_fit),4);
for j = 1:4
    y_spline = interp1(t,y_tmp(:,j),t_fit,'spline');
    y_fit(:,j) = y_spline;
end

figure(3)
plot(t_fit, y_fit(:,1),'k')
hold on
plot(t_fit, y_fit(:,2),'r')
hold on
plot(t_fit, y_fit(:,3),'g')
hold on
plot(t_fit, y_fit(:,4),'b')

%% save data
y = y_fit;
t = t_fit;
time_interval = 1/100;
save('data_estradiol','t','y','time_interval')