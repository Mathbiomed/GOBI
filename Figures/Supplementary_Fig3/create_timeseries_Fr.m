clc;
clear;
close all;

%% parameters
time_interval = 1/10;
period = 10;
num_component = 3;
num_data = 1;

%% create time-series with various initial value
% create initials
initials  =  rand([num_data,num_component]);

% simulate time-series
y_total = {};
tspan = linspace(0,period, period/time_interval+1);
for i = 1:length(initials(:,1))
    [t1, y1] = ode15s(@(t,y) Frzilator(t,y), tspan, initials(i,:));
    for j = 1:num_component
        y1(:,j) = y1(:,j) - min(y1(:,j));
        y1(:,j) = y1(:,j) / max(y1(:,j));
    end
    y_total{end+1} = y1;
end
t = t1;

figure(1) % various time-series
for i = 1:1
    y_tmp = cell2mat(y_total(i));
    plot(t/t(end), y_tmp(:,1), 'k')
    hold on
    plot(t/t(end), y_tmp(:,2), 'r')
    hold on
    plot(t/t(end), y_tmp(:,3), 'b')
    hold on
end
xlim([0,1])
xticks([0,1])
xticklabels({'0', 'period'})
yticklabels([])
xlabel('Time')
ylabel('Value of component')
set(gca,'fontsize',16)

% save data
filename = ['timeseries_',num2str(num_data)];
save(filename, 'y_total', 't', 'time_interval', 'period', 'num_component', 'num_data')

