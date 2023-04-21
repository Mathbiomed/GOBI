clc;
clear;
close all;

%% read raw data
data = readtable('data.xlsx');
data = table2array(data);

t = data(:,1);
t1 = t;
y = data(:,2:end);

t = t - t(1);
t_ori = t;

data = data(1:end,2:end);

for i = 1:length(data(1,:))
    data(:,i) = data(:,i) - min(data(:,i));
    data(:,i) = data(:,i) / max(data(:,i));
end

%% plot time series
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


%% parameter
window_size_ori = 0.9;
window_move_ori = 0.09;
thres_noise = 0;
num_component = 4;
time_interval = 1/500;

%% interpolate data
t_fit = linspace(0,2.45,2.45*500+1).';
y_fit = zeros(length(t_fit),num_component);
for i = 1:num_component
    y_fit(:,i) = interp1(t, y(:,i),t_fit,'spline');
end

%% cut and interpolate the data
y_total = {};
window_size = window_size_ori / time_interval;
window_move = window_move_ori / time_interval;

start = 1;
length_timeseries = length(y_fit(:,1));
while(1)
    if start + window_size-1 > length_timeseries
        break
    end
    
    % cut
    y_tmp = y_fit(start:start + window_size - 1,:);
    t = t_fit(1:window_size) / t_fit(window_size);
   
    % normalize
    for i = 1:num_component
        y_tmp(:,i) = y_tmp(:,i) - min(y_tmp(:,i));
        y_tmp(:,i) = y_tmp(:,i) / max(y_tmp(:,i));
    end
    y_total{end+1} = y_tmp; 
    start = start + window_move;
end
num_data = length(y_total);

%% plot interpolated data
idx = 1;
y_tmp = cell2mat(y_total(idx));

figure(2)
plot(t, y_tmp(:,1), 'b','Linewidth',1)
hold on
plot(t, y_tmp(:,2), 'r','Linewidth',1)

%% Save data

save('data_final_estradiol','t','y_total','time_interval','num_data','num_component')