clc;
clear;
close all;

%% load raw data
load('raw_data_repressilator')

%% parameter
window_size_ori = 40;
window_move_ori = 4;
thres_noise = 0;
num_component = 3;
time_interval = 1/4;

%% plot raw data
figure(1)
plot(t, y(:,1), 'bo-','Linewidth',3,'markersize',7,'markerfacecolor','b')
hold on
plot(t, y(:,2), 'ro-','Linewidth',3,'markersize',7,'markerfacecolor','r')
hold on
plot(t, y(:,3), 'go-','Linewidth',3,'markersize',7,'markerfacecolor','r')

xlim([0,420])
xticks([0,420])
set(gca,'fontsize',16)

%% interpolate data
t_fit = linspace(0,t(end),length(t)/time_interval+1).';
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
hold on
plot(t, y_tmp(:,3), 'g','Linewidth',1)
%% Save data

save('data_final_repressilator','t','y_total','time_interval','num_data','num_component')