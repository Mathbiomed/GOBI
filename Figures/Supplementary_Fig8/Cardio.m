clc;
clear;
close all;

%% load data
load('data_cardio.mat')
t = linspace(0,1031,1032).';
y_raw = [cardio,no2,so2,o3,rspar];
num_component = 5;

%% moving average
win_avg = 7;

y_movavg = y_raw;
for i = 1:1
    y_tmp = y_raw(:,i);
    y_movavg(:,i) = movmean(y_tmp, win_avg);
end

%% interpolation
time_interval = 1/2;

t_smooth = linspace(0, 1031, 1031/time_interval+1).';
y_smooth = zeros(length(t_smooth),5);
y_tmp = y_movavg;
for j = 1:5
    y_fit_tmp = interp1(t,y_tmp(:,j),t_smooth);
    y_smooth(:,j) = y_fit_tmp;
end

%% cut data
start_index = 1;
end_index = 1032;

start_index = (start_index-1)/time_interval+1;
end_index = end_index/time_interval-1;

y_smooth = y_smooth(start_index:end_index,:);
t_smooth = t_smooth(start_index:end_index);

window_size_ori = 365;
window_move = 30;
thres_noise = 0;

window_size = window_size_ori / time_interval;
window_move = window_move / time_interval;

y_total = {};
t_total = t_smooth(start_index:start_index+window_size-1);

start = 1;
length_timeseries = length(y_smooth(:,1));

while(1)
    if start + window_size > length_timeseries
        break
    end
    y_tmp = y_smooth(start:start + window_size - 1,:);
    for i = 1:length(y_tmp(1,:))
        y_tmp(:,i) = y_tmp(:,i) - min(y_tmp(:,i));
        y_tmp(:,i) = y_tmp(:,i) / max(y_tmp(:,i)) - min(y_tmp(:,i));
    end
    
    y_total{end+1} = y_tmp; 
    start = start + window_move;
end
num_data = length(y_total);

%% fourier fitting
y_total_fit = {};
for i = 1:num_data
    y_tmp = cell2mat(y_total(i));
    t_tmp = t_total;
    y_fit = y_tmp;
    
    for j = 1:num_component
        fitting = fit(t_tmp,y_tmp(:,j),'fourier4');
        w = fitting.w;
        fouriers = [
            ones(1,length(t_tmp));
            cos(w*t_tmp.');
            sin(w*t_tmp.');
            cos(2*w*t_tmp.');
            sin(2*w*t_tmp.');
            cos(3*w*t_tmp.');
            sin(3*w*t_tmp.');
            cos(4*w*t_tmp.');
            sin(4*w*t_tmp.')];
        coeffs = coeffvalues(fitting);
        y_tmp_fit = coeffs(1:end-1) * fouriers;
        y_fit(:,j) = y_tmp_fit;
    end
    y_total_fit{end+1} = y_fit;
end

%% calculate residual
length_data = length(y_total);
residual_total = zeros(length_data, num_component);
for i = 1:length_data
     y_fit = cell2mat(y_total_fit(i));
     y_noise = cell2mat(y_total(i));
     for j = 1:num_component
        residual_tmp = mean((y_fit(:,j) - y_noise(:,j)).^2);
        residual_total(i,j) = residual_tmp;
    end
end
mean(mean(residual_total))

