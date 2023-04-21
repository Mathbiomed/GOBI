clc;
clear;
close all;

%% load data
load('data_prey_predator.mat')
num_component = length(y(1,:));

%% parameter
% general
window_size_ori = 5;
window_move = 5;
thres_noise = 0;

%% process parameter
window_size = round(window_size_ori / time_interval);
window_move = round(window_move / time_interval);

%% cut data
disp('cut_data...')
y_total = {};
start = 1;
length_timeseries = length(y(:,1));
while(1)
    if start + window_size > length_timeseries
        break
    end
    y_tmp = y(start:start + window_size - 1,:);
    for i = 1:num_component
        y_tmp(:,i) = y_tmp(:,i) - min(y_tmp(:,i));
        y_tmp(:,i) = y_tmp(:,i) / max(y_tmp(:,i));
    end
    y_total{end+1} = y_tmp; 
    start = start + window_move;
end
length_data = length(y_total);

%% Fourier fitting
y_fit_total = {};
t_fit = linspace(0, window_size_ori, round(window_size_ori/time_interval)+1).';
t_fit = t_fit(1:end-1);

for i = 1:length_data
    y_fit = zeros(length(t_fit),num_component);
    for j = 1:num_component
        y_tmp = cell2mat(y_total(i));
        y1 = fit(t_fit,y_tmp(:,j),'fourier4');
        w = y1.w;
        fouriers = [
            ones(1,length(t_fit));
            cos(w*t_fit.');
            sin(w*t_fit.');
            cos(2*w*t_fit.');
            sin(2*w*t_fit.');
            cos(3*w*t_fit.');
            sin(3*w*t_fit.');
            cos(4*w*t_fit.');
            sin(4*w*t_fit.')];
        coeffs = coeffvalues(y1);
        y2 = coeffs(1:end-1) * fouriers;
        y_fit(:,j) = y2;
    end
    y_fit_total{end+1} = y_fit; 
end

%% calculate residual
residual_total = zeros(length_data, num_component);
for i = 1:length_data
     y_fit = cell2mat(y_fit_total(i));
     y_noise = cell2mat(y_total(i));
     for j = 1:num_component
        residual_tmp = mean((y_fit(:,j) - y_noise(:,j)).^2);
        residual_total(i,j) = residual_tmp;
    end
end
mean(mean(residual_total))

