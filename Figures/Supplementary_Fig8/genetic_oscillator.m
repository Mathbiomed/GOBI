clc;
clear;
close all;

num_component = 2;

%% load data
load('data_genetic_oscillator.mat')
%num_data = length(data_total(1,1,:));
num_data = 8;
num_component = 2;
thres_noise = 1e-5;

%% for each data, calculate period
period_list = [];
data_spline = zeros(101,2,num_data);
for i = 1:num_data
    y_tmp = reshape(data_total(:,:,i),[num_component, length(time)]);
    y_tmp = y_tmp .';
    
    t_int = linspace(0,1,101).';
    y_int = zeros(length(t_int),num_component);
    
    for j = 1:2
        y_int(:,j) = spline(time, y_tmp(:,j), t_int);
        y_int(:,j) = y_int(:,j) - min(y_int(:,j));
        y_int(:,j) = y_int(:,j) / max(y_int(:,j));
    end 
    data_spline(:,:,i) = y_int;
    
    fs = 100;
    [autocor1,lags1] = xcorr(y_int(:,1), fs,'coeff');
    [pks1,locs1] = findpeaks(autocor1);
    [autocor2,lags2] = xcorr(y_int(:,2), fs,'coeff');
    [pks2,locs2] = findpeaks(autocor2);
    
    if ~isempty(find(locs1 > 101))
        tmp1 = locs1(find(locs1 > 101));
    else
        tmp1 = [125];
    end
    
    if ~isempty(find(locs2 > 101))
        tmp2 = locs2(find(locs2 > 101));
    else
        tmp2 = [125];
    end
    
    per1 = (tmp1(1) - fs - 1) / fs;
    per2 = (tmp2(1) - fs - 1) / fs;
    %period_list = [period_list ; [per1,per2]];
        
%     figure(i)
%     plot(lags/fs,autocor)
%     xlabel('Lag (days)')
%     ylabel('Autocorrelation')
%     axis([-1 1 -0.4 1.1])
     
%     figure(i)
%     plot(t_int, y_int(:,1), 'r')
%     hold on
%     plot(t_int, y_int(:,2), 'b')
    
end

period_list = [0.5;0.3;0.15;0.15;0.15;1;1;1];

%% cut data using period
y_total = {};
t_total = {};
time_interval = 1/100;
for i = 1:8
    window_size_ori = period_list(i);
    window_move = window_size_ori;
    window_size = round(window_size_ori / time_interval);
    window_move = round(window_move / time_interval);
    
    y = reshape(data_spline(:,:,i),[101,2]);
    t_tmp = t_int(1:window_size);
    %% cut data
    start = 1;
    length_timeseries = length(y(:,1));
    while(1)
        if start + window_size > length_timeseries
            break
        end
        y_tmp = y(start:start + window_size - 1,:);
        y_total{end+1} = y_tmp; 
        t_total{end+1} = t_tmp;
        start = start + window_move;
    end
end
length_data = length(y_total);
t = time;

%% Fourier fitting
y_fit_total = {};

for i = 1:length_data
    y_tmp = cell2mat(y_total(i));
    t_fit = [0:time_interval:length(y_tmp(:,1))*time_interval];
    t_fit = t_fit(1:end-1).';
    y_fit = zeros(length(t_fit),num_component);
    for j = 1:num_component
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