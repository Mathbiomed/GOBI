clc;
clear;
close all;

%% load raw data
load('data')

%% parameters 1. interpolation

% interpolation method
% method = 1: linear interpolation
% method = 2: spline interpolation
% method = 3: fourier interpolation, 
%             In this case, users have to choose the order of fourier method (1~8)

method = 2;
num_fourier = 6;

%% parameters 2. cut the time series data
num_component = length(y(1,:));        % number of component in this system
window_size_ori = 40;     % For oscillatory data, 1 period is recommended
overlapping_ratio = 0.9;  % overlapping ratio of moving window technique

% choose sampling rate for interpolation. 
% window_size_ori/time_interval becomes number of time points per cutted data
% window_size_ori/time_interval = 100 is recommended
% window_size_ori/time_interval is high (low) make the inference accurate (less accurate) and slow (fast)
time_interval = 1/4;       

%% process parameters
window_size = window_size_ori / time_interval;
window_move_ori = ceil(window_size_ori * (1-overlapping_ratio));
window_move = window_move_ori / time_interval;


%% plot raw data
figure(1)
for i = 1:num_component
    plot(t, y(:,i))
    hold on
end
xlim([0,t(end)])
xticks([])
ylim([0,1])
yticks([])

%% interpolate data
t_fit = linspace(0,t(end),length(t)/time_interval+1).';
y_fit = zeros(length(t_fit),num_component);

if method == 1 % linear
    for i = 1:num_component
        y_fit(:,i) = interp1(t, y(:,i),t_fit,'linear');
    end
elseif method == 2 % spline
    for i = 1:num_component
        y_fit(:,i) = interp1(t, y(:,i),t_fit,'spline');
    end
elseif method == 3 % fourier
    for i = 1:num_component
        y_fit(:,i) = interp1(t, y(:,i),t_fit,'linear');
    end
    option = ['fourier',num2str(num_fourier)];
    for i = 1:num_component
        fitting = fit(t_fit,y_fit(:,i),option);
        w = fitting.w;
        fouriers = ones(1,length(t_fit));
        for j = 1:num_fourier
            fouriers = [fouriers; cos(j*w*t_fit.')];
            fouriers = [fouriers; sin(j*w*t_fit.')];
        end
        coeffs = coeffvalues(fitting);
        y_fit(:,i) = coeffs(1:end-1) * fouriers;
    end
end

%% plot interpolated data
figure(2)
for i = 1:num_component
    plot(t_fit, y_fit(:,i))
    hold on
end
xlim([0,t_fit(end)])

xlim([0,t_fit(end)])
xticks([])
ylim([0,1.02])
yticks([])
%% cut and interpolate the data
y_total = {};
start = 1; % it can be changed if users do not want to use the beginning of time series
length_timeseries = length(y_fit(:,1)); % it can be changed if users do not want to use the end of time series
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

%% plot cutted data
figure(3)
for j = 1:num_data
    y_tmp = cell2mat(y_total(j));
    for i = 1:num_component
        plot(t, y_tmp(:,i))
        hold on
    end
end
xlim([0,t(end)])

xlim([0,t(end)])
xticks([])
ylim([0,1.02])
yticks([])
%% Save data
time_interval = 1/window_size;
save('data_cut','t','y','y_total','time_interval','num_data','num_component')