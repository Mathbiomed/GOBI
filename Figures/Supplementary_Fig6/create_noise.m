clc;
clear;
close all;

%% Load timeseries data
load('SFL_timeseries.mat')

%% give noise
noise_list = [5:5:20];
for noise_percent = noise_list
    disp(noise_percent)
    y_total_noise = {};
    for i = 1:length(y_total)
        y_tmp = cell2mat(y_total(i));
        y_tmp_noise = y_tmp;
        for j = 1:3
            noise = normrnd(0, noise_percent/100,[length(t),1]);
            y_tmp_noise(:,j) = y_tmp(:,j) .* (1+noise);
            %y_tmp_noise(:,j) = y_tmp(:,j) + noise;
        end
        y_total_noise{end+1} = y_tmp_noise;
    end
    filename = ['SFL_timeseries_noise_',num2str(noise_percent)];
    save(filename, 'y_total_noise', 't', 'time_interval','noise_percent')
end


%% Load timeseries data
load('CFL_timeseries.mat')

%% give noise
noise_list = [5:5:20];
for noise_percent = noise_list
    disp(noise_percent)
    y_total_noise = {};
    for i = 1:length(y_total)
        y_tmp = cell2mat(y_total(i));
        y_tmp_noise = y_tmp;
        for j = 1:3
            noise = normrnd(0, noise_percent/100,[length(t),1]);
            y_tmp_noise(:,j) = y_tmp(:,j) .* (1+noise);
            %y_tmp_noise(:,j) = y_tmp(:,j) + noise;
        end
        y_total_noise{end+1} = y_tmp_noise;
    end
    filename = ['CFL_timeseries_noise_',num2str(noise_percent)];
    save(filename, 'y_total_noise', 't', 'time_interval','noise_percent')
end