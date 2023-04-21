clc;
clear;
close all;

%% parameters
%trial = 1;
noise_list = [2:2:30];
noise_list = [0];
%% Load timeseries data
filename = 'KF_timeseries';
load(filename)

%% give multiplicative noise
for noise_percent = noise_list
    disp(noise_percent)
    y_total_noise = {};
    for i = 1:length(y_total)
        y_tmp = cell2mat(y_total(i));
        noise = normrnd(0, noise_percent/100,[length(t),1]);
        y_tmp_noise = y_tmp.*(noise+1);
        y_total_noise{end+1} = y_tmp_noise;
    end
    filename = ['KF_timeseries_noise_',num2str(noise_percent)];
    save(filename, 'y_total_noise', 't', 'time_interval','noise_percent','num_component','num_data','period')
end