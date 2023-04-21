clc;
clear;
close all;

trial_list = [0:9];
system = 'cAMP';
%system = 'Fr';
%system = 'Gb';
%system = 'GW';
%system = 'KF';

for trial = trial_list
    %% parameters
    noise_list = [2:2:20];

    %% Load timeseries data
    filename = ['./Data_original/',system,'_timeseries_Trial',num2str(trial)];
    load(filename)

    %% give multiplicative noise
    for noise_percent = noise_list
        disp(noise_percent)
        y_total_noise = {};
        for i = 1:length(y_total)
            y_tmp = cell2mat(y_total(i));
            noise = normrnd(0, noise_percent/100,[length(t),num_component]);
            y_tmp_noise = y_tmp.*(noise+1);
            y_total_noise{end+1} = y_tmp_noise;
        end
        filename = ['./Data_multiplicative/',system,'_timeseries_multiplicative_',num2str(noise_percent),'_Trial',num2str(trial)];
        save(filename, 'y_total_noise', 't', 'time_interval','noise_percent','num_component','num_data','period')
    end
end