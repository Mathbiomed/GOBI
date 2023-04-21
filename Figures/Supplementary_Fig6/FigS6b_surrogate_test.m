clc;
clear;
close all;
addpath('../GOBI')
%% parameter
noise_list = [10,15,20];
num_boot = 500;
thres_noise = 0;

%% bootstrapping
%noise_list = 10;
for noise_percent = noise_list
    disp(noise_percent)
    % load time-series
    filename = ['CFL_timeseries_fit_',num2str(noise_percent)];
    load(filename)
    
    % paramters
    num_data = length(y_total);
    S_boot = zeros(num_data, num_boot);
    L_boot = zeros(num_data, num_boot);
    
    % calculate score with shuffled time-series
    for i = 1:num_data
        % import time-series
        y_tmp = cell2mat(y_total(i));
        y_fix = y_tmp(:,2);
        y_shu = y_tmp(:,1);
        y_tar = y_tmp(:,3);
        
        for j = 1:num_boot
            % shuffle
            y_shu_tmp = y_shu(randperm(length(y_shu)));
            [score_1, score_2, t_1, t_2] = ion_prime_dim2(y_fix, y_shu_tmp, y_tar, t, time_interval);
            
            % score of type 2
            loca_plus = find(score_1 > thres_noise);
            loca_minus = find(score_1 < -thres_noise);
            if isempty(loca_plus) && isempty(loca_minus)
                s1 = 1;
            else
                s1 = (sum(score_1(loca_plus)) + sum(score_1(loca_minus)))/ (abs(sum(score_1(loca_plus))) + abs(sum(score_1(loca_minus))));
            end
            l1 = (length(loca_minus) + length(loca_plus)) / (length(t_1)*length(t_2)/2);
            
            % save at S_boot and L_boot
            S_boot(i,j) = s1;
            L_boot(i,j) = l1;
        end
    end
    filename = ['CFL_bootstrapping_noise_',num2str(noise_percent)];
    save(filename, 'S_boot', 'L_boot')
end


for noise_percent = noise_list
    disp(noise_percent)
    % load time-series
    filename = ['SFL_timeseries_fit_',num2str(noise_percent)];
    load(filename)
    
    % paramters
    num_data = length(y_total);
    S_boot = zeros(num_data, num_boot);
    L_boot = zeros(num_data, num_boot);
    
    % calculate score with shuffled time-series
    for i = 1:num_data
        % import time-series
        y_tmp = cell2mat(y_total(i));
        y_fix = y_tmp(:,2);
        y_shu = y_tmp(:,1);
        y_tar = y_tmp(:,3);
        
        for j = 1:num_boot
            % shuffle
            y_shu_tmp = y_shu(randperm(length(y_shu)));
            [score_1, score_2, t_1, t_2] = ion_prime_dim2(y_fix, y_shu_tmp, y_tar, t, time_interval);
            
            % score of type 2
            loca_plus = find(score_1 > thres_noise);
            loca_minus = find(score_1 < -thres_noise);
            if isempty(loca_plus) && isempty(loca_minus)
                s1 = 1;
            else
                s1 = (sum(score_1(loca_plus)) + sum(score_1(loca_minus)))/ (abs(sum(score_1(loca_plus))) + abs(sum(score_1(loca_minus))));
            end
            l1 = (length(loca_minus) + length(loca_plus)) / (length(t_1)*length(t_2)/2);
            
            % save at S_boot and L_boot
            S_boot(i,j) = s1;
            L_boot(i,j) = l1;
        end
    end
    filename = ['SFL_bootstrapping_noise_',num2str(noise_percent)];
    save(filename, 'S_boot', 'L_boot')
end
