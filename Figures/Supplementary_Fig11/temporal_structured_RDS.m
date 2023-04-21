clc;
clear;
close all;

addpath('../GOBI') 

%% load data
load('data_temporal.mat')

%% parameter
size_window = 100;
move_window = 50;
num_window = (range / time_interval - size_window) / move_window + 1;
thres_noise = 1e-3;

S_total = zeros(num_data, num_window,2);
L_total = zeros(num_data, num_window,2);

%% compute RDS
for i = 1:num_data
     disp(i)
     % import time series
     y_tmp = cell2mat(y_total(i));
     
    % cut using moving window technique
    for j = 1:num_window
        % import single time series
        y_target_tmp = y_tmp(1 + move_window*(j-1):size_window + move_window*(j-1),:);
        y_target = y_target_tmp;
        t_target = tspan(1 + move_window*(j-1):size_window + move_window*(j-1));
        
        %normalize the time series
        for k = 1:3
            y_target(:,k) = y_target_tmp(:,k) - min(y_target_tmp(:,k));
            y_target(:,k) = y_target(:,k) / max(y_target(:,k));
        end
        
        % compute RDS for the regulation from X to Y
        [score_list, t_1, t_2] = RDS_ns_dim1(y_target(:,1), y_target(:,2), t_target, time_interval);
        for k = 1:2
            score = reshape(score_list(:,:,k),[length(t_target),length(t_target)]);
            loca_plus = find(score > thres_noise);
            loca_minus = find(score < -thres_noise);
            if isempty(loca_plus) && isempty(loca_minus)
                s = 1;    
            else
                s = (sum(score(loca_plus)) + sum(score(loca_minus)))/ (abs(sum(score(loca_plus))) + abs(sum(score(loca_minus))));
            end
            l = (length(loca_minus) + length(loca_plus)) / (length(t_1)*length(t_2)/2);
            if l < 0.1
                s = 1;
            end
            S_total(i,j,k) = s;
            L_total(i,j,k) = l;
        end        
    end
end

%% plot
figure(1) % figure of parameter
plot(t_pre, alpha_pre, 'r', 'LineWidth',2)
hold on
plot(t_pre, beta_pre, 'b', 'LineWidth',2)
hold on
plot(t_pre, gamma_pre, 'k', 'LineWidth',2)

xlim([0,40])
ylim([-1.2,1.2])
xticks([0:10:40])
yticks([-1,0,1])
xticklabels([])
yticklabels([])

figure(2)
idx = 1;
y_ori = cell2mat(y_total(i));
plot(tspan, y_ori(:,1), 'k')
hold on
plot(tspan, y_ori(:,2), 'r')
hold on
plot(tspan, y_ori(:,3), 'b')
hold on

figure(3)
S_pos = reshape(S_total(:,:,1),[num_data, num_window]);
S_neg = reshape(S_total(:,:,2),[num_data, num_window]);
window_list = [1:num_window];

errorbar(window_list, mean(S_pos), std(S_pos), 'LineWidth',2)
hold on
errorbar(window_list, mean(S_neg), std(S_neg), 'LineWidth',2)

xlim([0,40])
ylim([-1.2,1.2])
xticks([0:10:40])
yticks([-1,0,1])
xticklabels([])
yticklabels([])
