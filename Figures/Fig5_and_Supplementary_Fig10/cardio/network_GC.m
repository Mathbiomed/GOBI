clc;
clear;
close all;

%% load data
filename = ['data_cardio'];
load(filename)

y = y(32:end,:);
%y = y(1:720,:);
%% test GC
% parameters
num_component = 5;
sig_level = 0.05;
max_lag = 20;

% variables
cause_list = zeros(num_component,1);
F_list = zeros(num_component,1);
critic_list = zeros(num_component,1 );

%test
for i = 1:num_component
    for j = 1:1
        if i == j
            continue
        end

        [F_tmp, c_tmp] = granger_cause(y(:,j), y(:,i), sig_level, max_lag);
        F_list(i,j) = F_tmp;
        critic_list(i,j) = c_tmp;
        if F_tmp > c_tmp
            cause_list(i,j) = cause_list(i,j) + 1;
        end
    end
end

cause_list
