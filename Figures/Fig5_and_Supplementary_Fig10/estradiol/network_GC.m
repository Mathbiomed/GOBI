clc;
clear;
close all;

%% load data
filename = ['data_estradiol'];
load(filename)

%% test GC
% parameters
num_component = 4;
sig_level = 0.05;
max_lag = 25;

% variables
cause_list = zeros(num_component);
F_list = zeros(num_component);
critic_list = zeros(num_component);

%test
for i = 1:num_component
    for j = 1:num_component
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
