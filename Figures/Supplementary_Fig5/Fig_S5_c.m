clc;
clear;
close all;
addpath('../GOBI') 

%% parameter
system_list = {'KF','Fr','GW','rep','Gb','cAMP','SFL'};
num_system = length(system_list);
noise_list = [0:2:20];
acceptance = 0;
trial_list = [0:1];
% choose the threshold for regulation-detection score
% (i)   thres_S = 0.7
% (ii)  thres_S = 0.9
% (iii) thres_S = 0.99

thres_S = 0.9; 
tar_noise = [10];

c = {'k','r','b'};
%% threshold for S, L which maximize DB index

S_total_all = zeros(6,2,1000);
L_total_all = zeros(6,2,1000);

% i = 1 represents Kim-Forger model 
for i = 1
    for j = 1:length(tar_noise)
   
        for trial = trial_list
            % load data
            system_name = string(system_list(i));
            filename = append('./',system_name,'/',system_name,'_results_dim1_',num2str(noise_list(tar_noise(j))),'_Trial',num2str(trial));
            load(filename)
            S_total_all(:,:,trial*num_data+1:(trial+1)*num_data) = S_total;
            L_total_all(:,:,trial*num_data+1:(trial+1)*num_data) = L_total;
        end
    end
end
filename = append('true_network_',system_name);
load(filename)

for i = 1
    for j = 1:length(tar_noise)
   
        true_TRS = [];
        false_TRS = [];
        for k = 1:num_pair
            for l = 1:num_type
                S_tmp = reshape(S_total_all(k,l,:),[1000,1]);
                L_tmp = reshape(L_total_all(k,l,:),[1000,1]);

                L_processed = L_threshold(L_tmp, 0.1);
                S_processed = S_threshold(S_tmp, thres_S);

                S_processed = S_processed .* L_processed;
                TRS_tmp = sum(S_processed) / sum(L_processed);

                st = component_list_dim1(k,1);
                ed = component_list_dim1(k,2);
                if true_network(st,ed) == 1
                    if l == 1
                        true_TRS = [true_TRS ; TRS_tmp];
                    else
                        false_TRS = [false_TRS ; TRS_tmp];
                    end
                elseif true_network(st,ed) == -1
                    if l == 1
                        false_TRS = [false_TRS ; TRS_tmp];
                    else
                        true_TRS = [true_TRS ; TRS_tmp];
                    end
                else
                    false_TRS = [false_TRS ; TRS_tmp];
                end
            end
        end
    end 
end

true_TRS = flip(sort(true_TRS));
false_TRS = flip(sort(false_TRS));

figure(1)

scatter([1:length(true_TRS)], true_TRS,100, 'filled','b')
hold on
scatter([length(true_TRS)+1:length(true_TRS)+length(false_TRS)], false_TRS,100, 'filled','r')
xlim([0,13])
ylim([0,1])
xticks([1:12])
yticks([0,1])
xticklabels([]);
yticklabels([]);

