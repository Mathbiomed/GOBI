clc;
clear;
close all;
addpath('../../GOBI')
%% parameter
noise_list = [0:2:20];
thres_S_list = [0:0.01:1];
thres_L_list = [0.01:0.01:0.1];
thres_TRS_list = [0:0.01:1];
epsilon = 1e-7;
true_network = [
    [0, 1,0,0];
    [0, 0,1,0];
    [0, 0,0,1];
    [-1,0,0,0]];
filename = ['GW_results_dim1_0_Trial0'];
load(filename)
trial_list = [0:9];

%% create true index
true_index = zeros(num_pair, num_type);
for i = 1:num_pair
    st = component_list_dim1(i,1);
    ed = component_list_dim1(i,2);
    if true_network(st,ed) == 1
        true_index(i,1) = 1;
    elseif true_network(st,ed) == -1
        true_index(i,2) = 1;
    end
end

%% calculate F1
F1_list = zeros(length(thres_S_list),length(thres_L_list),length(thres_TRS_list),length(noise_list));

for trial = trial_list
for noise_idx = 1:length(noise_list)
    disp(noise_idx)
    % load data
    filename = ['GW_results_dim1_',num2str(noise_list(noise_idx)),'_Trial',num2str(trial)];
    load(filename)
    
    for thres_S_idx = 1:length(thres_S_list)
        for thres_L_idx = 1:length(thres_L_list)
            % calculate TRS
            for i = 1:num_pair
                for j = 1:num_type
                    S_tmp = reshape(S_total(i,j,:), [num_data,1]);
                    L_tmp = reshape(L_total(i,j,:), [num_data,1]);
                    
                    S_processed = S_threshold(S_tmp, thres_S_list(thres_S_idx));
                    L_processed = L_threshold(L_tmp, thres_L_list(thres_L_idx));
                    S_processed = S_processed .* L_processed;
                    if sum(L_processed) ~= 0
                        TRS(i,j) = sum(S_processed) / sum(L_processed);
                    end
                end
            end
            % calculate F1 for each threshold
            for thres_TRS_idx = 1:length(thres_TRS_list)
                TP = 0;
                TN = 0;
                FP = 0;
                FN = 0;
                for i = 1:num_pair
                    for j = 1:num_type
                        if true_index(i,j) == 1
                            if TRS(i,j) >= thres_TRS_list(thres_TRS_idx)
                                TP = TP + 1;
                            else
                                FN = FN + 1;
                            end
                        else
                            if TRS(i,j) >= thres_TRS_list(thres_TRS_idx)
                                FP = FP + 1;
                            else
                                TN = TN + 1;
                            end
                        end
                    end
                end
                Precision = TP / (TP + FP + epsilon);
                Recall = TP / (TP + FN + epsilon);
                F2 = 5/(1/(Precision+epsilon) + 4/(Recall+epsilon));
                F2_list(thres_S_idx,thres_L_idx,thres_TRS_idx,noise_idx) = F2;
            end
            
            
        end
    end
end

filename = ['F2_total_GW_Trial',num2str(trial)];
save(filename, 'F2_list', 'thres_S_list', 'thres_L_list','thres_TRS_list')
end