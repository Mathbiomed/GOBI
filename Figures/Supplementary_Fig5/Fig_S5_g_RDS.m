clc;
clear;
close all;
addpath('../GOBI') 
%% parameter
system_list = {'KF','Fr','GW','Gb','cAMP'};
num_system = length(system_list);
dimension = 1;
noise_list = [0:2:20];
trial_list = [0:9];
thres_L = 0.1;

noise_list = [0:10:20];

%% merge score of true and false regulations
score_total = [];
group_total = [];
for i = 1:5
    % load true network
    system_name = string(system_list(i));
    filename = append('true_network_',system_name);
    load(filename)
    
    % compute true index
    num_component = length(true_network(:,1));
    component = [1:num_component];
    component_list_dim1_tmp = nchoosek(component, 2);
    component_list_dim1 = [];
    for j = 1:length(component_list_dim1_tmp(:,1))
        component_list_dim1 = [component_list_dim1 ; [component_list_dim1_tmp(j,1), component_list_dim1_tmp(j,2)]];
        component_list_dim1 = [component_list_dim1 ; [component_list_dim1_tmp(j,2), component_list_dim1_tmp(j,1)]];
    end
    num_pair = length(component_list_dim1(:,1));
    num_type = 2^dimension;
    idx_true = zeros(num_pair, num_type);
    for j = 1:num_pair
        st = component_list_dim1(j,1);
        ed = component_list_dim1(j,2);
        if true_network(st,ed) == 1
            idx_true(j,1) = 1;
        elseif true_network(st,ed) == -1
            idx_true(j,2) = 1;
        end
    end
    
    % merge
    for noise_percent = noise_list
        score_true = [];
        score_false = [];
        for trial = trial_list
            % load data
            filename = append('./',system_name,'/',system_name,'_results_dim1_',num2str(noise_percent),'_Trial',num2str(trial));
            load(filename)
            
            % for each pair and type, merge
            for j = 1:num_pair
                for k = 1:num_type
                    S_tmp = reshape(S_total(j,k,:), [num_data,1]);
                    L_tmp = reshape(L_total(j,k,:), [num_data,1]);
                    
                    L_processed = L_threshold(L_tmp, thres_L);
                    S_processed = S_tmp .* L_processed;
                    S_processed(S_processed == 0) = nan;
                    if idx_true(j,k) == 1
                        score_true = [score_true ; S_processed];
                    else
                        score_false = [score_false ; S_processed];
                    end
                end
            end
        end
        
        group_idx_true = round(noise_percent/2*18) + i;
        score_total = [score_total ; score_true];
        group_total = [group_total ; group_idx_true * ones(length(score_true),1)];
        
        group_idx_false = round(noise_percent/2*18) + i + 6;
        score_total = [score_total ; score_false];
        group_total = [group_total ; group_idx_false * ones(length(score_false),1)];
        
        for j = 13:18
            group_idx_tmp = round(noise_percent/2*18) + j;
            score_total = [score_total ; 0];
            group_total = [group_total ; group_idx_tmp];
        end

    end
end

figure(1)
boxplot(score_total , group_total,'Colors','rgbkcm', 'Width', 0.9, 'OutlierSize', 2,'Symbol', '')
set(gca,'XTick',[6.5:18:186.5],'XTickLabel',{' '})
hold on

plot([6.5,42.5],[0.95,0.85],'r')
xlim([-1,49])
ylim([-1,1])
yticks([-1,-0.5])
yticklabels([])