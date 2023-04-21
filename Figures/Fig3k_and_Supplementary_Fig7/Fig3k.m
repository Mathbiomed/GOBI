clc;
clear;
close all;
addpath('../GOBI') 

%% parameters
trial_list = [0:9];
noise_list = [2:2:20];
%noise_list = [2];

dimension = 2;
thres_noise = 1e-7;

system_list = {'KF', 'Fr', 'GW', 'Gb', 'cAMP'};
num_component_list = [3,3,4,5,7];

true_network = [
    [0, 0, 0, 0, 1,1,0];
    [-1,0,-1, 0, 0,0,0];
    [0, 0, 0,-1, 0,0,0];
    [0, 0, 0, 0,-1,0,0];
    [0, 1, 0, 0, 0,0,0];
    [0, 0, 0, 0, 0,0,1];
    [1, 0, 1, 0, 0,0,0]];

%F1_list = zeros(length(noise_list), length(trial_list));
F2_list_total = [];
for system_idx = 1:length(system_list)
    system = char(system_list(system_idx));
    num_component = num_component_list(system_idx);
    disp(system)
    F2_list = [];
    for trial = trial_list
        F2_list_tmp = [];
        for noise_percent = noise_list
            regulation_network_dim1 = zeros(num_component,num_component);
            regulation_network_dim2 = zeros(num_component,num_component);

            %% load results
            filename = ['./TRS_dim2_multiplicative','/',system,'_TRS_dim2_',num2str(noise_percent),'_Trial',num2str(trial)];
            load(filename)
            filename = ['./TRS_dim1_multiplicative','/',system,'_TRS_dim1_',num2str(noise_percent),'_Trial',num2str(trial)];
            load(filename)
            filename = ['./Delta_dim2_multiplicative','/',system,'_Delta_dim2_',num2str(noise_percent),'_Trial',num2str(trial)];
            load(filename)
            filename = ['./Surrogate_dim2_multiplicative','/',system,'_surrogate_dim2_',num2str(noise_percent),'_Trial',num2str(trial)];
            load(filename)

            %% 1D regulation  
            TRS_1D_result = zeros(num_component,1);
            for i = 1:num_component
                % indexing every cause of target
                index_list = [];
                for j = 1:length(component_list_dim1(:,1))
                    if component_list_dim1(j,2) == i
                        index_list = [index_list, j];
                    end
                end

                % chech inferred 1D regulation
                regulation_idx_list = [];
                TRS_value_list = [];
                for j = index_list
                    for k = 1:length(regulation_1dim(1,:))
                        if regulation_1dim(j,k) == 1
                            regulation_idx_list = [regulation_idx_list;[j,k]];
                            TRS_value_list = [TRS_value_list;TRS_dim1(j,k)];
                        end
                    end
                end

                if ~isempty(regulation_idx_list)
                    % choose maximum TRS if there are multiple regulation
                    [val,loc] = max(TRS_value_list);
                    TRS_1D_result(i) = val;
                    % save at regulation_list
                    pair_idx = regulation_idx_list(loc,1);
                    type_idx = regulation_idx_list(loc,2);

                    st = component_list_dim1(pair_idx,1);
                    ed = component_list_dim1(pair_idx,2);
                    if type_idx == 1
                        regulation_network_dim1(st,ed) = 1;
                    else
                        regulation_network_dim1(st,ed) = -1;
                    end
                end
            end
            %regulation_network_dim1

            %% 2D regulation       
            result_boot_2dim = zeros(length(component_list_dim2(:,1)),4);
            if isempty(candidate_2dim)
            else 
            for i = 1:length(candidate_2dim(:,1))
                result_boot_2dim(candidate_2dim(i,1),candidate_2dim(i,2)) = 1;
            end
            end
            TRS_2D_result = zeros(num_component,1);
            for i = 1:num_component
                % indexing every cause of target
                index_list = [];
                for j = 1:length(component_list_dim2(:,1))
                    if component_list_dim2(j,3) == i
                        index_list = [index_list, j];
                    end
                end

                % check inferred 2D regulation
                regulation_idx_list = [];
                TRS_value_list = [];
                for j = index_list
                    for k = 1:length(delta_2dim(1,:))
                        if result_boot_2dim(j,k) == 1
                            regulation_idx_list = [regulation_idx_list;[j,k]];
                            TRS_value_list = [TRS_value_list;TRS_dim2(j,k)];
                        end
                    end
                end

                if ~isempty(regulation_idx_list)
                    % choose maximum TRS if there are multiple regulation
                    [val,loc] = max(TRS_value_list);
                    TRS_2D_result(i) = val;
                    % save at regulation_list
                    pair_idx = regulation_idx_list(loc,1);
                    type_idx = regulation_idx_list(loc,2);
                    if val > TRS_1D_result(i)
                        st1 = component_list_dim2(pair_idx,1);
                        st2 = component_list_dim2(pair_idx,2);
                        ed = component_list_dim2(pair_idx,3);
                        if type_idx == 1
                            regulation_network_dim2(st1,ed) = 1;
                            regulation_network_dim2(st2,ed) = 1;
                        elseif type_idx == 2
                            regulation_network_dim2(st1,ed) = 1;
                            regulation_network_dim2(st2,ed) = -1;
                        elseif type_idx == 3
                            regulation_network_dim2(st1,ed) = -1;
                            regulation_network_dim2(st2,ed) = 1;
                        else
                            regulation_network_dim2(st1,ed) = -1;
                            regulation_network_dim2(st2,ed) = -1;
                        end
                    else
                        regulation_network_dim2(:,i) = regulation_network_dim1(:,i);

                    end
                else
                    regulation_network_dim2(:,i) = regulation_network_dim1(:,i);
                end
            end
            %regulation_network_dim2
            
            if system_idx == 1
                true_network = [
                    [0, 1, 0];
                    [0, 0, 1];
                    [-1,0, 0]];
            elseif system_idx == 2
                true_network = [
                    [0, 1, 0];
                    [0, 0, 1];
                    [-1,0, 0]];
            elseif system_idx == 3
                true_network = [
                    [0, 1, 0, 0];
                    [0, 0, 1, 0];
                    [0, 0, 0, 1]
                    [-1,0, 0, 0]];
            elseif system_idx == 4
                true_network = [
                    [0, 1,0,0,0];
                    [0, 0,1,0,0];
                    [0, 1,0,1,0];
                    [0, 0,1,0,0];
                    [-1,0,0,1,0]];
            elseif system_idx == 5
                true_network = [
                    [0, 0, 0, 0, 1,1,0];
                    [-1,0,-1, 0, 0,0,0];
                    [0, 0, 0,-1, 0,0,0];
                    [0, 0, 0, 0,-1,0,0];
                    [0, 1, 0, 0, 0,0,0];
                    [0, 0, 0, 0, 0,0,1];
                    [1, 0, 1, 0, 0,0,0]];
            end
            F2_score_tmp = cal_F2(regulation_network_dim2, true_network);
            F2_list_tmp = [F2_list_tmp,F2_score_tmp];
        end
        F2_list = [F2_list;F2_list_tmp];
    end
    F2_list_total = [F2_list_total; [1,mean(F2_list)]];
end

x = [0:2:20];
for i = 1:length(system_list)
    plot(x, F2_list_total(i,:),'linewidth',2)
    hold on
end
legend(system_list,'Location','southeast')
ylim([0,1])
xlim([0,20])
xticks([0:2:20])
yticks([0:0.2:1])
xticklabels([])
yticklabels([])