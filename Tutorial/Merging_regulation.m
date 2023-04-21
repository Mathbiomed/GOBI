clc;
clear;
close all;

%% import 1D, 2D results
load('data_with_options')
load('TRS_dim1')
load('TRS_dim2')
load('Delta_dim2')
load('Surrogate_dim2')

%% 1D regulation  
TRS_1D_result = zeros(num_component,1);
for i = 1:num_component
    % indexing all the causes of each target
    index_list = [];
    for j = 1:length(component_list_dim1(:,1))
        if component_list_dim1(j,2) == i
            index_list = [index_list, j];
        end
    end

    % check inferred 1D regulation
    regulation_idx_list = [];
    TRS_value_list = [];
    for j = index_list
        for k = 1:length(regulation_1dim(1,:))
            if regulation_1dim(j,k) == 1
                regulation_idx_list = [regulation_idx_list;[j,k]];
                TRS_value_list = [TRS_value_list;TRS_total_dim1(j,k)];
            end
        end
    end

    if ~isempty(regulation_idx_list)
        if length(regulation_idx_list(:,1)) == 1
            %[val,loc] = max(TRS_value_list);
            %TRS_1D_result(i) = val;
            % save at regulation_list
            pair_idx = regulation_idx_list(1,1);
            type_idx = regulation_idx_list(1,2);

            st = component_list_dim1(pair_idx,1);
            ed = component_list_dim1(pair_idx,2);
            if type_idx == 1
                regulation_network_dim1(st,ed) = 1;
            else
                regulation_network_dim1(st,ed) = -1;
            end
        end
    end
end

regulation_network_dim1

%% 2D regulation  
regulation_network_dim2 = zeros(num_component);
num_candidate_boot = length(boot_candidate_list(:,1));

% find the index for candidates of 2D regulations
candidate = [];
for i = 1:num_candidate_boot
    idx1 = find(component_list_dim2(:,1) == boot_candidate_list(i,1));
    idx2 = find(component_list_dim2(:,2) == boot_candidate_list(i,2));
    idx3 = find(component_list_dim2(:,3) == boot_candidate_list(i,3));
    idx = intersect(intersect(idx1, idx2), idx3);
    
    candidate = [candidate; [idx, boot_type_list(i)]];
end

%% make basis network to find the potential indirect effect (merging every inferred 2D candidates with inferred 1D regulations)
regulation_2D = regulation_network_dim1;
for i = 1:num_candidate_boot
    indirect_idx = [0,0];
    tar_pair = candidate(i,1);
    tar_type = candidate(i,2);
    if tar_type == 1
        type_idx = [1,1];
    elseif tar_type == 2
        type_idx = [1,-1];
    elseif tar_type == 3
        type_idx = [-1,1];
    else
        type_idx = [-1,-1];
    end

    cause_1 = component_list_dim2(tar_pair,1);
    cause_2 = component_list_dim2(tar_pair,2);
    target  = component_list_dim2(tar_pair,3);

    regulation_2D(cause_1, target) = type_idx(1);
    regulation_2D(cause_2, target) = type_idx(2);
end

%% check each candidate whether indirect or not
indirect_list = [];
for i = 1:num_candidate_boot
    indirect_idx = [0,0];
    tar_pair = candidate(i,1);
    tar_type = candidate(i,2);
    
    if tar_type == 1
        type_idx = [1,1];
    elseif tar_type == 2
        type_idx = [1,-1];
    elseif tar_type == 3
        type_idx = [-1,1];
    else
        type_idx = [-1,-1];
    end
    cause_1 = component_list_dim2(tar_pair,1);
    cause_2 = component_list_dim2(tar_pair,2);
    target  = component_list_dim2(tar_pair,3);
    
    % check first causal variable
    result_indirect_1 = [];
    for j = 1:num_component
        if regulation_2D(cause_1,j) ~= 0 && j ~= target
            % [ori, prev, curr, tar, type, results, network]
            result_indirect_1 = [result_indirect_1;isIndirect(cause_1, cause_1,j,target,regulation_2D(cause_1,j),result_indirect_1,regulation_2D)];
        end
    end
    if ~isempty(find(result_indirect_1 == type_idx(1)))
        indirect_idx(1) = 1;
    end

    % check second causal variable
    result_indirect_2 = [];
    for j = 1:num_component
        if regulation_2D(cause_2,j) ~= 0 && j ~= target
            % [ori, prev, curr, tar, type, results, network]
            result_indirect_2 = [result_indirect_2;isIndirect(cause_2, cause_2,j,target,regulation_2D(cause_2,j),result_indirect_2,regulation_2D)];
        end
    end
    if ~isempty(find(result_indirect_2 == type_idx(2)))
        indirect_idx(2) = 1;
    end
    indirect_list = [indirect_list;indirect_idx];
end

% it is possible to consider all the regulations as potential indirect regulations
%indirect_list = ones(num_candidate_boot, 2);

%% infer the regulation using p_surrogate
% using the results of surrogate test and potential indirect regulations
boot_2dim = ones(num_candidate_boot,2);
for i = 1:num_candidate_boot
    for k = 1:2
        if indirect_list(i,k) == 1 && surrogate_list(i,k) > surrogate_list(i,k+2)
            boot_2dim(i,k) = 0;
        end
    end
end

% merge the results (nan means indirect regulation)
for i = 1:num_candidate_boot
    if boot_2dim(i,1) == 1
        st = boot_candidate_list(i,1);
        ed = boot_candidate_list(i,3);
        if ~isnan(regulation_network_dim2(st,ed))
            if boot_type_list(i) == 1 || boot_type_list(i) == 2
                regulation_network_dim2(st,ed) = 1;
            else
                regulation_network_dim2(st,ed) = -1;
            end
        end
    else
        st = boot_candidate_list(i,1);
        ed = boot_candidate_list(i,3);
        regulation_network_dim2(st,ed) = nan;
        
    end
    if boot_2dim(i,2) == 1
        st = boot_candidate_list(i,2);
        ed = boot_candidate_list(i,3);
        if ~isnan(regulation_network_dim2(st,ed))
            if boot_type_list(i) == 1 || boot_type_list(i) == 3
                regulation_network_dim2(st,ed) = 1;
            else
                regulation_network_dim2(st,ed) = -1;
            end
        end
    else
        st = boot_candidate_list(i,2);
        ed = boot_candidate_list(i,3);
        regulation_network_dim2(st,ed) = nan;
    end
end
for i = 1:num_component
    for j = 1:num_component
        if isnan(regulation_network_dim2(i,j))
            regulation_network_dim2(i,j) = 0;
        end
    end
end

regulation_network_dim2