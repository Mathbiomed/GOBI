clc;
clear;
close all;

%% load data
system = 'IFL';
filename = [system, '_TRS_dim1'];
load(filename)
filename = [system, '_TRS_dim2'];
load(filename)
filename = [system, '_Delta_dim2'];
load(filename)
filename = [system, '_Surrogate_dim2'];
load(filename)

num_component = 3;

%% 1D regulation
regulation_network_dim1 = zeros(num_component);
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
                TRS_value_list = [TRS_value_list;TRS_total_dim1(j,k)];
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
regulation_network_dim2 = zeros(num_component);
num_candidate_boot = length(boot_candidate_list(:,1));
boot_2dim = ones(num_candidate_boot,2);

for i = 1:num_candidate_boot
    for k = 1:2
        if surrogate_list(i,k) > surrogate_list(i,k+2)
            boot_2dim(i,k) = 0;
        end
    end
end

indirect_list = ones(num_candidate_boot, 2);
candidate_2dim = [];

for i = 1:num_candidate_boot
    if indirect_list(i,1) == 1 && boot_2dim(i,1) == 0
        continue
    elseif indirect_list(i,2) == 1 && boot_2dim(i,2) == 0
        continue
    else
        idx = intersect(find(component_list_dim2(:,1) == boot_candidate_list(i,1)), find(component_list_dim2(:,2) == boot_candidate_list(i,2)));
        candidate_2dim = [candidate_2dim;[idx, boot_type_list(i)]];
    end
end

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
        for k = 1:4
            if result_boot_2dim(j,k) == 1
                regulation_idx_list = [regulation_idx_list;[j,k]];
                TRS_value_list = [TRS_value_list;TRS_total_dim2(j,k)];
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

regulation_network_dim2
