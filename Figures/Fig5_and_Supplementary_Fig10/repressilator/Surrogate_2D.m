clc;
clear;
close all;
addpath('../../GOBI') 

%% parameters
dimension = 2;
thres_noise = 1e-7;
p_thres = 0.001;
t_L = 0.01;
thres_noise = 1e-5;

%% load data
type = 'full';
%type = 'half';
%type = 'quarter';

filename = ['RDS_dim2_',type];
load(filename)
component_list_dim2 = component_list;
filename = ['TRS_1D_',type];
load(filename)
component_list_dim1 = component_list;
filename = ['Delta_dim2_',type];
load(filename)
filename = ['data_merge_cut_',type];
load(filename)

num_boot = 200;
parpool threads;
    
%% 1D regulation       
regulation_network_dim1 = zeros(num_component,num_component);
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
%% find regulation pair which pass the criteria of TRS (idx_pair, idx_type)
candidate = [];
for i = 1:num_pair
    for j = 1:num_type
        if delta_2dim(i,j) == 1
            candidate = [candidate;[i,j]];
        end
    end
end
%candidate = [[2,1];[20,1];[29,1]];
if isempty(candidate)
    num_candidate = 0;
else
    num_candidate = length(candidate(:,1));
end

%% check each candidate whether indirect or not
indirect_list = [];
for i = 1:num_candidate
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

    regulation_2D = regulation_network_dim1;
    regulation_2D(cause_1, target) = type_idx(1);
    regulation_2D(cause_2, target) = type_idx(2);

%         % first component
%         result_indirect_1 = [];
%         for j = 1:num_component
%             if regulation_2D(cause_1,j) ~= 0 && j ~= target
%                 % [ori, prev, curr, tar, type, results, network]
%                 result_indirect_1 = [result_indirect_1;isIndirect(cause_1, cause_1,j,target,regulation_2D(cause_1,j),result_indirect_1,regulation_2D)];
%             end
%         end
%         if ~isempty(find(result_indirect_1 == type_idx(1)))
%             indirect_idx(1) = 1;
%         end
% 
%         % second component
%         result_indirect_2 = [];
%         for j = 1:num_component
%             if regulation_2D(cause_2,j) ~= 0 && j ~= target
%                 % [ori, prev, curr, tar, type, results, network]
%                 result_indirect_2 = [result_indirect_2;isIndirect(cause_2, cause_2,j,target,regulation_2D(cause_2,j),result_indirect_2,regulation_2D)];
%             end
%         end
%         if ~isempty(find(result_indirect_2 == type_idx(2)))
%             indirect_idx(2) = 1;
%         end
    %indirect_list = [indirect_list;indirect_idx];
    indirect_list = [indirect_list;[1,1]];
end

%% Bootstrapping
S_boot_total = zeros(num_boot,dimension,num_data,num_candidate);
parfor i = 1:num_candidate
    disp(num2str(i/num_candidate))
    for j = 1:2
        if indirect_list(i,j) == 1
            % indexing [fix_index,shuffling_index,target_index]
            fix_idx = component_list_dim2(candidate(i,1),3-j);
            shu_idx = component_list_dim2(candidate(i,1),j);
            tar_idx = component_list_dim2(candidate(i,1),3);
            type_idx = candidate(i,2);

            % import score and region
            S_ori_tmp = reshape(S_total(candidate(i,1),candidate(i,2),:),[num_data,1]);
            L_ori_tmp = reshape(L_total(candidate(i,1),candidate(i,2),:),[num_data,1]);

            L_ori_tmp = L_threshold(L_ori_tmp, t_L);
            S_ori_tmp = S_ori_tmp .* L_ori_tmp;

            % bootstrap only when score is not empty
            % for each data
            for k = 1:num_data
                %send(dq, j)
                if S_ori_tmp(k) ~= 0
                    y_tmp = cell2mat(y_total(k));
                    ts_fix = y_tmp(:,fix_idx);
                    ts_shu = y_tmp(:,shu_idx);
                    ts_tar = y_tmp(:,tar_idx);

                    % calculate score for bootstrapping
                    S_boot_tmp = zeros(num_boot,1);
                    for m = 1:num_boot
                        ts_shuffled = ts_shu(randperm(length(ts_shu)));
                        if j == 2
                            [score_list, t_1, t_2] = RDS_ns_dim2(ts_fix, ts_shuffled, ts_tar, t, time_interval);
                        else
                            [score_list, t_1, t_2] = RDS_ns_dim2(ts_shuffled, ts_fix, ts_tar, t, time_interval);
                        end
                        score = reshape(score_list(:,:,type_idx),[length(t),length(t)]);
                        loca_plus = find(score > thres_noise);
                        loca_minus = find(score < -thres_noise);
                        if isempty(loca_plus) && isempty(loca_minus)
                            s = nan;
                        else
                            s = (sum(score(loca_plus)) + sum(score(loca_minus)))/ (abs(sum(score(loca_plus))) + abs(sum(score(loca_minus))));
                        end

                        l = (length(loca_minus) + length(loca_plus)) / (length(t_1)*length(t_2)/2);

                        if l > t_L
                            S_boot_tmp(m) = s;
                        else
                            S_boot_tmp(m) = nan;
                        end
                    end

                    S_boot_total(:,j,k,i) = S_boot_tmp;
                end
            end
        end
    end
end

%% calculte p-value using z-test and Fisher's method
boot_2dim = indirect_list;
p_total = [];
for i = 1:num_candidate
    for k = 1:2
        if indirect_list(i,k) == 1
            p_list = [];
            r_list = [];

            %import score
            S_ori_tmp = reshape(S_total(candidate(i,1),candidate(i,2),:),[num_data,1]);
            L_ori_tmp = reshape(L_total(candidate(i,1),candidate(i,2),:),[num_data,1]);
            L_ori_tmp = L_threshold(L_ori_tmp, t_L);
            S_ori_tmp = S_ori_tmp .* L_ori_tmp;
            % ztest
            for j = 1:num_data
                if S_ori_tmp(j) ~= 0
                    S_boot_tmp = reshape(S_boot_total(:,k,j,i),[num_boot,1]);

                    [h,p] = ztest(S_ori_tmp(j), mean(S_boot_tmp),std(S_boot_tmp),'Tail','right');

                    p_list = [p_list;p];

                    r = (length(find(S_boot_tmp > S_ori_tmp(j)))+1) / length(S_boot_tmp);
                    r_list = [r_list;r];
                end
            end
            % Fisher's method
            num_p = length(nonzeros(rmmissing(p_list)));
            sum_p = nansum(-2*log(nonzeros(p_list)));

            %sum_p = nansum(-2*log(r_list));
            %num_p = [length(r_list(:,1)),length(r_list(:,2))];

            Fisher_thres = -2*log(p_thres)*num_p;
            if sum_p < Fisher_thres
                boot_2dim(i,k) = 0;
            end
            p_total = [p_total ; [sum_p,num_p,Fisher_thres]];
        end
    end
end

%% remove indirect regulation from candidate
candidate_2dim = [];
for i = 1:num_candidate
    if indirect_list(i,1) == 1 && boot_2dim(i,1) == 0
        continue
    elseif indirect_list(i,2) == 1 && boot_2dim(i,2) == 0
        continue
    else
        candidate_2dim = [candidate_2dim;candidate(i,:)];
    end
end
candidate_2dim

filename = ['Surrogate_dim2_',type];
save(filename, 'boot_2dim', 'S_boot_total','candidate_2dim','indirect_list','component_list_dim2','num_data','num_type','num_pair','dimension','S_boot_total')


delete(gcp('nocreate'))