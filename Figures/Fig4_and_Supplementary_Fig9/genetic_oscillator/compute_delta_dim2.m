clc;
clear;
close all;

%% parameter
dimension = 2;
noise_level = 5;

thres_S = 0.9 - 0.005 * noise_level;
thres_L = 0.05;
thres_TRS = 0.9 - 0.01 * noise_level;

thres_delta = 0;
thres_noise = 1e-5;
thres_boot = 0.001;
%% load data
load('RDS_dim2')

%% compute TRS
num_type = 2^dimension;
TRS_total = zeros(num_pair, num_type);
for i = 1:num_pair
    for j = 1:num_type
        S_tmp = reshape(S_total_list(i,j,:),[num_data,1]);
        L_tmp = reshape(L_total_list(i,j,:),[num_data,1]);
        
        S_processed = S_threshold(S_tmp, thres_S);
        L_processed = L_threshold(L_tmp, thres_L);
        
        if sum(L_processed) == 0
            TRS_tmp = 0;
        else
            TRS_tmp = sum(S_processed .* L_processed) / sum(L_processed);
        end
        TRS_total(i,j) = TRS_tmp;
    end
end

%% change to matrix
TRS_list = zeros(num_component, num_type * num_component);
cause_list = nchoosek([1:num_component], 2);
cause_list = [[1,2]];
for i = 1:1
    for j = 1:num_component
%         if ismember(j, cause_list(i,:))
%             TRS_list(i,(j-1)*num_type+1:num_type*j) = nan;
%             continue
%         end
        
        % find index of TRS_total corresponding to TRS_list
        st1 = cause_list(1);
        st2 = cause_list(2);
        ed = j;
        idx_st1 = find(component_list(:,1) == st1);
        idx_st2 = find(component_list(:,2) == st2);
        idx_ed = find(component_list(:,3) == ed);
        idx = intersect(intersect(idx_st1,idx_st2),idx_ed);
        
        TRS_tmp = TRS_total(idx,:);
        TRS_list(i,(j-1)*num_type+1:num_type*j) = TRS_tmp;
    end
end

%% find candidate for delta test
delta_candidate_list = [];
type_list = [];
for i = 1:1
    for j = 1:num_component
%         if ismember(j, cause_list(i,:))
%             continue
%         end
        
        for k = 1:num_type
            if TRS_list(i,(j-1)*num_type + k) > thres_TRS
                delta_candidate_list = [delta_candidate_list ; [cause_list,j]];
                type_list = [type_list ; k];
            end
        end
     end
end
num_candidate_delta = length(delta_candidate_list(:,1));

%% calculate delta
delta_list = zeros(num_candidate_delta, 2);
cor_type = [
    [2,3];
    [1,4];
    [4,1];
    [3,2]];
for i = 1:num_candidate_delta
    st1 = delta_candidate_list(i,1);
    st2 = delta_candidate_list(i,2);
    ed = delta_candidate_list(i,3);
    idx_st1 = find(component_list(:,1) == st1);
    idx_st2 = find(component_list(:,2) == st2);
    idx_ed = find(component_list(:,3) == ed);
    idx = intersect(intersect(idx_st1, idx_st2),idx_ed);
    
    type_tmp = type_list(i);
    S_tmp_ori = reshape(S_total_list(idx,type_tmp,:), [num_data,1]);
    S_tmp_1 = reshape(S_total_list(idx,cor_type(type_tmp,1),:), [num_data,1]);
    S_tmp_2 = reshape(S_total_list(idx,cor_type(type_tmp,2),:), [num_data,1]);
    L_ori = reshape(L_total_list(idx,type_tmp,:), [num_data,1]);
    L_tmp_1 = reshape(L_total_list(idx,cor_type(type_tmp,1),:), [num_data,1]);
    L_tmp_2 = reshape(L_total_list(idx,cor_type(type_tmp,2),:), [num_data,1]);
    
    L_processed_ori = L_threshold(L_ori,thres_L);
    L_processed_1 = L_threshold(L_tmp_1,thres_L);
    L_processed_2 = L_threshold(L_tmp_2,thres_L);
    S_processed_ori = S_tmp_ori .* L_processed_ori;
    S_processed_1 = S_tmp_1 .* L_processed_1;
    S_processed_2 = S_tmp_2 .* L_processed_2;
    
    
    S_processed_1(find(S_processed_1 == 0)) = NaN;
    S_processed_2(find(S_processed_2 == 0)) = NaN;
    S_processed_ori(find(S_processed_ori == 0)) = NaN;
    %S_processed_2(find(isnan(S_processed_2) == 1)) = -1;
    %S_processed_ori(find(isnan(S_processed_ori) == 1)) = -1;
    
    delta_1 = S_processed_ori - S_processed_1;
    delta_2 = S_processed_ori - S_processed_2;
    
    %delta_1 = S_tmp_ori - S_tmp_1;
    %delta_2 = S_tmp_ori - S_tmp_2;
    
    p1 = 0;
    p2 = 0;
    %delta_1(find(isnan(delta_1) == 1)) = 1;
    %delta_2(find(isnan(delta_2) == 1)) = 1;
    if ~isempty(rmmissing(delta_1))
        [p1,h1,stats1] = signrank(rmmissing(delta_1),-1e-2,'tail','right');
        %p1 = length(find(rmmissing(delta_1) < -1e-3)) / length(rmmissing(delta_1));
    end
    if ~isempty(rmmissing(delta_2))
        [p2,h2,stats2] = signrank(rmmissing(delta_2),-1e-2,'tail','right');
        %p2 = length(find(rmmissing(delta_2) < -1e-3)) / length(rmmissing(delta_2));
    end
    delta_list(i,:) = [p1,p2];
end
