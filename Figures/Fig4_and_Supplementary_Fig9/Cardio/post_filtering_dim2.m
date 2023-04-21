clc;
clear;
close all;

%% parameter
dimension = 2;
noise_level = 30;

thres_S = 0.9 - 0.005 * noise_level;
thres_L = 0.05;
thres_TRS = 0.9 - 0.01 * noise_level;

thres_delta = 0;
thres_noise = 0;
thres_boot = 0.001;
%% load data

load('Cardio_score_dim2')
load('data_cardio')


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
%TRS_list = zeros(num_component, num_type * num_component);
TRS_list = zeros(num_component, num_type * 1);
%cause_list = nchoosek([1:num_component], 2);
cause_list = nchoosek([2:num_component], 2);
for i = 1:length(cause_list)
    %for j = 1:num_component
    for j = 1
        if ismember(j, cause_list(i,:))
            TRS_list(i,(j-1)*num_type+1:num_type*j) = nan;
            continue
        end
        
        % find index of TRS_total corresponding to TRS_list
        st1 = cause_list(i,1);
        st2 = cause_list(i,2);
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
for i = 1:length(cause_list)
    %for j = 1:num_component
    for j = 1
        if ismember(j, cause_list(i,:))
            continue
        end
        
        for k = 1:num_type
            if TRS_list(i,(j-1)*num_type + k) >= thres_TRS
                delta_candidate_list = [delta_candidate_list ; [cause_list(i,:),j]];
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
%for i = 6
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
    %S_processed_ori(find(S_processed_ori == 0)) = NaN;
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
        %delta_1(find(isnan(delta_1) == 1)) = 1;
        %delta_2(find(isnan(delta_2) == 1)) = 1;
        %[p1,h1,stats1] = signrank(rmmissing(delta_1),-1e-3,'tail','right');
        
        p1 = length(find(rmmissing(delta_1) < -1e-3)) / length(rmmissing(delta_1));
    end
    if ~isempty(rmmissing(delta_2))
        %delta_1(find(isnan(delta_1) == 1)) = 1;
        %delta_2(find(isnan(delta_2) == 1)) = 1;
        %[p2,h2,stats2] = signrank(rmmissing(delta_2),-1e-3,'tail','right');
        p2 = length(find(rmmissing(delta_2) < -1e-3)) / length(rmmissing(delta_2));
    end
    delta_list(i,:) = [p1,p2];
end

%% find candidate for surrogate time-series analysis
boot_candidate_list = [];
boot_type_list = [];
for i = 1:num_candidate_delta
    if delta_list(i,1) <= thres_delta && delta_list(i,2) <= thres_delta
        boot_candidate_list = [boot_candidate_list ; delta_candidate_list(i,:)];
        boot_type_list = [boot_type_list; type_list(i)];
    end
end
num_candidate_boot = length(boot_candidate_list(:,1));

%% start parallel pool
parpool threads;
clear completedJobs;
dq = parallel.pool.DataQueue;
wb = waitbar(0,'Processing');
N = num_data * num_candidate_boot;
Listener = afterEach(dq, @(varargin) waitbar((completedJobs/N),wb,sprintf('Completed: %d', completedJobs(1))));


%% bootstrapping
num_boot = 100;
surrogate_list = [];
for i = 1:num_candidate_boot
    disp(i)
    st1 = boot_candidate_list(i,1);
    st2 = boot_candidate_list(i,2);
    ed = boot_candidate_list(i,3);
    type_tmp = boot_type_list(i);
    p_total = [];
    parfor j = 1:num_data
        y_tmp = cell2mat(y_total(j));    
        C1 = y_tmp(:,st1);
        C2 = y_tmp(:,st2);
        T = y_tmp(:,ed);
        t_target = t(1:length(y_tmp(:,1)));
        
        boot_tmp = [];
        for k = 1:num_boot
            C1_shuffled = C1(randperm(length(C1)));
            [score_list, t_1, t_2] = RDS_ns_dim2(C1_shuffled, C2, T, t_target, time_interval);
           
            
            score_tmp = reshape(score_list(:,:,type_tmp),[length(t_1(:,1)),length(t_1(1,:))]);
            
            loca_plus = find(score_tmp > thres_noise);
            loca_minus = find(score_tmp < -thres_noise);
            if isempty(loca_plus) && isempty(loca_minus)
                s_tmp_1 = 1;
            else
                s_tmp_1 = (sum(score_tmp(loca_plus)) + sum(score_tmp(loca_minus)))/ (abs(sum(score_tmp(loca_plus))) + abs(sum(score_tmp(loca_minus))));
            end
            
            C2_shuffled = C2(randperm(length(C2)));
            [score_list, t_1, t_2] = RDS_ns_dim2(C1, C2_shuffled, T, t_target, time_interval);
            
            
            score_tmp = reshape(score_list(:,:,type_tmp),[length(t_1(:,1)),length(t_1(1,:))]);
            
            loca_plus = find(score_tmp > thres_noise);
            loca_minus = find(score_tmp < -thres_noise);
            if isempty(loca_plus) && isempty(loca_minus)
                s_tmp_2 = 1;
            else
                s_tmp_2 = (sum(score_tmp(loca_plus)) + sum(score_tmp(loca_minus)))/ (abs(sum(score_tmp(loca_plus))) + abs(sum(score_tmp(loca_minus))));
            end
            
            boot_tmp = [boot_tmp ;[s_tmp_1,s_tmp_2]];
        end
        [score_list, t_1, t_2] = RDS_ns_dim2(C1, C2, T, t_target, time_interval);
       
        score_tmp = reshape(score_list(:,:,type_tmp),[length(t_1(:,1)),length(t_1(1,:))]);

        loca_plus = find(score_tmp > thres_noise);
        loca_minus = find(score_tmp < -thres_noise);
        if isempty(loca_plus) && isempty(loca_minus)
            s_ori = -1;
        else
            s_ori = (sum(score_tmp(loca_plus)) + sum(score_tmp(loca_minus)))/ (abs(sum(score_tmp(loca_plus))) + abs(sum(score_tmp(loca_minus))));
        end
        l_ori = (length(loca_minus) + length(loca_plus)) / (length(t_1)*length(t_2)/2);
        %if l_ori > thres_L        
            [h1,p1] = ztest(s_ori, mean(boot_tmp(:,1)),std(boot_tmp(:,1)),'Tail','right');
            [h1_pre,p1_pre] = ztest(1+1e-2, mean(boot_tmp(:,1)),std(boot_tmp(:,1)),'Tail','right');
            [h2,p2] = ztest(s_ori, mean(boot_tmp(:,2)),std(boot_tmp(:,2)),'Tail','right');
            [h2_pre,p2_pre] = ztest(1+1e-2, mean(boot_tmp(:,2)),std(boot_tmp(:,2)),'Tail','right');
            %p_total = [p_total ; [p1-p1_pre,p2-p2_pre]];
            p_total = [p_total ; [p1,p2]];
        %end
    end
    if ~isempty(p_total)
        p_total(p_total == nan) = 0.5;
        p_tmp_1 = nonzeros(rmmissing(p_total(:,1)));
        
        p_tmp_2 = nonzeros(rmmissing(p_total(:,2)));
        
        fisher_tmp_1 = 2* sum(-log(p_tmp_1));
        fisher_tmp_2 = 2* sum(-log(p_tmp_2));
        num_p_1 = length(p_tmp_1);
        num_p_2 = length(p_tmp_2);
        fisher_thres_1 = chi2cdf(-2*log(thres_boot)*num_p_1, 2*num_p_1, 'upper');
        fisher_thres_2 = chi2cdf(-2*log(thres_boot)*num_p_2, 2*num_p_2, 'upper');
        fisher_tmp = [chi2cdf(fisher_tmp_1, 2*num_p_1, 'upper'), chi2cdf(fisher_tmp_2, 2*num_p_2, 'upper'),fisher_thres_1,fisher_thres_2];
    else
        fisher_tmp = [0,0,0,0];
    end
    
    surrogate_list = [surrogate_list; [fisher_tmp]];
end
fisher_thres = chi2cdf(-2*log(thres_boot)*num_data, 2*num_data, 'upper')


filename = ['post_filtering_dim2'];
save(filename, 'boot_candidate_list','boot_type_list','num_data','surrogate_list', 'delta_list', 'type_list', 'delta_candidate_list')



delete(gcp('nocreate'))

%% function for parallel pool
function j = completedJobs(varargin)
    persistent n
    if isempty(n)
        n = 0;
    end
    if numel(varargin) ~=0
    else
        n = n+1;
    end
    j=n;
end

