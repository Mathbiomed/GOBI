clc;
clear;
close all;

%% parameter
dimension = 4;
noise_level = 30;

thres_S = 0.9 - 0.005 * noise_level;
thres_L = 0.005;
thres_TRS = 0.9 - 0.01 * noise_level;

thres_delta = 0.01;
thres_noise = 0;
thres_boot = 0.001;
%% load data

load('Cardio_score_dim4')
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

%% find candidate for delta test
delta_candidate_list = [];
type_list = [];
for i = 1:num_pair
        
    for j = 1:num_type
        if TRS_total(i, j) >= thres_TRS
            delta_candidate_list = [delta_candidate_list ; component_list(i,1:end-1)];
            type_list = [type_list ; j];
        end
    end
end
num_candidate_delta = length(delta_candidate_list(:,1));

%% calculate delta
delta_list = zeros(num_candidate_delta, 4);
cor_type = [
    [2,3,5,9];
    [1,4,6,10];
    [4,1,7,11];
    [3,2,8,12];
    [6,7,1,13];
    [5,8,2,14];
    [8,5,3,15];
    [7,6,4,16];
    [10,11,13,1];
    [9,12,14,2];
    [12,9,15,3];
    [11,10,16,4];
    [14,15,9,5];
    [13,16,10,6];
    [16,13,11,7];
    [15,14,12,8]];
for i = 1:num_candidate_delta
%for i = 6
    st1 = delta_candidate_list(i,1);
    st2 = delta_candidate_list(i,2);
    st3 = delta_candidate_list(i,3);
    st4 = delta_candidate_list(i,4);
    idx_st1 = find(component_list(:,1) == st1);
    idx_st2 = find(component_list(:,2) == st2);
    idx_st3 = find(component_list(:,3) == st3);
    idx_st4 = find(component_list(:,4) == st4);
    idx = intersect(intersect(intersect(idx_st1, idx_st2),idx_st3),idx_st4);
    
    type_tmp = type_list(i);
    S_tmp_ori = reshape(S_total_list(idx,type_tmp,:), [num_data,1]);
    S_tmp_1 = reshape(S_total_list(idx,cor_type(type_tmp,1),:), [num_data,1]);
    S_tmp_2 = reshape(S_total_list(idx,cor_type(type_tmp,2),:), [num_data,1]);
    S_tmp_3 = reshape(S_total_list(idx,cor_type(type_tmp,3),:), [num_data,1]);
    S_tmp_4 = reshape(S_total_list(idx,cor_type(type_tmp,4),:), [num_data,1]);
    
    L_ori = reshape(L_total_list(idx,type_tmp,:), [num_data,1]);
    L_tmp_1 = reshape(L_total_list(idx,cor_type(type_tmp,1),:), [num_data,1]);
    L_tmp_2 = reshape(L_total_list(idx,cor_type(type_tmp,2),:), [num_data,1]);
    L_tmp_3 = reshape(L_total_list(idx,cor_type(type_tmp,3),:), [num_data,1]);
    L_tmp_4 = reshape(L_total_list(idx,cor_type(type_tmp,4),:), [num_data,1]);
    
    L_processed_ori = L_threshold(L_ori,thres_L);
    L_processed_1 = L_threshold(L_tmp_1,thres_L);
    L_processed_2 = L_threshold(L_tmp_2,thres_L);
    L_processed_3 = L_threshold(L_tmp_3,thres_L);
    L_processed_4 = L_threshold(L_tmp_4,thres_L);
    
    S_processed_ori = S_tmp_ori .* L_processed_ori;
    S_processed_1 = S_tmp_1 .* L_processed_1;
    S_processed_2 = S_tmp_2 .* L_processed_2;
    S_processed_3 = S_tmp_3 .* L_processed_3;
    S_processed_4 = S_tmp_4 .* L_processed_4;
    
    
    S_processed_1(find(S_processed_1 == 0)) = NaN;
    S_processed_2(find(S_processed_2 == 0)) = NaN;
    S_processed_3(find(S_processed_3 == 0)) = NaN;
    S_processed_4(find(S_processed_4 == 0)) = NaN;
    %S_processed_ori(find(S_processed_ori == 0)) = NaN;
    %S_processed_2(find(isnan(S_processed_2) == 1)) = -1;
    %S_processed_ori(find(isnan(S_processed_ori) == 1)) = -1;
    
    delta_1 = S_processed_ori - S_processed_1;
    delta_2 = S_processed_ori - S_processed_2;
    delta_3 = S_processed_ori - S_processed_3;
    delta_4 = S_processed_ori - S_processed_4;
    
    %delta_1 = S_tmp_ori - S_tmp_1;
    %delta_2 = S_tmp_ori - S_tmp_2;
    
    p1 = NaN;
    p2 = NaN;
    p3 = NaN;
    p4 = NaN;
    %delta_1(find(isnan(delta_1) == 1)) = 1;
    %delta_2(find(isnan(delta_2) == 1)) = 1;
    if ~isempty(rmmissing(delta_1))
        %delta_1(find(isnan(delta_1) == 1)) = 1;
        %[p1,h1,stats1] = signrank(rmmissing(delta_1),-1e-3,'tail','right');
        
        p1 = length(find(rmmissing(delta_1) < -1e-3)) / length(rmmissing(delta_1));
    end
    
    if ~isempty(rmmissing(delta_2))
        %delta_2(find(isnan(delta_2) == 1)) = 1;
        %[p2,h2,stats2] = signrank(rmmissing(delta_2),-1e-3,'tail','right');
        p2 = length(find(rmmissing(delta_2) < -1e-3)) / length(rmmissing(delta_2));
    end
    
    if ~isempty(rmmissing(delta_3))
        %delta_3(find(isnan(delta_3) == 1)) = 1;
        %[p3,h3,stats3] = signrank(rmmissing(delta_3),-1e-3,'tail','right');
        p3 = length(find(rmmissing(delta_3) < -1e-3)) / length(rmmissing(delta_3));
    end
    
    if ~isempty(rmmissing(delta_4))
        %delta_3(find(isnan(delta_3) == 1)) = 1;
        %[p3,h3,stats3] = signrank(rmmissing(delta_3),-1e-3,'tail','right');
        p4 = length(find(rmmissing(delta_4) < -1e-3)) / length(rmmissing(delta_4));
    end
    
    delta_list(i,:) = [p1,p2,p3,p4];
end

%% find candidate for surrogate time-series analysis
boot_candidate_list = [];
boot_type_list = [];
for i = 1:num_candidate_delta
    if delta_list(i,1) <= thres_delta && delta_list(i,2) <= thres_delta && delta_list(i,3) < thres_delta
        boot_candidate_list = [boot_candidate_list ; delta_candidate_list(i,:)];
        boot_type_list = [boot_type_list; type_list(i)];
    end
end
num_candidate_boot = length(boot_candidate_list(:,1));

% %% bootstrapping
% num_boot = 100;
% fisher_total = [];
% for i = 1:num_candidate_boot
% %for i = 5
%     disp(i)
%     st1 = boot_candidate_list(i,1);
%     st2 = boot_candidate_list(i,2);
%     st3 = boot_candidate_list(i,3);
%     ed = 1;
%     type_tmp = boot_type_list(i);
%     p_total = [];
%     for j = 1:num_data
%         y_tmp = cell2mat(y_total(j));    
%         C1 = y_tmp(:,st1);
%         C2 = y_tmp(:,st2);
%         C3 = y_tmp(:,st3);
%         T = y_tmp(:,ed);
%         t_target = t(1:length(y_tmp(:,1)));
%         
%         boot_tmp = [];
%         for k = 1:num_boot
%             C1_shuffled = C1(randperm(length(C1)));
%             [score_1, score_2, score_3, score_4, score_5, score_6, score_7, score_8, t_1, t_2] = ion_prime_dim3_ns(C1_shuffled, C2, C3, T, t_target, time_interval);
%             score_total = zeros(length(t_1(:,1)),length(t_1(1,:)),4);
%             score_total(:,:,1) = score_1;
%             score_total(:,:,2) = score_2;
%             score_total(:,:,3) = score_3;
%             score_total(:,:,4) = score_4;
%             score_total(:,:,5) = score_5;
%             score_total(:,:,6) = score_6;
%             score_total(:,:,7) = score_7;
%             score_total(:,:,8) = score_8;
%             
%             score_tmp = reshape(score_total(:,:,type_tmp),[length(t_1(:,1)),length(t_1(1,:))]);
%             
%             loca_plus = find(score_tmp > thres_noise);
%             loca_minus = find(score_tmp < -thres_noise);
%             if isempty(loca_plus) && isempty(loca_minus)
%                 s_tmp_1 = 1;
%             else
%                 s_tmp_1 = (sum(score_tmp(loca_plus)) + sum(score_tmp(loca_minus)))/ (abs(sum(score_tmp(loca_plus))) + abs(sum(score_tmp(loca_minus))));
%             end
%             
%             C2_shuffled = C2(randperm(length(C2)));
%             [score_1, score_2, score_3, score_4, score_5, score_6, score_7, score_8, t_1, t_2] = ion_prime_dim3_ns(C1, C2_shuffled, C3, T, t_target, time_interval);
%             score_total = zeros(length(t_1(:,1)),length(t_1(1,:)),4);
%             score_total(:,:,1) = score_1;
%             score_total(:,:,2) = score_2;
%             score_total(:,:,3) = score_3;
%             score_total(:,:,4) = score_4;
%             score_total(:,:,5) = score_5;
%             score_total(:,:,6) = score_6;
%             score_total(:,:,7) = score_7;
%             score_total(:,:,8) = score_8;
%             
%             score_tmp = reshape(score_total(:,:,type_tmp),[length(t_1(:,1)),length(t_1(1,:))]);
%             
%             loca_plus = find(score_tmp > thres_noise);
%             loca_minus = find(score_tmp < -thres_noise);
%             if isempty(loca_plus) && isempty(loca_minus)
%                 s_tmp_2 = 1;
%             else
%                 s_tmp_2 = (sum(score_tmp(loca_plus)) + sum(score_tmp(loca_minus)))/ (abs(sum(score_tmp(loca_plus))) + abs(sum(score_tmp(loca_minus))));
%             end
%             
%             C3_shuffled = C3(randperm(length(C3)));
%             [score_1, score_2, score_3, score_4, score_5, score_6, score_7, score_8, t_1, t_2] = ion_prime_dim3_ns(C1, C2, C3_shuffled, T, t_target, time_interval);
%             score_total = zeros(length(t_1(:,1)),length(t_1(1,:)),4);
%             score_total(:,:,1) = score_1;
%             score_total(:,:,2) = score_2;
%             score_total(:,:,3) = score_3;
%             score_total(:,:,4) = score_4;
%             score_total(:,:,5) = score_5;
%             score_total(:,:,6) = score_6;
%             score_total(:,:,7) = score_7;
%             score_total(:,:,8) = score_8;
%             
%             score_tmp = reshape(score_total(:,:,type_tmp),[length(t_1(:,1)),length(t_1(1,:))]);
%             
%             loca_plus = find(score_tmp > thres_noise);
%             loca_minus = find(score_tmp < -thres_noise);
%             if isempty(loca_plus) && isempty(loca_minus)
%                 s_tmp_3 = 1;
%             else
%                 s_tmp_3 = (sum(score_tmp(loca_plus)) + sum(score_tmp(loca_minus)))/ (abs(sum(score_tmp(loca_plus))) + abs(sum(score_tmp(loca_minus))));
%             end
%             
%             boot_tmp = [boot_tmp ;[s_tmp_1,s_tmp_2,s_tmp_3]];
%         end
%         
%         [score_1, score_2, score_3, score_4,score_5, score_6, score_7, score_8, t_1, t_2] = ion_prime_dim3_ns(C1, C2, C3, T, t_target, time_interval);
%         score_total = zeros(length(t_1(:,1)),length(t_1(1,:)),4);
%         score_total(:,:,1) = score_1;
%         score_total(:,:,2) = score_2;
%         score_total(:,:,3) = score_3;
%         score_total(:,:,4) = score_4;
%         score_total(:,:,5) = score_5;
%         score_total(:,:,6) = score_6;
%         score_total(:,:,7) = score_7;
%         score_total(:,:,8) = score_8;
%         
%         score_tmp = reshape(score_total(:,:,type_tmp),[length(t_1(:,1)),length(t_1(1,:))]);
% 
%         loca_plus = find(score_tmp > thres_noise);
%         loca_minus = find(score_tmp < -thres_noise);
%         if isempty(loca_plus) && isempty(loca_minus)
%             s_ori = -1;
%         else
%             s_ori = (sum(score_tmp(loca_plus)) + sum(score_tmp(loca_minus)))/ (abs(sum(score_tmp(loca_plus))) + abs(sum(score_tmp(loca_minus))));
%         end
%         l_ori = (length(loca_minus) + length(loca_plus)) / (length(t_1)*length(t_2)/2);
%         %if l_ori > thres_L        
%             [h1,p1] = ztest(s_ori, mean(boot_tmp(:,1)),std(boot_tmp(:,1)),'Tail','right');
%             %[h1_pre,p1_pre] = ztest(1+1e-2, mean(boot_tmp(:,1)),std(boot_tmp(:,1)),'Tail','right');
%             [h2,p2] = ztest(s_ori, mean(boot_tmp(:,2)),std(boot_tmp(:,2)),'Tail','right');
%             %[h2_pre,p2_pre] = ztest(1+1e-2, mean(boot_tmp(:,2)),std(boot_tmp(:,2)),'Tail','right');
%             [h3,p3] = ztest(s_ori, mean(boot_tmp(:,3)),std(boot_tmp(:,3)),'Tail','right');
%             
%             %p_total = [p_total ; [p1-p1_pre,p2-p2_pre]];
%             p_total = [p_total ; [p1,p2,p3]];
%         %end
%     end
%     if ~isempty(p_total)
%         %p_total(p_total == NaN) = 0.5;
%         p_tmp_1 = nonzeros(rmmissing(p_total(:,1)));
%         p_tmp_2 = nonzeros(rmmissing(p_total(:,2)));
%         p_tmp_3 = nonzeros(rmmissing(p_total(:,3)));
%         
%         fisher_tmp_1 = 2* sum(-log(p_tmp_1));
%         fisher_tmp_2 = 2* sum(-log(p_tmp_2));
%         fisher_tmp_3 = 2* sum(-log(p_tmp_3));
%         
%         num_p_1 = length(p_tmp_1);
%         num_p_2 = length(p_tmp_2);
%         num_p_3 = length(p_tmp_3);
%         
%         fisher_thres_1 = chi2cdf(-2*log(thres_boot)*num_p_1, 2*num_p_1, 'upper');
%         fisher_thres_2 = chi2cdf(-2*log(thres_boot)*num_p_2, 2*num_p_2, 'upper');
%         fisher_thres_3 = chi2cdf(-2*log(thres_boot)*num_p_3, 2*num_p_3, 'upper');
%         
%         fisher_tmp = [chi2cdf(fisher_tmp_1, 2*num_p_1, 'upper'), chi2cdf(fisher_tmp_2, 2*num_p_2, 'upper'), chi2cdf(fisher_tmp_3, 2*num_p_3, 'upper'), fisher_thres_1,fisher_thres_2,fisher_thres_3];
%     else
%         fisher_tmp = [0,0,0,0];
%     end
%     
%     fisher_total = [fisher_total; [fisher_tmp]];
% end
% fisher_thres = chi2cdf(-2*log(thres_boot)*num_data, 2*num_data, 'upper')

