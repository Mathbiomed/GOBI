clc;
clear;
close all;
addpath('GOBI')

%% load data
load('data_with_options')
load('RDS_dim2')
load('Delta_dim2')

%% parameters
if thres_R == 0
    thres_R = 0.05;
end
thres_noise = 0;

%% find candidate for surrogate test (find the regulation passing the delta test)
boot_candidate_list = [];
boot_type_list = [];

for i = 1:num_candidate_delta
    if isnan(delta_list(i,1))
        delta_list(i,1) = 0;
    end
    if isnan(delta_list(i,2))
        delta_list(i,2) = 0;
    end
    if delta_list(i,1) <= p_delta && delta_list(i,2) <= p_delta
        boot_candidate_list = [boot_candidate_list ; delta_candidate_list(i,:)];
        boot_type_list = [boot_type_list; delta_type_list(i)];
    end
end
num_candidate_boot = length(boot_candidate_list(:,1)); % candidate for the surrogate test (passing the delta test)

%% surrogate test
num_boot = 100;
surrogate_list = [];
for i = 1:num_candidate_boot
    disp(i)
    st1 = boot_candidate_list(i,1); % first causal variable
    st2 = boot_candidate_list(i,2); % second causal variable
    ed = boot_candidate_list(i,3);  % target variable
    type_tmp = boot_type_list(i);   
    p_total = [];
    for j = 1:num_data
        y_tmp = cell2mat(y_total(j));    
        C1 = y_tmp(:,st1); % import time series of the first variable
        C2 = y_tmp(:,st2); % import time series of the second variable
        T = y_tmp(:,ed);   % import time series of the target variable
        t_target = t(1:length(y_tmp(:,1)));
        
        boot_tmp = [];
        for k = 1:num_boot
            % bootstrapping of regulation-detection score with shuffled time series of cause 1
            C1_shuffled = C1(randperm(length(C1))); % shuffle the time seires of the first causal variable
            
            % compute the regulation-detection score using shuffled time series
            if type_self == -1
                [score_list, t_1, t_2] = RDS_ns_dim2(C1_shuffled, C2, T, t_target, time_interval);
            elseif type_self == 1
                [score_list, t_1, t_2] = RDS_ps_dim2(C1_shuffled, C2, T, t_target, time_interval);
            else
                [score_list, t_1, t_2] = RDS_dim2(C1_shuffled, C2, T, t_target, time_interval);
            end
            
            score_tmp = reshape(score_list(:,:,type_tmp),[length(t_1(:,1)),length(t_1(1,:))]);
            
            loca_plus = find(score_tmp > thres_noise);
            loca_minus = find(score_tmp < -thres_noise);
            if isempty(loca_plus) && isempty(loca_minus)
                s_tmp_1 = 1;
            else
                s_tmp_1 = (sum(score_tmp(loca_plus)) + sum(score_tmp(loca_minus)))/ (abs(sum(score_tmp(loca_plus))) + abs(sum(score_tmp(loca_minus))));
            end
            
            % bootstrapping of regulation-detection score with shuffled time series of cause 1
            C2_shuffled = C2(randperm(length(C2))); % shuffle the time seires of the second causal variable
            
            % compute the regulation-detection score using shuffled time series
            if type_self == -1
                [score_list, t_1, t_2] = RDS_ns_dim2(C1, C2_shuffled, T, t_target, time_interval);
            elseif type_self == 1
                [score_list, t_1, t_2] = RDS_ps_dim2(C1, C2_shuffled, T, t_target, time_interval);
            else
                [score_list, t_1, t_2] = RDS_dim2(C1, C2_shuffled, T, t_target, time_interval);
            end
            
            score_tmp = reshape(score_list(:,:,type_tmp),[length(t_1(:,1)),length(t_1(1,:))]);
            
            loca_plus = find(score_tmp > thres_noise);
            loca_minus = find(score_tmp < -thres_noise);
            if isempty(loca_plus) && isempty(loca_minus)
                s_tmp_2 = 1;
            else
                s_tmp_2 = (sum(score_tmp(loca_plus)) + sum(score_tmp(loca_minus)))/ (abs(sum(score_tmp(loca_plus))) + abs(sum(score_tmp(loca_minus))));
            end
            
            boot_tmp = [boot_tmp ;[s_tmp_1,s_tmp_2]]; % save the RDS from shuffled time series
        end
        
        % compute the regulation-detection score with original time series
        if type_self == -1
            [score_list, t_1, t_2] = RDS_ns_dim2(C1, C2, T, t_target, time_interval);
        elseif type_self == 1
            [score_list, t_1, t_2] = RDS_ps_dim2(C1, C2, T, t_target, time_interval);
        else
            [score_list, t_1, t_2] = RDS_dim2(C1, C2, T, t_target, time_interval);
        end
        
        score_tmp = reshape(score_list(:,:,type_tmp),[length(t_1(:,1)),length(t_1(1,:))]);

        loca_plus = find(score_tmp > thres_noise);
        loca_minus = find(score_tmp < -thres_noise);
        if isempty(loca_plus) && isempty(loca_minus)
            s_ori = 1;
        else
            s_ori = (sum(score_tmp(loca_plus)) + sum(score_tmp(loca_minus)))/ (abs(sum(score_tmp(loca_plus))) + abs(sum(score_tmp(loca_minus))));
        end
        r_ori = (length(loca_minus) + length(loca_plus)) / (length(t_1)*length(t_2)/2);
     
        % using one sided Z test to compute p value (z score)
        [h1,p1] = ztest(s_ori, mean(boot_tmp(:,1)),std(boot_tmp(:,1)),'Tail','right');
        [h2,p2] = ztest(s_ori, mean(boot_tmp(:,2)),std(boot_tmp(:,2)),'Tail','right');
        p_total = [p_total ; [p1,p2]];
    end
    if ~isempty(p_total)
        p_tmp_1 = nonzeros(rmmissing(p_total(:,1)));
        p_tmp_2 = nonzeros(rmmissing(p_total(:,2)));
        
        % combine p-value using Fisher's method
        fisher_tmp_1 = 2* sum(-log(p_tmp_1));
        fisher_tmp_2 = 2* sum(-log(p_tmp_2));
        num_p_1 = length(p_tmp_1);
        num_p_2 = length(p_tmp_2);
        
        % combine p_surrogate to define the threshold for the surrogate test
        fisher_thres_1 = chi2cdf(-2*log(p_surrogate)*num_p_1, 2*num_p_1, 'upper'); % threshold of combined p-value for cause1
        fisher_thres_2 = chi2cdf(-2*log(p_surrogate)*num_p_2, 2*num_p_2, 'upper'); % threshold of combined p-value for cause2
        fisher_tmp = [chi2cdf(fisher_tmp_1, 2*num_p_1, 'upper'), chi2cdf(fisher_tmp_2, 2*num_p_2, 'upper'),fisher_thres_1,fisher_thres_2];
    else
        fisher_tmp = [0,0,0,0];
    end
    surrogate_list = [surrogate_list; [fisher_tmp]];
end

% save the results of surrogate test
filename = ['Surrogate_dim2'];
save(filename, 'boot_candidate_list','boot_type_list','component_list_dim2','num_data','surrogate_list', 'delta_list', 'delta_type_list', 'delta_candidate_list')

