clc;
clear;
close all;

%% parameters
trial_list = [0:9];
noise_list = [2:2:20];
%noise_list = [2];

dimension = 2;
thres_noise = 1e-7;
system_list = {'cAMP','Fr','Gb','GW','KF'};

noise_type_list = {'additive','blue','brown','pink','purple','dynamical','multiplicative'};
noise_type_list = {'additive','blue','pink','multiplicative'};
%noise_type_list = {'additive'};

system_list = {'cAMP'};

p_thres = 0.01;
delta_pair = [
    [2,3];
    [1,4];
    [1,4];
    [2,3]];
t_L = 0.01;

for noise_type_idx = 1:length(noise_type_list)
    noise_type = char(noise_type_list(noise_type_idx));
    disp(noise_type)
    for system_idx = 1:length(system_list)
        system = char(system_list(system_idx));
        disp(system)
        for trial = trial_list
            for noise_percent = noise_list
                %% load data
                filename = ['./RDS_dim2_',noise_type,'/',system,'_results_dim2_',num2str(noise_percent),'_Trial',num2str(trial)];
                load(filename)
                filename = ['./TRS_dim2_',noise_type,'/',system,'_TRS_dim2_',num2str(noise_percent),'_Trial',num2str(trial)];
                load(filename)
                %% find regulation pair which pass the criteria of TRS (idx_pair, idx_type)
                candidate = [];
                for i = 1:num_pair
                    for j = 1:num_type
                        if regulation_2dim(i,j) == 1
                            candidate = [candidate;[i,j]];
                        end
                    end
                end
                if isempty(candidate)
                    num_candidate = 0;
                else
                    num_candidate = length(candidate(:,1));
                end
                
                %% create delta
                delta_list = zeros(num_data,dimension,num_candidate);
                for i = 1:num_candidate
                    % import score and region
                    S_ori_tmp = reshape(S_total(candidate(i,1),candidate(i,2),:),[num_data,1]);
                    S_delta1_tmp = reshape(S_total(candidate(i,1),delta_pair(candidate(i,2),1),:),[num_data,1]);
                    S_delta2_tmp = reshape(S_total(candidate(i,1),delta_pair(candidate(i,2),2),:),[num_data,1]);

                    L_ori_tmp = reshape(L_total(candidate(i,1),candidate(i,2),:),[num_data,1]);
                    L_delta1_tmp = reshape(L_total(candidate(i,1),delta_pair(candidate(i,2),1),:),[num_data,1]);
                    L_delta2_tmp = reshape(L_total(candidate(i,1),delta_pair(candidate(i,2),2),:),[num_data,1]);

                    % remove empty region
                    L_ori_tmp = L_threshold(L_ori_tmp, t_L);
                    L_delta1_tmp = L_threshold(L_delta1_tmp, t_L);
                    L_delta2_tmp = L_threshold(L_delta2_tmp, t_L);

                    S_ori_tmp = S_ori_tmp .* L_ori_tmp;
                    S_delta1_tmp = S_delta1_tmp .* L_delta1_tmp;
                    S_delta2_tmp = S_delta2_tmp .* L_delta2_tmp;

                    S_ori_tmp(S_ori_tmp == 0) = NaN;
                    S_delta1_tmp(S_delta1_tmp == 0) = NaN;
                    S_delta2_tmp(S_delta2_tmp == 0) = NaN;

                    % calculate delta
                    delta_list_tmp = [S_ori_tmp - S_delta1_tmp, S_ori_tmp - S_delta2_tmp];
                    delta_list(:,:,i) = delta_list_tmp;
                end

                %% calculte p-value using Wilcoxon sign rank
                p_list = [];
                for i = 1:num_candidate
                    delta_list_tmp = reshape(delta_list(:,:,i),[num_data,dimension]);
                    p1 = 0;
                    p2 = 0;
                    delta_tmp1 = rmmissing(delta_list_tmp(:,1));
                    delta_tmp2 = rmmissing(delta_list_tmp(:,2));
                    if ~isempty(delta_tmp1)
                        [p1,h1,stats1] = signrank(delta_list_tmp(:,1),-1e-3,'tail','right');
                    end
                    if ~isempty(delta_tmp2)
                        [p2,h2,stats2] = signrank(delta_list_tmp(:,2),-1e-3,'tail','right');
                    end
                    p_list = [p_list;[p1,p2]];
                end
                %p_list
                %% infer 2D regulation using the p-value
                delta_2dim = regulation_2dim;

                for i = 1:num_candidate
                    if p_list(i,1) > p_thres || p_list(i,2) > p_thres
                        delta_2dim(candidate(i,1),candidate(i,2)) = 0;
                    end
                end
                %[delta_2dim(3,1),delta_2dim(20,1),delta_2dim(29,1)]
                
                filename = ['./Delta_dim2_',noise_type,'/',system,'_Delta_dim2_',num2str(noise_percent),'_Trial',num2str(trial)];

                save(filename, 'delta_2dim', 'delta_list','component_list_dim2','num_data','num_type','num_pair','dimension')

            end
        end
    end
end
