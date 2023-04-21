clc;
clear;
close all;
addpath('../../GOBI')
%% start parallel pool
parpool threads;
clear completedJobs;
dq = parallel.pool.DataQueue;
wb = waitbar(0,'Processing');
N = 10*100*11;
Listener = afterEach(dq, @(varargin) waitbar((completedJobs/N),wb,sprintf('Completed: %d', completedJobs(1))));

trial_list = [0:9];
for trial = trial_list
    %% parameters
    %trial = 1;
    dimension = 1;
    thres_noise = 1e-5;
    load('KF_timeseries_Trial1') % for import parameters
    noise_list = [0:2:20];

    %% Create pairs of component for one dimensional regulation
    component = [1:num_component];
    component_list_dim1_tmp = nchoosek(component, 2);
    component_list_dim1 = [];
    for i = 1:length(component_list_dim1_tmp(:,1))
        component_list_dim1 = [component_list_dim1 ; [component_list_dim1_tmp(i,1), component_list_dim1_tmp(i,2)]];
        component_list_dim1 = [component_list_dim1 ; [component_list_dim1_tmp(i,2), component_list_dim1_tmp(i,1)]];
    end
    num_pair = length(component_list_dim1(:,1));
    num_type = 2^dimension;

    %% Start inference framework
    for noise_percent = noise_list
        disp(noise_percent)

        %% load data
        filename = ['KF_timeseries_fit_',num2str(noise_percent),'_Trial',num2str(trial)];
        load(filename)

        %% Calculate regulation-detection-score (S) & size of the regulation-detection region (L)
        S_total = zeros(length(component_list_dim1(:,1)),2^dimension,num_data);
        L_total = zeros(length(component_list_dim1(:,1)),2^dimension,num_data);

        parfor j = 1:length(y_total)
            send(dq, j)

            y_target = cell2mat(y_total(j));
            S = zeros(num_pair,num_type);
            L = zeros(num_pair,num_type);

            for i = 1:num_pair
                % calculate regulation detection function
                st = component_list_dim1(i,1);
                ed = component_list_dim1(i,2);
                [score_list, t_1, t_2] = RDS_ns_dim1(y_target(:,st), y_target(:,ed), t, time_interval);

                % calculate S & R
                s_tmp = zeros(1,num_type);
                l_tmp = zeros(1,num_type);
                for k = 1:num_type
                    score = reshape(score_list(:,:,k),[length(t),length(t)]);
                    loca_plus = find(score > thres_noise);
                    loca_minus = find(score < -thres_noise);
                    if isempty(loca_plus) && isempty(loca_minus)
                        s = 1;
                    else
                        s = (sum(score(loca_plus)) + sum(score(loca_minus)))/ (abs(sum(score(loca_plus))) + abs(sum(score(loca_minus))));
                    end

                    l = (length(loca_minus) + length(loca_plus)) / (length(t_1)*length(t_2)/2);

                    s_tmp(k) = s;
                    l_tmp(k) = l;
                end

                % save as list
                S(i,:) = s_tmp;
                L(i,:) = l_tmp;
            end

            S_total(:,:,j) = S;
            L_total(:,:,j) = L;
        end

        %% Save results
        filename = ['KF_results_dim1_',num2str(noise_percent),'_Trial',num2str(trial)];
        save(filename, 'S_total', 'L_total','component_list_dim1','num_data','num_type','num_pair','dimension')
    end
end
%% end parallel pool
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