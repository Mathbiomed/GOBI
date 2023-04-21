clc;
clear;
close all;
addpath('GOBI') 

%% load data
load('data_with_options.mat')

%% parameters
thres_noise = 0; % threshold for regulation-detection function, 0 is default
dimension = 1;   % dimension of the framework

%% all the pairs of 1D regulation (cause, target)

component = [1:num_component];
component_list_dim1_tmp = nchoosek(component, 1);
component_list_dim1 = [];
for i = 1:length(component_list_dim1_tmp)
    for j = 1:num_component
        if i == j
            if isnan(type_self)
                component_list_dim1 = [component_list_dim1 ; [component_list_dim1_tmp(i), j]];
            else
                continue
            end
        else
            component_list_dim1 = [component_list_dim1 ; [component_list_dim1_tmp(i), j]];
        end
    end
end

num_pair = length(component_list_dim1(:,1));
num_type = 2.^dimension;

%% start parallel pool
parpool threads;
clear completedJobs;
dq = parallel.pool.DataQueue;
wb = waitbar(0,'Processing');
N = num_data;
Listener = afterEach(dq, @(varargin) waitbar((completedJobs/N),wb,sprintf('Completed: %d', completedJobs(1))));

%% from all data, calculate regulation detection score & region
disp('calculate regulation detection score...')
S_total_list = zeros(num_pair,num_type,num_data); % save regulation-detection score for all data
L_total_list = zeros(num_pair,num_type,num_data); % save regulation-detection region for all data
parfor i = 1:num_data
    send(dq,i)
    y_target = cell2mat(y_total(i));
    t_target = t;
    
    S_total = zeros(num_pair,num_type); % save regulation-detection score for each data
    L_total = zeros(num_pair,num_type); % save regulation-detection region for eacg data
    
    for j = 1:length(component_list_dim1(:,1))
        st = component_list_dim1(j,1); % index for cause
        ed = component_list_dim1(j,2); % index for target
        
        % compute regulation-detection function
        if type_self == -1
            [score_list, t_1, t_2] = RDS_ns_dim1(y_target(:,st), y_target(:,ed), t_target, time_interval);
        elseif type_self == 1
            [score_list, t_1, t_2] = RDS_ps_dim1(y_target(:,st), y_target(:,ed), t_target, time_interval);
        else
            [score_list, t_1, t_2] = RDS_dim1(y_target(:,st), y_target(:,ed), t_target, time_interval);
        end
        
        % compute regulation-detection score
        for k = 1:num_type
            score = reshape(score_list(:,:,k),[length(t_target),length(t_target)]);
            loca_plus = find(score > thres_noise);
            loca_minus = find(score < -thres_noise);
            if isempty(loca_plus) && isempty(loca_minus)
                s = 1;    
            else
                s = (sum(score(loca_plus)) + sum(score(loca_minus)))/ (abs(sum(score(loca_plus))) + abs(sum(score(loca_minus))));
            end
            l = (length(loca_minus) + length(loca_plus)) / (length(t_1)*length(t_2)/2);
            S_total(j,k) = s;
            L_total(j,k) = l;
        end
    end
    S_total_list(:,:,i) = S_total;
    L_total_list(:,:,i) = L_total;    
end

filename = 'RDS_dim1';
save(filename, 'S_total_list', 'L_total_list','component_list_dim1', 'num_component', 'num_pair','num_type','num_data')

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

