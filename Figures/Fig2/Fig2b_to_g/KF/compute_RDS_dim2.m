clc;
clear;
close all;

addpath('../../../GOBI')
addpath('../../../Timeseries_Data/in_silico/KF')

%% parameter
thres_noise = 1e-5;
num_component = 3;
dimension = 2;

%% load data
filename = 'KF_timeseries_Trial0';
load(filename)

%% create pair
component = [1:num_component];
component_list_dim2_tmp = nchoosek(component, dimension);
component_list_dim2 = [];
for i = 1:length(component_list_dim2_tmp(:,1))
    for j = 1:length(component)
        if ~ismember(component(j), component_list_dim2_tmp(i,:))
            component_list_dim2 = [component_list_dim2 ; [component_list_dim2_tmp(i,:),component(j)]];
        end
    end
end
num_pair = length(component_list_dim2(:,1));
num_type = 2^dimension;

%% start parallel pool
parpool threads;
clear completedJobs;
dq = parallel.pool.DataQueue;
wb = waitbar(0,'Processing');
N = num_data;
Listener = afterEach(dq, @(varargin) waitbar((completedJobs/N),wb,sprintf('Completed: %d', completedJobs(1))));


%% calculate score
component_list = component_list_dim2;
S_total = zeros(num_pair,num_type,num_data);
L_total = zeros(num_pair,num_type,num_data);

parfor i = 1:num_data
    send(dq, i)
    y_target = cell2mat(y_total(i));
    t_target = t;

    S_tmp = zeros(num_pair,num_type);
    L_tmp = zeros(num_pair,num_type);
    for j = 1:num_pair
        st1 = component_list(j,1);
        st2 = component_list(j,2);
        ed = component_list(j,3);

        [score_list, t_1, t_2] = RDS_ns_dim2(y_target(:,st1), y_target(:,st2), y_target(:,ed), t_target, time_interval);
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
            S_tmp(j,k) = s;
            L_tmp(j,k) = l;
        end

    end
    S_total(:,:,i) = S_tmp;
    L_total(:,:,i) = L_tmp;  
end

filename = 'KF_RDS_dim2';
save(filename,'S_total','L_total','component_list','num_data','component_list')

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