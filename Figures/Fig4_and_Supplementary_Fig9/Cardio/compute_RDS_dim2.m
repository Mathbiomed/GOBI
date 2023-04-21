clc;
clear;
close all;

addpath('../../GOBI')

%% import data
load('data_cardio.mat')

%% parameter from data
num_component = 5;
length_data = length(y_total);
thres_noise = 1e-5;

%% all the pairs of components for dim1
component = [2:num_component];
component_list_dim2_tmp = nchoosek(component, 2);
component_list_dim2 = [];
for i = 1:length(component_list_dim2_tmp(:,1))
    for j = 1
        if ismember(j, component_list_dim2_tmp(i,:))
            continue
        end
        component_list_dim2 = [component_list_dim2 ; [component_list_dim2_tmp(i,:), j]];
    end
end

num_pair = length(component_list_dim2);

%% for every data, calculate regulation detection score
disp('calculate regulation detection score...')
S_total_list = zeros(num_pair,4,length_data);
L_total_list = zeros(num_pair,4,length_data);
for i = 1:length_data
    y_target = cell2mat(y_total(i));
    t_target = t(1:length(y_target(:,1)));
    
    S_total = zeros(num_pair,4);
    L_total = zeros(num_pair,4);
    for j = 1:length(component_list_dim2(:,1))
        st1 = component_list_dim2(j,1);
        st2 = component_list_dim2(j,2);
        ed = component_list_dim2(j,3);
        [score_list, t_1, t_2] = RDS_ns_dim2(y_target(:,st1), y_target(:,st2), y_target(:,ed), t_target, time_interval);
        
        for k = 1:4
            score = reshape(score_list(:,:,k),[length(t),length(t)]);
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
    %save s,r,l at the list
    S_total_list(:,:,i) = S_total;
    L_total_list(:,:,i) = L_total;    
end

component_list = component_list_dim2;
num_data = length_data;
filename = 'RDS_dim2';
save(filename, 'S_total_list', 'L_total_list','component_list', 'num_component', 'num_pair', 'num_data')

