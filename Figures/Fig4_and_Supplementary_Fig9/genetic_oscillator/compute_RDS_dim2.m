clc;
clear;
close all;
addpath('../../GOBI')

%% import data
load('data_final_gene.mat')

%% parameter
dimension = 2;
num_component = 2;
thres_noise = 1e-5;

%% all the pairs of components for dim1
component_list_dim2 = [[1,2,1];[1,2,2]];

num_pair = length(component_list_dim2(:,1));
num_data = length(y_total);

%% for every data, calculate regulation detection score
disp('calculate regulation detection score...')
S_total_list = zeros(num_pair,4,num_data);
L_total_list = zeros(num_pair,4,num_data);
for i = 1:num_data
    y_target = cell2mat(y_total(i));
    t_target = linspace(0,1,1/time_interval+1).';
    
    S_total = zeros(num_pair,4);
    L_total = zeros(num_pair,4);
    for j = 1:length(component_list_dim2(:,1))
        st1 = component_list_dim2(j,1);
        st2 = component_list_dim2(j,2);
        ed = component_list_dim2(j,3);
        [score_list, t_1, t_2] = RDS_dim2(y_target(:,st1), y_target(:,st2), y_target(:,ed), t_target, time_interval);
        
        for k = 1:4
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
    %save s,r,l at the list
    S_total_list(:,:,i) = S_total;
    L_total_list(:,:,i) = L_total;    
end

component_list = component_list_dim2;
filename = 'RDS_dim2';
save(filename, 'S_total_list', 'L_total_list','component_list', 'num_component', 'num_pair', 'num_data')
