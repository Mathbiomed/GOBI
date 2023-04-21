clc;
clear;
close all;

addpath('../../GOBI')
%% load data
load('data_final_prey_predator.mat')

%% parameters
thres_noise = 0;
%% all the pairs of components for dim1
component_list_dim1= [[1,2];[2,1]];
num_pair = length(component_list_dim1(:,1));

%% for every data, calculate regulation detection score
disp('calculate regulation detection score...')
S_total_list = zeros(num_pair,2,num_data);
L_total_list = zeros(num_pair,2,num_data);
for i = 1:num_data
    y_target = cell2mat(y_total(i));
    t_target = t;
    
    S_total = zeros(num_pair,2);
    L_total = zeros(num_pair,2);
    for j = 1:num_pair
        st = component_list_dim1(j,1);
        ed = component_list_dim1(j,2);
        [score_list, t_1, t_2] = RDS_dim1(y_target(:,st), y_target(:,ed), t_target, time_interval);
        
       for k = 1:2
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

component_list = component_list_dim1;
filename = 'RDS_dim1';
save(filename, 'S_total_list', 'L_total_list','component_list', 'num_component', 'num_pair', 'num_data')

