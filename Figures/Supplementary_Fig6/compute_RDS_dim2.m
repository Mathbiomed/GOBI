clc;
clear;
close all;
addpath('../GOBI')
%% parameter
period = 10;
num_component = 3;
thres_noise = 1e-5;
dimension = 2;
num_data = 100;
noise_list = [0:5:20];

%% create pair
component = [1:num_component];
component_list_tmp = nchoosek(component, dimension+1);
component_list = [];
for i = 1:length(component_list_tmp(:,1))
    component_list = [component_list ; [component_list_tmp(i,1), component_list_tmp(i,2), component_list_tmp(i,3)]];
    component_list = [component_list ; [component_list_tmp(i,3), component_list_tmp(i,1), component_list_tmp(i,2)]];
    component_list = [component_list ; [component_list_tmp(i,2), component_list_tmp(i,3), component_list_tmp(i,1)]];
end

%% compute RDS
for l = 1:length(noise_list)
    noise_percent = noise_list(l);
    %% load data
    filename = ['SFL_timeseries_fit_',num2str(noise_percent)];
    load(filename)
    
    %% calculate score
    S_total = zeros(length(component_list(:,1)),2^dimension,num_data);
    L_total = zeros(length(component_list(:,1)),2^dimension,num_data);

    for i = 1:num_data
        y_target = cell2mat(y_total(i));
        t_target = t;

        for j = 1:length(component_list(:,1))
            st1 = component_list(j,1);
            st2 = component_list(j,2);
            ed = component_list(j,3);
            
            [score_list, t_1, t_2] = RDS_dim2(y_target(:,st1),y_target(:,st2), y_target(:,ed), t_target, time_interval);
            for k = 1:2^dimension
                score = reshape(score_list(:,:,k),[length(t_target),length(t_target)]);
                loca_plus = find(score > thres_noise);
                loca_minus = find(score < -thres_noise);
                if isempty(loca_plus) && isempty(loca_minus)
                    s = 1;    
                else
                    s = (sum(score(loca_plus)) + sum(score(loca_minus)))/ (abs(sum(score(loca_plus))) + abs(sum(score(loca_minus))));
                end
                l = (length(loca_minus) + length(loca_plus)) / (length(t_1)*length(t_2)/2);

                S_total(j,k,i) = s;
                L_total(j,k,i) = l;
            end
        end
    end
    filename = ['SFL_result_dim2_',num2str(noise_percent)];
    save(filename,'S_total','L_total','component_list','dimension','num_component','period','num_data')
end

for l = 1:length(noise_list)
    noise_percent = noise_list(l);
    %% load data
    filename = ['CFL_timeseries_fit_',num2str(noise_percent)];
    load(filename)
    
    %% calculate score
    S_total = zeros(length(component_list(:,1)),2^dimension,num_data);
    L_total = zeros(length(component_list(:,1)),2^dimension,num_data);

    for i = 1:num_data
        y_target = cell2mat(y_total(i));
        t_target = t;

        for j = 1:length(component_list(:,1))
            st1 = component_list(j,1);
            st2 = component_list(j,2);
            ed = component_list(j,3);

            [score_list, t_1, t_2] = RDS_dim2(y_target(:,st1),y_target(:,st2), y_target(:,ed), t_target, time_interval);
            for k = 1:2^dimension
                score = reshape(score_list(:,:,k),[length(t_target),length(t_target)]);
                loca_plus = find(score > thres_noise);
                loca_minus = find(score < -thres_noise);
                if isempty(loca_plus) && isempty(loca_minus)
                    s = 1;    
                else
                    s = (sum(score(loca_plus)) + sum(score(loca_minus)))/ (abs(sum(score(loca_plus))) + abs(sum(score(loca_minus))));
                end
                l = (length(loca_minus) + length(loca_plus)) / (length(t_1)*length(t_2)/2);

                S_total(j,k,i) = s;
                L_total(j,k,i) = l;
            end
        end
    end
    filename = ['CFL_result_dim2_',num2str(noise_percent)];
    save(filename,'S_total','L_total','component_list','dimension','num_component','period','num_data')
end