clc;
clear;
close all;

%% parameter
period = 10;
num_component = 3;
thres_noise = 0;
dimension = 2;
length_data = 100;
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
    %filename = ['IFL_timeseries_fit_',num2str(noise_percent)];
    filename = ['CFL_timeseries_fit_',num2str(noise_percent)];
    %filename = ['SFL_timeseries_fit_',num2str(noise_percent)];
    load(filename)
    
    %% calculate score
    S_total = zeros(length(component_list(:,1)),2^dimension,length_data);
    L_total = zeros(length(component_list(:,1)),2^dimension,length_data);

    for i = 1:length_data

        y_target = cell2mat(y_total(i));
        t_target = t;

        S_tmp = zeros(length(component_list(:,1)),2^dimension);
        L_tmp = zeros(length(component_list(:,1)),2^dimension);
        for j = 1:length(component_list(:,1))
            
            % calculate regulation detection function
            st1 = component_list(j,1);
            st2 = component_list(j,2);
            ed = component_list(j,3);
            
            [score_list, t_1, t_2] = RDS_dim2(y_target(:,st1), y_target(:,st2), y_target(:,ed), t, time_interval);
        
            % calculate S & R
            s_tmp = zeros(1,4);
            l_tmp = zeros(1,4);
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

                s_tmp(k) = s;
                l_tmp(k) = l;
            end
            S_tmp(j,:) = s_tmp;
            L_tmp(j,:) = l_tmp;
        end
        
        S_total(:,:,i) = S_tmp;
        L_total(:,:,i) = L_tmp;
        
    end
    %filename = ['IFL_result_dim2_',num2str(noise_percent)];
    filename = ['CFL_result_dim2_',num2str(noise_percent)];
    %filename = ['SFL_result_dim2_',num2str(noise_percent)];
    save(filename,'S_total','L_total','component_list','dimension','num_component','period','length_data')
end