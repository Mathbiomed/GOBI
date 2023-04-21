clc;
clear;
close all;

%% parameter
system_list = {'KF','Fr','GW','Gb','cAMP'};
num_component_list = [3,3,4,5,7];
dimension = 3; % dimension of framework
S_thres = 0.999;
L_thres = 0.01;
num_type = 2^dimension;
max_name_length = [3,1,3,3,9,6];

%% process
% i indicates the system
% i = 1: Kim-Forger model
% i = 2: Frzilator
% i = 3: Goodwin oscillator
% i = 4: Goldbeter model
% i = 5: cAMP oscillator

for i = 5
    % load data
    system_name = string(system_list(i));
    filename = append('./',system_name,'/',system_name,'_RDS_dim',num2str(dimension));
    load(filename)
    filename = append('components_',system_name);
    load(filename)
    
    % variables
    num_component = num_component_list(i);
    num_pair = length(component_list(:,1));
    
    
    % create the list of causes
    causes = [];
    components_index = [1:num_component];
    
    for j = 1:num_pair
        cause_pair = '(';
        for k = 1:dimension+1
            cause_tmp = string(components(component_list(j,k)));
            cause_tmp_length = strlength(cause_tmp);
            for l = cause_tmp_length+1:max_name_length(i)
                cause_tmp = append(' ',cause_tmp);
            end
            cause_pair = append(cause_pair,cause_tmp);
            if k ~= dimension+1
                cause_pair = append(cause_pair,',');
            end
        end
        cause_pair = append(cause_pair, ')');
        causes = [causes ; cause_pair];
    end
    
    total_list = zeros(num_pair, num_type); % results of regulation-detection scores
    for j = 1:num_pair
        for k = 1:num_type
            S_tmp = reshape(S_total(j,k,:),[num_data,1]);
            L_tmp = reshape(L_total(j,k,:),[num_data,1]);
            
            L_processed = L_threshold(L_tmp, L_thres);
            S_processed = S_threshold(S_tmp, S_thres);
            S_processed = S_processed .* L_processed;
            
            result_tmp = sum(S_processed) / sum(L_processed);
            total_list(j,k) = ceil(result_tmp*100)/100;
        end
    end
end


