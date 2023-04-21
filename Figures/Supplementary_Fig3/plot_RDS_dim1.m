clc;
clear;
close all;
addpath('../GOBI') 

%% load data
num_data = 10;
num_component = 3;

filename = ['RDS_dim1_',num2str(num_data)];

load(filename)

%% process
% row: cause / column: target
acti_list = zeros(num_component);
repp_list = zeros(num_component);
for i = 1:num_pair
    for j = 1:num_type
        S = reshape(S_total(i,j,:),[num_data,1]);
        L = reshape(L_total(i,j,:),[num_data,1]);
        
        L_processed = L_threshold(L,0.0001);
        S_processed = S_threshold(S,0.9999);
        
        st = component_list_dim1(i,1);
        ed = component_list_dim1(i,2);
        
        if j ==1
            acti_list(st,ed) = sum(S .* L_processed) / length(find(L_processed == 1));
        else
            repp_list(st,ed) = sum(S .* L_processed) / length(find(L_processed == 1));
        end
    end
end

for i = 1:num_component
    for j = 1:num_component
        if acti_list(i,j) > 0.9999
            acti_list(i,j) = 1;
        else
            acti_list(i,j) = 0;
        end
        if repp_list(i,j) > 0.9999
            repp_list(i,j) = 1;
        else
            repp_list(i,j) = 0;
        end
    end
end

acti_list % matrix of inferred positive regulations 
repp_list % matrix of inferred negative regulations