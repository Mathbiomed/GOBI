clc;
clear;
close all;
addpath('../GOBI') 

%% load data
num_data = 5;
num_component = 3;

filename = ['RDS_dim2_wo_self_',num2str(num_data)];
load(filename)

%% process
score_list = zeros([num_pair, num_type]);

for i = 1:num_pair
    for j = 1:num_type
        S = reshape(S_total(i,j,:),[num_data,1]);
        L = reshape(L_total(i,j,:),[num_data,1]);
        
        L_processed = L_threshold(L,0.0001);
        S_processed = S_threshold(S,0.9999);
        
        score_list(i,j) = sum(S .* L_processed) / length(find(L_processed == 1));
      
    end
end

for i = 1:num_pair
    for j = 1:num_type
        st1 = component_list_dim2(i,1);
        st2 = component_list_dim2(i,2);
        ed = component_list_dim2(i,3);
        
        if score_list(i,j) > 0.9999
            disp([st1, st2, ed, j]) % [cause1, cause2, target, type]
        end
    end
end
