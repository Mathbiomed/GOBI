clc;
clear;
close all;

%% load data
type = 'full';
%type = 'half';
%type = 'quarter';

filename = ['RDS_dim1_',type];
load(filename)

noise_percent = 8;

%% choose threshold using guide
t_S = 0.9 - 0.005 * noise_percent;
t_L = 0.05;
t_TRS = 0.9 - 0.01 * noise_percent;

%% Calculate Total Regulation Score (TRS)
TRS_dim1 = zeros(num_pair,2^dimension);
for i = 1:num_pair
    for j = 1:num_type
        S_tmp = reshape(S_total(i,j,:),[num_data,1]);
        L_tmp = reshape(L_total(i,j,:),[num_data,1]);
        S_processed = S_threshold(S_tmp, t_S);
        L_processed = L_threshold(L_tmp, t_L);
        N = sum(L_processed);
        TRS_dim1(i,j) = sum(S_processed.*L_processed)/N;
    end
end
regulation_1dim = zeros(num_pair,2^dimension);
for i = 1:num_pair
    for j = 1:num_type
        if TRS_dim1(i,j) >= t_TRS
            regulation_1dim(i,j) = 1;
        end
    end
end
regulation_1dim
    
filename = ['TRS_1D_',type];
save(filename, 'regulation_1dim','TRS_dim1', 'component_list','num_type','num_pair','dimension')
                
