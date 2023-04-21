clc;
clear;
close all;
addpath('../GOBI') 

%% parameter
dimension = 1;
thres_L = 0.01;
true_index = 2;
false_index = 5;
%% load data
filename = ['IFL_result_dim1_0'];
load(filename)
S_total_0 = S_total;
L_total_0 = L_total;

filename = ['IFL_result_dim1_5'];
load(filename)
S_total_1 = S_total;
L_total_1 = L_total;

filename = ['IFL_result_dim1_10'];
load(filename)
S_total_2 = S_total;
L_total_2 = L_total;

filename = ['IFL_result_dim1_15'];
load(filename)
S_total_3 = S_total;
L_total_3 = L_total;

filename = ['IFL_result_dim1_20'];
load(filename)
S_total_4 = S_total;
L_total_4 = L_total;

num_pair = length(component_list(:,1));
num_type = 2^dimension;


%% process
x = [1:num_pair * num_type];
S_0 = [];
S_1 = [];
S_2 = [];
S_3 = [];
S_4 = [];
for i = 1:num_pair
    for j = 1:num_type
        S_tmp_0 = reshape(S_total_0(i,j,:),[num_data,1]);
        L_tmp_0 = reshape(L_total_0(i,j,:),[num_data,1]);
        
        S_tmp_1 = reshape(S_total_1(i,j,:),[num_data,1]);
        L_tmp_1 = reshape(L_total_1(i,j,:),[num_data,1]);
        
        S_tmp_2 = reshape(S_total_2(i,j,:),[num_data,1]);
        L_tmp_2 = reshape(L_total_2(i,j,:),[num_data,1]);
        
        S_tmp_3 = reshape(S_total_3(i,j,:),[num_data,1]);
        L_tmp_3 = reshape(L_total_3(i,j,:),[num_data,1]);
        
        S_tmp_4 = reshape(S_total_4(i,j,:),[num_data,1]);
        L_tmp_4 = reshape(L_total_4(i,j,:),[num_data,1]);
        
        L_processed_0 = L_threshold(L_tmp_0, thres_L);
        S_tmp_0 = S_threshold(S_tmp_0, 0.9);
        S_processed_0 = S_tmp_0 .* L_processed_0;
        TRS_0 = sum(S_processed_0) / sum(L_processed_0);
        
        S_0 = [S_0 ; TRS_0];
        
        L_processed_1 = L_threshold(L_tmp_1, thres_L);
        S_tmp_1 = S_threshold(S_tmp_1, 0.875);
        S_processed_1 = S_tmp_1 .* L_processed_1;
        TRS_1 = sum(S_processed_1) / sum(L_processed_1);
        S_1 = [S_1 ; TRS_1];
        
        L_processed_2 = L_threshold(L_tmp_2, thres_L);
        S_tmp_2 = S_threshold(S_tmp_2, 0.85);
        S_processed_2 = S_tmp_2 .* L_processed_2;
        TRS_2 = sum(S_processed_2) / sum(L_processed_2);
        S_2 = [S_2 ; TRS_2];
        
        L_processed_3 = L_threshold(L_tmp_3, thres_L);
        S_tmp_3 = S_threshold(S_tmp_3, 0.825);
        S_processed_3 = S_tmp_3 .* L_processed_3;
        TRS_3 = sum(S_processed_3) / sum(L_processed_3);
        S_3 = [S_3 ; TRS_3];
        
        L_processed_4 = L_threshold(L_tmp_4, thres_L);
        S_tmp_4 = S_threshold(S_tmp_4, 0.8);
        S_processed_4 = S_tmp_4 .* L_processed_4;
        TRS_4 = sum(S_processed_4) / sum(L_processed_4);
        S_4 = [S_4 ; TRS_4];
    end
end
S_plot = [S_0,S_1,S_2,S_3,S_4];

figure(1)
x = [0:5:20];
plot(x,S_plot(1,:),'Color',[1,0,0,0.5],'Marker','o', 'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','m')
hold on
plot(x,S_plot(2,:),'b','Marker','o', 'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','b')
hold on
plot(x,S_plot(5,:),'r','Marker','o', 'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','r')
hold on
plot(x,S_plot(6,:),'Color',[1,0,0,0.5],'Marker','o', 'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','m')
hold on
plot(x,S_plot(9,:),'Color',[1,0,0,0.5],'Marker','o', 'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','m')
hold on
plot(x,S_plot(10,:),'Color',[1,0,0,0.5],'Marker','o', 'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','m')
hold on
plot(x,S_plot(11,:),'Color',[1,0,0,0.5],'Marker','o', 'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','m')
hold on
plot(x,S_plot(12,:),'Color',[1,0,0,0.5],'Marker','o', 'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','m')
hold on
xlim([-2.5,22.5])
ylim([-0.1,1])
xticks([0:5:20])
yticks([0,1])
xticklabels([0,5,10,15,20])
xlabel('Multiplicative Noise (%)')
ylabel('TRS')
set(gca, 'FontSize',14)
