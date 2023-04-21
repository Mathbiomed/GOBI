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
S_plot_0 = [];
S_plot_1 = [];
S_plot_2 = [];
S_plot_3 = [];
S_plot_4 = [];
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
        S_processed_0 = S_tmp_0 .* L_processed_0;
        S_processed_0(S_processed_0 == 0) = NaN;
        S_plot_0 = [S_plot_0 , S_processed_0];
        
        L_processed_1 = L_threshold(L_tmp_1, thres_L);
        S_processed_1 = S_tmp_1 .* L_processed_1;
        S_processed_1(S_processed_1 == 0) = NaN;
        S_plot_1 = [S_plot_1 , S_processed_1];
        
        L_processed_2 = L_threshold(L_tmp_2, thres_L);
        S_processed_2 = S_tmp_2 .* L_processed_2;
        S_processed_2(S_processed_2 == 0) = NaN;
        S_plot_2 = [S_plot_2 , S_processed_2];
        
        L_processed_3 = L_threshold(L_tmp_3, thres_L);
        S_processed_3 = S_tmp_3 .* L_processed_3;
        S_processed_3(S_processed_3 == 0) = NaN;
        S_plot_3 = [S_plot_3 , S_processed_3];
        
        L_processed_4 = L_threshold(L_tmp_4, thres_L);
        S_processed_4 = S_tmp_4 .* L_processed_4;
        S_processed_4(S_processed_4 == 0) = NaN;
        S_plot_4 = [S_plot_4 , S_processed_4];
    end
end

S_plot_total = [
    S_plot_0(:,false_index);
    S_plot_0(:,true_index);
    S_plot_1(:,false_index);
    S_plot_1(:,true_index);
    S_plot_2(:,false_index);
    S_plot_2(:,true_index);
    S_plot_3(:,false_index);
    S_plot_3(:,true_index);
    S_plot_4(:,false_index);
    S_plot_4(:,true_index);];
index_tf = [
    ones(num_data,1);
    zeros(num_data,1);
    ones(num_data,1);
    zeros(num_data,1);
    ones(num_data,1);
    zeros(num_data,1);
    ones(num_data,1);
    zeros(num_data,1);
    ones(num_data,1);
    zeros(num_data,1)];
index_noise = [
    ones(num_data,1) * 0;
    ones(num_data,1) * 0;
    ones(num_data,1) * 2;
    ones(num_data,1) * 2;
    ones(num_data,1) * 4;
    ones(num_data,1) * 4;
    ones(num_data,1) * 6;
    ones(num_data,1) * 6;
    ones(num_data,1) * 8;
    ones(num_data,1) * 8];

plot_total = array2table([index_noise,S_plot_total, index_tf]);
plot_total.Properties.VariableNames = {'noise' 'score' 'tf'};
%% plot

figure(1)
b = boxchart(plot_total.noise, plot_total.score, 'GroupByColor', plot_total.tf);
b(1).BoxWidth = 1;
b(1).BoxFaceColor = [0,0,1];
b(1).LineWidth = 1;
b(1).BoxFaceAlpha = 0.1;
b(1).MarkerStyle = 'none';

b(2).BoxWidth = 1;
b(2).BoxFaceColor = [1,0,0];
b(2).LineWidth = 1;
b(2).BoxFaceAlpha = 0.1;
b(2).MarkerStyle = 'none';

xlim([-1,9])
ylim([-1,1])
xticks([0,2,4,6,8])
yticks([-1,0,1])
xticklabels([0,5,10,15,20])
xlabel('Multiplicative Noise (%)')
ylabel('S_{\sigma}')
set(gca, 'FontSize',14)


