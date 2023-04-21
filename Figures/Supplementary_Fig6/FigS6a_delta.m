clc;
clear;
close all;
addpath('../GOBI')
%% parameter
thres_L = 0.01;
true_index = 2;
false_index = 5;

delta_2 = [];
delta_3 = [];
for i = 2:3
    for j = [0:5:20]
        if i == 2
            filename = ['CFL','_result_dim2_',num2str(j)];
        else
            filename = ['SFL','_result_dim2_',num2str(j)];
        end
        load(filename)

        %% process
        S_tmp = reshape(S_total(1,1,:),[num_data,1]);
        L_tmp = reshape(L_total(1,1,:),[num_data,1]);
        S_tmp_minus = reshape(S_total(1,3,:),[num_data,1]);
        L_tmp_minus = reshape(L_total(1,3,:),[num_data,1]);
        L_processed = L_threshold(L_tmp, thres_L);
        L_processed_minus = L_threshold(L_tmp_minus, thres_L);
        S_processed = S_tmp .* L_processed;
        S_processed_minus = S_tmp_minus .* L_processed_minus;
        S_processed(S_processed == 0) = NaN;
        S_processed_minus(S_processed_minus == 0) = NaN;
        delta = S_processed - S_processed_minus;
        if i == 2
            delta_2 = [delta_2 , delta];
        else
            delta_3 = [delta_3 , delta];
        end
    end
end

% compute p value using Wilcoxon signed rank test
% idx = 1~5 indicate the index for noise level 0~20
idx = 4;
[p,h,stats] = signrank(delta_2(:,idx),0,'tail','left')
[p,h,stats] = signrank(delta_3(:,idx),0,'tail','left')


delta_total = [
    
    delta_3(:,1);
    delta_3(:,2);
    delta_3(:,3);
    delta_3(:,4);
    delta_3(:,5);
    delta_2(:,1);
    delta_2(:,2);
    delta_2(:,3);
    delta_2(:,4);
    delta_2(:,5)];

index_case = [
    ones(num_data,1);
    ones(num_data,1);
    ones(num_data,1);
    ones(num_data,1);
    ones(num_data,1);
    zeros(num_data,1);
    zeros(num_data,1);
    zeros(num_data,1);
    zeros(num_data,1);
    zeros(num_data,1);];
index_noise = [
    ones(num_data,1) * 0;
    ones(num_data,1) * 2;
    ones(num_data,1) * 4;
    ones(num_data,1) * 6;
    ones(num_data,1) * 8;
    ones(num_data,1) * 0;
    ones(num_data,1) * 2;
    ones(num_data,1) * 4;
    ones(num_data,1) * 6;
    ones(num_data,1) * 8;];

plot_total = array2table([index_noise,delta_total, index_case]);
plot_total.Properties.VariableNames = {'noise' 'delta' 'case'};
%% plot

figure(1)
plot([-1,15],[0,0],'Color',[0.5,0.5,0.5],'LineWidth',1)
hold on
b = boxchart(plot_total.noise, plot_total.delta, 'GroupByColor', plot_total.case);
b(1).BoxWidth = 1;
b(1).BoxFaceColor = [0,0,1];
b(1).LineWidth = 1;
b(1).BoxFaceAlpha = 0.1;
b(1).MarkerStyle = 'none';
%b(1).WhiskerLineStyle = 'none';

b(2).BoxWidth = 1;
b(2).BoxFaceColor = [1,0,0];
b(2).LineWidth = 1;
b(2).BoxFaceAlpha = 0.1;
b(2).MarkerStyle = 'none';
%b(2).WhiskerLineStyle = 'none';
%legend
xlim([-1,9])
ylim([-0.06,0.03])
xticks([0,2,4,6,8])
yticks([-0.06,0,0.03])
xticklabels([0,5,10,15,20])
set(gca, 'FontSize',14)
xlabel('Multiplicative Noise (%)')
ylabel('\Delta')
        