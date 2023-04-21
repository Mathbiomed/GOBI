clc;
clear;
close all;
addpath('../GOBI')
%% parameter
noise_list = [10,15,20];
num_boot = 500;
num_data = 100;
%% calculate p-value
p_list_CFL = zeros(length(noise_list),num_data);
p_list_SFL = zeros(length(noise_list),num_data);

for noise_idx = 1:length(noise_list)
    
    noise_percent = noise_list(noise_idx);
    
    % load data
    filename = ['CFL_bootstrapping_noise_',num2str(noise_percent)];
    load(filename)
    S_boot_CFL = S_boot;
    L_boot_CFL = L_boot;
    
    filename = ['CFL_result_dim2_',num2str(noise_percent)];
    load(filename)
    S_ori_CFL_total = S_total;
    L_ori_CFL_total = L_total;
    
    filename = ['SFL_bootstrapping_noise_',num2str(noise_percent)];
    load(filename)
    S_boot_SFL = S_boot;
    L_boot_SFL = L_boot;
    
    filename = ['SFL_result_dim2_',num2str(noise_percent)];
    load(filename)
    S_ori_SFL_total = S_total;
    L_ori_SFL_total = L_total;
    
    % import original score
    S_ori_CFL = reshape(S_ori_CFL_total(1,3,:),[num_data,1]);
    L_ori_CFL = reshape(L_ori_CFL_total(1,3,:),[num_data,1]);
    S_ori_SFL = reshape(S_ori_SFL_total(1,3,:),[num_data,1]);
    L_ori_SFL = reshape(L_ori_SFL_total(1,3,:),[num_data,1]);
    
    % cal p
    for i = 1:num_data
        [h,p_CFL] = ztest(S_ori_CFL(i), mean(S_boot_CFL(i,:)),std(S_boot_CFL(i,:)),'Tail','right');
        [h,p_SFL] = ztest(S_ori_SFL(i), mean(S_boot_SFL(i,:)),std(S_boot_SFL(i,:)),'Tail','right');
        
        p_list_CFL(noise_idx,i) = p_CFL;
        p_list_SFL(noise_idx,i) = p_SFL;
    end
end
%% calculate combined p-value
Fisher_list_CFL = [];
combined_p_CFL = [];
Fisher_list_SFL = [];
combined_p_SFL = [];

for noise_idx = 1:length(noise_list)
    sum_p_CFL = -2*nansum(log(p_list_CFL(noise_idx,:)));
    num_p_CFL = length(find(~isnan(p_list_CFL(noise_idx,:))));
    Fisher_list_CFL = [Fisher_list_CFL ; [sum_p_CFL, num_p_CFL]];
    combined_p_CFL = [combined_p_CFL; chi2cdf(sum_p_CFL, 2*num_p_CFL,'upper')]
    
    sum_p_SFL = -2*nansum(log(p_list_SFL(noise_idx,:)));
    num_p_SFL = length(find(~isnan(p_list_SFL(noise_idx,:))));
    Fisher_list_SFL = [Fisher_list_SFL ; [sum_p_SFL, num_p_SFL]];
    combined_p_SFL = [combined_p_SFL; chi2cdf(sum_p_SFL, 2*num_p_SFL,'upper')]
end

a = -2*log(0.001)*num_data;
p_thres = chi2cdf(a,2*num_data,'upper');

% combined p-value
figure(1)
plot(noise_list, -log(combined_p_CFL), '-ob')
hold on
plot(noise_list, -log(combined_p_SFL), '-or')
hold on
plot(noise_list, -log(p_thres)*ones(length(noise_list),1), '-og')

xlim([9,21])
xticks([10,15,20])
ylim([250,500])
yticks([250:50:500])



