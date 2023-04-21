clc;
clear;
close all;
addpath('../GOBI') 

%% load data=
noise_percent = 20;
data_index_2 = 38; % index of data for CFL
data_index_3 = 73; % index of data for SFL
num_data = 100;
num_boot = 100;
thres_noise = 0;
%% load timeseries
filename = ['CFL_timeseries_fit_',num2str(noise_percent)];
load(filename)
y_case2 = y_total;

filename = ['SFL_timeseries_fit_',num2str(noise_percent)];
load(filename)
y_case3 = y_total;

%% load score
filename = ['CFL_result_dim2_',num2str(noise_percent)];
load(filename)
S_case2 = S_total;
L_case2 = L_total;

S_tmp_plus  = reshape(S_case2(1,1,:),[num_data,1]);
S_tmp_minus = reshape(S_case2(1,3,:),[num_data,1]);
delta_2 = S_tmp_plus - S_tmp_minus; % compute delta for CFL

filename = ['SFL_result_dim2_',num2str(noise_percent)];
load(filename)
S_case3 = S_total;
L_case3 = L_total;

S_tmp_plus  = reshape(S_case3(1,1,:),[num_data,1]);
S_tmp_minus = reshape(S_case3(1,3,:),[num_data,1]);
delta_3 = S_tmp_plus - S_tmp_minus; % compute delta for SFL

delta_total = [delta_2 , delta_3];
%% bootstrap
y2 = cell2mat(y_case2(data_index_2));
y3 = cell2mat(y_case3(data_index_3));

% shuffle the time series
Z_shuffled_2 = [];
Z_shuffled_3 = [];
for i = 1:num_boot
    Z_tmp = y2(randperm(length(y2(:,1))),1);
    Z_shuffled_2 = [Z_shuffled_2,Z_tmp];
    Z_tmp = y3(randperm(length(y3(:,1))),1);
    Z_shuffled_3 = [Z_shuffled_3,Z_tmp];
end

% compute RDS with shuffled time series
boot2 = [];
boot3 = [];
for i = 1:num_boot
    % CFL
    [score_list, t_1, t_2] = RDS_dim2(Z_shuffled_2(:,i), y2(:,2), y2(:,3), t, time_interval);
    
    score_1 = reshape(score_list(:,:,3),[length(t),length(t)]);
    loca_plus = find(score_1 > thres_noise);
    loca_minus = find(score_1 < -thres_noise);
    if isempty(loca_plus) && isempty(loca_minus)
        s1 = 1;
    else
        s1 = (sum(score_1(loca_plus)) + sum(score_1(loca_minus)))/ (abs(sum(score_1(loca_plus))) + abs(sum(score_1(loca_minus))));
    end
    l1 = (length(loca_minus) + length(loca_plus)) / (length(t_1)*length(t_2)/2);
    
    boot2 = [boot2;[s1,l1]];
    
    % SFL
    [score_list, t_1, t_2] = RDS_dim2(Z_shuffled_3(:,i), y3(:,2), y3(:,3), t, time_interval);
    score_1 = reshape(score_list(:,:,3),[length(t),length(t)]);
    loca_plus = find(score_1 > thres_noise);
    loca_minus = find(score_1 < -thres_noise);
    if isempty(loca_plus) && isempty(loca_minus)
        s1 = 1;
    else
        s1 = (sum(score_1(loca_plus)) + sum(score_1(loca_minus)))/ (abs(sum(score_1(loca_plus))) + abs(sum(score_1(loca_minus))));
    end
    l1 = (length(loca_minus) + length(loca_plus)) / (length(t_1)*length(t_2)/2);
    
    boot3 = [boot3;[s1,l1]];
end

score_2 = S_case2(1,3,data_index_2);
score_3 = S_case3(1,3,data_index_3);

score_22 = S_case2(1,1,data_index_2);
score_33 = S_case3(1,1,data_index_3);

S_2 = reshape(S_case2(1,3,:),[num_data,1]);
S_3 = reshape(S_case3(1,3,:),[num_data,1]);
target_S = [S_2 , S_3];

%% plot histogram
figure(4) % CFL
x = [0.92:.000001:1];
pd = fitdist(boot2(:,1),'Normal');
y = normpdf(x,pd.mu,pd.sigma);
histogram(boot2(:,1),10, 'Normalization','pdf')
hold on
plot(x,y)
hold on
scatter(score_2,0,100,'*','b')
hold on
scatter(score_22,0,100,'o','b')
xlim([0.92,1])
yticks([])
%xticks([0.9,1])

figure(5) % SFL
x = [0.92:.000001:1];
pd = fitdist(boot3(:,1),'Normal');
y = normpdf(x,pd.mu,pd.sigma);
plot(x,y)
hold on
histogram(boot3(:,1),10, 'Normalization','pdf')
hold on
scatter(score_3,0,100,'*','b')
hold on
scatter(score_33,0,100,'o','b')
xlim([0.92,1])
yticks([])
%xticks([0.9,1])