clc;
clear;
close all;

%% parameters
noise_percent = 20;
data_index_2 = 38; % index of data for CFL
data_index_3 = 73; % index of data for SFL

%% load timeseries
filename = ['CFL_timeseries_noise_',num2str(noise_percent)];
load(filename)
y_case2 = y_total_noise;

filename = ['SFL_timeseries_noise_',num2str(noise_percent)];
load(filename)
y_case3 = y_total_noise;

%% load interpolated timeseries
filename = ['CFL_timeseries_fit_',num2str(noise_percent)];
load(filename)
y_fit = y_total;

%% plot time series
y_CFL = cell2mat(y_case2(data_index_2));
y_SFL = cell2mat(y_case3(data_index_3));

% normalize
for i = 1:3
    y_CFL(:,i) = y_CFL(:,i) - min(y_CFL(:,i));
    y_CFL(:,i) = y_CFL(:,i) / max(y_CFL(:,i));
    y_SFL(:,i) = y_SFL(:,i) - min(y_SFL(:,i));
    y_SFL(:,i) = y_SFL(:,i) / max(y_SFL(:,i));
end

figure(1) % Fig 3f
plot(t, y_CFL(:,1), 'k')
hold on
plot(t, y_CFL(:,2), 'r')
hold on
plot(t, y_CFL(:,3), 'b')

xlim([0,10])
xticks([0,10])
xticklabels([])
ylim([0,1])
yticks([0,1])
yticklabels([])


figure(2) % Fig 3g
plot(t, y_SFL(:,1), 'k')
hold on
plot(t, y_SFL(:,2), 'r')
hold on
plot(t, y_SFL(:,3), 'b')

xlim([0,10])
xticks([0,10])
xticklabels([])
ylim([0,1])
yticks([0,1])
yticklabels([])

%% shuffle the time series
y_CFL_fit = cell2mat(y_fit(data_index_2));

%normalize
for i = 1:3
    y_CFL_fit(:,i) = y_CFL_fit(:,i) - min(y_CFL_fit(:,i));
    y_CFL_fit(:,i) = y_CFL_fit(:,i) / max(y_CFL_fit(:,i));
end

figure(3) % Fig 3h left
plot(t, y_CFL_fit(:,1),'k')
xlim([0,10])
xticks([0,10])
xticklabels([])
ylim([0,1])
yticks([0,1])
yticklabels([])

A_shuffle_1 = y_CFL_fit(randperm(length(y_CFL_fit(:,1))),1);
A_shuffle_2 = y_CFL_fit(randperm(length(y_CFL_fit(:,1))),1);

figure(4) % Fig 3h top-right
plot(t, A_shuffle_1,'k')
xlim([0,10])
xticks([0,10])
xticklabels([])
ylim([0,1])
yticks([0,1])
yticklabels([])

figure(5) % Fig 3h top-left
plot(t, A_shuffle_2,'k')
xlim([0,10])
xticks([0,10])
xticklabels([])
ylim([0,1])
yticks([0,1])
yticklabels([])
