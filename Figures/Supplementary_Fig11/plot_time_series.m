clc;
clear;
close all;

load('data_temporal.mat')

%% plot
figure(1)
idx = 8;
y_ori = cell2mat(y_total(idx));
plot(tspan, y_ori(:,1), 'k', 'LineWidth',2)
hold on
plot(tspan, y_ori(:,2), 'r','LineWidth',2)
hold on
xlim([0,100])
ylim([0,1])
xticks([0:25:100])
yticks([0:0.5:1])
xticklabels([])
yticklabels([])


figure(2) % figure of parameter
plot(t_pre(1:1000), alpha_pre(1:1000), 'r', 'LineWidth',2)
hold on
plot(t_pre(500:1500), beta_pre(500:1500), 'b', 'LineWidth',2)
hold on
plot(t_pre(1501:2000), gamma_pre(1501:2000), 'k', 'LineWidth',2)

xlim([0,100])
ylim([0,1])
xticks([0:25:100])
yticks([0,1])
xticklabels([])
yticklabels([])