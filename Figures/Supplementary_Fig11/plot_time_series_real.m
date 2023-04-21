clc;
clear;
close all;

%% load data
load('data_genetic_oscillator.mat')

%% parameters
num_data = 8;
num_component = 2;

%% spline fit
t_fit = linspace(0,1,201).';
time_interval = 1/200;
data_spline = zeros(length(t_fit),num_component,num_data);
for i = 1:num_data
    y_tmp = reshape(data_total(:,:,i),[num_component, length(time)]);
    y_tmp = y_tmp .';
    
    
    y_int = zeros(length(t_fit),num_component);
    for j = 1:2
        y_int(:,j) = spline(time, y_tmp(:,j), t_fit);
    end 
    data_spline(:,:,i) = y_int;
end
%% plot
figure(1)
idx = 1;
y_tmp = reshape(data_spline(:,:,idx), [length(t_fit),num_component]);
plot(t_fit, y_tmp(:,1), 'r')
hold on
plot(t_fit, y_tmp(:,2), 'b')
xlim([0,1])
ylim([0,800])
xticks([0,1])
yticks([0,800])

figure(2)
idx = 2;
y_tmp = reshape(data_spline(:,:,idx), [length(t_fit),num_component]);
plot(t_fit, y_tmp(:,1), 'r')
hold on
plot(t_fit, y_tmp(:,2), 'b')
xlim([0,1])
ylim([0,300])
xticks([0,1])
yticks([0,300])

figure(3)
idx = 3;
y_tmp = reshape(data_spline(:,:,idx), [length(t_fit),num_component]);
plot(t_fit, y_tmp(:,1), 'r')
hold on
plot(t_fit, y_tmp(:,2), 'b')
xlim([0,1])
ylim([0,200])
xticks([0,1])
yticks([0,200])

figure(4)
idx = 4;
y_tmp = reshape(data_spline(:,:,idx), [length(t_fit),num_component]);
plot(t_fit, y_tmp(:,1), 'r')
hold on
plot(t_fit, y_tmp(:,2), 'b')
xlim([0,1])
ylim([0,200])
xticks([0,1])
yticks([0,200])

figure(5)

scatter([1,2,3], [0.5796, 0.5306, 0.0667])
xlim([0.5,3.5])
ylim([0,0.6])
xticks([1:3])
yticks([0,0.6])
axis square

