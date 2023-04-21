clc;
clear;
close all;

%% import data
load('data_estradiol_cut.mat')

%% parameter from data
num_component = 4;
length_data = length(y_total);

%% norm
y_total_norm = {};
for i = 1:length_data
    y_tmp = cell2mat(y_total(i));
    for j = 1:num_component
        y_tmp(:,j) = y_tmp(:,j) - min(y_tmp(:,j));
        y_tmp(:,j) = y_tmp(:,j) / max(y_tmp(:,j));
    end
    y_total_norm{end+1} = y_tmp;
end
y_total = y_total_norm;

%% Fourier fitting
t_fit = t(1:length(y_tmp(:,1)));
y_fit_total = {};
for i = 1:length_data
    y_fit = zeros(length(t_fit),num_component);
    for j = 1:num_component
        y_tmp = cell2mat(y_total(i));
        y1 = fit(t_fit,y_tmp(:,j),'fourier4');
        w = y1.w;
        fouriers = [
            ones(1,length(t_fit));
            cos(w*t_fit.');
            sin(w*t_fit.');
            cos(2*w*t_fit.');
            sin(2*w*t_fit.');
            cos(3*w*t_fit.');
            sin(3*w*t_fit.');
            cos(4*w*t_fit.');
            sin(4*w*t_fit.')];
        coeffs = coeffvalues(y1);
        y2 = coeffs(1:end-1) * fouriers;
        y_fit(:,j) = y2;
    end
    y_fit_total{end+1} = y_fit; 
end
%% calculate residual
residual_total = zeros(length_data, num_component);
for i = 1:length_data
     y_fit = cell2mat(y_fit_total(i));
     y_noise = cell2mat(y_total(i));
     for j = 1:num_component
        residual_tmp = mean((y_fit(:,j) - y_noise(:,j)).^2);
        residual_total(i,j) = residual_tmp;
    end
end
mean(mean(residual_total))

