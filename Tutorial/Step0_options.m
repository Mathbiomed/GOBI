clc;
clear;
close all;
%% import data
load('data_cut')

%% options 1. choose thresholds and critical values

% defaults
noise_level = 8; % compute noise level using residuals
thres_L = 0;     % threshold for regulation-detection region
                 % thres_L = 0 represents default (decreases as dimension increases)
                 
thres_S = 0.9 - 0.005 * noise_level; % threshold for regulation-detection score
thres_TRS = 0.9 - 0.01 * noise_level;% threshold for total regulation score

p_delta = 0.01; % critical value for delta test
p_surrogate = 0.001; % critical value for surrogate test

% users can adjust the threshold values
% thres_L = 0;
% thres_S = 1;
% thres_TRS = 1;

%% options 2. Choose the types of self-regulations

% type_self = nan: infer without any assumptions of self regulations
% type_self = -1 : negative self regulations
% type_self = 0  : no self regulations
% type_self = 1  : positive self regulations

type_self = -1;

%% options 3. maximum dimension of framework

% default
if isnan(type_self)
    max_D = num_component;
else
    max_D = num_component - 1;
end

% users can adjust the max_D values
% max_D = 2;

filename = ['data_with_options'];
save(filename, 't','y','y_total','time_interval','num_data','num_component','thres_L','thres_S','thres_TRS','p_delta','p_surrogate','type_self','max_D')


