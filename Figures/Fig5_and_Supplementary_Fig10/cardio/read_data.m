clc;
clear;
close all;

%% load data
load('hk_data_v1.mat')
t = linspace(0,1031,1032).';
y = [cardio,no2,so2,o3,rspar];

save('data_cardio','t','y')



