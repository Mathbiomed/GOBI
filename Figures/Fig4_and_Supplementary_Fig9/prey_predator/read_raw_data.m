clc;
clear;
close all;

%% Read data file
data = readtable('raw_data_prey_predator.xlsx');
t = table2array(data(:,1));
P = table2array(data(:,2)); % Paramecium (prey)
D = table2array(data(:,3)); % Didinium (predator)

y = [P,D];
save('raw_data_prey_predator','t','y')