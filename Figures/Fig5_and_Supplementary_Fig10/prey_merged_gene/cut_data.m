clc;
clear;
close all;

%% load data
%type = 'full';
%type = 'half';
type = 'quarter';

filename = ['data_merge_',type];
load(filename)

%% cut the data
%parameters
window_size = 101;
window_move = 10;
length_timeseries = length(t);
num_component = 4;

% variables
y_cut_total = {};
t_cut = t(1:window_size) / t(window_size);

% cut
start = 1;
while(1)
    if start + window_size-1 > length_timeseries
        break
    end
    y_cut = y(start:start + window_size - 1,:);
    
    for i = 1:4
        y_cut(:,i) = y_cut(:,i) - min(y_cut(:,i));
        y_cut(:,i) = y_cut(:,i) / max(y_cut(:,i));
    end
    
    y_cut_total{end+1} = y_cut; 
    start = start + window_move;
end

%% save data
y_total = y_cut_total;
t = t_cut;
num_data = length(y_total);
time_interval = 1/100;

filename = ['data_merge_cut_',type];
save(filename, 'y_total', 't', 'num_component', 'num_data', 'time_interval')



