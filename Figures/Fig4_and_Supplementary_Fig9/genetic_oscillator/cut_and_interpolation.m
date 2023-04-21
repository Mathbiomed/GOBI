clc;
clear;
close all;

period_list = [40,20,15,10,10,50,15,15];

num_component = 2;
y_total = {};
for i = 1:8
    %% load data
    filename = ['data_gene_raw_',num2str(i)];
    load(filename)
    y = y';
    
    %% plot raw data
    figure(1)
    if i <= 4
    subplot(4,2,i)
    plot([7:length(t)], y(7:end,1),'b')
    hold on
    plot([7:length(t)], y(7:end,2),'r')
    end
    if i > 4
    subplot(4,2,i)
    plot([7:70], y(7:70,1),'b')
    hold on
    plot([7:70], y(7:70,2),'r')
    end
    %% cut and interpolate the data
    window_size = period_list(i);
    window_move = ceil(period_list(i));

    start = 7;
    length_timeseries = length(y(:,1));
    if i > 4
        length_timeseries = 70;
    end
    while(1)
        if start + window_size-1 > length_timeseries
            break
        end
        y_tmp = y(start:start + window_size - 1,:);
        
        % spline fit
        time_interval = 1/100;
        t_tmp = t(1:window_size);
        t_tmp = t_tmp / t_tmp(end);
        t_fit = linspace(0,1,1/time_interval+1).';
        
        y_fit = zeros(length(t_fit),2);
        y_fit(:,1) = interp1(t_tmp, y_tmp(:,1), t_fit, 'spline');
        y_fit(:,2) = interp1(t_tmp, y_tmp(:,2), t_fit, 'spline');
        
        % normalize
        y_fit(:,1) = y_fit(:,1) - min(y_fit(:,1));
        y_fit(:,1) = y_fit(:,1) / max(y_fit(:,1));
        y_fit(:,2) = y_fit(:,2) - min(y_fit(:,2));
        y_fit(:,2) = y_fit(:,2) / max(y_fit(:,2));
        
        % save
        y_total{end+1} = y_fit; 

        start = start + window_move;
    end
    
end
num_data = length(y_total);

%% Save data

save('data_final_gene','y_total','num_data','num_component','time_interval')