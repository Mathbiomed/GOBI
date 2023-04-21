clc;
clear;
close all;

%% parameter
system_list = {'KF','Fr','GW','Gb','cAMP'};
%system_list = {'Fr'};
num_system = length(system_list);
noise_list = [0:2:30];

%% calculate residual
residual_list = zeros(num_system, length(noise_list));
for i = 1:num_system
%for i = 1
    for j = 1:length(noise_list)
        %load data
        system_name = string(system_list(i));
        filename = append('./',system_name,'/',system_name,'_timeseries_fit_',num2str(noise_list(j)));
        load(filename)
        filename = append('./',system_name,'/',system_name,'_timeseries_noise_',num2str(noise_list(j)));
        load(filename)
        
        % for each data, each component, calculate residual
        residual_total = zeros(num_data, num_component);
        for k = 1:num_data
            y_fit = cell2mat(y_total(k));
            y_noise = cell2mat(y_total_noise(k));
            for l = 1:num_component
                residual_tmp = mean((y_fit(:,l) - y_noise(:,l)).^2);
                residual_total(k,l) = mean(residual_tmp);
            end
        end
        
        % make average for all the data
        residual_list(i,j) = mean(mean(residual_total));
    end
end

figure(1)
for i = 1:num_system
    plot(noise_list, residual_list(i,:))
    hold on
end
xlim([0,30])
ylim([0,0.03])
xticks([0:10:30])
yticks([0,0.01,0.02,0.03])
ax = gca;
ax.FontSize = 16; 
legend(system_list)

figure(2)
plot(noise_list, mean(residual_list))
hold on
ax = gca;
ax.FontSize = 16; 
