clc;
clear;
close all;

%% parameters
noise_list = [22:2:30];
noise_list = [0];
%% fitting the noisy time-series using fourier fit
for noise_percent = noise_list
    disp(noise_percent)

    % load noisy time-series data
    filename = ['Gb_timeseries_noise_',num2str(noise_percent)];
    load(filename)

    % fit
    y_total = {};
    t_fit = linspace(0, period, period/time_interval+1).';
    for i = 1:num_data
        y_tmp = cell2mat(y_total_noise(i));
        y_fit = zeros(length(t_fit),num_component);
        for j = 1:num_component
            y1 = fit(t,y_tmp(:,j),'fourier4');
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
        y_total{end+1} = y_fit;
    end

    y_noise = cell2mat(y_total_noise(i));
    y = cell2mat(y_total(end));
    plot(t, y_noise(:,1),'k')
    hold on
    plot(t, y(:,1), 'r')

    filename = ['Gb_timeseries_fit_',num2str(noise_percent)];
    save(filename, 'y_total', 't', 'time_interval','noise_percent','num_component','num_data','period')
end