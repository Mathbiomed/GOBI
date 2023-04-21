clc;
clear;
close all;

noise_list = [10];

for noise_percent = noise_list
    disp(noise_percent)
    filename = ['SFL_timeseries_noise_',num2str(noise_percent)];
    load(filename)
    y_total = {};
    period = t(end);
    t_fit = linspace(0, period, period/time_interval+1).';
    for i = 1:length(y_total_noise)
    %for i = 1
        y_tmp = cell2mat(y_total_noise(i));
        y_fit = zeros(length(t_fit),2);
        for j = 1:3
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
    filename = ['SFL_timeseries_fit_',num2str(noise_percent)];
    save(filename, 'y_total', 't', 'time_interval')
end
% 
% for noise_percent = noise_list
%     disp(noise_percent)
%     filename = ['CFL_timeseries_noise_',num2str(noise_percent)];
%     load(filename)
%     y_total = {};
%     period = t(end);
%     t_fit = linspace(0, period, period/time_interval+1).';
%     for i = 1:length(y_total_noise)
%         y_tmp = cell2mat(y_total_noise(i));
%         y_fit = zeros(length(t_fit),2);
%         for j = 1:3
%             y1 = fit(t,y_tmp(:,j),'fourier4');
%             w = y1.w;
%             fouriers = [
%                 ones(1,length(t_fit));
%                 cos(w*t_fit.');
%                 sin(w*t_fit.');
%                 cos(2*w*t_fit.');
%                 sin(2*w*t_fit.');
%                 cos(3*w*t_fit.');
%                 sin(3*w*t_fit.');
%                 cos(4*w*t_fit.');
%                 sin(4*w*t_fit.')];
%             coeffs = coeffvalues(y1);
%             y2 = coeffs(1:end-1) * fouriers;
%             y_fit(:,j) = y2;
%         end
%         y_total{end+1} = y_fit;
%     end
%     y_noise = cell2mat(y_total_noise(i));
%     filename = ['CFL_timeseries_fit_',num2str(noise_percent)];
%     save(filename, 'y_total', 't', 'time_interval')
% end

% for noise_percent = noise_list
%     disp(noise_percent)
%     filename = ['IFL_timeseries_noise_',num2str(noise_percent)];
%     load(filename)
%     y_total = {};
%     period = t(end);
%     t_fit = linspace(0, period, period/time_interval+1).';
%     for i = 1:length(y_total_noise)
%         y_tmp = cell2mat(y_total_noise(i));
%         y_fit = zeros(length(t_fit),2);
%         for j = 1:3
%             y1 = fit(t,y_tmp(:,j),'fourier4');
%             w = y1.w;
%             fouriers = [
%                 ones(1,length(t_fit));
%                 cos(w*t_fit.');
%                 sin(w*t_fit.');
%                 cos(2*w*t_fit.');
%                 sin(2*w*t_fit.');
%                 cos(3*w*t_fit.');
%                 sin(3*w*t_fit.');
%                 cos(4*w*t_fit.');
%                 sin(4*w*t_fit.')];
%             coeffs = coeffvalues(y1);
%             y2 = coeffs(1:end-1) * fouriers;
%             y_fit(:,j) = y2;
%         end
%         y_total{end+1} = y_fit;
%     end
%     y_noise = cell2mat(y_total_noise(i));
%     filename = ['IFL_timeseries_fit_',num2str(noise_percent)];
%     save(filename, 'y_total', 't', 'time_interval')
% end