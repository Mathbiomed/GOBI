clc;
clear;
close all;

%system_list = {'cAMP','Fr','Gb','GW','KF'};
%noise_type_list = {'additive','blue','brown','pink','purple','dynamical','multiplicative'};

system_list = {'cAMP'};
noise_type_list = {'additive','blue','pink','multiplicative'};
%noise_type_list = {'dynamical'};
trial_list = [0:9];
noise_list = [2:2:20];

idx = 0;
figure(1)
for noise_type_idx = 1:length(noise_type_list)
    noise_type = char(noise_type_list(noise_type_idx));
    disp(noise_type)
    for system_idx = 1:length(system_list)
        system = char(system_list(system_idx));
        disp(system)
        
        idx = idx + 1;
        
        subplot(2,2,idx)
        for trial = trial_list
            %% fitting the noisy time-series using fourier fit
            for noise_percent = noise_list
                
                % load noisy time-series data
                filename = ['./Data_',noise_type,'/',system,'_timeseries_',noise_type, '_',num2str(noise_percent),'_Trial',num2str(trial)];
                load(filename)

                % fit
                y_total = {};
                t_fit = linspace(0, period, period/time_interval+1).';
                for i = 1:num_data
                    y_tmp = cell2mat(y_total_noise(i));
                    y_tmp = double(y_tmp);
                    y_fit = zeros(length(t_fit),num_component);
                    for j = 1:num_component
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
                    y_total{end+1} = y_fit;
                end
                
                y_noise = cell2mat(y_total_noise(i));
                y = cell2mat(y_total(end));
                plot(t, y_noise(:,1),'k')
                hold on
                plot(t, y(:,1), 'r')
                hold on

                filename = ['./Data_',noise_type,'_fit/',system,'_timeseries_fit_',num2str(noise_percent),'_Trial',num2str(trial)];
                save(filename, 'y_total', 't', 'time_interval','noise_percent','num_component','num_data','period')
            end
        end
    end
end