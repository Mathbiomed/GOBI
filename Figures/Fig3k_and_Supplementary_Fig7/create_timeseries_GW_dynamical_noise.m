clc;
clear;
close all;

trial_list = [0:9];
noise_list = [2:2:20];

time_interval = 1/15;
num_component = 4;
num_data = 100;
range_extention = 0.5;
parpool threads;

for trial = trial_list
    disp(trial)
    for noise_percent = noise_list
        
        %% Find the range of original time-series
        % simulate in large time scale
        initial = ones(1,num_component);
        period = 100;
        tspan = linspace(0,period, period/time_interval+1);

        [t_ori, y_ori] = ode45(@(t,y) Goodwin(t,y), tspan, initial);

        figure(1) % original time-series
        plot(t_ori, y_ori(:,1), 'k')
        hold on
        plot(t_ori, y_ori(:,2), 'r')
        hold on
        plot(t_ori, y_ori(:,3), 'b')
        hold on
        plot(t_ori, y_ori(:,4), 'g')
        hold on

        % find period using peaks
        [~,locs]=findpeaks(y_ori(:,1));
        period = (locs(end) - locs(end-1)) * time_interval;
        %period = ceil(period)/2;
        period = 20;
        % find range of time-series
        st_period = locs(end-1);
        ed_period = locs(end);
        range = zeros(num_component,2);
        mean_range = zeros(num_component,1);
        for i = 1:num_component
            range(i,1) = min(y_ori(st_period:ed_period,i)) * range_extention;
            range(i,2) = max(y_ori(st_period:ed_period,i)) * (1+range_extention);
            mean_range(i,1) = mean(y_ori(st_period:ed_period,i));
        end

        %% create time-series with various initial value
        % create initials
        initials  =  rand([num_data,num_component])*1;
        for i = 1:num_component
            initials(:,i) = initials(:,i) * (range(i,2) - range(i,1)) + range(i,1);

        end

        % simulate time-series
        y_total_noise = {1,100};
        tspan = linspace(0,period, period/time_interval+1);
        parfor i = 1:length(initials(:,1))
            [t1, y1] = ode45(@(t,y) Goodwin_noise(t,y,noise_percent/100,mean_range(1),mean_range(2),mean_range(3),mean_range(4)), tspan, initials(i,:));
            for j = 1:num_component
                y1(:,j) = y1(:,j) - min(y1(:,j));
                y1(:,j) = y1(:,j) / max(y1(:,j));
            end
            y_total_noise{i} = y1;
        end
        t = tspan;

        figure(2) % various time-series
        for i = 1:1
            y_tmp = cell2mat(y_total_noise(i));
            plot(t/t(end), y_tmp(:,1), 'k')
            hold on
            plot(t/t(end), y_tmp(:,2), 'r')
            hold on
            plot(t/t(end), y_tmp(:,3), 'b')
            hold on
            plot(t/t(end), y_tmp(:,4), 'g')
            hold on
        end
        xlim([0,1])
        xticks([0,1])
        xticklabels({'0', 'period'})
        yticklabels([])
        xlabel('Time')
        ylabel('Value of component')
        set(gca,'fontsize',16)

        % save data
        filename = ['./Data_dynamical/GW_timeseries_dynamical_',num2str(noise_percent),'_Trial',num2str(trial)];
        save(filename, 'y_total_noise', 't', 'time_interval','noise_percent','num_component','num_data','period')
    end
end
delete(gcp('nocreate'))

