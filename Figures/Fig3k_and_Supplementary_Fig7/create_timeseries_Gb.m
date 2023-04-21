clc;
clear;
close all;

trial_list = [0:9];
%trial_list = 1;
for trial = trial_list
    %% parameters
    time_interval = 1/10;
    num_component = 5;
    num_data = 100;
    range_extention = 0.5;
    %trial = 1;

    %% Find the range of original time-series
    % simulate in large time scale
    initial = ones(1,num_component);
    period = 100;
    tspan = linspace(0,period, period/time_interval+1);

    [t_ori, y_ori] = ode45(@(t,y) Goldbeter(t,y), tspan, initial);

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
    %period = ceil(period);
    period = 10;
    % find range of time-series
    st_period = locs(end-1);
    ed_period = locs(end);
    range = zeros(num_component,2);
    for i = 1:num_component
        range(i,1) = min(y_ori(st_period:ed_period,i)) * range_extention;
        range(i,2) = max(y_ori(st_period:ed_period,i)) * (1+range_extention);
    end

    %% create time-series with various initial value
    % create initials
    initials  =  rand([num_data,num_component])*1;
    for i = 1:num_component
        %initials(:,i) = initials(:,i) * (range(i,2) - range(i,1)) + range(i,1);
        initials(:,i) = initials(:,i) * 2;

    end

    % simulate time-series
    y_total = {};
    tspan = linspace(0,period, period/time_interval+1);
    for i = 1:length(initials(:,1))
        [t1, y1] = ode45(@(t,y) Goldbeter(t,y), tspan, initials(i,:));
        for j = 1:num_component
            y1(:,j) = y1(:,j) - min(y1(:,j));
            y1(:,j) = y1(:,j) / max(y1(:,j));
        end
        y_total{end+1} = y1;
    end
    t = t1;

    figure(2) % various time-series
    for i = 1:20
        y_tmp = cell2mat(y_total(i));
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
    filename = ['./Data_original/Gb_timeseries_Trial',num2str(trial)];
    %filename = ['Gb_timeseries_large'];
    save(filename, 'y_total', 't', 'time_interval', 'period', 'num_component', 'num_data')
end


