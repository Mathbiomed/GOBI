clc;
clear;
close all;

addpath('../GOBI')

global idx

%% parameters
lag_list = [0:2:50]; % time lag for extended regulation-detection score
data_list = [0:0.1:2]; % data with different initial conditions
data_list = data_list(1:end-1);
num_data = length(data_list);

num_lag = length(lag_list);
S_total = zeros(20, num_lag);
y_total = {};

for i = 1:num_data
        
    idx = data_list(i);
    
    %% solve DDE
    tspan = linspace(0,20,201);
    lag = [5,3]; % parameter for DDE
    options = ddeset('RelTol',1e-4,'AbsTol',1e-4);
    dde_sol = dde23(@dde_fun,lag,@dde_hist,tspan,options);
    
    % plot time series simulated from DDE
    figure(1)
    plot(dde_sol.x, dde_sol.y(1,:),'b')
    hold on
    plot(dde_sol.x, dde_sol.y(2,:),'r')
    xlabel('t')

    %% interpolate time series
    y_tmp = dde_sol.y;
    y_tmp = y_tmp.';
    t_tmp = dde_sol.x;
    t_tmp = t_tmp.';

    y_fit = [];
    for j = 1:2
        y_tmp_fit = spline(t_tmp, y_tmp(:,j), tspan);
        y_fit = [y_fit ; y_tmp_fit];
    end
    y_fit = y_fit.';
    
    y_total{end+1} = y_fit;
    
    % plot interpolated time series which were simulated from DDE
    figure(2)
    plot(tspan, y_fit(:,1), 'b')
    hold on
    plot(tspan, y_fit(:,2), 'r')
    hold on


    %% compute extended RDS with time lag
    
    thres_noise = 0;
    S_tmp = zeros(num_lag,2);
    L_tmp = zeros(num_lag,2);

    for j = 1:num_lag
        lag_tmp = lag_list(j);
        cause  = y_fit(1:end-lag_tmp,1);
        target = y_fit(1+lag_tmp:end,2);
        t = tspan(40:end-lag_tmp);
        time_interval = 1/20;

        [score_list, t_1, t_2] = RDS_ns_dim1(cause, target, t, time_interval);
        for k = 1:2
            score = reshape(score_list(:,:,k),[length(t),length(t)]);
            loca_plus = find(score > thres_noise);
            loca_minus = find(score < -thres_noise);
            if isempty(loca_plus) && isempty(loca_minus)
                s = 1;    
            else
                s = (sum(score(loca_plus)) + sum(score(loca_minus)))/ (abs(sum(score(loca_plus))) + abs(sum(score(loca_minus))));
            end
            l = (length(loca_minus) + length(loca_plus)) / (length(t_1)*length(t_2)/2);
            S_tmp(j,k) = s;
            L_tmp(j,k) = l;
        end
    end

    S_total(i,:) = S_tmp(:,2).'; 
end

figure(3) % plot mean and std of extended RDS
S_mean = mean(S_total);
S_std = std(S_total);
errorbar(lag_list, S_mean, S_std,"-o","MarkerSize",5)

xlim([0,50])
xticks([0:10:50])
xticklabels([0:1:5])
ylim([-1.5,1.5])
yticks([-1:0.5:1])

%% DDE
function dydt = dde_fun(t,y,Z)
alpha_x = 1; alpha_y = 1;
n_x = 6; n_y = 6;
beta_x = 1; beta_y = 1;


X_lag1 = Z(:,1);
X_lag2 = Z(:,2);

dydt = zeros(2,1);
dydt(1) = alpha_x / [1 + X_lag1(1).^n_x] - beta_x * y(1);
dydt(2) = alpha_y / [1 + X_lag2(1).^n_y] - beta_y * y(2);

end

%% initial conditions
function H = dde_hist(t)

global idx;
H = ones(2,1);
H(1) = sin(t + pi * idx)+1;

end