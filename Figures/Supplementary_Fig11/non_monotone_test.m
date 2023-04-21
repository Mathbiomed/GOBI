clc;
clear;
close all;

addpath('../GOBI') 

%% parameter
range = 60; % range of time-series
time_interval = 1/20; % 1/time_interval represents sampling rate
range_interval = 20;
num_data = 10;

tspan = linspace(0,range, range/time_interval+1);
length_interval = range_interval / time_interval;

%% parameter of X (source)
time_interval_pre = 1/4;
xt_pre = linspace(0,range, range/time_interval_pre+1);
xt = linspace(0,range, range/time_interval+1);

%% parameter for ODE
alpha = 1;
beta = 0;
gamma = 1;
delta = 1;

%% create time series
y_total = {};
for i = 1:num_data
    %% create source (X)
    x_pre = rand([length(xt_pre),1])*4-2;
    x_source = spline(xt_pre,x_pre,xt);
    
    %% initial condition
    y_initial = rand([1,2])*4-2; % initial condition for Y and Z
    
    %% simulate ODE
    y_ori = zeros(length(tspan),3);
    y_ori(:,1) = x_source;

    % first ode: positive regulation from X to Y and Z
    [t1, y1] = ode45(@(t,y) XY_pos(t,y,xt(1:length_interval+1),x_source(1:length_interval+1),alpha,delta), tspan(1:length_interval+1), y_initial);
    y_ori(1:length_interval+1,2:3) = y1;
    
    % second ode: non-monotone regulation from X to Y / positive regulation
    % from X to Z
    [t2, y2] = ode45(@(t,y) XY_mix(t,y,xt(length_interval+1:2*length_interval+1),x_source(length_interval+1:2*length_interval+1),beta,gamma,delta), tspan(length_interval+1:2*length_interval+1), y1(end,:));
    y_ori(length_interval+1:2*length_interval+1,2:3) = y2;

    % third ode: positive, negative regulation from X to Y and Z,
    % respectively
    [t3, y3] = ode45(@(t,y) XY_neg(t,y,xt(length_interval*2+1:3*length_interval+1),x_source(length_interval*2+1:3*length_interval+1),alpha,delta), tspan(length_interval*2+1:3*length_interval+1), y2(end,:));
    y_ori(length_interval*2+1:3*length_interval+1,2:3) = y3;
    
    % save time series at the variable 'y_total'
    y_total{end+1} = y_ori;
    if i == 1
    figure(1)
    plot(tspan, y_ori(:,1), 'k')
    hold on
    plot(tspan, y_ori(:,2), 'r')
    hold on
    plot(tspan, y_ori(:,3), 'b')
    end
end

%% compute RDS
% parameters for moving window
size_window = 100;
move_window = 100;
num_window = (range_interval / time_interval - size_window) / move_window + 1;
thres_noise = 1e-5;

I_total_XY = zeros(num_data, num_window, size_window, size_window,2);
I_total_YZ = zeros(num_data, num_window, size_window, size_window,2);
X_trans = zeros(num_data, num_window, size_window);
Y_trans = zeros(num_data, num_window, size_window);

for i = 1:num_data
    disp(i)
    % import time series
    y_tmp = cell2mat(y_total(i));
    
    % cut using moving window technique
    for j = 1:num_window
        % import single time series
        y_target_tmp = y_tmp(1 + 400 + move_window*(j-1) : 400 + size_window + move_window*(j-1),:);
        y_target = y_target_tmp;
        t_target = tspan(1 + 400 + move_window*(j-1) : 400 + size_window + move_window*(j-1));
        
        %normalize the time series
%         for k = 1:3
%             y_target(:,k) = y_target_tmp(:,k) - min(y_target_tmp(:,k));
%             y_target(:,k) = y_target(:,k) / max(y_target(:,k));
%         end
        
        % compute RDf for the regulation from X to Y
        [score_list, t_1, t_2] = RDS_ns_dim1(y_target(:,1), y_target(:,2), t_target, time_interval);
        for k = 1:2
            score = reshape(score_list(:,:,k),[length(t_target),length(t_target)]);
            I_total_XY(i,j,:,:,k) = score;            
        end
        
        [score_list, t_1, t_2] = RDS_ns_dim1(y_target(:,2), y_target(:,3), t_target, time_interval);
        for k = 1:2
            score = reshape(score_list(:,:,k),[length(t_target),length(t_target)]);
            I_total_YZ(i,j,:,:,k) = score;            
        end
        
        % save the value of X and Y
        X_trans(i,j,:) = y_target(:,1);
        Y_trans(i,j,:) = y_target(:,2);
    end
end

%% transform (t,t*) into X(t),X(t*) and Y(t),Y(t*)
trans_total_XY = cell(num_data, num_window,2);
trans_total_YZ = cell(num_data, num_window,2);

for i = 1:num_data
    for j = 1:num_window
        score_XY_pos = reshape(I_total_XY(i,j,:,:,1), [size_window, size_window]);
        score_XY_neg = reshape(I_total_XY(i,j,:,:,2), [size_window, size_window]);
        score_YZ_pos = reshape(I_total_YZ(i,j,:,:,1), [size_window, size_window]);
        score_YZ_neg = reshape(I_total_YZ(i,j,:,:,2), [size_window, size_window]);
        
        X_var = reshape(X_trans(i,j,:), [size_window,1]);
        Y_var = reshape(Y_trans(i,j,:), [size_window,1]);
        
        trans_XY_pos = [];
        trans_YZ_pos = [];
        trans_XY_neg = [];
        trans_YZ_neg = [];
        
        for m = 1:size_window
            for n = 1:size_window
                if abs(score_XY_pos(m,n)) > thres_noise
                    trans_XY_pos = [trans_XY_pos ; [X_var(m), X_var(n), Y_var(m), Y_var(n), score_XY_pos(m,n)]];
                end
                if abs(score_YZ_pos(m,n)) > thres_noise
                    trans_YZ_pos = [trans_YZ_pos ; [Y_var(m), Y_var(n), score_YZ_pos(m,n)]];
                end
                if abs(score_XY_neg(m,n)) > thres_noise
                    trans_XY_neg = [trans_XY_neg ; [X_var(m), X_var(n), Y_var(m), Y_var(n), score_XY_neg(m,n)]];
                end
                if abs(score_YZ_neg(m,n)) > thres_noise
                    trans_YZ_neg = [trans_YZ_neg ; [Y_var(m), Y_var(n), score_YZ_neg(m,n)]];
                end
            end
        end
        
        trans_total_XY{i,j,1} = trans_XY_pos;
        trans_total_YZ{i,j,1} = trans_YZ_pos;
        trans_total_XY{i,j,2} = trans_XY_neg;
        trans_total_YZ{i,j,2} = trans_YZ_neg;
    end
end

%% plot
thres_noise = 1e-3;
figure(2)
for i = 1:num_data
    for j = 1:num_window
    %for j = 1:1
        trans_XY_pos = cell2mat(trans_total_XY(i,j,1));
        trans_XY_neg = cell2mat(trans_total_XY(i,j,2));
        
%         if isempty(trans_XY)
%             continue
%         end
%         
%         if isempty(find(trans_XY(:,3) > 0)) || isempty(find(trans_XY(:,3) < 0))
%             continue
%         end

        loca_pos = find(trans_XY_pos(:,5) > 0);
        loca_neg = find(trans_XY_pos(:,5) < 0);

        scatter(trans_XY_pos(loca_pos,1),trans_XY_pos(loca_pos,2),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.005,'MarkerEdgeAlpha',.005)
        hold on
        scatter(trans_XY_pos(loca_neg,1),trans_XY_pos(loca_neg,2),'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.005,'MarkerEdgeAlpha',.005)
        hold on
        
        loca_pos = find(trans_XY_neg(:,5) > 0);
        loca_neg = find(trans_XY_neg(:,5) < 0);

        scatter(trans_XY_neg(loca_pos,1),trans_XY_neg(loca_pos,2),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.005,'MarkerEdgeAlpha',.005)
        hold on
        scatter(trans_XY_neg(loca_neg,1),trans_XY_neg(loca_neg,2),'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.005,'MarkerEdgeAlpha',.005)
        hold on
    end
end
plot([-1,1],[0,0],'--k', 'LineWidth',2)
hold on
plot([0,0],[-1,1],'--k', 'LineWidth',2)
hold on
plot([-1,1],[1,-1],'--k', 'LineWidth',2)

% xlim([-0.5,0.5])
% ylim([-0.5,0.5])

% figure(4)
% for i = 1:num_data
%     for j = 1:num_window
%     %for j = 1:1
%         trans_XY_pos = cell2mat(trans_total_XY(i,j,1));
%         trans_XY_neg = cell2mat(trans_total_XY(i,j,2));
% 
%         loca_pos = find(trans_XY_pos(:,3) - trans_XY_pos(:,4) > -thres_noise);
%         loca_neg = find(trans_XY_pos(:,3) - trans_XY_pos(:,4) < -thres_noise);
% 
%         scatter(trans_XY_pos(loca_pos,1),trans_XY_pos(loca_pos,2),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)
%         hold on
%         %scatter(trans_XY_pos(loca_neg,1),trans_XY_pos(loca_neg,2),'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.005,'MarkerEdgeAlpha',.005)
%         %hold on
%         
%         loca_pos = find(trans_XY_neg(:,3) - trans_XY_neg(:,4) > -thres_noise);
%         loca_neg = find(trans_XY_neg(:,3) - trans_XY_neg(:,4) < -thres_noise);
% 
%         scatter(trans_XY_neg(loca_pos,1),trans_XY_neg(loca_pos,2),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)
%         hold on
%         %scatter(trans_XY_neg(loca_neg,1),trans_XY_neg(loca_neg,2),'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.005,'MarkerEdgeAlpha',.005)
%         %hold on
%     end
% end
% 
% figure(3)
% for i = 1:5
%     for j = 1:num_window
%         
%         trans_YZ_pos = cell2mat(trans_total_YZ(i,j,1));
%         trans_YZ_neg = cell2mat(trans_total_YZ(i,j,2));
%         
% %         if isempty(trans_YZ)
% %             continue
% %         end
%         
%         loca_pos = find(trans_YZ_pos(:,3) > 0);
%         loca_neg = find(trans_YZ_pos(:,3) < 0);
% 
%         scatter(trans_YZ_pos(loca_pos,1),trans_YZ_pos(loca_pos,2),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.005,'MarkerEdgeAlpha',.005)
%         hold on
%         scatter(trans_YZ_pos(loca_neg,1),trans_YZ_pos(loca_neg,2),'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.005,'MarkerEdgeAlpha',.005)
%         hold on
%         
%         loca_pos = find(trans_YZ_neg(:,3) > 0);
%         loca_neg = find(trans_YZ_neg(:,3) < 0);
% 
%         scatter(trans_YZ_neg(loca_pos,1),trans_YZ_neg(loca_pos,2),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',.005,'MarkerEdgeAlpha',.005)
%         hold on
%         scatter(trans_YZ_neg(loca_neg,1),trans_YZ_neg(loca_neg,2),'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',.005,'MarkerEdgeAlpha',.005)
%         hold on
%     end
% end
