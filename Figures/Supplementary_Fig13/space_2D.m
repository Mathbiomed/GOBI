clc;
clear;
close all;
addpath('../GOBI') 

%% load data
load('2D_timeseries.mat')

%% parameter
thres_noise = 1e-7;
size_window = length(t);
I_total_XY = zeros(num_data, size_window, size_window, 2);
X_trans = zeros(num_data, size_window);
Y_trans = zeros(num_data, size_window);

%% compute RDF
disp('compute RDF')
for i = 1:num_data
    % import time series
    y_target = cell2mat(y_total(i));
    t_target = t;
    % compute RDf for the regulation from X to Y
    [score_list, t_1, t_2] = RDS_dim2(y_target(:,2), y_target(:,3), y_target(:,1), t_target, time_interval);
    for k = 1:4
        score = reshape(score_list(:,:,k),[length(t_target),length(t_target)]);
        I_total_XY(i,:,:,k) = score;            
    end

    % save the value of X and Y
    X_trans(i,:) = y_target(:,2);
    Y_trans(i,:) = y_target(:,3);
end

%% transform (t,t*) into X(t),X(t*)
disp('transformation')
trans_total_XY = cell(num_data,1);

for i = 1:num_data
    score_XY = reshape(I_total_XY(i,:,:,1), [size_window, size_window]);
    
    X_var = reshape(X_trans(i,:), [size_window,1]);
    Y_var = reshape(Y_trans(i,:), [size_window,1]);

    trans_XY = [];
    
    for m = 1:size_window
        for n = 1:size_window
            if abs(score_XY(m,n)) > thres_noise
                trans_XY = [trans_XY ; [t(m), t(n), score_XY(m,n)]];
            end
        end
    end        
    trans_total_XY{i,1} = trans_XY;
end

%% plot
disp('plot')
alpha = 0.05;

figure(1)
num_data_used = 25;
for i = 1:num_data_used
    
    trans_XY = cell2mat(trans_total_XY(i,1));

    if isempty(trans_XY)
        continue
    end

    loca_pos = find(trans_XY(:,3) > 0);
    loca_neg = find(trans_XY(:,3) < 0);

    scatter(trans_XY(loca_pos,1),trans_XY(loca_pos,2),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha)
    hold on
    scatter(trans_XY(loca_neg,1),trans_XY(loca_neg,2),'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha)
    hold on

    xlim([0,1])
    ylim([0,1])
    xticks([0,1])
    yticks([0,1])
    set(gca,'XColor', 'none','YColor','none')
    axis square
end
filename = ['figure_2D_',num2str(num_data_used),'.png'];
saveas(gcf,filename)