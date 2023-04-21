clc;
clear;
close all;

addpath('../GOBI') 

%% load data
load('data_temporal.mat')

%% parameter
st = 250;
size_window = 100;
move_window = 500;
num_window = ceil((range / time_interval - size_window) / move_window + 1)-1;
thres_noise = 1e-7;

num_data = 50;
I_total_XY = zeros(num_data, num_window, size_window, size_window,2);
X_trans = zeros(num_data, num_window, size_window);

%% compute RDF
disp('compute RDF')
for i = 1:num_data
    
    % import time series
    y_tmp = cell2mat(y_total(i));
    
    % cut using moving window technique
    for j = 1:num_window
        % import single time series
        y_target_tmp = y_tmp(st+1 + move_window*(j-1) : st+size_window + move_window*(j-1),:);
        y_target = y_target_tmp;
        t_target = tspan(st+1 + move_window*(j-1) : st+size_window + move_window*(j-1));
        
        % compute RDf for the regulation from X to Y
        [score_list, t_1, t_2] = RDS_ns_dim1(y_target(:,1), y_target(:,2), t_target, time_interval);
        for k = 1:2
            score = reshape(score_list(:,:,k),[length(t_target),length(t_target)]);
            I_total_XY(i,j,:,:,k) = score;            
        end
        
        % save the value of X and Y
        X_trans(i,j,:) = y_target(:,1);
    end
end

%% transform (t,t*) into X(t),X(t*)
disp('transformation')
trans_total_XY = cell(num_data, num_window,2);

for i = 1:num_data
    for j = 1:num_window
        score_XY_pos = reshape(I_total_XY(i,j,:,:,1), [size_window, size_window]);
        score_XY_neg = reshape(I_total_XY(i,j,:,:,2), [size_window, size_window]);
        
        X_var = reshape(X_trans(i,j,:), [size_window,1]);
        
        trans_XY_pos = [];
        trans_XY_neg = [];
        
        for m = 1:size_window
            for n = 1:size_window
                if abs(score_XY_pos(m,n)) > thres_noise
                    trans_XY_pos = [trans_XY_pos ; [X_var(m), X_var(n), score_XY_pos(m,n)]];
                end
                if abs(score_XY_neg(m,n)) > thres_noise
                    trans_XY_neg = [trans_XY_neg ; [X_var(m), X_var(n), score_XY_neg(m,n)]];
                end
            end
        end        
        trans_total_XY{i,j,1} = trans_XY_pos;
        trans_total_XY{i,j,2} = trans_XY_neg;
    end
end

%% plot
disp('plot')

for j = 1:num_window
    for i = 1:num_data
    
        trans_XY_pos = cell2mat(trans_total_XY(i,j,1));
        trans_XY_neg = cell2mat(trans_total_XY(i,j,2));
        
        if isempty(trans_XY_pos)
            continue
        end
        
        if isempty(find(trans_XY_pos(:,3) > 0)) || isempty(find(trans_XY_pos(:,3) < 0))
            continue
        end
        if isempty(trans_XY_neg)
            continue
        end
        
        if isempty(find(trans_XY_neg(:,3) > 0)) || isempty(find(trans_XY_neg(:,3) < 0))
            continue
        end
        
        loca_pos = find(trans_XY_pos(:,3) > 0);
        loca_neg = find(trans_XY_pos(:,3) < 0);
        
        if j == 1
            alpha = 0.01;
        elseif j == 3
            alpha = 0.02;
        else
            alpha = 0.005;
        end
        
        figure(j)
        scatter(trans_XY_pos(loca_pos,1),trans_XY_pos(loca_pos,2),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha)
        hold on
        scatter(trans_XY_pos(loca_neg,1),trans_XY_pos(loca_neg,2),'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha)
        hold on
        
        loca_pos = find(trans_XY_neg(:,3) > 0);
        loca_neg = find(trans_XY_neg(:,3) < 0);

        scatter(trans_XY_neg(loca_pos,1),trans_XY_neg(loca_pos,2),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha)
        hold on
        scatter(trans_XY_neg(loca_neg,1),trans_XY_neg(loca_neg,2),'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha)
        hold on
        %plot([0,1],[0,1],'--k', 'LineWidth',2)

        xlim([0,1])
        ylim([0,1])
        xticks([0,1])
        yticks([0,1])
        set(gca,'XColor', 'none','YColor','none')
        axis square
    end
    filename = ['figure',num2str(j),'.png'];
    saveas(gcf,filename)
end
