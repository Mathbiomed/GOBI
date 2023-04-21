clc;
clear;
close all;

addpath('../GOBI') 

%% load data
load('data_temporal.mat')

%% parameter
size_window = 100;
move_window = 50;
num_window = (range / time_interval - size_window) / move_window + 1;
thres_noise = 1e-7;

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
        y_target_tmp = y_tmp(1 + move_window*(j-1) : size_window + move_window*(j-1),:);
        y_target = y_target_tmp;
        t_target = tspan(1 + move_window*(j-1) : size_window + move_window*(j-1));
        
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
disp('transform into (t,t^*) into X(t),X(t^*)')

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

%% probability distribution
disp('compute pdf')
num_bin = 4;
plus_count = zeros(num_window,2*num_bin.^2);
minus_count = zeros(num_window,2*num_bin.^2);
for i = 1:num_data
    for j = 1:num_window
        trans_XY_pos = cell2mat(trans_total_XY(i,j,1));
        for k = 1:length(trans_XY_pos(:,1))
            X_var1 = trans_XY_pos(k,1);
            X_var2 = trans_XY_pos(k,2);
            score_var = trans_XY_pos(k,3);
            
            m = fix(X_var1/(1/num_bin));
            n = fix(X_var2/(1/num_bin));
         
            if m == num_bin
                m = m-1;
            end
            if n == num_bin
                n = n-1;
            end
            idx = 2*num_bin*n + 2*m;
            
            if X_var1-m/num_bin > X_var2-n/num_bin
                idx = idx+2;
            else
                idx = idx+1;
            end
            
            if score_var > 0
                plus_count(j,idx) = plus_count(j,idx)+1;
            else
                minus_count(j,idx) = minus_count(j,idx)+1;
            end
        end
        
        trans_XY_neg = cell2mat(trans_total_XY(i,j,2));
        for k = 1:length(trans_XY_neg(:,1))
            X_var1 = trans_XY_neg(k,1);
            X_var2 = trans_XY_neg(k,2);
            score_var = trans_XY_neg(k,3);
            
            m = fix(X_var1/(1/num_bin));
            n = fix(X_var2/(1/num_bin));
            
            if m == num_bin
                m = m-1;
            end
            if n == num_bin
                n = n-1;
            end
            idx = 2*num_bin*n + 2*m;
            
            if X_var1-m/num_bin > X_var2-n/num_bin
                idx = idx+2;
            else
                idx = idx+1;
            end
            
            if score_var > 0
                plus_count(j,idx) = plus_count(j,idx)+1;
            else
                minus_count(j,idx) = minus_count(j,idx)+1;
            end
        end
    end
end

for i = 1:num_window
    plus_count(i,:) = plus_count(i,:) / sum(plus_count(i,:)) + eps;
    minus_count(i,:) = minus_count(i,:) / sum(minus_count(i,:)) + eps;
end

%% variational information
disp('compute KL divergence')

KL_list = zeros(num_window,1);
X = [1:2*num_bin^2].';
for i = 1:num_window
    KL_list(i) = kldiv(X,plus_count(i,:).', minus_count(i,:).','js');
end

figure(1)
plot([1:num_window].', KL_list, 'r','LineWidth',2)
xlim([0,40])
xticks([0:10:40])
ylim([0,0.75])
yticks([0,0.25,0.5,0.75])

%% MESH
figure(3)
x_axis = [-2/3, -1/3, 1/3, 2/3, 4/3, 5/3, 7/3, 8/3, 10/3, 11/3, 13/3, 14/3];
x_axis = [x_axis, x_axis, x_axis, x_axis, x_axis, x_axis];
y_axis = [2/3, 1/3, 2/3, 1/3, 2/3, 1/3, 2/3, 1/3, 2/3, 1/3, 2/3, 1/3];
y_axis = [y_axis-1, y_axis, y_axis+1, y_axis+2, y_axis+3, y_axis+4];

plus_count = plus_count(16,:);
minus_count = minus_count(16,:);
plus_count = [zeros(1,12), [0,0,plus_count(1:8),0,0], [0,0,plus_count(9:16),0,0], [0,0,plus_count(17:24),0,0], [0,0,plus_count(25:32),0,0], zeros(1,12)];
minus_count = [zeros(1,12), [0,0,minus_count(1:8),0,0], [0,0,minus_count(9:16),0,0], [0,0,minus_count(17:24),0,0], [0,0,minus_count(25:32),0,0], zeros(1,12)];

[xq,yq] = meshgrid(-1:1/20:5, -1:1/20:5);
vq = griddata(x_axis,y_axis,minus_count,xq,yq);  %(x,y,v) being your original data for plotting points
mesh(xq,yq,vq, 'EdgeColor','None','FaceColor', 'flat')
xlim([0,4])
ylim([0,4])
xticks([0,4])
yticks([0,4])
colorbar
set(gca,'XColor', 'none','YColor','none')
axis square