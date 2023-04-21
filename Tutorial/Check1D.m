clc;
clear;
close all;
addpath('GOBI')

%% load data
load('data_with_options')

%% parameter
st = 2; % choose causal variable
ed = 1; % choose target variable
type = 2;

thres_noise = 0;
size_window = length(t);
X_trans = zeros(num_data, size_window);

%% compute RDF
disp('compute RDF')

I_total = zeros(num_data, size_window, size_window);

for i = 1:num_data
    % import time series
    y_target = cell2mat(y_total(i));
    t_target = t;
    
    % compute RDF for the regulation
    if type_self == -1
        [score_list, t_1, t_2] = RDS_ns_dim1(y_target(:,st), y_target(:,ed), t_target, time_interval);
    elseif type_self == 1
        [score_list, t_1, t_2] = RDS_ps_dim1(y_target(:,st), y_target(:,ed), t_target, time_interval);
    else
        [score_list, t_1, t_2] = RDS_dim1(y_target(:,st), y_target(:,ed), t_target, time_interval);
    end
    
    score = reshape(score_list(:,:,type),[length(t_target),length(t_target)]);
    I_total(i,:,:) = score;
end

%% linearlize into (t, t*, RDF(t,t*))
disp('transformation')
trans_total = cell(num_data,1);

for i = 1:num_data
    rdf_tmp = reshape(I_total(i,:,:), [size_window, size_window]);
    
    trans_tmp = [];

    for m = 1:size_window
        for n = 1:size_window
            if abs(rdf_tmp(m,n)) > thres_noise
                trans_tmp = [trans_tmp ; [t(m), t(n), rdf_tmp(m,n)]];
            end
        end
    end        
    trans_total{i,1} = trans_tmp;
end

%% plot
disp('plot')
alpha = 0.05;

total_space = zeros(length(t));
figure(1)
for i = 1:num_data
    
    trans_tmp = cell2mat(trans_total(i,1));

    if isempty(trans_tmp)
        continue
    end

    loca_pos = find(trans_tmp(:,3) > 0);
    loca_neg = find(trans_tmp(:,3) < 0);
    
    color_tmp = linspace(0,1,101).';
    color_list = [color_tmp, color_tmp, ones(101,1)];
    
    scatter(trans_tmp(loca_pos,1),trans_tmp(loca_pos,2),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha)
    hold on
    scatter(trans_tmp(loca_neg,1),trans_tmp(loca_neg,2),'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha)
    hold on
    
    for i = 1:length(loca_pos) 
        m = find(t == trans_tmp(loca_pos(i),1));
        n = find(t == trans_tmp(loca_pos(i),2));
        total_space(m,n) = 1;
    end
    
    xlim([0,1])
    ylim([0,1])
    xticks([0,1])
    yticks([0,1])
    set(gca,'XColor', 'none','YColor','none')
    axis square
end
filename = ['Check_1D_3.png'];
%saveas(gcf,filename)


%% warning if the space for RDF is insufficient
filled = length(find(total_space == 1));
filled_ratio = filled / (length(t).^2);

if filled_ratio < 0.8
    disp('data is insufficient for confident inference')
else
end