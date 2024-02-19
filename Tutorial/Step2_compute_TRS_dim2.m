clc;
clear;
close all;
addpath('GOBI')

%% load data
load('data_with_options')
load('RDS_dim2')

%% parameter
if thres_R == 0
    thres_R = 0.05;
end
dimension = 2;
%% compute TRS for each pair of 1D regulation
TRS_total = zeros(num_pair, num_type);
for i = 1:num_pair
    for j = 1:num_type
        S_tmp = reshape(S_total_list(i,j,1:end),[num_data,1]); % import regulation-detection score
        R_tmp = reshape(R_total_list(i,j,1:end),[num_data,1]); % import regulation-detection region
        
        S_processed = S_threshold(S_tmp, thres_S); % test whether S > S^thres
        R_processed = R_threshold(R_tmp, thres_R); % test whether R < R_thres
        
        if sum(R_processed) == 0
            TRS_tmp = nan;
        else
            TRS_tmp = sum(S_processed .* R_processed) / sum(R_processed); % compute TRS
        end
        TRS_total(i,j) = TRS_tmp;
    end
end

%% infer 2D regulation using the criteria TRS > TRS^thres
regulation_2dim = zeros(num_pair,2^dimension);
for i = 1:num_pair
    for j = 1:num_type
        if TRS_total(i,j) >= thres_TRS
            regulation_2dim(i,j) = 1;
        end
    end
end

%% save TRS & inferred regulation
TRS_total_dim2 = TRS_total;

filename = ['TRS_dim2'];
save(filename, 'regulation_2dim','TRS_total_dim2','component_list_dim2')

%% plot heatmap

c_nan = [0.9 0.9 0.9];
font_s = 14;
tmp = 1 - linspace(0, 1, 501)';
cmap_score = [[tmp],[tmp],[ones(501,1)]];


h = heatmap(TRS_total);

TRS_range = [0,1];
h.FontName = 'Arial';
h.Colormap = cmap_score;
h.ColorLimits = TRS_range;
h.FontSize = font_s;
h.MissingDataColor = c_nan;
h.CellLabelColor = 'none';
s = struct(h);
cbh = s.Colorbar;
set(cbh, 'YTick', TRS_range)