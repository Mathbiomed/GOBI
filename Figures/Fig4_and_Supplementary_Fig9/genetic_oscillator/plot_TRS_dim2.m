clc;
clear;
close all;
addpath('../../GOBI')

%% parameter
dimension = 2;
noise_level = 8; 

thres_S = 0.9 - 0.005 * noise_level;
thres_L = 0.05;
%% load data
load('RDS_dim2')

%% compute TRS
num_type = 2^dimension;
TRS_total = zeros(num_pair, num_type);
for i = 1:num_pair
    for j = 1:num_type
        S_tmp = reshape(S_total_list(i,j,1:end),[num_data,1]);
        L_tmp = reshape(L_total_list(i,j,1:end),[num_data,1]);
        
        S_processed = S_threshold(S_tmp, thres_S);
        L_processed = L_threshold(L_tmp, thres_L);
        
        if sum(L_processed) == 0
            TRS_tmp = 0;
        else
            TRS_tmp = sum(S_processed .* L_processed) / sum(L_processed);
        end
        TRS_total(i,j) = TRS_tmp;
    end
end

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
