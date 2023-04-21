clc;
clear;
close all;
addpath('../GOBI') 
%% parameter
system_list = {'KF','Fr','GW','Gb','cAMP'};
num_system = length(system_list);
noise_list = [0:2:20];
acceptance = 0.05;
trial_list = [0:9];
thres_S_list = [0:0.01:1];
thres_L_list = [0.01:0.01:0.1];
thres_TRS_list = [0:0.01:1];
noise_list = [10];
%trial_list = [0];


c_nan = [0.9 0.9 0.9];

font_s = 14;
tmp = linspace(0, 1, 501)';
cmap_score = [[ones(500,1);1-tmp],[tmp(1:end-1);1-tmp],[tmp; ones(500,1)]];
%cmap_value = [[ones(500,1);1-tmp],ones(1001,1),[tmp; ones(500,1)]];
cmap_value = [[ones(500,1);1-tmp],[tmp*0.3+0.7; (1-tmp(1:end-1))*0.6+0.4],[tmp(1:end-1);1-tmp]];



%% threshold for S, L which maximize F2 score
loc_list_S = [];
loc_list_L = [];
loc_list_TRS = [];

idx_list = [];
for i = 1:5
    for j = 1:length(noise_list)
        % find the maximum F score 
        max_F_list_tmp = zeros(length(thres_S_list),length(thres_L_list),length(trial_list));
        
        for trial_idx = 1:length(trial_list)
            % load data
            system_name = string(system_list(i));
            
            filename = append('./',system_name,'/F2_total_',system_name,'_Trial',num2str(trial_list(trial_idx)));
            load(filename)
            
            F_score = F2_list;
            
            % merge F score for each pair of thres_L and thres_S
            for x1 = 1:length(thres_S_list)
                for x2 = 1:length(thres_L_list)
                    F_tmp = reshape(F_score(x1,x2,:,6),[length(thres_TRS_list),1]);
                    max_F_list_tmp(x1,x2,trial_idx) = max(F_tmp);
                end
            end
        end
        max_F_list = zeros(length(thres_S_list)-1,length(thres_L_list));
        for x1 = 1:length(thres_S_list)-1
            for x2 = 1:length(thres_L_list)
                max_F_list(x1,x2) = mean(max_F_list_tmp(x1,x2,:));
            end
        end
        
        figure(i)
        h = heatmap(flipud(max_F_list));

       
        h.GridVisible = 'off';
        h.FontName = 'Arial';
        h.Colormap = cmap_score;
        if i == 1
            h.ColorLimits = [0.9 1];
        elseif i == 2
            h.ColorLimits = [0.8 1];
        elseif i == 3
            h.ColorLimits = [0.8 1];
        elseif i == 4
            h.ColorLimits = [0.7 0.8];
        elseif i == 5
            h.ColorLimits = [0.9 1];
        end
        h.FontSize = font_s;
        h.XLabel = 'thres_L';
        h.YLabel = 'thres_S';

    end
 
end