clc;
clear;
close all;
addpath('../GOBI') 
%% parameter
system_list = {'KF','Fr','GW','Gb','cAMP'};
num_system = length(system_list);
noise_list = [0:2:20];
acceptance = 0;
trial_list = [0:9];

tar_noise = [1,6,11];
c = {'k','b','r'};

%% find optimal threshold for S, L which maximize F2 score using exponential fitting
tmp = 10; % fix threshold for regulation-detection region (0.1)
idx_list = [];

for i = 1
    for j = 1:length(tar_noise)
        
        loc_list_S = [];
        loc_list_L = [];
        loc_list_TRS = [];

        % find threshold which maximize F2 score    
        for trial = trial_list
            % load data
            system_name = string(system_list(i));
            filename = append('./',system_name,'/F2_total_',system_name,'_Trial',num2str(trial));
            load(filename)
            
            F_tmp = F2_list(:,tmp,:,tar_noise(j));
            F_max = max(max(F_tmp));
            F_min = min(min(F_tmp));

            F1_thres = F_max - acceptance * (F_max - F_min);
            loc1 = [];
            loc2 = [];
            for x1 = 1:length(thres_S_list)
                for x2 = 1:length(thres_TRS_list)
                    if F_tmp(x1,x2) >= F1_thres
                        loc1 = [loc1;x1];
                        loc2 = [loc2;x2];
                    end
                end
            end

            
            loc_list_S = [loc_list_S, thres_S_list(loc1)];
            loc_list_TRS = [loc_list_TRS, thres_TRS_list(loc2)];
            
        end
        % exponential fitting
        x = loc_list_S;
        y = loc_list_TRS;
        f = @(a,x) a(1) .* exp(a(2) .* x) + a(3);
        %f = @(a,x) a(1) .* (1 - x.^a(2)).^(0.5) + a(3);
        options = optimset('MaxFunEvals',1e15, 'MaxIter',1e15 ,'TolFun',1e-15,'TolX',1e-20);
        A = fminsearch(@(a) norm(y-f(a,x)), [-1;1;2],options);
        %A = fminsearch(@(a) norm(y-f(a,x)), [1;2;0],options);
        
        figure(j)
        scatter(x, y,'filled','MarkerFaceColor',string(c(j)),'MarkerFaceAlpha',0.1)
        xlim([0,1])
        ylim([0,1])
        xticks([0,0.5,1])
        yticks([0,0.5,1])
        xticklabels([])
        yticklabels([])
        
        figure(4)
        domain_x = [0:0.01:1];
        scatter(x, y,125,'s','filled','MarkerFaceColor',string(c(j)),'MarkerFaceAlpha',0.02)
        hold on
        patchline(domain_x, A(1) .* exp(A(2) .* domain_x) + A(3),'linestyle','-','edgecolor', string(c(j)), 'Linewidth',2,'edgealpha',0.8)
        hold on
        xlim([0.5,1])
        ylim([0.5,1])
        xticks([0.5,1])
        yticks([0.5,1])
        xticklabels([])
        yticklabels([])
        
    end
end
