clc;
clear;
close all;
addpath('../GOBI') 
%% parameter
system_list = {'KF','Fr','GW','Gb','cAMP'};
num_system = length(system_list);
noise_list = [0:2:20];
acceptance = 0.05; % use threshold which make 95% of maximum F2 score
trial_list = [0:9];

%% threshold for S, L which maximize DB index
loc_list_S = [];
loc_list_L = [];
loc_list_TRS = [];

tmp = 10; % fix threshold for regulation-detection region (0.1)

thres_S = [];
thres_TRS = [];
for i = 1:5
    disp(i)
    thres_S_tmp = [];
    thres_TRS_tmp = [];
    for j = 1:length(noise_list)
        loc_list_S = [];
        loc_list_TRS = [];
        for trial = trial_list
            % load data
            system_name = string(system_list(i));
            
            filename = append('./',system_name,'/F2_total_',system_name,'_Trial',num2str(trial));
            
            load(filename)

            % compute the maximum of F score 
            F_tmp = F2_list(:,tmp,:,j);
            F_max = max(max(F_tmp));
            F_min = min(min(F_tmp));

            F_thres = F_max - acceptance * (F_max - F_min);
            
            % find threshold which maximize F score
            loc1 = [];
            loc2 = [];
            for x1 = 1:length(thres_S_list)
                for x2 = 1:length(thres_TRS_list)
                    if F_tmp(x1,x2) >= F_thres
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
        options = optimset('MaxFunEvals',1e15, 'MaxIter',1e15 ,'TolFun',1e-20,'TolX',1e-20);
        A = fminsearch(@(a) norm(y-f(a,x)), [-1;1;2],options);
       

        domain_x = [0:0.01:1];
        range_y = A(1) .* exp(A(2) .* domain_x) + A(3);
        distance = sqrt((range_y - 1).^2 + (domain_x - 1).^2);
        [val, loc] = min(distance);
        thres_S_tmp = [thres_S_tmp ; domain_x(median(loc))];
        thres_TRS_tmp = [thres_TRS_tmp ; range_y(median(loc))];
    end
    
    thres_S = [thres_S, thres_S_tmp];
    thres_TRS = [thres_TRS, thres_TRS_tmp];
    
end
% 
figure(13)
for i = 1:5
    plot(noise_list, thres_S(:,i),'-o','DisplayName',string(system_list(i)))
    hold on
end
plot(noise_list, 0.9 - noise_list * 0.005,'r', 'LineWidth',2)
    
xlim([-1,21])
ylim([0.5,1])
xticks([0:2:20])
yticks([0.5,0.75,1])
%xticklabels([])
%yticklabels([])
legend('Location','southwest')
title('threshold for score')

figure(14)
for i = 1:5
    plot(noise_list, thres_TRS(:,i),'-o','DisplayName',string(system_list(i)))
    hold on
end
plot(noise_list, 0.9 - noise_list * 0.01, 'r', 'LineWidth',2)

xlim([-1,21])
ylim([0.5,1])
xticks([0:2:20])
yticks([0.5,0.75,1])
%xticklabels([])
%yticklabels([])
legend('Location','southwest')
title('threshold for TRS')