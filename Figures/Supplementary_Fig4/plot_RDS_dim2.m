clc;
clear;
close all;
addpath('../GOBI') 

%% load data
filename = ['RDS_dim2'];
load(filename)

%% compute TRS

figure(1)
S_box = [];
for i = 1:num_pair
    for j = 1:num_type
        S = reshape(S_total(i,j,:),[num_data,1]);
        L = reshape(L_total(i,j,:),[num_data,1]);
        
        L_processed = L_threshold(L, 0.01);
        S = S .* L_processed;
        S(find(S == 0)) = nan;
        
        S_box = [S_box, S];
    end
end

boxplot(S_box)
ylim([-1,1])
yticks([-1,0,1])
ylabel('Regulation detection score')

xticklabels({'A->A & B->A', 'A->A & B-|A', 'A-|A & B->A', 'A-|A & B-|A','A->B & B->B', 'A->B & B-|B', 'A-|B & B->B', 'A-|B & B-|B'})
xlabel('regulations')

title('2D regulation')

set(gca,'FontSize',15)

save('fig_2D', 'S_box')