clc;
clear;
close all;
addpath('../GOBI') 

%% load data
filename = ['RDS_dim1'];
load(filename)

%% compute TRS

figure(1)
S_box = [];
for i = 1:num_pair
    for j = 1:num_type
        S = reshape(S_total(i,j,:),[num_data,1]);
        S_box = [S_box, S];
    end
end
boxplot(S_box)

ylim([-1,1])
yticks([-1,0,1])
ylabel('Regulation detection score')

xticklabels({'A->A','A-|A', 'A->B','A-|B', 'B->A','B-|A', 'B->B','B-|B'})
xlabel('regulations')

title('1D regulation')

set(gca,'FontSize',15)

save('fig_1D', 'S_box')