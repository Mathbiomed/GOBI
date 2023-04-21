clc;
clear;
close all;

%% choose the data
idx = 4;

%% load data
filename_1 = ['data_gene_raw_', num2str(idx)];
load(filename_1)
y_gene = y.';
t_gene = t / t(end);

filename_2 = ['data_prey_predator_raw'];
load(filename_2)
y_prey = [y(:,2),y(:,1)];
t_prey = t / t(end);

%% change sampling rate
y_gene_1 = y_gene;
y_gene_2 = y_gene(1:2:end,:);
y_gene_3 = y_gene(1:4:end,:);

t_gene_1 = t_gene;
t_gene_2 = t_gene(1:2:end);
t_gene_3 = t_gene(1:4:end);

y_prey_1 = y_prey;
y_prey_2 = y_prey(1:2:end,:);
y_prey_3 = y_prey(1:4:end,:);

t_prey_1 = t_prey;
t_prey_2 = t_prey(1:2:end);
t_prey_3 = t_prey(1:4:end);

t_gene_1 = t_gene_1 / t_gene_1(end);
t_gene_2 = t_gene_2 / t_gene_2(end);
t_gene_3 = t_gene_3 / t_gene_3(end);

t_prey_1 = t_prey_1 / t_prey_1(end);
t_prey_2 = t_prey_2 / t_prey_2(end);
t_prey_3 = t_prey_3 / t_prey_3(end);
%% interpolation
time_interval = 1/700;
t_fit = linspace(0,1,1/time_interval+1).';

y_gene_1_int = zeros(length(t_fit),2);
y_gene_2_int = zeros(length(t_fit),2);
y_gene_3_int = zeros(length(t_fit),2);

y_prey_1_int = zeros(length(t_fit),2);
y_prey_2_int = zeros(length(t_fit),2);
y_prey_3_int = zeros(length(t_fit),2);

for i = 1:2
    y_gene_1_int(:,i) = interp1(t_gene_1, y_gene_1(:,i),t_fit,'spline');
    y_gene_2_int(:,i) = interp1(t_gene_2, y_gene_2(:,i),t_fit,'spline');
    y_gene_3_int(:,i) = interp1(t_gene_3, y_gene_3(:,i),t_fit,'spline');
    
    y_prey_1_int(:,i) = interp1(t_prey_1, y_prey_1(:,i),t_fit,'spline');
    y_prey_2_int(:,i) = interp1(t_prey_2, y_prey_2(:,i),t_fit,'spline');
    y_prey_3_int(:,i) = interp1(t_prey_3, y_prey_3(:,i),t_fit,'spline');
end

for i = 1:2
    y_gene_1_int(:,i) = y_gene_1_int(:,i) - min(y_gene_1_int(:,i));
    y_gene_1_int(:,i) = y_gene_1_int(:,i) / max(y_gene_1_int(:,i));
    y_gene_2_int(:,i) = y_gene_2_int(:,i) - min(y_gene_2_int(:,i));
    y_gene_2_int(:,i) = y_gene_2_int(:,i) / max(y_gene_2_int(:,i));
    y_gene_3_int(:,i) = y_gene_3_int(:,i) - min(y_gene_3_int(:,i));
    y_gene_3_int(:,i) = y_gene_3_int(:,i) / max(y_gene_3_int(:,i));
    
    y_prey_1_int(:,i) = y_prey_1_int(:,i) - min(y_prey_1_int(:,i));
    y_prey_1_int(:,i) = y_prey_1_int(:,i) / max(y_prey_1_int(:,i));
    y_prey_2_int(:,i) = y_prey_2_int(:,i) - min(y_prey_2_int(:,i));
    y_prey_2_int(:,i) = y_prey_2_int(:,i) / max(y_prey_2_int(:,i));
    y_prey_3_int(:,i) = y_prey_3_int(:,i) - min(y_prey_3_int(:,i));
    y_prey_3_int(:,i) = y_prey_3_int(:,i) / max(y_prey_3_int(:,i));
end

%% merge data
t = t_fit;
% 
% y = [y_gene_1_int,y_prey_1_int];
% save('data_merge_full','y','t')
% 
% y = [y_gene_2_int,y_prey_2_int];
% save('data_merge_half','y','t')
% 
% y = [y_gene_3_int,y_prey_3_int];
% save('data_merge_quarter','y','t')


%% plot
figure(1)
yyaxis right
plot(t_gene_1, y_gene_1(:,1) , 'r')
hold on
plot(t_gene_1, y_gene_1(:,2), 'b')
ylim([0,200])
yticks([0,200])
set(gca,'YColor','k');

yyaxis left
plot(t_prey_1, y_prey_1(:,1), 'm')
hold on
plot(t_prey_1, y_prey_1(:,2), 'g')
ylim([0,350])
yticks([0,350])
set(gca,'YColor','k');

xlim([0,1])
xticks([0,1])

figure(2)
subplot(3,1,1)

plot(t_gene_1, y_gene_1(:,1) , '-or', 'MarkerFaceColor', 'r', 'Markersize',5)
hold on
plot(t_gene_1, y_gene_1(:,2), '-ob', 'MarkerFaceColor', 'b', 'Markersize',5)
hold on
plot(t_prey_1, y_prey_1(:,1)* 200 / 350, '-om', 'MarkerFaceColor', 'm', 'Markersize',5)
hold on
plot(t_prey_1, y_prey_1(:,2)* 200 / 350, '-og', 'MarkerFaceColor', 'g', 'Markersize',5)
ylim([0,200])
yticks([0,200])
set(gca,'YColor','k');

xlim([0,1])
xticks([0,1])

subplot(3,1,2)

plot(t_gene_2, y_gene_2(:,1) , '-or', 'MarkerFaceColor', 'r', 'Markersize',5)
hold on
plot(t_gene_2, y_gene_2(:,2), '-ob', 'MarkerFaceColor', 'b', 'Markersize',5)
hold on
plot(t_prey_2, y_prey_2(:,1)* 200 / 350, '-om', 'MarkerFaceColor', 'm', 'Markersize',5)
hold on
plot(t_prey_2, y_prey_2(:,2)* 200 / 350, '-og', 'MarkerFaceColor', 'g', 'Markersize',5)

ylim([0,200])
yticks([0,200])
set(gca,'YColor','k');

xlim([0,1])
xticks([0,1])

subplot(3,1,3)
plot(t_gene_3, y_gene_3(:,1) , '-or', 'MarkerFaceColor', 'r', 'Markersize',5)
hold on
plot(t_gene_3, y_gene_3(:,2), '-ob', 'MarkerFaceColor', 'b', 'Markersize',5)
hold on
plot(t_prey_3, y_prey_3(:,1) * 200 / 350, '-om', 'MarkerFaceColor', 'm', 'Markersize',5)
hold on
plot(t_prey_3, y_prey_3(:,2) * 200 / 350, '-og', 'MarkerFaceColor', 'g', 'Markersize',5)

ylim([0,200])
yticks([0,200])
set(gca,'YColor','k');

xlim([0,1])
xticks([0,1])
