clc;
clear;
close all;

%% load data
filename = ['data_repressilator_raw'];
load(filename)

t = t / t(end);

%% change sampling rate
y_1 = y;
y_2 = y(1:2:end,:);
y_3 = y(1:4:end,:);

t_1 = t;
t_2 = t(1:2:end);
t_3 = t(1:4:end);

%% plot
figure(1)
subplot(3,1,1)
plot(t_1, y_1(:,1), '-or', 'MarkerFaceColor', 'r', 'Markersize',5')
hold on
plot(t_1, y_1(:,2), '-ob', 'MarkerFaceColor', 'b', 'Markersize',5')
hold on
plot(t_1, y_1(:,3), '-og', 'MarkerFaceColor', 'g', 'Markersize',5')

ylim([0,1])
yticks([0,1])
xlim([0,1])
xticks([0,1])

subplot(3,1,2)
plot(t_2, y_2(:,1), '-or', 'MarkerFaceColor', 'r', 'Markersize',5')
hold on
plot(t_2, y_2(:,2), '-ob', 'MarkerFaceColor', 'b', 'Markersize',5')
hold on
plot(t_2, y_2(:,3), '-og', 'MarkerFaceColor', 'g', 'Markersize',5')

ylim([0,1])
yticks([0,1])
xlim([0,1])
xticks([0,1])

subplot(3,1,3)
plot(t_3, y_3(:,1), '-or', 'MarkerFaceColor', 'r', 'Markersize',5')
hold on
plot(t_3, y_3(:,2), '-ob', 'MarkerFaceColor', 'b', 'Markersize',5')
hold on
plot(t_3, y_3(:,3), '-og', 'MarkerFaceColor', 'g', 'Markersize',5')

ylim([0,1])
yticks([0,1])
xlim([0,1])
xticks([0,1])

%% interpolation
time_interval = 1/200;
t_fit = linspace(0,1,1/time_interval+1).';

y_1_int = zeros(length(t_fit),2);
y_2_int = zeros(length(t_fit),2);
y_3_int = zeros(length(t_fit),2);

for i = 1:3
    y_1_int(:,i) = interp1(t_1, y_1(:,i),t_fit,'spline');
    y_2_int(:,i) = interp1(t_2, y_2(:,i),t_fit,'spline');
    y_3_int(:,i) = interp1(t_3, y_3(:,i),t_fit,'spline');
end

for i = 1:3
    y_1_int(:,i) = y_1_int(:,i) - min(y_1_int(:,i));
    y_1_int(:,i) = y_1_int(:,i) / max(y_1_int(:,i));
    y_2_int(:,i) = y_2_int(:,i) - min(y_2_int(:,i));
    y_2_int(:,i) = y_2_int(:,i) / max(y_2_int(:,i));
    y_3_int(:,i) = y_3_int(:,i) - min(y_3_int(:,i));
    y_3_int(:,i) = y_3_int(:,i) / max(y_3_int(:,i));
end

%% merge data
t = t_fit;

y = y_1_int;
save('data_merge_full','y','t')

y = y_2_int;
save('data_merge_half','y','t')

y = y_3_int;
save('data_merge_quarter','y','t')


%% plot
figure(2)
subplot(3,1,1)
plot(t, y_1_int(:,1), 'r')
hold on
plot(t, y_1_int(:,2), 'b')
hold on
plot(t, y_1_int(:,3), 'm')

subplot(3,1,2)
plot(t, y_2_int(:,1), 'r')
hold on
plot(t, y_2_int(:,2), 'b')
hold on
plot(t, y_2_int(:,3), 'm')

subplot(3,1,3)
plot(t, y_3_int(:,1), 'r')
hold on
plot(t, y_3_int(:,2), 'b')
hold on
plot(t, y_3_int(:,3), 'm')
