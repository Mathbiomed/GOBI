%% Step 4. Plot the result using 

data_B = readmatrix('data(b).csv');
result_data = readmatrix('simul_result.csv');
mean_a = result_data(1,:);
std_a = result_data(2,:);

t=0:0.01:50; % time span

figure(4)
x1 = [t fliplr(t)];
x2 = [ mean_a+std_a fliplr(mean_a-std_a)];
figure(1)
plot(t,data_B, '-k','LineWidth', 2);
hold on;
plot(t,mean_a, '-b','LineWidth', 2);
fill(x1,x2,'b','FaceAlpha',.3);
hold off;
title('time series of component A (input)')

