clc;
clear;
close all;


%% initial 
time_interval = 1/10;
period = 10;
tspan = linspace(0,period, period/time_interval+1);
initial = [0.4,0.24,0.7,0.7,0.73,0.44,0.47];
idx = 4;

[t_ori, y_ori] = ode45(@(t,y) cAMP(t,y), tspan, initial);
y_ori = y_ori(:,idx);
min_ori = min(y_ori);
y_ori = y_ori - min_ori;
max_ori = max(y_ori);
y_ori = y_ori / max_ori;

%% simulate dynamical noise
initial_tmp = ones(1,7);
period = 100;
tspan_tmp = linspace(0,period, period/time_interval+1);

[t_ori_tmp, y_ori_tmp] = ode45(@(t,y) cAMP(t,y), tspan_tmp, initial_tmp);

% find period using peaks
[~,locs]=findpeaks(y_ori_tmp(:,1));
st_period = locs(end-1);
ed_period = locs(end);
mean_range = zeros(7,1);
for i = 1:7
    mean_range(i,1) = mean(y_ori_tmp(st_period:ed_period,i));
end

[t_dyn, y_dyn] = ode45(@(t,y) cAMP_noise(t,y,0.2,mean_range(1),mean_range(2),mean_range(3),mean_range(4),mean_range(5),mean_range(6),mean_range(7)), tspan, initial);
y_dyn = y_dyn(:,idx);
y_dyn = y_dyn - min_ori;
y_dyn = y_dyn / max_ori;

%% create noise
noise_white = 0.2 * randn([101,1]);
cn = dsp.ColoredNoise('blue',101,1,OutputDataType='single');
noise_blue = 0.2 * cn();
cn = dsp.ColoredNoise('pink',101,1,OutputDataType='single');
noise_pink = 0.2 * cn();

%% create noisy time series
y_multi = y_ori .* (1+noise_white);
y_add = y_ori + noise_white;
y_pink = y_ori .* (1+noise_pink);
y_blue = y_ori .* (1+noise_blue);

%% plot
figure(1)
plot(t_ori, y_ori, 'k')
hold on
plot(t_ori, y_multi, 'r')
xlim([0,10])
xticks([0,10])
xticklabels([])
ylim([-0.2,1.2])
yticks([0,1])
yticklabels([])

figure(2)
plot(t_ori, y_ori, 'k')
hold on
plot(t_ori, y_add, 'r')
xlim([0,10])
xticks([0,10])
xticklabels([])
ylim([-0.2,1.2])
yticks([0,1])
yticklabels([])

figure(3)
plot(t_ori, y_ori, 'k')
hold on
plot(t_ori, y_pink, 'r')
xlim([0,10])
xticks([0,10])
xticklabels([])
ylim([-0.2,1.2])
yticks([0,1])
yticklabels([])

figure(4)
plot(t_ori, y_ori, 'k')
hold on
plot(t_ori, y_blue, 'r')
xlim([0,10])
xticks([0,10])
xticklabels([])
ylim([-0.2,1.2])
yticks([0,1])
yticklabels([])

figure(5)
plot(t_ori, y_ori, 'k')
hold on
plot(t_ori, y_dyn, 'r')
xlim([0,10])
xticks([0,10])
xticklabels([])
ylim([-0.5,1.5])
yticks([0,1])
yticklabels([])

t = t_ori;
save('data_example', 't_ori', 'y_ori', 'y_multi','y_add','y_pink','y_blue','y_dyn','initial')