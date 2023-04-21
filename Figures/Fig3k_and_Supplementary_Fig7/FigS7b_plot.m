clc;
clear;
close all;

%% load data
load('data_example')

%% fit
y_tmp = fit(t_ori,y_multi,'fourier4');
w = y_tmp.w;
fouriers = [
    ones(1,length(t_ori));
    cos(w*t_ori.');
    sin(w*t_ori.');
    cos(2*w*t_ori.');
    sin(2*w*t_ori.');
    cos(3*w*t_ori.');
    sin(3*w*t_ori.');
    cos(4*w*t_ori.');
    sin(4*w*t_ori.')];
coeffs = coeffvalues(y_tmp);
y_multi_fit = coeffs(1:end-1) * fouriers;


y_tmp = fit(t_ori,y_add,'fourier4');
w = y_tmp.w;
fouriers = [
    ones(1,length(t_ori));
    cos(w*t_ori.');
    sin(w*t_ori.');
    cos(2*w*t_ori.');
    sin(2*w*t_ori.');
    cos(3*w*t_ori.');
    sin(3*w*t_ori.');
    cos(4*w*t_ori.');
    sin(4*w*t_ori.')];
coeffs = coeffvalues(y_tmp);
y_add_fit = coeffs(1:end-1) * fouriers;

y_blue = double(y_blue);
y_tmp = fit(t_ori,y_blue,'fourier4');
w = y_tmp.w;
fouriers = [
    ones(1,length(t_ori));
    cos(w*t_ori.');
    sin(w*t_ori.');
    cos(2*w*t_ori.');
    sin(2*w*t_ori.');
    cos(3*w*t_ori.');
    sin(3*w*t_ori.');
    cos(4*w*t_ori.');
    sin(4*w*t_ori.')];
coeffs = coeffvalues(y_tmp);
y_blue_fit = coeffs(1:end-1) * fouriers;

y_pink = double(y_pink);
y_tmp = fit(t_ori,y_pink,'fourier4');
w = y_tmp.w;
fouriers = [
    ones(1,length(t_ori));
    cos(w*t_ori.');
    sin(w*t_ori.');
    cos(2*w*t_ori.');
    sin(2*w*t_ori.');
    cos(3*w*t_ori.');
    sin(3*w*t_ori.');
    cos(4*w*t_ori.');
    sin(4*w*t_ori.')];
coeffs = coeffvalues(y_tmp);
y_pink_fit = coeffs(1:end-1) * fouriers;

y_tmp = fit(t_ori,y_dyn,'fourier4');
w = y_tmp.w;
fouriers = [
    ones(1,length(t_ori));
    cos(w*t_ori.');
    sin(w*t_ori.');
    cos(2*w*t_ori.');
    sin(2*w*t_ori.');
    cos(3*w*t_ori.');
    sin(3*w*t_ori.');
    cos(4*w*t_ori.');
    sin(4*w*t_ori.')];
coeffs = coeffvalues(y_tmp);
y_dyn_fit = coeffs(1:end-1) * fouriers;

%% plot
figure(1)
plot(t_ori, y_ori, 'k')
hold on
plot(t_ori, y_multi, 'r')
hold on
plot(t_ori, y_multi_fit, 'b')

xlim([0,10])
xticks([0,10])
xticklabels([])
ylim([-0.5,1.5])
yticks([-0.5:0.5:1.5])
yticklabels([])

figure(2)
plot(t_ori, y_ori, 'k')
hold on
plot(t_ori, y_add, 'r')
hold on
plot(t_ori, y_add_fit, 'b')

xlim([0,10])
xticks([0,10])
xticklabels([])
ylim([-0.5,1.5])
yticks([-0.5:0.5:1.5])
yticklabels([])

figure(3)
plot(t_ori, y_ori, 'k')
hold on
plot(t_ori, y_pink, 'r')
hold on
plot(t_ori, y_pink_fit, 'b')

xlim([0,10])
xticks([0,10])
xticklabels([])
ylim([-0.5,1.5])
yticks([-0.5:0.5:1.5])
yticklabels([])

figure(4)
plot(t_ori, y_ori, 'k')
hold on
plot(t_ori, y_blue, 'r')
hold on
plot(t_ori, y_blue_fit, 'b')

xlim([0,10])
xticks([0,10])
xticklabels([])
ylim([-0.5,1.5])
yticks([-0.5:0.5:1.5])
yticklabels([])

figure(5)
plot(t_ori, y_ori, 'k')
hold on
plot(t_ori, y_dyn, 'r')
hold on
plot(t_ori, y_dyn_fit, 'b')

xlim([0,10])
xticks([0,10])
xticklabels([])
ylim([-0.5,1.5])
yticks([-0.5:0.5:1.5])
yticklabels([])