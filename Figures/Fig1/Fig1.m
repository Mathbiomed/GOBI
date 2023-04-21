clc;
clear;
close all;

addpath('../GOBI') 

%% parameter
range = 1; % range of time-series
time_interval = 1/100; % 1/time_interval represents sampling rate

%% parameter for figure
c_nan = [0.9 0.9 0.9];
font_s = 14;
tmp = linspace(0, 1, 501)';
cmap_score = [[ones(500,1);1-tmp],[tmp(1:end-1);1-tmp],[tmp; ones(500,1)]];
cmap_value = [[ones(500,1);1-tmp],[tmp*0.3+0.7; (1-tmp(1:end-1))*0.6+0.4],[tmp(1:end-1);1-tmp]];

%% create source (Z = cos(2*pi*t), Y = sin(2*pi*t))
zt_pre = linspace(0,range, range/time_interval+1);
zt = linspace(0,range, range/time_interval+1);
yt = zt;
z = cos(2*pi*zt);
yy = sin(2*pi*yt);

%% Fig. 1a: time-series Z, X, X' with positive 1D regulation from Z to X
% create time-series
tspan = linspace(0,range, range/time_interval+1);
[t0, x0] = ode45(@(t,y) sample_a(t,y,zt,z), tspan, [0]);
ts = [z.', x0];

% normalize the time series Z, X
ts_norm = ts;
for i = 1:2
    ts_tmp = ts(:,i);
    ts_tmp = ts_tmp / max(abs(ts_tmp));
    ts_norm(:,i) = ts_tmp;
end

% create X'
X_prime = gradient(ts(:,2), time_interval);
X_prime = X_prime / max(abs(X_prime)); % normalize X'

% plot Fig. 1a
figure(1)
plot(t0, ts_norm(:,1), 'k','linewidth',2)
hold on
plot(t0, ts_norm(:,2), 'r','linewidth',2)
hold on
plot(t0, X_prime(:,1), 'b','linewidth',2)

xlim([0,range])
ylim([-1,1])
xticks([0,range])
yticks([-1,0,1])
xticklabels({'0', '1'})
set(gca, 'FontSize', font_s)
xlabel('Norm. time');
title('Fig. 1a')

%% Fig. 1b. regulation-detection function with positive 1D regulation from Z to X
% compute regulation-detection function
t_target = t0;
y_target = ts_norm;
[score_list, t_1, t_2] = RDS_dim1(y_target(:,1), y_target(:,2), t_target, time_interval);
score_target = reshape(score_list(:,:,1), [length(t0),length(t0)]); 
score_target(score_target == 0) = NaN;
score_target = score_target / max(max(abs(score_target))); %normalize regulation-detection function

%plot Fig. 1b
figure(2)
h = heatmap(flipud(score_target));
timelabel = string(t0);
timelabel(mod(t0, range) ~= 0) = '';
timelabel(end) = '1';
h.XDisplayLabels = timelabel;
h.YDisplayLabels = flipud(timelabel);
h.GridVisible = 'off';
h.FontName = 'Arial';
h.Colormap = cmap_score;
h.ColorLimits = [-1 1];
h.FontSize = font_s;
h.MissingDataColor = c_nan;
h.XLabel = 't';
h.YLabel = 't*';

s = struct(h);
s.XAxis.TickLabelRotation = 0;
s.YAxis.TickLabelRotation = 90;
cbh = s.Colorbar;
set(cbh, 'YTick', [-1,0,1])
title('Fig. 1b')
%% Fig. 1c. regulation-detection function with negative 1D regulation from Z to X
% create time-series
tspan = linspace(0,range, range/time_interval+1);
[t0, x0] = ode45(@(t,y) sample_c(t,y,zt,z), tspan, [0]);
ts = [z.', x0];
% normalize the time series Z, X
ts_norm = ts;
for i = 1:2
    ts_tmp = ts(:,i);
    ts_tmp = ts_tmp / max(abs(ts_tmp));
    ts_norm(:,i) = ts_tmp;
end
% compute regulation-detection function
t_target = t0;
y_target = ts_norm;
[score_list, t_1, t_2] = RDS_dim1(y_target(:,1), y_target(:,2), t_target, time_interval);
score_target = reshape(score_list(:,:,2), [length(t0),length(t0)]); 
score_target(score_target == 0) = NaN;
score_target = score_target / max(max(abs(score_target))); %normalize regulation-detection function

%plot Fig. 1c
figure(3)
h = heatmap(flipud(score_target));
timelabel = string(t0);
timelabel(mod(t0, range) ~= 0) = '';
timelabel(end) = '1';
h.XDisplayLabels = timelabel;
h.YDisplayLabels = flipud(timelabel);
h.GridVisible = 'off';
h.FontName = 'Arial';
h.Colormap = cmap_score;
h.ColorLimits = [-1 1];
h.FontSize = font_s;
h.MissingDataColor = c_nan;
h.XLabel = 't';
h.YLabel = 't*';
s = struct(h);
s.XAxis.TickLabelRotation = 0;
s.YAxis.TickLabelRotation = 90;
cbh = s.Colorbar;
set(cbh, 'YTick', [-1,0,1])
title('Fig. 1c')
%% Fig. 1d-f. regulation-detection function with 2D regulation from Z,Y to X. Specifically, both Z and Y positively regulate X
% create time-series
tspan = linspace(0,range, range/time_interval+1);
[t0, x0] = ode45(@(t,y) sample_d(t,y,zt,z,yt,yy), tspan, [0]);
ts = [z.',yy.', x0];
% normalize the time series Z, X
ts_norm = ts;
for i = 1:3
    ts_tmp = ts(:,i);
    ts_tmp = ts_tmp / max(abs(ts_tmp));
    ts_norm(:,i) = ts_tmp;
end
% compute regulation-detection function
t_target = t0;
y_target = ts_norm;
[score_list, t_1, t_2] = RDS_dim2(y_target(:,1), y_target(:,2), y_target(:,3), t_target, time_interval);
score_target = reshape(score_list(:,:,1), [length(t0),length(t0)]); 
score_target(score_target == 0) = NaN;
score_target = score_target / max(max(abs(score_target))); %normalize regulation-detection function

%plot Fig. 1e
figure(4)
h = heatmap(flipud(score_target));
timelabel = string(t0);
timelabel(mod(t0, range) ~= 0) = '';
timelabel(end) = '1';
h.XDisplayLabels = timelabel;
h.YDisplayLabels = flipud(timelabel);
h.GridVisible = 'off';
h.FontName = 'Arial';
h.Colormap = cmap_score;
h.ColorLimits = [-1 1];
h.FontSize = font_s;
h.MissingDataColor = c_nan;
h.XLabel = 't';
h.YLabel = 't*';
s = struct(h);
s.XAxis.TickLabelRotation = 0;
s.YAxis.TickLabelRotation = 90;
cbh = s.Colorbar;
set(cbh, 'YTick', [-1,0,1])
title('Fig. 1e')

score_target = reshape(score_list(:,:,2), [length(t0),length(t0)]); 
score_target(score_target == 0) = NaN;
score_target = score_target / max(max(abs(score_target))); %normalize regulation-detection function

%plot Fig. 1f
figure(5)
h = heatmap(flipud(score_target));
timelabel = string(t0);
timelabel(mod(t0, range) ~= 0) = '';
timelabel(end) = '1';
h.XDisplayLabels = timelabel;
h.YDisplayLabels = flipud(timelabel);
h.GridVisible = 'off';
h.FontName = 'Arial';
h.Colormap = cmap_score;
h.ColorLimits = [-1 1];
h.FontSize = font_s;
h.MissingDataColor = c_nan;
h.XLabel = 't';
h.YLabel = 't*';
s = struct(h);
s.XAxis.TickLabelRotation = 0;
s.YAxis.TickLabelRotation = 90;
cbh = s.Colorbar;
set(cbh, 'YTick', [-1,0,1])
title('Fig. 1f')
%% Fig. 1g-i. regulation-detection function with 2D regulation from Z,Y to X. Specifically, Z (Y) positively (negatively) regulate X, respectively.
% create time-series
tspan = linspace(0,range, range/time_interval+1);
[t0, x0] = ode45(@(t,y) sample_g(t,y,zt,z,yt,yy), tspan, [0]);
ts = [z.',yy.', x0];
% normalize the time series Z, X
ts_norm = ts;
for i = 1:3
    ts_tmp = ts(:,i);
    ts_tmp = ts_tmp / max(abs(ts_tmp));
    ts_norm(:,i) = ts_tmp;
end
% compute regulation-detection function
t_target = t0;
y_target = ts_norm;
[score_list, t_1, t_2] = RDS_dim2(y_target(:,1), y_target(:,2), y_target(:,3), t_target, time_interval);
score_target = reshape(score_list(:,:,1), [length(t0),length(t0)]); 
score_target(score_target == 0) = NaN;
score_target = score_target / max(max(abs(score_target))); %normalize regulation-detection function

%plot Fig. 1h
figure(6)
h = heatmap(flipud(score_target));
timelabel = string(t0);
timelabel(mod(t0, range) ~= 0) = '';
timelabel(end) = '1';
h.XDisplayLabels = timelabel;
h.YDisplayLabels = flipud(timelabel);
h.GridVisible = 'off';
h.FontName = 'Arial';
h.Colormap = cmap_score;
h.ColorLimits = [-1 1];
h.FontSize = font_s;
h.MissingDataColor = c_nan;
h.XLabel = 't';
h.YLabel = 't*';
s = struct(h);
s.XAxis.TickLabelRotation = 0;
s.YAxis.TickLabelRotation = 90;
cbh = s.Colorbar;
set(cbh, 'YTick', [-1,0,1])
title('Fig. 1h')

score_target = reshape(score_list(:,:,2), [length(t0),length(t0)]); 
score_target(score_target == 0) = NaN;
score_target = score_target / max(max(abs(score_target))); %normalize regulation-detection function

%plot Fig. 1i
figure(7)
h = heatmap(flipud(score_target));
timelabel = string(t0);
timelabel(mod(t0, range) ~= 0) = '';
timelabel(end) = '1';
h.XDisplayLabels = timelabel;
h.YDisplayLabels = flipud(timelabel);
h.GridVisible = 'off';
h.FontName = 'Arial';
h.Colormap = cmap_score;
h.ColorLimits = [-1 1];
h.FontSize = font_s;
h.MissingDataColor = c_nan;
h.XLabel = 't';
h.YLabel = 't*';
s = struct(h);
s.XAxis.TickLabelRotation = 0;
s.YAxis.TickLabelRotation = 90;
cbh = s.Colorbar;
set(cbh, 'YTick', [-1,0,1])
title('Fig. 1i')
%% Fig. 1j-l. regulation-detection function with positive 1D regulation from Z to X and Y does not affect X
% create time-series
tspan = linspace(0,range, range/time_interval+1);
[t0, x0] = ode45(@(t,y) sample_j(t,y,zt,z,yt,yy), tspan, [0]);
ts = [z.',yy.', x0];
% normalize the time series Z, X
ts_norm = ts;
for i = 1:3
    ts_tmp = ts(:,i);
    ts_tmp = ts_tmp / max(abs(ts_tmp));
    ts_norm(:,i) = ts_tmp;
end
% compute regulation-detection function
t_target = t0;
y_target = ts_norm;
[score_list, t_1, t_2] = RDS_dim2(y_target(:,1), y_target(:,2), y_target(:,3), t_target, time_interval);
score_target = reshape(score_list(:,:,1), [length(t0),length(t0)]); 
score_target(score_target == 0) = NaN;
score_target = score_target / max(max(abs(score_target))); %normalize regulation-detection function

%plot Fig. 1k
figure(8)
h = heatmap(flipud(score_target));
timelabel = string(t0);
timelabel(mod(t0, range) ~= 0) = '';
timelabel(end) = '1';
h.XDisplayLabels = timelabel;
h.YDisplayLabels = flipud(timelabel);
h.GridVisible = 'off';
h.FontName = 'Arial';
h.Colormap = cmap_score;
h.ColorLimits = [-1 1];
h.FontSize = font_s;
h.MissingDataColor = c_nan;
h.XLabel = 't';
h.YLabel = 't*';
s = struct(h);
s.XAxis.TickLabelRotation = 0;
s.YAxis.TickLabelRotation = 90;
cbh = s.Colorbar;
set(cbh, 'YTick', [-1,0,1])
title('Fig. 1k')

score_target = reshape(score_list(:,:,2), [length(t0),length(t0)]); 
score_target(score_target == 0) = NaN;
score_target = score_target / max(max(abs(score_target))); %normalize regulation-detection function

%plot Fig. 1l
figure(9)
h = heatmap(flipud(score_target));
timelabel = string(t0);
timelabel(mod(t0, range) ~= 0) = '';
timelabel(end) = '1';
h.XDisplayLabels = timelabel;
h.YDisplayLabels = flipud(timelabel);
h.GridVisible = 'off';
h.FontName = 'Arial';
h.Colormap = cmap_score;
h.ColorLimits = [-1 1];
h.FontSize = font_s;
h.MissingDataColor = c_nan;
h.XLabel = 't';
h.YLabel = 't*';
s = struct(h);
s.XAxis.TickLabelRotation = 0;
s.YAxis.TickLabelRotation = 90;
cbh = s.Colorbar;
set(cbh, 'YTick', [-1,0,1])
title('Fig. 1l')
