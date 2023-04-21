clc;
clear;
close all;
addpath('../GOBI')

%% parameter
period = 2*pi;

time_interval = 1/60;
num_component = 3;
thres_noise = 0;

%% create source
time_interval_pre = 2;
zt = linspace(0,period, period/time_interval+1);
z = cos(zt);
yt = zt;
yy = sin(zt);
%% create time series
y_0 =  zeros([num_component-2,1]);
tspan = linspace(0,period, period/time_interval+1);
[t0, y0] = ode15s(@(t,y) sample_2dim_pn(t,y,zt,z,yt,yy), tspan, y_0);
y0 = [z.', yy.', y0];

%% normalize time series
y_norm = y0;
for i = 1:num_component
    y_tmp = y0(:,i);
    y_tmp = y_tmp / max(abs(y_tmp));
    y_norm(:,i) = y_tmp;
end

%% create 2 dim score 1
st1_dim2 = 1;
st2_dim2 = 2;
ed_dim2 = 3;
t_target = t0;
y_target = y_norm;

[score_list, t_1, t_2] = RDS_dim2(y_target(:,st1_dim2), y_target(:,st2_dim2), y_target(:,ed_dim2), t_target, time_interval);

score_dim2_1 = reshape(score_list(:,:,1),[length(t0), length(t0)]);
score_dim2_2 = reshape(score_list(:,:,2),[length(t0), length(t0)]);
score_dim2_3 = reshape(score_list(:,:,3),[length(t0), length(t0)]);
score_dim2_4 = reshape(score_list(:,:,4),[length(t0), length(t0)]);

score_dim2_1(find(score_dim2_1 == 0)) = nan;
score_dim2_2(find(score_dim2_2 == 0)) = nan;
score_dim2_3(find(score_dim2_3 == 0)) = nan;
score_dim2_4(find(score_dim2_4 == 0)) = nan;

score_dim2_1 = score_dim2_1 / max(max(abs(score_dim2_1)));
score_dim2_2 = score_dim2_2 / max(max(abs(score_dim2_2)));
score_dim2_3 = score_dim2_3 / max(max(abs(score_dim2_3)));
score_dim2_4 = score_dim2_4 / max(max(abs(score_dim2_4)));


%% figure
c_nan = [0.9 0.9 0.9];
timelabel = string(t0);
timelabel(mod(t0, period) ~= 0) = '';
timelabel(end) = '1';
font_s = 14;
tmp = linspace(0, 1, 501)';
cmap_score = [[ones(500,1);1-tmp],[tmp(1:end-1);1-tmp],[tmp; ones(500,1)]];
%cmap_value = [[ones(500,1);1-tmp],ones(1001,1),[tmp; ones(500,1)]];
cmap_value = [[ones(500,1);1-tmp],[tmp*0.3+0.7; (1-tmp(1:end-1))*0.6+0.4],[tmp(1:end-1);1-tmp]];

figure(1) % heatmap of 2dim score pp 
h = heatmap(flipud(score_dim2_1));

timelabel = string(t0);
timelabel(mod(t0, period) ~= 0) = '';
timelabel(end) = '{\tau}';
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

figure(2) % heatmap of 2dim score pn
h = heatmap(flipud(score_dim2_2));

timelabel = string(t0);
timelabel(mod(t0, period) ~= 0) = '';
timelabel(end) = '{\tau}';
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

figure(3) % heatmap of 2dim score np
h = heatmap(flipud(score_dim2_3));

timelabel = string(t0);
timelabel(mod(t0, period) ~= 0) = '';
timelabel(end) = '{\tau}';
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

figure(4) % heatmap of 2dim score nn
h = heatmap(flipud(score_dim2_4));

timelabel = string(t0);
timelabel(mod(t0, period) ~= 0) = '';
timelabel(end) = '{\tau}';
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







