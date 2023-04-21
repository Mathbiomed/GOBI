clc;
clear;
close all;

addpath('../GOBI') 

%% load data
load('data_genetic_oscillator')

%% parameters
num_data = 8;
num_component = 2;

%% spline fit
t_fit = linspace(0,1,101).';
time_interval = 1/100;
data_spline = zeros(length(t_fit),num_component,num_data);
for i = 1:num_data
    y_tmp = reshape(data_total(:,:,i),[num_component, length(time)]);
    y_tmp = y_tmp .';
    
    y_int = zeros(length(t_fit),num_component);
    for j = 1:2
        y_int(:,j) = spline(time, y_tmp(:,j), t_fit);
    end 
    data_spline(:,:,i) = y_int;
end


%% compute RDF sigma^28 -> exp(TetR)*sigma^28
thres_noise = 0;
I_total_XY = zeros(num_data,length(t_fit),length(t_fit), 2);
Y_trans = zeros(num_data, length(t_fit));

for i = 1:num_data
    % import time series
    y_tmp = reshape(data_spline(:,:,i),[length(t_fit), num_component]);

    % import single time series
    cause_1 = y_tmp(:,1);
    cause_2 = y_tmp(:,2);
    target = exp(y_tmp(:,1)) .* y_tmp(:,2);
    %target = -1./y_tmp(:,1) + 1./ y_tmp(:,2);
    
    % compute RDf for the regulation from X to Y
    [score_list, t_1, t_2] = RDS_dim2(cause_1, cause_2, target, t_fit, time_interval);
    type_list = [3,4];
    for k = 1:2
        score = reshape(score_list(:,:,type_list(k)),[length(t_fit),length(t_fit)]);
        I_total_XY(i,:,:,k) = score;            
    end

    % save the value of Y
    Y_trans(i,:) = cause_2;
end

%% transform (t,t*) into X(t),X(t*)
trans_total_XY = cell(num_data,2);

for i = 1:num_data
    score_XY_pos = reshape(I_total_XY(i,:,:,1), [length(t_fit), length(t_fit)]);
    score_XY_neg = reshape(I_total_XY(i,:,:,2), [length(t_fit), length(t_fit)]);

    Y_var = Y_trans(i,:);
    
    trans_XY_pos = [];
    trans_XY_neg = [];

    for m = 1:length(t_fit)
        for n = 1:length(t_fit)
            if abs(score_XY_pos(m,n)) > thres_noise
                trans_XY_pos = [trans_XY_pos ; [Y_var(m), Y_var(n), score_XY_pos(m,n)]];
            end
            if abs(score_XY_neg(m,n)) > thres_noise
                trans_XY_neg = [trans_XY_neg ; [Y_var(m), Y_var(n), score_XY_neg(m,n)]];
            end
        end
    end        
    trans_total_XY{i,1} = trans_XY_pos;
    trans_total_XY{i,2} = trans_XY_neg;
end

%% plot

figure(2)
for i = 1:4
    trans_XY_pos = cell2mat(trans_total_XY(i,1));
    trans_XY_neg = cell2mat(trans_total_XY(i,2));

    loca_pos = find(trans_XY_pos(:,3) > 0);
    loca_neg = find(trans_XY_pos(:,3) < 0);

    alpha = 0.02;

    %subplot(4,2,i)
    scatter(trans_XY_pos(loca_pos,1),trans_XY_pos(loca_pos,2),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha)
    hold on
    scatter(trans_XY_pos(loca_neg,1),trans_XY_pos(loca_neg,2),'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha)
    hold on

    loca_pos = find(trans_XY_neg(:,3) > 0);
    loca_neg = find(trans_XY_neg(:,3) < 0);

    scatter(trans_XY_neg(loca_pos,1),trans_XY_neg(loca_pos,2),'MarkerFaceColor','b','MarkerEdgeColor','b','MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha)
    hold on
    scatter(trans_XY_neg(loca_neg,1),trans_XY_neg(loca_neg,2),'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha)
    hold on
    xlim([0,180])
    ylim([0,180])
    xticks([0,180])
    yticks([0,180])
    set(gca,'XColor', 'none','YColor','none')
    axis square
end

%% probability distribution
disp('compute pdf')
num_bin = 4;
plus_count = zeros(1,2*num_bin.^2);
minus_count = zeros(1,2*num_bin.^2);
for i = 1:4
    trans_XY_pos = cell2mat(trans_total_XY(i,1));
    for k = 1:length(trans_XY_pos(:,1))
        X_var1 = trans_XY_pos(k,1);
        X_var2 = trans_XY_pos(k,2);
        score_var = trans_XY_pos(k,3);

        m = fix(X_var1/(180/num_bin));
        n = fix(X_var2/(180/num_bin));

        if m == num_bin
            m = m-1;
        end
        if n == num_bin
            n = n-1;
        end
        idx = 2*num_bin*n + 2*m;

        if X_var1-m*180/num_bin > X_var2-n*180/num_bin
            idx = idx+2;
        else
            idx = idx+1;
        end

        if score_var > 0
            plus_count(idx) = plus_count(idx)+1;
        else
            minus_count(idx) = minus_count(idx)+1;
        end
    end

    trans_XY_neg = cell2mat(trans_total_XY(i,2));
    for k = 1:length(trans_XY_neg(:,1))
        X_var1 = trans_XY_neg(k,1);
        X_var2 = trans_XY_neg(k,2);
        score_var = trans_XY_neg(k,3);

        m = fix(X_var1/(180/num_bin));
        n = fix(X_var2/(180/num_bin));

        if m == num_bin
            m = m-1;
        end
        if n == num_bin
            n = n-1;
        end
        idx = 2*num_bin*n + 2*m;

        if X_var1-m*180/num_bin > X_var2-n*180/num_bin
            idx = idx+2;
        else
            idx = idx+1;
        end

        if score_var > 0
            plus_count(idx) = plus_count(idx)+1;
        else
            minus_count(idx) = minus_count(idx)+1;
        end
    end
    
end


plus_count = plus_count / sum(plus_count) + eps;
minus_count = minus_count/ sum(minus_count) + eps;

%% variational information
disp('compute KL divergence')

X = [1:2*num_bin^2].';
KL = kldiv(X,plus_count.', minus_count.','js')

%% mesh
figure(3)
x_axis = [-2/3, -1/3, 1/3, 2/3, 4/3, 5/3, 7/3, 8/3, 10/3, 11/3, 13/3, 14/3];
x_axis = [x_axis, x_axis, x_axis, x_axis, x_axis, x_axis];
y_axis = [2/3, 1/3, 2/3, 1/3, 2/3, 1/3, 2/3, 1/3, 2/3, 1/3, 2/3, 1/3];
y_axis = [y_axis-1, y_axis, y_axis+1, y_axis+2, y_axis+3, y_axis+4];

%[X,Y] = meshgrid(x_axis,y_axis);
%Z = [plus_count(1:8);plus_count(9:16);plus_count(17:24);plus_count(25:32)];
%contour(X,Y,Z)
plus_count = [zeros(1,12), [0,0,plus_count(1:8),0,0], [0,0,plus_count(9:16),0,0], [0,0,plus_count(17:24),0,0], [0,0,plus_count(25:32),0,0], zeros(1,12)];
minus_count = [zeros(1,12), [0,0,minus_count(1:8),0,0], [0,0,minus_count(9:16),0,0], [0,0,minus_count(17:24),0,0], [0,0,minus_count(25:32),0,0], zeros(1,12)];

[xq,yq] = meshgrid(-1:1/20:5, -1:1/20:5);
vq = griddata(x_axis,y_axis,minus_count,xq,yq);  %(x,y,v) being your original data for plotting points
T = delaunay(xq,yq);
mesh(xq,yq,vq, 'EdgeColor','None','FaceColor', 'flat')
%trimesh(T, xq,yq,vq, 'EdgeColor','None','FaceColor', 'flat')
xlim([0,4])
ylim([0,4])
xticks([0,4])
yticks([0,4])
xlabel('x')
ylabel('y')
%set(gca,'XColor', 'none','YColor','none')
axis square
%hold on
%plot3(x_axis,y_axis,plus_count,'o')

