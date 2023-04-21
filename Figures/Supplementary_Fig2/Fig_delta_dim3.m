clc;
clear;
close all;

%% parameter
dimension = 3;
thres_L = 0.01;

%% load data
filename = 'sample_result_dim3';
load(filename)

num_pair = length(component_list(:,1));
num_type = 2^dimension;

%% tick
delta_index = [
    [1,3,7];
    [1,5,7];
    [1,7,6];
    
    [2,1,5];
    [2,1,3];
    [2,1,2];
    
    [2,2,6];
    [2,2,4];
    [2,1,2];
    
    [3,1,5];
    [3,5,7];
    [3,5,6];
    
    [3,2,6];
    [3,6,8];
    [3,5,6];
    
    [3,3,7];
    [3,5,7];
    [3,7,8];
    
    [3,4,8];
    [3,6,8];
    [3,7,8]];

tick_total_x = {
    '\Delta_{B^{-}C^{+}}^{D} (A)';
    '\Delta_{A^{-}C^{+}}^{D} (B)';
    '\Delta_{A^{-}B^{-}}^{D} (C)';
    
    '\Delta_{B^{+}D^{+}}^{C} (A)';
    '\Delta_{A^{+}D^{+}}^{C} (B)';
    '\Delta_{A^{+}B^{+}}^{C} (D)';
    
    '\Delta_{B^{+}D^{-}}^{C} (A)';
    '\Delta_{A^{+}D^{-}}^{C} (B)';
    '\Delta_{A^{+}B^{+}}^{C} (D)';
    
    '\Delta_{C^{+}D^{+}}^{B} (A)';
    '\Delta_{A^{-}D^{+}}^{B} (C)';
    '\Delta_{A^{-}C^{+}}^{B} (D)';
    
    '\Delta_{C^{+}D^{-}}^{B} (A)';
    '\Delta_{A^{-}D^{-}}^{B} (C)';
    '\Delta_{A^{-}C^{+}}^{B} (D)';
    
    '\Delta_{C^{-}D^{+}}^{B} (A)';
    '\Delta_{A^{-}D^{+}}^{B} (C)';
    '\Delta_{A^{-}C^{-}}^{B} (D)';
   
    '\Delta_{C^{-}D^{-}}^{B} (A)';
    '\Delta_{A^{-}D^{-}}^{B} (C)';
    '\Delta_{A^{-}C^{-}}^{B} (D)'};

%% process
x = [1:num_pair * num_type];
delta_plot = [];
for i = 1:length(delta_index(:,1))
    S_tmp_1 = reshape(S_total(delta_index(i,1),delta_index(i,2),:),[num_data,1]);
    L_tmp_1 = reshape(L_total(delta_index(i,1),delta_index(i,2),:),[num_data,1]);

    L_processed_1 = L_threshold(L_tmp_1, thres_L);
    S_processed_1 = S_tmp_1 .* L_processed_1;
    S_processed_1(S_processed_1 == 0) = NaN;
    
    S_tmp_2 = reshape(S_total(delta_index(i,1),delta_index(i,3),:),[num_data,1]);
    L_tmp_2 = reshape(L_total(delta_index(i,1),delta_index(i,3),:),[num_data,1]);

    L_processed_2 = L_threshold(L_tmp_2, thres_L);
    S_processed_2 = S_tmp_2 .* L_processed_2;
    S_processed_2(S_processed_2 == 0) = NaN;
    
    delta_plot = [delta_plot, S_processed_1 - S_processed_2];
end

%% plot
figure(1)
b = boxchart(delta_plot);

b.BoxWidth = 0.5;
b.BoxFaceColor = [0,0,0];
b.LineWidth = 2;
b.BoxFaceAlpha = 0.1;
b.MarkerStyle = 'none';
yticks([-2,0,2])
xticklabels(tick_total_x)
xlabel('\sigma')
ylabel('\Delta')
set(gca, 'FontSize',14)





