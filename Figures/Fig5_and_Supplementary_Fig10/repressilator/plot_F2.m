clc;
clear;
close all;

network_GC = {
    [[0,1,0];
     [1,0,1];
     [0,0,0]],
    
    [[0,1,1];
     [0,0,1];
     [0,1,0]],
    
    [[0,1,0];
     [0,0,1];
     [0,0,0]]};

network_CCM = {
    [[0,1,1];
     [0,0,1];
     [1,1,0]],
    
    [[0,1,1];
     [1,0,1];
     [1,1,0]],
    
    [[0,1,1];
     [0,0,1];
     [1,1,0]]};
 
 network_PCM = {
    [[0,1,0];
     [1,0,0];
     [0,1,0]],
    
    [[0,1,0];
     [1,0,1];
     [0,1,0]],
    
    [[0,1,0];
     [1,0,1];
     [0,1,0]]};
 
 network_GOBI = {
    [[0,0,-1];
     [-1,0,0];
     [0,-1,0]],
    
    [[0,0,-1];
     [-1,0,0];
     [0,-1,0]],
    
    [[0,0,0];
     [-1,0,0];
     [0,-1,0]]};
 
 true_network = [
     [0,0,-1];
     [-1,0,0];
     [0,-1,0]];
 
 F2_list = zeros(3,4);
 for i = 1:3
     GC_tmp = cell2mat(network_GC(i));
     CCM_tmp = cell2mat(network_CCM(i));
     PCM_tmp = cell2mat(network_PCM(i));
     GOBI_tmp = cell2mat(network_GOBI(i));
     
     F2_list(i,1) = cal_F2(GC_tmp, abs(true_network));
     F2_list(i,2) = cal_F2(CCM_tmp, abs(true_network));
     F2_list(i,3) = cal_F2(PCM_tmp, abs(true_network));
     F2_list(i,4) = cal_F2(GOBI_tmp, true_network);
 end
 
  sampling_list = [1,1/2,1/4];
 plot(sampling_list, F2_list(:,1), '-or','MarkerFaceColor', 'r', 'Markersize',10')
 hold on
 plot(sampling_list, F2_list(:,2), '-og','MarkerFaceColor', 'g', 'Markersize',10')
 hold on
 plot(sampling_list, F2_list(:,3), '-ob','MarkerFaceColor', 'b', 'Markersize',10')
 hold on
 plot(sampling_list, F2_list(:,4), '-om','MarkerFaceColor', 'm', 'Markersize',10')
 
 xlim([1/4,1])
 xticks([1/4,1/2,1])
 ylim([0,1])
 yticks([0,1])
 xlabel('fraction of data')
 ylabel('F_2 score')
 set(gca,'FontSize',12)
 
 legend({'GC', 'CCM', 'PCM', 'GOBI'}, 'Location','southeast')
 