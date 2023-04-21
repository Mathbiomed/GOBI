%% function used in Simulated Annealing 
%% adjust the parameter so that they are within the range

function [parav1, parav2, parav3] = adjust_range(parav1, parav2, parav3)

% parav1: a,b,c
% parav2: K, K_2
% parav3: n, m


parav1(parav1>5)=5;
parav1(parav1<0)=0;
parav2(parav2>20)=20;
parav2(parav2<0)=0;
parav3(parav3>5)=5;
parav3(parav3<1)=1;

end