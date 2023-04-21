clc;
clear;
close all;

%filename = ['PCM_result'];
filename = ['PCM_result_2y'];

load(filename)

PPC = PPC.';

PPC_total = PPC(2:end,1);

idx_k = kmeans(PPC_total,2,'MaxIter',10);

PPC_1 = PPC_total(find(idx_k == 1));
PPC_2 = PPC_total(find(idx_k == 2));

if mean(PPC_1) > mean(PPC_2)
    thres = (max(PPC_2) + min(PPC_1)) / 2;
else
    thres = (max(PPC_1) + min(PPC_2)) / 2;
end

thres

causal_network = zeros([4,1]);
for j = 1:4
    if PPC(j+1,1) > thres
        causal_network(j,1) = 1;
    end
end
causal_network
