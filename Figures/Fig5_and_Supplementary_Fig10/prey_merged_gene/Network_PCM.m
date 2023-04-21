clc;
clear;
close all;

%type = 'full';
type = 'half';
%type = 'quarter';

filename = ['PCM_result_',type];
load(filename)

PPC = PPC.';
PPC_total = [];
for i = 1:4
    for j = 1:4
        if i == j 
            continue
        end
        PPC_total = [PPC_total; PPC(i,j)];
    end
end

figure(1)
scatter(PPC_total, [1:length(PPC_total)])

idx_k = kmeans(PPC_total,2,'MaxIter',10000,'Display','final','Replicates',100);

PPC_1 = PPC_total(find(idx_k == 1));
PPC_2 = PPC_total(find(idx_k == 2));

if mean(PPC_1) > mean(PPC_2)
    thres = (max(PPC_2) + min(PPC_1)) / 2;
else
    thres = (max(PPC_1) + min(PPC_2)) / 2;
end

thres
%thres = 0.5684

causal_network = zeros(4);
for i = 1:4
    for j = 1:4
        if PPC(i,j) > thres
            causal_network(i,j) = 1;
        end
    end
end
causal_network
