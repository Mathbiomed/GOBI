function [L_tmp] = L_threshold(L, thres)
L_tmp = L;
for i = 1:length(L)
    if L(i) < thres
        L_tmp(i) = 0;
    else
        L_tmp(i) = 1;
        %L_tmp(i,j) = L_tmp(i,j);
    end
end
end

