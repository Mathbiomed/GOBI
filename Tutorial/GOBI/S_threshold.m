function [S_tmp] = S_threshold(S, thres)

S_tmp = S;
for i = 1:length(S)
    if S(i) > thres
        S_tmp(i) = 1;
    else
        S_tmp(i) = 0;
    end
end

end

