function [score_list, t_1, t_2] = RDS_ns_dim1(X, Y, t, time_interval)

f = gradient(Y, time_interval);
score_list = zeros(length(t),length(t),2);
t_1 = zeros(length(t),length(t));
t_2 = zeros(length(t),length(t));

for t_ori = 1:length(t)
    for t_prime = 1:length(t)
        t_1(t_ori, t_prime) = t(t_ori);
        t_2(t_ori, t_prime) = t(t_prime);
        
        X_d = X(t_ori) - X(t_prime);
        Y_d = Y(t_ori) - Y(t_prime);
        f_d = f(t_ori) - f(t_prime);
        score_tmp = X_d * Y_d * f_d;
        
        if     X_d >= 0 && Y_d < 0         %type 1
            score_list(t_ori, t_prime,1) = -score_tmp;
        elseif X_d < 0 && Y_d <  0         %type 2
            score_list(t_ori, t_prime,2) = score_tmp;
        end
    end
end

end

