function [score_list, t_1, t_2] = RDS_ps_dim4(X, Y, Z, W, K, t, time_interval)

f = gradient(K, time_interval);
score_list = zeros(length(t),length(t),16);
t_1 = zeros(length(t),length(t));
t_2 = zeros(length(t),length(t));

for t_ori = 1:length(t)
    for t_prime = 1:length(t)
        t_1(t_ori, t_prime) = t(t_ori);
        t_2(t_ori, t_prime) = t(t_prime);
        
        X_d = X(t_ori) - X(t_prime);
        Y_d = Y(t_ori) - Y(t_prime);
        Z_d = Z(t_ori) - Z(t_prime);
        W_d = X(t_ori) - W(t_prime);
        K_d = K(t_ori) - K(t_prime);
        f_d = f(t_ori) - f(t_prime);
        score_tmp = X_d * Y_d * Z_d * W_d * K_d * f_d;
        
        if     X_d >= 0 && Y_d >= 0 && Z_d >= 0 && W_d >= 0 && K_d >=  0         %type 1
            score_list(t_ori, t_prime, 1) = score_tmp;
        elseif X_d >= 0 && Y_d >= 0 && Z_d >= 0 && W_d < 0 && K_d >=  0     %type 2
            score_list(t_ori, t_prime, 2) = -score_tmp;
        elseif X_d >= 0 && Y_d >= 0 && Z_d < 0 && W_d >= 0 && K_d >=  0     %type 3
            score_list(t_ori, t_prime, 3) = -score_tmp;
        elseif X_d >= 0 && Y_d >= 0 && Z_d < 0 && W_d < 0 && K_d >=  0      %type 4
            score_list(t_ori, t_prime, 4) = score_tmp;
        elseif X_d >= 0 && Y_d < 0 && Z_d >= 0 && W_d >= 0 && K_d >=  0         %type 1
            score_list(t_ori, t_prime, 5) = -score_tmp;
        elseif X_d >= 0 && Y_d < 0 && Z_d >= 0 && W_d < 0 && K_d >=  0     %type 2
            score_list(t_ori, t_prime, 6) = score_tmp;
        elseif X_d >= 0 && Y_d < 0 && Z_d < 0 && W_d >= 0 && K_d >=  0     %type 3
            score_list(t_ori, t_prime, 7) = score_tmp;
        elseif X_d >=  0 && Y_d < 0 && Z_d < 0 && W_d < 0 && K_d >=  0      %type 4
            score_list(t_ori, t_prime, 8) = -score_tmp;
        elseif X_d < 0 && Y_d >= 0 && Z_d >= 0 && W_d >= 0 && K_d >=  0         %type 1
            score_list(t_ori, t_prime, 9) = -score_tmp;
        elseif X_d < 0 && Y_d >= 0 && Z_d >= 0 && W_d < 0 && K_d >=  0     %type 2
            score_list(t_ori, t_prime, 10) = score_tmp;
        elseif X_d <  0 && Y_d >= 0 && Z_d < 0 && W_d >= 0 && K_d >=  0     %type 3
            score_list(t_ori, t_prime, 11) = score_tmp;
        elseif X_d < 0 && Y_d >= 0 && Z_d < 0 && W_d < 0 && K_d >=  0      %type 4
            score_list(t_ori, t_prime, 12) = -score_tmp;
        elseif X_d < 0 && Y_d < 0 && Z_d >= 0 && W_d >= 0 && K_d >=  0         %type 1
            score_list(t_ori, t_prime, 13) = score_tmp;
        elseif X_d < 0 && Y_d < 0 && Z_d >= 0 && W_d < 0 && K_d >= 0     %type 2
            score_list(t_ori, t_prime, 14) = -score_tmp;
        elseif X_d <  0 && Y_d < 0 && Z_d < 0 && W_d >= 0 && K_d >=  0     %type 3
            score_list(t_ori, t_prime, 15) = -score_tmp;
        elseif X_d <  0 && Y_d < 0 && Z_d < 0 && W_d < 0 && K_d >= 0      %type 4
            score_list(t_ori, t_prime, 16) = score_tmp;
        end
        
    end
end
end

