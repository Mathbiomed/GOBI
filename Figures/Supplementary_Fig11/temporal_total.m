function dydt = temporal_total(t,y,t_pre, alpha_pre, beta_pre, gamma_pre, X_pre, Z_pre)

alpha = interp1(t_pre,alpha_pre,t);
beta = interp1(t_pre,beta_pre,t);
gamma = interp1(t_pre,gamma_pre,t);
X = interp1(t_pre,X_pre,t);
Z = interp1(t_pre,Z_pre,t);
n = 2;
Kd_1 = 0.1;
Kd_2 = 0.5;

% [Y]
dydt = zeros(1,1);
dydt(1) = alpha .* (X.^n) / (Kd_1 + X.^n) + 0.4*beta ./ (Kd_2 + X.^n) + 2*gamma .* (Z.^3) - y(1);

end