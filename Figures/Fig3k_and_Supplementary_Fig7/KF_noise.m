function dydt = KF_noise(t,y, noise_level, a,b,c)
% parameter
k_1 = 1;   k_2 = 0.16; k_3 = 0.29; k_4 = 0.3;
A_t = 0.6; K_d = 1e-5;

f_P = (A_t - K_d - y(3)+sqrt((A_t - K_d - y(3))^2 + 4*A_t*K_d))/(2*A_t);
% Y= [M, PC, P]
dydt=zeros(3,1);
dydt(1) = k_1*f_P - k_2*y(1) + noise_level * randn() * a;
dydt(2) = k_1*y(1) - k_3*y(2) + noise_level * randn() * b;
dydt(3) = k_1*y(2) - k_4*y(3) + noise_level * randn() * c;