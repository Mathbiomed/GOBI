function dydt = Frzilator(t,y)
% parameter
phi = 0.08;
k_c = 4; k_e = 4; 
d_f = 1; d_c = 2; d_e = 2;

% Y= [f, c, e]
dydt=zeros(3,1);
dydt(1) = phi*((1-y(1))/(0.01 + (1-y(1)))) - d_f*(y(1)/(0.005+y(1)))*y(3);
dydt(2) = k_c*((1-y(2))/(0.005 + (1-y(2))))*y(1) - d_c*(y(2)/(0.005+y(2)));
dydt(3) = k_e*((1-y(3))/(0.005 + (1-y(3))))*y(2) - d_e*(y(3)/(0.005+y(3)));
end