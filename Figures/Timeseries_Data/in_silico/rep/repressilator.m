function dydt = repressilator(t,y)

alpha = 7;
alpha_0 = 0;
beta = 3;
n = 10;

%    [1     , 2     , 3     , 4     , 5   , 6   ]
% Y= [m_lacl, p_lacl, m_tetR, p_tetR, m_cl, p_cl] 

dydt=zeros(6,1);
dydt(1) = -y(1) + alpha/(1 + (y(6).^n)) + alpha_0;
dydt(2) = beta*(y(1) - y(2));
dydt(3) = -y(3) + alpha/(1 + (y(2).^n)) + alpha_0;
dydt(4) = beta*(y(3) - y(4));
dydt(5) = -y(5) + alpha/(1 + (y(4).^n)) + alpha_0;
dydt(6) = beta*(y(5) - y(6));

end
