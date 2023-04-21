function dydt = cAMP_noise(t,y, noise_level,a,b,c,d,e,f,g)

% Computational modelling suggests dynamic interactions between Ca 2+, IP3
% and G protein-coupled modules are key to robust Dictyostelium aggregation

k_1  = 2.0; k_2  = 0.9; k_3  = 2.5; k_4  = 1.5; k_5  = 0.6;
k_6  = 0.8; k_7  = 1.0; k_8  = 1.3; k_9  = 0.3; k_10 = 0.8;
k_11 = 0.7; k_12 = 4.9; k_13 = 23;  k_14 = 4.5;


%    [1  , 2  , 3   , 4   , 5    , 6    , 7   ]
% Y= [ACA, PKA, ERK2, RegA, cAMPi, cAMPe, CAR1] 

dydt=zeros(7,1);
dydt(1) = k_1 * y(7) - k_2 * y(1) * y(2) + noise_level * randn() * a;
dydt(2) = k_3 * y(5) - k_4 * y(2) + noise_level * randn() * b;
dydt(3) = k_5 * y(7) - k_6 * y(2) * y(3) + noise_level * randn() * c;
dydt(4) = k_7 - k_8 * y(3) * y(4) + noise_level * randn() * d;
dydt(5) = k_9 * y(1) - k_10 * y(4) * y(5) + noise_level * randn() * e;
dydt(6) = k_11 * y(1) - k_12 * y(6) + noise_level * randn() * f;
dydt(7) = k_13 * y(6) - k_14 * y(7) + noise_level * randn() * g;

end
