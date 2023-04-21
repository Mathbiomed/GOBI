function dydt = sample(t,y,zt,z)

Z = interp1(zt,z,t);

dydt = zeros(2,1);
dydt(1) = -Z ;
dydt(2) = y(1) + 2*Z.^3;
end