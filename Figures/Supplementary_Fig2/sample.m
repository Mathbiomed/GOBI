function dydt = sample(t,y,zt,z)

Z = interp1(zt,z,t);

dydt = zeros(3,1);
dydt(1) = -Z.^3 ;
dydt(2) = y(1).^3 + Z.^5;
dydt(3) = -1.5*y(1).^3 - y(1).^5 + y(2).^3 - 0.5*Z.^5;
end