function dydt = CFL(t,y,zt,z)

Z = interp1(zt,z,t);

dydt = zeros(2,1);
dydt(1) = -Z;
dydt(2) = y(1) - Z.^3 - 0.2*Z ;

end