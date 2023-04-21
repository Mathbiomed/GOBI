function dydt = SFL(t,y,zt,z)

Z = interp1(zt,z,t);

dydt = zeros(2,1);
dydt(1) = -Z;
dydt(2) = 0.5 + y(1);

end