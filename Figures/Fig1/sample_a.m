function dydt = sample_a(t,y,zt,z)

Z = interp1(zt,z,t);

%x y
dydt = zeros(1,1);
dydt(1) = Z.^3 ;

end