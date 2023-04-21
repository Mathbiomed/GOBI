function dydt = example_3D(t,y,bt,b,ct,c,dt,d)

B = interp1(bt,b,t);
C = interp1(ct,c,t);
D = interp1(dt,d,t);

dydt = zeros(1,1);
%dydt(1) = B - C - y(1);
dydt(1) = B + C + D;
end