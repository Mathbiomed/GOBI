function dydt = example_2D(t,y,bt,b,ct,c)

B = interp1(bt,b,t);
C = interp1(ct,c,t);

dydt = zeros(1,1);
%dydt(1) = B - C - y(1);
dydt(1) = B + C;
end