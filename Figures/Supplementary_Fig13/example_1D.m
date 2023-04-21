function dydt = example_1D(t,y,bt,b)

B = interp1(bt,b,t);

dydt = zeros(1,1);
dydt(1) = B;
end