function dydt = Goodwin(t,y)

dydt = zeros(4,1);
dydt(1) = 1 / (1 + y(4)^10) - 0.4 * y(1);
dydt(2) = y(1) - 0.4 * y(2);
dydt(3) = y(2) - 0.4 * y(3);
dydt(4) = y(3) - 0.4 * y(4);

end