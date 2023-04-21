function dydt = Goodwin_noise(t,y, noise_level, a, b, c, d)

dydt = zeros(4,1);
dydt(1) = 1 / (1 + y(4)^10) - 0.4 * y(1) + noise_level * randn() * a;
dydt(2) = y(1) - 0.4 * y(2) + noise_level * randn() * b;
dydt(3) = y(2) - 0.4 * y(3) + noise_level * randn() * c;
dydt(4) = y(3) - 0.4 * y(4) + noise_level * randn() * d;

end