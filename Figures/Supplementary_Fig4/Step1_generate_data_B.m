%% Step 1. Generating time series of B using ODE model

t=0:0.01:50; % time span



[t,x] = ode45(@ODE_function, t, [1]); % solve ODE for B using ODE 'ODE_function' with initial condition B(0)=1

figure(2)
subplot(2,1,1) % time series of A
plot(t,x2, '-r');
title('time series of component A (input)')

subplot(2,1,2) % time series of B
plot(t,x, '-b');
title('time series of component B')


csvwrite('data(B).csv', x) % save the data to csv file 
