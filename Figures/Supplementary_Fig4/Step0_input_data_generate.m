%% Step 0. Generating time series of input signal A

rng(10); % for the reproducibility
t=0:0.01:50; % time span


x1=normrnd(1,0.5,1,51); % pick 51 points from normal distribution N(1, 0.5)

x2=spline(0:1:50, x1, t); % use spline to smoothly connect the points

figure(1) % plot the input dat
plot(0:1:50,x1,'o', t, x2)
title('input signal (C)') 

csvwrite('input(A).csv', x2) % save the data to csv file 
