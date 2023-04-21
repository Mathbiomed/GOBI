%% ODE model 
function dxdt = func(t,x)

dxdt = (1/3*input(t)^3-input(t)^2+input(t))+1-x(1);


    function [outputArg1] = input(inputArg1) % input data
        input_d = readmatrix('input(A).csv');
        outputArg1=makima(0:0.01:50,input_d,inputArg1);
    end
end
