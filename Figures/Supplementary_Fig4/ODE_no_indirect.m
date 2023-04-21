%% function used in Simulated Annealing 
%% ODE function assuming hill type form of regulation

function dC = ODE_no_indirect(t,C) % change interpolation to pchip

global a b c K K2 n m

%  ODE function

dC= a * input(t)^n / (K + input(t)^n)+ b * K2 / (K2 + input(t)^m) - c* C;

    function [outputArg1] = input(inputArg1)
        input_d = readmatrix('input(A).csv');
        outputArg1=makima(0:0.01:50,input_d,inputArg1);
    end

end

