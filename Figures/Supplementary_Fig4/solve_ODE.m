function [A_time] = solve_ODE(parav1, parav2, parav3)


global a b c ...
K K2 n m

dmp=num2cell(parav1); kmp=num2cell(parav2); nmp=num2cell(parav3);


[a b c]=deal(dmp{:});
[K K2]=deal(kmp{:});
[n m]=deal(nmp{:});



t=0:0.01:50;

[t,x] = ode45(@ODE_no_indirect, t, [1]);
A_time = x;
end