%% Step 3. Analyze the result using 

global a b c K K2 n m


data=readmatrix('data.csv');

a_time=[];

para =  readmatrix('para_full_wrong_mixed3.csv');

figure(3) % boxplot comparing alpha and beta
boxplot(para(1:100,1:2),'Labels',{'positive','negative'})

for i=1:size(para,1)
    i
	paraset = para(i,:);
    tmp=num2cell(paraset);
    [a b c K K2 n m] = deal(tmp{:});
    
    [A_time] = solve_ODE([a,b,c],[K, K2],[n, m]);
    
    
    a_time=[a_time; A_time'];
end


csvwrite('simul_result.csv', [mean(a_time);std(a_time)]);




