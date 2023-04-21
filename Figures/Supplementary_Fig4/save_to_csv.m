%% function used in Simulated Annealing 
%% save the result

function [result1,result2]=save_to_csv(ran_num, iteration, value, temp, result1,result2)
global a b c K K2 n m

            result1=[result1; [a b c K K2 n m]];
            result2=[result2; [value, temp, iteration]];
            csvwrite(['para_',num2str(ran_num),'.csv'],result1)
            csvwrite(['error_',num2str(ran_num),'.csv'],result2)  
            
end