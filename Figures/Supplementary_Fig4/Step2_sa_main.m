%% Step 0. Simulated Annealing (SA) code for parameter estimation

function Step2_sa_main(ran_num)

rng(ran_num);

global a b c K K2 n m % 7 parameters



result1=[]; % save the parameter
result2=[]; % save the error


data=readmatrix('data(B).csv');


for wb=1:200
    
    value=100000; % initial error value
    testmind=8;
    iteration=0;
    temp=1; % initial temperature is 1
    
    idran=0+(5-0)*rand(1,3); % a,b,c
    ikran=0+(20-0)*rand(1,2); % K K2
    inran = 1+(5-1)*rand(1,2); %n m
    

    [idran, ikran, inran] = adjust_range(idran, ikran, inran); % adjust the parameter range
    dir=idran; kir=ikran; nir=inran;
    
    maxp=0.4;
    hold=7; % threshold value for stopping the SA
    
    while iteration < 501 && value > hold

        
        % perturb the parameter
        thres1=max([-maxp -value/testmind]);
        thres2=min([maxp value/testmind]);
        
        perturb1=exp(thres1+(thres2-thres1)*rand(1,3));
        perturb2=exp(thres1+(thres2-thres1)*rand(1,2));
        perturb3=exp(thres1+(thres2-thres1)*rand(1,2));
        
        parav1=perturb1.*dir; parav2=perturb2.*kir; parav3=perturb3.*nir;
        
        
        [parav1, parav2, parav3] = adjust_range(parav1, parav2, parav3);
        
        [A_time] = solve_ODE(parav1, parav2, parav3); % solve the ODE using the parameter
        

        error = norm(data(:,1)-A_time); % calculate the error
        temp=0.9801*temp; % update the temperature
        iteration=iteration+1; % increase the iteration
        
        if error < value % new parameter fits the data better
            
            % update the value and parameters
            value=error; 
            paraset1=[a b c];
            paraset2=[K K2];
            paraset3=[n m];
            dir=paraset1; kir=paraset2; nir=paraset3;
            
            % display the value
            disp([value, temp, iteration/1000, a, b])
            
        elseif exp(-abs(value-error)/temp) > rand % new parameter is worse at fitting data, but due to high temperature, we update
     
            % update the value and parameters
            value=error;            
            paraset1=[a b c];
            paraset2=[K K2];
            paraset3=[n m];
            
            dir=paraset1; kir=paraset2; nir=paraset3;
            
            % display the value
            disp([1,value, iteration/1000, a, b])
            
        else % new parameter is worse at fitting data and we do not update
            disp([2, value, a, b])
        end
    end
    
    [result1,result2] = save_to_csv(ran_num, iteration, value, temp, result1,result2); % after the simulation, save the data
    
end

end