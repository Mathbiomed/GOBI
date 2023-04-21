clc;
clear;
close all;

%% parameters
trial_list = [0:9];
noise_list = [2:2:20];
%noise_list = [2];

dimension = 1;
thres_noise = 1e-7;
%system_list = {'cAMP','Fr','Gb','GW','KF'};
%noise_type_list = {'additive','blue','brown','pink','purple','dynamical','multiplicative'};
noise_type_list = {'additive','blue','pink','multiplicative'};
%noise_type_list = {'additive'};

system_list = {'cAMP'};

for noise_type_idx = 1:length(noise_type_list)
    noise_type = char(noise_type_list(noise_type_idx));
    disp(noise_type)
    for system_idx = 1:length(system_list)
        system = char(system_list(system_idx));
        disp(system)
        for trial = trial_list
            for noise_percent = noise_list
                 %% load data
                filename = ['./RDS_dim1_',noise_type,'/',system,'_results_dim1_',num2str(noise_percent),'_Trial',num2str(trial)];
                load(filename)

                %% choose threshold using guide
                t_S = 0.9 - 0.005 * noise_percent;
                t_L = 0.05;
                t_TRS = 0.9 - 0.01 * noise_percent;

                %% Calculate Total Regulation Score (TRS)
                TRS_dim1 = zeros(num_pair,2^dimension);
                for i = 1:num_pair
                    for j = 1:num_type
                        S_tmp = reshape(S_total(i,j,:),[num_data,1]);
                        L_tmp = reshape(L_total(i,j,:),[num_data,1]);
                        S_processed = S_threshold(S_tmp, t_S);
                        L_processed = L_threshold(L_tmp, t_L);
                        N = sum(L_processed);
                        TRS_dim1(i,j) = sum(S_processed.*L_processed)/N;
                    end
                end

                %% infer 1D regulation using the criteria TRS > TRS^thres
                regulation_1dim = zeros(num_pair,2^dimension);
                for i = 1:num_pair
                    for j = 1:num_type
                        if TRS_dim1(i,j) >= t_TRS
                            regulation_1dim(i,j) = 1;
                        end
                    end
                end
                %% Save results
                filename = ['./TRS_dim1_',noise_type,'/',system,'_TRS_dim1_',num2str(noise_percent),'_Trial',num2str(trial)];

                save(filename, 'regulation_1dim','TRS_dim1', 'component_list_dim1','num_data','num_type','num_pair','dimension')
                %[regulation_1dim(1,1),regulation_1dim(7,1),regulation_1dim(11,1),regulation_1dim(6,2)]
                [TRS_dim1(2,2), TRS_dim1(4,2), TRS_dim1(6,2)]
            end
        end
    end
end
