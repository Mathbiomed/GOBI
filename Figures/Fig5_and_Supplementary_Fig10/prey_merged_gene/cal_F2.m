function [score] = cal_F1(pred_net,true_net)

num_component = length(pred_net(:,1));

FN = 0; TP = 0; FP = 0; TN = 0;

for i = 1:num_component
    for j = 1:num_component
        if i ~= j
            pred_val = pred_net(i,j);
            true_val = true_net(i,j);

            if pred_val == 0
                if true_val == 0
                    TN = TN + 1;
                else
                    FN = FN + 1;
                end
            else
                if true_val == pred_val
                    TP = TP + 1;
                elseif true_val == 0
                    FP = FP + 1;
                else
                    FP = FP + 1;
                end
            end
        end
    end
end
Precision = TP / (TP + FP);
Recall = TP / (TP + FN);

score = 5/(4/Precision + 1/Recall);

end

