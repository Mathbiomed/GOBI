function [results] = isIndirect(ori_node,prev_node,current_node,target_node,prev_type,prev_results,network)
    
results = [];
%% if current node is same with target, terminate
if current_node == target_node
    results = [results;prev_type];
%% else, search the path 1 step
else
    for i = 1:length(network(1,:))    
        if network(current_node,i) ~= 0 && i ~= prev_node && i ~= ori_node
            type_tmp = network(current_node,i);
            results = [results ; isIndirect(ori_node,current_node, i,target_node,prev_type * type_tmp,prev_results,network)];
        end
    end
end
end

