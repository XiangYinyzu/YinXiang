function [E_consume_sum, task_allocated_node, E_consume_allocated_node] = task_allocation_in_cluster(mu, num_cluster, node_cluster, L, E_elec, Efs, node_resource, task_resource, dist_head_node, allocated_task_cluster, C_residual_E)
%
% Record of revisions:
%   Data               Programmer            Description of language
%  ======            =============          =========================
% 20/09/2017          Xiang Yin                    Original
% 

node_in_allocated_cluster = find(node_cluster == allocated_task_cluster);   % find out the nodes which are in the allocated cluster
task_required_resource = find(task_resource == 1);     % remove "0" from task resource
task_allocated_node = zeros(1, length(task_required_resource));
E_consume_allocated_node = zeros(1, length(task_required_resource));
for i = 1 : length(task_required_resource)
    node_provide_required_resource_index = find(task_required_resource(i) == node_resource(node_in_allocated_cluster));  % nodes index in allocated cluster that provide required resource i
    node_provide_required_resource = node_in_allocated_cluster(node_provide_required_resource_index);       % find the nodes in the cluster that can provide the required resource i
    for j = 1 : length(node_provide_required_resource_index)
        if mod(node_provide_required_resource(j), num_cluster) == 0
            dist(i, j) = dist_head_node(allocated_task_cluster, (node_provide_required_resource(j) / num_cluster));
        else
            dist(i, j) = dist_head_node(allocated_task_cluster, fix(node_provide_required_resource(j) / num_cluster) + 1);
        end
        E_consume(i, j) = (L * E_elec + L * Efs * dist(i, j)^2);
        %p(i, j) = 1 / (times(node_provide_required_resource(j)) ^ alpha * (1 + E_consume(i, j) ^ belta));         % the probability of allocating subtask i to the node which can provide jth resource of the task
    end
    temp1 = 0;
    temp2 = 0;
    temp3 = 0;
    temp4 = 0;
    for j = 1 : length(node_provide_required_resource_index)
        temp1 = temp1 + (1 / E_consume(i, j));
    end
    temp2 = sum(C_residual_E(node_provide_required_resource_index));
    for j = 1 : length(node_provide_required_resource_index)
        temp3 = 1 / E_consume(i, j);
        temp4 = C_residual_E(j);
        p(i, j) = (1 - mu) * (temp3 / temp1)  + mu * (temp4 / temp2);
    end
    
    [max_p, max_p_index] = max(p(i, :));          % find out the maximum probability and the correlated index
    task_allocated_node(i) = node_provide_required_resource(max_p_index);
    E_consume_allocated_node(i) = E_consume(i, max_p_index);
    %task_allocated_node = [task_allocated_node, node_provide_required_resource(max_p_index)];
    %E_consume_allocated_node = [E_consume_allocated_node, E_consume(i, max_p_index)];     % record the energy consumption of the allocated node
end

E_consume_sum = sum(E_consume_allocated_node);  % the total energy consumption for executing the task

end

