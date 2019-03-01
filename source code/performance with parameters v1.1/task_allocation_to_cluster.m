function [WSN_E_consum, E_consume_allocated_node_result, allocated_task_cluster, task_allocated_node_result] = task_allocation_to_cluster(mu, num_node, num_task, num_cluster, task_cluster, node_cluster, L, E_elec, Efs, node_resource, task_resource, dist_head_node, ini_E)

%
% Record of revisions:
%   Data               Programmer            Description of language
%  ======            =============          =========================
% 25/09/2017          Xiang Yin                    Original
% 

allocated_task_cluster = zeros(1, num_task);
task_cluster_backup = task_cluster;                     
for i = 1 : num_task
    cluster_index = find(task_cluster_backup(i, :) == 1);     % find out the cluster that can execute task i
    allocated_task_cluster(i) = cluster_index(1);         
    task_cluster_backup(:, cluster_index(1)) = 0;         
end
stable_flag = 0;

% initialize parameters for dynamic programming
delta = 1;
sigma = 0.00001;

max_iteration = 100;       % the max iteration time
while stable_flag == 0
    v_task = zeros(1, num_task + 1);    % initialize state value function
    iteration = 0;
    task_allocated_node_result = [];
    E_consume_allocated_node_result = [];
    % policy evaluation
    while (delta > sigma) && (iteration <= max_iteration)
        delta = 0;
        %temp_times = zeros(1, num_node);
        for i = 1 : num_task
            v = v_task(i);
            [E_consume_sum, task_allocated_node, E_consume_allocated_node] = task_allocation_in_cluster(mu, num_cluster, node_cluster, L, E_elec, Efs, node_resource, task_resource(i, :), dist_head_node, allocated_task_cluster(i), ini_E)
            v_task(i) = E_consume_sum + v_task(i + 1);
            delta = max(delta, abs(v - v_task(i)));
        end
        iteration = iteration + 1;
    end

    % policy improvement
    task_cluster_backup = task_cluster;
    b = allocated_task_cluster;
    task_E_consum = zeros(1, num_task);
    for i = 1 : num_task
        task_cluster_index = find(task_cluster_backup(i, :) == 1);        % find out the clusters that can complete the task
        task_allocated_node_temp = [];
        E_consume_allocated_node_temp = [];
        for j = 1 : length(task_cluster_index)
            [E_consume_sum, task_allocated_node, E_consume_allocated_node] = task_allocation_in_cluster(mu, num_cluster, node_cluster, L, E_elec, Efs, node_resource, task_resource(i, :), dist_head_node, task_cluster_index(j), ini_E)
            v_task_all(i, j) = E_consume_sum + v_task(i + 1);
            task_allocated_node_temp(j, :) = task_allocated_node;    % record the temporary values of "task_allocated_node"
            E_consume_allocated_node_temp(j, :) = E_consume_allocated_node;  % the temporary values of energy consumed by each node
        end
        v_task_all_temp = v_task_all(i, :);
        v_task_all_temp(v_task_all_temp == 0) = [];       % delete element "0" from "v_task_all"
        [min_value, min_value_index] = min(v_task_all_temp);    % find out index of the allocated cluster that has the maximum value
        allocated_task_cluster_updated(i) = task_cluster_index(min_value_index);             % update the task allocated cluster
        task_cluster_backup(:, task_cluster_index(min_value_index)) = 0;             % if the cluster has been allocated with task, it can not be associated with any other task
        task_E_consum(i) = min_value;        % energy consumption for completing task i
        
        task_allocated_node_result = [task_allocated_node_result, task_allocated_node_temp(min_value_index, :)];    
        E_consume_allocated_node_result = [E_consume_allocated_node_result, E_consume_allocated_node_temp(min_value_index, :)];     
    end
    allocated_task_cluster = allocated_task_cluster_updated;
    WSN_E_consum = sum(task_E_consum);       

    if all(b == allocated_task_cluster_updated)
        stable_flag = 1;
    end
end

end

