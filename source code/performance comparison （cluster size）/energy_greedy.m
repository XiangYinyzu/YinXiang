function [EG_WSN_E_consume, EG_residual_E, EG_allocated_task_cluster, unallocated_task_num, EG_task_allocated_node] = energy_greedy(radius, num_task,  num_cluster, node_cluster, EG_residual_E, L, E_elec, Efs, node_resource, task_resource, task_cluster, dist_head_node)
% Record of revisions:
%   Data               Programmer            Description of language
%  ======            =============          =========================
% 29/11/2017          Xiang Yin                    Original
% 

EG_task_allocated_node = [];

unallocated_task_num = 0;
task_cluster_backup = task_cluster;
EG_WSN_E_consume = 0; 
EG_allocated_task_cluster = zeros(1, num_task);

for i = 1 : num_task
    cluster_index = find(task_cluster_backup(i, :) == 1);     % find out the cluster that can execute task i
%     temp = randperm(length(cluster_index));
%     if isempty(temp)
%         unallocated_task_num = unallocated_task_num + 1;
%         continue;
%     end
    if isempty(cluster_index)
        unallocated_task_num = unallocated_task_num + 1;
        continue;
    end
    EG_allocated_task_cluster(i) = cluster_index(randi([1, length(cluster_index)]));
    task_cluster_backup(:, cluster_index(randi([1, length(cluster_index)]))) = 0;
    %EG_allocated_task_cluster(i) = cluster_index(temp(1));         % assign task i to the first cluster that can complete it, this is an initial allocation
    %task_cluster_backup(:, cluster_index(temp(1))) = 0;         % each cluster can only execute one task in each application
    
    node_in_allocated_cluster = find(node_cluster == EG_allocated_task_cluster(i));   % find out the nodes which are in the allocated cluster
    task_required_resource = find(task_resource(i, :) == 1);
    for j = 1 : length(task_required_resource)
        node_provide_required_resource_index = find(task_required_resource(j) == node_resource(node_in_allocated_cluster));  % nodes index in allocated cluster that provide required resource i
        node_provide_required_resource = node_in_allocated_cluster(node_provide_required_resource_index);    % find the nodes in the cluster that can provide the required resource j
        dist = zeros(1, length(node_provide_required_resource_index));
        % compute the distance between the cluster head and the nodes that can provide resource j
        for k = 1 : length(node_provide_required_resource_index)
            if mod(node_provide_required_resource(k), num_cluster) == 0
                dist(k) = dist_head_node(EG_allocated_task_cluster(i), (node_provide_required_resource(k) / num_cluster));
            else
                dist(k) = dist_head_node(EG_allocated_task_cluster(i), fix(node_provide_required_resource(k) / num_cluster) + 1);
            end
            if dist(k) == 0
                dist(k) = 100000;         
            end
        end

        if k == 1
            dist(k) = radius / 2;
        end
        [sorted_dist, index] = sort(dist);     % sort the distance in ascending order
        % allocate sub-task to the node which will consume the least energy
        EG_task_allocated_node(i, j) = node_provide_required_resource(index(1));
        EG_energy_consume = (L * E_elec + L * Efs * dist(index(1))^2);
        EG_residual_E(node_provide_required_resource(index(1))) = EG_residual_E(node_provide_required_resource(index(1))) - EG_energy_consume;         % update the residual of the node
        %        for k = 1 : length(node_provide_required_resource_index)
        %            if (EG_residual_E(node_provide_required_resource(index(k))) > 0)        % the residual energy of the node with the minimum energy consumption should be positive
        %                EG_task_allocated_node(i, j) = node_provide_required_resource(index(k));
        %                EG_energy_consume = (L * E_elec + L * Efs * dist(index(k))^2);
        %                EG_residual_E(node_provide_required_resource(index(k))) = EG_residual_E(node_provide_required_resource(index(k))) - EG_energy_consume;         % update the residual of the node
        %                break;
        %            else
        %                EG_energy_consume = 0;
        %            end
        %       end
        EG_WSN_E_consume = EG_WSN_E_consume + EG_energy_consume;
    end
end

end
