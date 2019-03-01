function energy_consumption_comparison
%
%
% Record of revisions:
%   Data               Programmer            Description of language
%  ======            =============          =========================
% 18/09/2017          Xiang Yin                    Original
% 

% alpha = 0.5 : 0.5 : 4;
% belta = 0.01 : 0.01 : 0.06;
mu = 0.1 : 0.1 : 0.9;

num_app = 300;          % number of application
num_task = 15;          % number of tasks in an application

%%%%%%%%%%%%%%%%%%%% parameters for WSN %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Field Dimensions - we account for a circular area 
radius = 500;             % the radius of the circular area

% number of resource type in the system
resource_type = 6;

% Coordinates of Sink node which is placed on origin
sink_x = 0;
sink_y = 0;

% Number of nodes
num_node = 400;
num_node_in_cluster = 20;

% Number of cluster
num_cluster = num_node / num_node_in_cluster;  % each cluster contains 10 sensor nodes on average

% Initialize coordinate of each sensor node, here, the sensor nodes are
% uniformly distributed in the area
node_cluster = zeros(1, num_node);
for i = 1 : num_node
    rand_radius = radius * rand;
    temp = mod(i, num_cluster);
    if temp == 0
        temp = num_cluster;
    end
    node_position_x(i) = rand_radius * cos((temp - 1) * (2 * pi / num_cluster) + (2 * pi / num_cluster) * rand);
    node_position_y(i) = rand_radius * sin((temp - 1) * (2 * pi / num_cluster) + (2 * pi / num_cluster) * rand);
    node_cluster(i) = temp;         % the cluster node i is located
end
figure(1);
rectangle('Position',[-radius,-radius,2 * radius,2 * radius],'Curvature',[1,1]),axis equal;
hold on;
plot(node_position_x, node_position_y, 'bo');

%%%%%%%%%%%% parameters for sensor nodes %%%%%%%%%%%%%%%
node_resource = zeros(1, num_node);
node_resource = randi([1, resource_type], 1, num_node);   % assign each node a specific resource type

%%%%%%%%%%%% parameters for sensor energy %%%%%%%%%%%%%%%

E_elec = 50 * 0.000000001;     % electroincs energy, which is the same for transimtter and receiver
Efs = 10 * 0.000000000001;      % amplifier energy for free space channel
% E_elec = 50 * 0.0000001;     % electroincs energy, which is the same for transimtter and receiver
% Efs = 10 * 0.0000000001;      % amplifier energy for free space channel
L = randi([40000, 60000]);       % bits of information for completing a sub-task

%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%

% in each cluster, randomly designate a node as the cluster head
cluster_head = zeros(1, num_cluster);
for i = 1 : num_cluster
    cluster_head(i) = i;
end
% compute the distance between each pair of cluster head and the sink node
dist_sink_head = zeros(1, num_cluster);
for i = 1 : num_cluster
    dist_sink_head(i) = sqrt((node_position_x(cluster_head(i))) ^ 2 + (node_position_y(cluster_head(i))) ^ 2);
end

% compute the distance between each node and the cluster head
dist_head_node = zeros(num_cluster, num_node_in_cluster);
for i = 1 : num_node
    cluster_index = node_cluster(i);
    index = fix(i / num_cluster) + 1;      % when i=num_cluster, a singular value appears, it is handled below
    dist_head_node(cluster_index, index) = sqrt((node_position_x(cluster_index) - node_position_x(i)) ^ 2 + (node_position_y(cluster_index) - node_position_y(i)) ^ 2);
end
% adjust the matrix of cluster head node distance
for i = 1 : num_node_in_cluster
    dist_head_node(num_cluster, i) = dist_head_node(num_cluster, i + 1);
end
dist_head_node(:, num_node_in_cluster + 1) = [];

% count the resources of each cluster
cluster_resource = zeros(num_cluster, resource_type);
for i = 1 : num_node
    temp = node_resource(i);
    cluster_resource(node_cluster(i), temp) = 1;
end

task_resource = cell(1, num_app);
task_cluster = cell(1, num_app);
for n = 1 : num_app

    task_resource_temp = zeros(num_task, resource_type);
    for i = 1 : num_task
        task_resource_type = randi([2, resource_type]);  % for a task, it demands a specific number of resources which is in the range [2, application_resource_type]
        temp = randperm(resource_type);
        task_resource_temp(i, temp(1 : task_resource_type)) = 1;
    end
    task_resource{1, n} = task_resource_temp;
    
    task_cluster_temp = zeros(num_task, num_cluster);
    for i = 1 : num_task
        for j = 1 : num_cluster
            if all(cluster_resource(j, :) >= task_resource_temp(i, :))  %the resources of the cluster cover the required resources of the task
                task_cluster_temp(i, j) = 1;
            end
        end
    end
    task_cluster{1, n} = task_cluster_temp;
end

WSN_E_consum_all = zeros(1, length(mu));
std_residual_E = zeros(1, length(mu));

% for k = 1 : length(alpha)
%     for m = 1 : length(belta)
%         %times = 0.5 * ones(1, num_node);    % initialize the times a node has been used
%         WSN_E_consum_app = zeros(1, num_app);
%         ini_E = 2 * ones(1, num_node);             % initial energy of each node
%         for n = 1 : num_app
%             % allocate tasks to clusters using Dynamic Programming
%             [WSN_E_consum, E_consume_allocated_node_result, allocated_task_cluster, task_allocated_node_result] = task_allocation_to_cluster(mu, num_node, num_task, num_cluster, task_cluster{1, n}, node_cluster, L, E_elec, Efs, node_resource, task_resource{1, n}, dist_head_node, ini_E, alpha(k), belta(m))
%             for i = 1 : length(task_allocated_node_result)
%                 %times(task_allocated_node_result(i)) = times(task_allocated_node_result(i)) + 1;     % update the times a node has been used
%                 ini_E(task_allocated_node_result(i)) = ini_E(task_allocated_node_result(i)) - E_consume_allocated_node_result(i);    % update the residual of each node
%             end
%             
%             E_consum_head_sink = zeros(1, num_task);
%             for i = 1 : num_task
%                 L_task = length(find(task_resource{1, n}(i, :) ~= 0));      % the number of required resources of task i
%                 E_consum_head_sink(i) = (L_task * E_elec + L_task * Efs * dist_sink_head(allocated_task_cluster(i))^2);   % energy consumption for sending data from cluster head of task i to sink node
%             end
%             
%             WSN_E_consum_app(n) = WSN_E_consum + sum(E_consum_head_sink);             % record the overall energy consumption in each application
%         end
%         WSN_E_consum_belta(k, m) = sum(WSN_E_consum_app);
%         
%         ini_E(1 : num_cluster) = [];         % remove the residual energy of cluster head when computing the std of residual energy
%         std_residual_E_belta(k, m) = std(ini_E);     % calculate the standard deviation of the residual energy of nodes
%     end
% end

% allocated_task_cluster_matrix = cell(1, length(mu));
% task_allocated_node_result_matrix = cell(1, length(mu));
for k = 1 : length(mu)
    WSN_E_consum_app = zeros(1, num_app);
    ini_E = 12 * ones(1, num_node);             % initial energy of each node
    for n = 1 : num_app
        % allocate tasks to clusters using Dynamic Programming
        [WSN_E_consum, E_consume_allocated_node_result, allocated_task_cluster, task_allocated_node_result] = task_allocation_to_cluster(mu(k), num_node, num_task, num_cluster, task_cluster{1, n}, node_cluster, L, E_elec, Efs, node_resource, task_resource{1, n}, dist_head_node, ini_E)
%         allocated_task_cluster_matrix{1, k}(n, :) = allocated_task_cluster;
%         task_allocated_node_result_matrix{1, k}(n, 1 : length(task_allocated_node_result)) = task_allocated_node_result;
        
        for i = 1 : length(task_allocated_node_result)
            %times(task_allocated_node_result(i)) = times(task_allocated_node_result(i)) + 1;     % update the times a node has been used
            ini_E(task_allocated_node_result(i)) = ini_E(task_allocated_node_result(i)) - E_consume_allocated_node_result(i) * 5;    % update the residual of each node
        end
        
        E_consum_head_sink = zeros(1, num_task);
        for i = 1 : num_task
            L_task = length(find(task_resource{1, n}(i, :) ~= 0));      % the number of required resources of task i
            E_consum_head_sink(i) = (L_task * E_elec + L_task * Efs * dist_sink_head(allocated_task_cluster(i))^2);   % energy consumption for sending data from cluster head of task i to sink node
        end
        
        WSN_E_consum_app(n) = WSN_E_consum + sum(E_consum_head_sink);             % record the overall energy consumption in each application
    end
    WSN_E_consum_all(k) = sum(WSN_E_consum_app);
    
    ini_E(1 : num_cluster) = [];         % remove the residual energy of cluster head when computing the std of residual energy
    E_std_residual(k) = std(ini_E);     % calculate the standard deviation of the residual energy of nodes
end

figure(2);
plot(mu, WSN_E_consum_all, 'ro', mu, WSN_E_consum_all, '--');
xlabel('\mu');
ylabel('energy consumption');

figure(3);
plot(mu, E_std_residual, 'bx', mu, E_std_residual, '--');
xlabel('\mu');
ylabel('standard deviation of residual energy');

end