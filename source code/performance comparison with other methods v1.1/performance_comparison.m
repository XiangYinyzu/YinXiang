function performance_comparison
%
% Record of revisions:
%   Data               Programmer            Description of language
%  ======            =============          =========================
% 18/09/2017          Xiang Yin                    Original
% 

mu = 0.2;
num_app = 350;          % number of application
num_task = 15;          % number of tasks in an application

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
num_cluster = num_node / num_node_in_cluster;  

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
% rectangle('Position',[-radius,-radius,2 * radius,2 * radius],'Curvature',[1,1]),axis equal;
% hold on;
% plot(node_position_x, node_position_y, 'bo');

node_resource = randi([1, resource_type], 1, num_node);   % assign each node a specific resource type

ini_E = 1 * ones(1, num_node);             % initial energy of each node
E_elec = 50 * 0.000000001;     % electroincs energy, which is the same for transimtter and receiver
Efs = 10 * 0.000000000001;      % amplifier energy for free space channel
L = randi([40000, 60000]);       % bits of information for completing a sub-task

w1 = 0.3;  % W1-W3 为fitness权重系数
w2 = 0.3;
maxgen = 30;
swarm_size = 10;
c1 = 2;
c2 = 2;
% Wmax = 0.9;
% Wmin = 0.4;
w = rand(1) * (0.4 - 0.9) + 0.9;

cluster_head = zeros(1, num_cluster);
for i = 1 : num_cluster
    cluster_head(i) = i;
end

dist_sink_node = sqrt(power(node_position_x, 2) + power(node_position_y, 2));

dist_sink_head = zeros(1, num_cluster);
for i = 1 : num_cluster
    dist_sink_head(i) = sqrt((node_position_x(cluster_head(i))) ^ 2 + (node_position_y(cluster_head(i))) ^ 2);
end

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

C_residual_E = ini_E;
EG_residual_E = ini_E;
PSO_residual_E = ini_E;

C_node_resource = node_resource;
EG_node_resource = node_resource;
PSO_node_resource = node_resource;

C_WSN_E_consume_app = zeros(1, num_app);
EG_WSN_E_consume_app = zeros(1, num_app);
PSO_WSN_E_consume_app = zeros(1, num_app);

std_C_residual_E = zeros(1, num_app);
std_EG_residual_E = zeros(1, num_app);
std_PSO_residual_E = zeros(1, num_app);

C_unallocated_task_num = zeros(1, num_app);             % the number of task that can not be assigned for the sake of energy depletion of some nodes
EG_unallocated_task_num = zeros(1, num_app);
PSO_unallocated_task_num = zeros(1, num_app);

C_live_node_num = zeros(1, num_app);
EG_live_node_num = zeros(1, num_app);
PSO_live_node_num = zeros(1, num_app);

C_energy_consumption_per_task = zeros(1, num_app);
EG_energy_consumption_per_task = zeros(1, num_app);
PSO_energy_consumption_per_task = zeros(1, num_app);

C_allocated_task_cluster_all = [];
EG_allocated_task_cluster_all = [];

for m = 1 : num_app
    C_num_task = num_task;
    EG_num_task = num_task;
    PSO_num_task = num_task;
    
    % a flag that indicates whether tasks can be allocated, it is "1" if
    % any tasks can be assigned
    C_task_allocation_flag = 1;
    EG_task_allocation_flag = 1;

    task_resource = zeros(num_task, resource_type);
    for i = 1 : num_task
        task_resource_type = randi([2, resource_type]);  % for a task, it demands a specific number of resources which is in the range [2, application_resource_type]
        temp = randperm(resource_type);
        task_resource(i, temp(1 : task_resource_type)) = 1;
    end
    C_task_resource = task_resource;
    EG_task_resource = task_resource;
    PSO_task_resource = task_resource;
    
%     if any(EG_residual_E < 0)
%         keyboard;
%     end
    
    C_node_resource(find(C_residual_E < 0)) = 0;
    EG_node_resource(find(EG_residual_E < 0)) = 0;
    PSO_node_resource(find(PSO_residual_E < 0)) = 0;
    
    C_cluster_resource = zeros(num_cluster, resource_type);
    EG_cluster_resource = zeros(num_cluster, resource_type);
    for i = 1 : num_node
        C_temp = C_node_resource(i);
        if C_temp ~= 0;
            C_cluster_resource(node_cluster(i), C_temp) = 1;
        end
    end
    
    for i = 1 : num_node
        EG_temp = EG_node_resource(i);
        if EG_temp ~= 0;
            EG_cluster_resource(node_cluster(i), EG_temp) = 1;
        end
    end
    
    C_task_cluster = zeros(num_task, num_cluster);
    for i = 1 : num_task
        for j = 1 : num_cluster
            if all(C_cluster_resource(j, :) >= task_resource(i, :))  % the resources of the cluster cover the required resources of the task
                C_task_cluster(i, j) = 1;
            end
        end
        if all(C_task_cluster(i, :) == 0)        % if no cluster can perform the task, the task is signed as an unallocated task
            C_unallocated_task_num(m) = C_unallocated_task_num(m) + 1;
            C_num_task = C_num_task - 1;
        end
    end
    % if a task can not be allocated, remove it and the corresponding
    % resource requirment
    for i = num_task : -1 : 1
        if all(C_task_cluster(i, :) == 0)
            C_task_resource(i, :) = [];                 % delete the ith task from current application
            C_task_cluster(i, :) = [];
        end
    end
    % if all tasks in current application can not be assigned, set the
    % corresponding flag to 0
    if isempty(C_task_cluster)
        C_task_allocation_flag = 0;
    end

    EG_task_cluster = zeros(num_task, num_cluster);
    for i = 1 : num_task
        for j = 1 : num_cluster
            if all(EG_cluster_resource(j, :) >= task_resource(i, :))  %the resources of the cluster cover the required resources of the task
                EG_task_cluster(i, j) = 1;
            end
        end
        if all(EG_task_cluster(i, :) == 0)
            EG_unallocated_task_num(m) = EG_unallocated_task_num(m) + 1;
            EG_num_task = EG_num_task - 1;
        end
    end
    for i = num_task : -1 : 1
        if all(EG_task_cluster(i, :) == 0)
            EG_task_resource(i, :) = [];                 % delete the ith task from current application
            EG_task_cluster(i, :) = [];
        end
    end
    if isempty(EG_task_cluster)
        EG_task_allocation_flag = 0;
    end
    
    % count the task that can not be allocated by PSO algorithm because the
    % depletion of energy
    PSO_unallocated_task_flag = zeros(1, num_task);         % indicate whether a task can be allocated successfully
    if any(PSO_node_resource < 0)
        for i = 1 : num_task
            for j = 1 : resource_type
                while task_resource(i, j) == 1
                    if ismember(j, PSO_node_resource) == 0
                        PSO_unallocated_task_num(m) = PSO_unallocated_task_num(m) + 1;
                        PSO_unallocated_task_flag(i) = 1;
                    end
                    break;
                end
            end
        end
        for i = num_task : -1 : 1
            if PSO_unallocated_task_flag(i) == 1
                PSO_task_resource(i, :) = [];
            end
        end
        PSO_num_task = num_task - length(find(PSO_unallocated_task_flag));
    end
    
    if C_task_allocation_flag
        [C_WSN_E_consume, E_consume_allocated_node_result, allocated_task_cluster, task_allocated_node_result, unallocated_task_num] = task_allocation_to_cluster(mu, C_num_task, num_cluster, C_cluster_resource, C_task_cluster, node_cluster, L, E_elec, Efs, C_node_resource, C_task_resource, dist_head_node, C_residual_E)
    end
    
    for i = 1 : length(task_allocated_node_result)
        %times(task_allocated_node_result(i)) = times(task_allocated_node_result(i)) + 1;     % update the times a node has been used
        C_residual_E(task_allocated_node_result(i)) = C_residual_E(task_allocated_node_result(i)) - E_consume_allocated_node_result(i);    % update the residual of each node
    end
    std_C_residual_E(m) = std(C_residual_E);                % calculate the standard deviation of the residual energy of nodes
    C_live_node_num(m) = num_node - length(find(C_residual_E <= 0));        % the live node in each application
    if num_task - C_unallocated_task_num(m) - unallocated_task_num == 0
        C_energy_consumption_per_task(m) = 0;
    else
        C_energy_consumption_per_task(m) = C_WSN_E_consume / (num_task - C_unallocated_task_num(m) - unallocated_task_num) %(num_task - C_unallocated_task_num(m) - unallocated_task_num);     % average energy consumption for one task in the application
    end
    C_unallocated_task_num(m) = C_unallocated_task_num(m) + unallocated_task_num; % add the task that can not be allocated in cluster
    
    C_allocated_task_cluster_all = [C_allocated_task_cluster_all, allocated_task_cluster];
    
    E_consum_head_sink = zeros(1, num_task);
    for i = 1 : length(allocated_task_cluster)
        L_task = length(find(task_resource(i, :) ~= 0));      % the number of required resources of task i
        E_consum_head_sink(i) = (L_task * E_elec + L_task * Efs * dist_sink_head(allocated_task_cluster(i))^2);   % energy consumption for sending data from cluster head of task i to sink node
    end
    C_WSN_E_consume_app(m) = C_WSN_E_consume + sum(E_consum_head_sink);             % record the overall energy consumption in each application
    
%     if m == 60
%         keyboard;
%     end      
    
    if EG_task_allocation_flag
        [EG_WSN_E_consume, EG_residual_E, EG_allocated_task_cluster, unallocated_task_num, EG_task_allocated_node] = energy_greedy(radius, EG_num_task, num_cluster, node_cluster, EG_residual_E, L, E_elec, Efs, EG_node_resource, EG_task_resource, EG_task_cluster, dist_head_node)   % assign tasks by energy greedy algorithm
    end
    
    testa = find(EG_node_resource == 0);
    testb = ismember(testa, EG_task_allocated_node);
    if any(testb == 1)
        keyboard;
    end
    
    EG_E_consum_head_sink = zeros(1, num_task);
    for i = 1 : EG_num_task
        if EG_allocated_task_cluster(i) == 0
            continue;
        else
            L_task = length(find(task_resource(i, :) ~= 0));      % the number of required resources of task i
            EG_E_consum_head_sink(i) = (L_task * E_elec + L_task * Efs * dist_sink_head(EG_allocated_task_cluster(i))^2);   % energy consumption for sending data from cluster head of task i to sink node
        end
    end
    EG_WSN_E_consume_app(m) = EG_WSN_E_consume + sum(EG_E_consum_head_sink);             % record the overall energy consumption in each application
    std_EG_residual_E(m) = std(EG_residual_E);                % calculate the standard deviation of the residual energy of nodes
    EG_live_node_num(m) = num_node - length(find(EG_residual_E <= 0));
    if num_task - EG_unallocated_task_num(m) - unallocated_task_num == 0
        EG_energy_consumption_per_task(m) = 0;
    else
        EG_energy_consumption_per_task(m) = EG_WSN_E_consume / (num_task - EG_unallocated_task_num(m) - unallocated_task_num);
    end
%     EG_AA(m) = unallocated_task_num;
%     EG_BB(m) = EG_WSN_E_consume;
    EG_unallocated_task_num(m) = EG_unallocated_task_num(m) + unallocated_task_num;
    
    EG_allocated_task_cluster_all = [EG_allocated_task_cluster_all, EG_allocated_task_cluster];
    
    % allocate tasks by PSO
     %按资源类型将节点归类
     R_type = cell(1, resource_type);
     for i = 1 : resource_type
         R_type{i} = find(node_resource == i);     % to a kind of resource, R_type records which nodes have that resource
     end
 
    T_cell = cell(1, PSO_num_task);
    for i = 1 : PSO_num_task
        PSO_task_resource_index = find(PSO_task_resource(i, :) == 1);
        T_cell{i} = PSO_task_resource_index;
    end
    
    [PSO_WSN_E_consume, PSO_task_allocated_node_result, unallocated_task_num] = pso(swarm_size, maxgen, w, w1, w2, c1, c2, L, E_elec, Efs, PSO_num_task, resource_type, dist_sink_node, T_cell, PSO_node_resource, PSO_task_resource, R_type)
    PSO_WSN_E_consume_app(m) = PSO_WSN_E_consume;
    % compute the residual energy of each node
    PSO_task_allocated_node_result(PSO_task_allocated_node_result == 0) = [];
    for i = 1 : length(PSO_task_allocated_node_result)
        temp = PSO_task_allocated_node_result(i);
        PSO_residual_E(temp) = PSO_residual_E(temp) - (L * E_elec + L * Efs * dist_sink_node(temp)^2);
    end
    
    std_PSO_residual_E(m) = std(PSO_residual_E);
    PSO_live_node_num(m) = num_node - length(find(PSO_residual_E <= 0));
    if num_task - PSO_unallocated_task_num(m) - (unallocated_task_num / swarm_size) == 0
        PSO_energy_consumption_per_task(m) = 0;
    else
        PSO_energy_consumption_per_task(m) = PSO_WSN_E_consume / (num_task - PSO_unallocated_task_num(m) - (unallocated_task_num / swarm_size));
    end
    PSO_unallocated_task_num(m) = PSO_unallocated_task_num(m) + (unallocated_task_num / swarm_size);
end

C_allocated_task_cluster_each = zeros(1, num_cluster);
EG_allocated_task_cluster_each = zeros(1, num_cluster);
for i = 1 : num_cluster
    C_allocated_task_cluster_each(i) = length(find(C_allocated_task_cluster_all == i));
    EG_allocated_task_cluster_each(i) = length(find(EG_allocated_task_cluster_all == i));
end

std_C_allocated_task_cluster_each = std(C_allocated_task_cluster_each);
std_EG_allocated_task_cluster_each = std(EG_allocated_task_cluster_each);

C_WSN_E_consume_all = zeros(1, num_app);
EG_WSN_E_consume_all = zeros(1, num_app);
PSO_WSN_E_consume_all = zeros(1, num_app);

C_unallocated_task_num_all = zeros(1, num_app);
EG_unallocated_task_num_all = zeros(1, num_app);
PSO_unallocated_task_num_all = zeros(1, num_app);

for i = 1 : num_app
    C_WSN_E_consume_all(i) = sum(C_WSN_E_consume_app(1 : i));
    EG_WSN_E_consume_all(i) = sum(EG_WSN_E_consume_app(1 : i));
    PSO_WSN_E_consume_all(i) = sum(PSO_WSN_E_consume_app(1 : i));
    
    C_unallocated_task_num_all(i) = sum(C_unallocated_task_num(1 : i));
    EG_unallocated_task_num_all(i) = sum(EG_unallocated_task_num(1 : i));
    PSO_unallocated_task_num_all(i) = sum(PSO_unallocated_task_num(1 : i));
end
app_index = 1 : num_app;
figure(1);
plot(app_index, C_WSN_E_consume_all, 'r--', app_index, EG_WSN_E_consume_all, 'b-.', app_index, PSO_WSN_E_consume_all, 'g-');
xlabel('application number');
ylabel('energy consumption');
legend('TP-CAA-HH', 'G-CAA-HH', 'BPSO-CAA-H');

figure(2);
plot(app_index, std_C_residual_E, 'r--', app_index, std_EG_residual_E, 'b-.', app_index, std_PSO_residual_E, 'g-')
xlabel('application number');
ylabel('standard deviation of residual energy');
legend('TP-CAA-HH', 'G-CAA-HH', 'BPSO-CAA-H');

figure(3)
plot(app_index, C_unallocated_task_num_all, 'r--', app_index, EG_unallocated_task_num_all, 'b-.', app_index, PSO_unallocated_task_num_all, 'g-')
xlabel('application number');
ylabel('number of unallocated task');
legend('TP-CAA-HH', 'G-CAA-HH', 'BPSO-CAA-H');

figure(4)
plot(app_index, C_unallocated_task_num, 'ro', app_index, EG_unallocated_task_num, 'bx', app_index, PSO_unallocated_task_num, 'gd')
xlabel('application number');
ylabel('number of unallocated task in an application');
legend('TP-CAA-HH', 'G-CAA-HH', 'BPSO-CAA-H');

figure(5)
plot(app_index, C_energy_consumption_per_task, 'r--', app_index, EG_energy_consumption_per_task, 'b-.')
xlabel('application number');
ylabel('energy consumption per task in an application');
legend('TP-CAA-HH', 'G-CAA-HH');

figure(6)
plot(app_index, C_live_node_num, 'r--', app_index, EG_live_node_num, 'b-.', app_index, PSO_live_node_num, 'g-')
xlabel('application number');
ylabel('live node number');
legend('TP-CAA-HH', 'G-CAA-HH', 'BPSO-CAA-H');

end