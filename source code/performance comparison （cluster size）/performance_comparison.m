function performance_comparison
%
% Record of revisions:
%   Data               Programmer            Description of language
%  ======            =============          =========================
% 18/09/2017          Xiang Yin                    Original
% 

mu = 0.2;
num_app = 200;          % number of application
num_task = 5;          % number of tasks in an application

radius = 2000;             % the radius of the circular area

% the number of average
average_count = 10;

% number of resource type in the system
resource_type = 6;

% Coordinates of Sink node which is placed on origin
sink_x = 0;
sink_y = 0;

% Number of nodes
%num_node = 400;
%num_node_in_cluster = 20;

% Number of cluster
% num_cluster = num_node / num_node_in_cluster;  % each cluster contains 10 sensor nodes on average
num_cluster = 15;

num_node_in_cluster_min = 16;
num_node_in_cluster_max = 30

w1 = 0.3;  % W1-W3 为fitness权重系数
w2 = 0.3;
maxgen = 30;
swarm_size = 10;
c1 = 2;
c2 = 2;
% Wmax = 0.9;
% Wmin = 0.4;
w = rand(1) * (0.4 - 0.9) + 0.9;

E_elec = 50 * 0.000000001;     % electroincs energy, which is the same for transimtter and receiver
Efs = 10 * 0.000000000001;      % amplifier energy for free space channel
L = randi([40000, 60000]);       % bits of information for completing a sub-task

task_resource = zeros(num_task, resource_type);
for i = 1 : num_task
    task_resource_type = randi([2, resource_type]);  % for a task, it demands a specific number of resources which is in the range [2, application_resource_type]
    temp = randperm(resource_type);
    task_resource(i, temp(1 : task_resource_type)) = 1;
end

num_node_in_cluster = num_node_in_cluster_min : 2 : num_node_in_cluster_max;

C_E_consume_node = zeros(average_count, length(num_node_in_cluster));
C_residual_E_node = zeros(average_count, length(num_node_in_cluster));
C_live_node_node = zeros(average_count, length(num_node_in_cluster));
C_unallocated_task_node = zeros(average_count, length(num_node_in_cluster));

EG_E_consume_node = zeros(average_count, length(num_node_in_cluster));
EG_residual_E_node = zeros(average_count, length(num_node_in_cluster));
EG_live_node_node = zeros(average_count, length(num_node_in_cluster));
EG_unallocated_task_node = zeros(average_count, length(num_node_in_cluster));

PSO_E_consume_node = zeros(average_count, length(num_node_in_cluster));
PSO_residual_E_node = zeros(average_count, length(num_node_in_cluster));
PSO_live_node_node = zeros(average_count, length(num_node_in_cluster));
PSO_unallocated_task_node = zeros(average_count, length(num_node_in_cluster));

for count = 1 : average_count
     num_node_base = num_node_in_cluster_min * num_cluster;
        node_cluster_base = zeros(1, num_node_base);
        for i = 1 : num_node_base
            rand_radius = radius * rand;
            temp = mod(i, num_cluster);
            if temp == 0
                temp = num_cluster;
            end
            node_position_x_base(i) = rand_radius * cos((temp - 1) * (2 * pi / num_cluster) + (2 * pi / num_cluster) * rand);
            node_position_y_base(i) = rand_radius * sin((temp - 1) * (2 * pi / num_cluster) + (2 * pi / num_cluster) * rand);
            node_cluster_base(i) = temp;         % the cluster node i is located
        end
        % rectangle('Position',[-radius,-radius,2 * radius,2 * radius],'Curvature',[1,1]),axis equal;
        % hold on;
        % plot(node_position_x, node_position_y, 'bo');
        
       
        node_resource_base = randi([1, resource_type], 1, num_node_base);   % assign each node a specific resource type
        
        cluster_head = zeros(1, num_cluster);
        for i = 1 : num_cluster
            cluster_head(i) = i;
        end
        
    for t = 1 : length(num_node_in_cluster)
        num_node_new = num_node_in_cluster(t) * num_cluster - num_node_base;
        num_node = num_node_base + num_node_new;
        % Initialize coordinate of each sensor node, here, the sensor nodes are
        % uniformly distributed in the area
        node_cluster = zeros(1, num_node_new);
        for i = 1 : num_node_new
            rand_radius = radius * rand;
            temp = mod(i, num_cluster);
            if temp == 0
                temp = num_cluster;
            end
            node_position_x_new(i) = rand_radius * cos((temp - 1) * (2 * pi / num_cluster) + (2 * pi / num_cluster) * rand);
            node_position_y_new(i) = rand_radius * sin((temp - 1) * (2 * pi / num_cluster) + (2 * pi / num_cluster) * rand);
            node_cluster_new(i) = temp;         % the cluster node i is located
        end
        if num_node_new == 0
            node_position_x_new = [];
            node_position_y_new = [];
            node_cluster_new = [];
        end
        node_position_x = [node_position_x_base, node_position_x_new];
        node_position_y = [node_position_y_base, node_position_y_new];
        node_cluster = [node_cluster_base, node_cluster_new];
        % rectangle('Position',[-radius,-radius,2 * radius,2 * radius],'Curvature',[1,1]),axis equal;
        % hold on;
        % plot(node_position_x, node_position_y, 'bo');
        
        node_resource_new = randi([1, resource_type], 1, num_node_new);   % assign each node a specific resource type
        node_resource = [node_resource_base, node_resource_new];
        
        ini_E = 1 * ones(1, num_node);             % initial energy of each node
        
        dist_sink_node = sqrt(power(node_position_x, 2) + power(node_position_y, 2));
        
        dist_sink_head = zeros(1, num_cluster);
        for i = 1 : num_cluster
            dist_sink_head(i) = sqrt((node_position_x(cluster_head(i))) ^ 2 + (node_position_y(cluster_head(i))) ^ 2);
        end
        
        dist_head_node = zeros(num_cluster, num_node_in_cluster(t));
        for i = 1 : num_node
            cluster_index = node_cluster(i);
            index = fix(i / num_cluster) + 1;      % when i=num_cluster, a singular value appears, it is handled below
            dist_head_node(cluster_index, index) = sqrt((node_position_x(cluster_index) - node_position_x(i)) ^ 2 + (node_position_y(cluster_index) - node_position_y(i)) ^ 2);
        end
        for i = 1 : num_node_in_cluster(t)
            dist_head_node(num_cluster, i) = dist_head_node(num_cluster, i + 1);
        end
        dist_head_node(:, num_node_in_cluster(t) + 1) = [];
        
        %times = 0.5 * ones(1, num_node);    % initialize the times a node has been used
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
            
            C_task_allocation_flag = 1;
            EG_task_allocation_flag = 1;
            

            C_task_resource = task_resource;
            EG_task_resource = task_resource;
            PSO_task_resource = task_resource;
            
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
            for i = num_task : -1 : 1
                if all(C_task_cluster(i, :) == 0)
                    C_task_resource(i, :) = [];                 % delete the ith task from current application
                    C_task_cluster(i, :) = [];
                end
            end
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
            end
            if C_task_allocation_flag == 0
                std_C_residual_E(m) = std_C_residual_E(m - 1);
                C_live_node_num(m) = C_live_node_num(m - 1);
                C_unallocated_task_num(m) = 0;
                C_WSN_E_consume_app(m) = 0;
                C_unallocated_task_num(m) = num_task;
            end
            
            if EG_task_allocation_flag
                [EG_WSN_E_consume, EG_residual_E, EG_allocated_task_cluster, unallocated_task_num, EG_task_allocated_node] = energy_greedy(radius, EG_num_task, num_cluster, node_cluster, EG_residual_E, L, E_elec, Efs, EG_node_resource, EG_task_resource, EG_task_cluster, dist_head_node)   % assign tasks by energy greedy algorithm
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
                EG_unallocated_task_num(m) = EG_unallocated_task_num(m) + unallocated_task_num;
                
                EG_allocated_task_cluster_all = [EG_allocated_task_cluster_all, EG_allocated_task_cluster];
            end
            if EG_task_allocation_flag == 0
                EG_WSN_E_consume_app(m) = 0;
                std_EG_residual_E(m) = std_EG_residual_E(m - 1);
                EG_live_node_num(m) = EG_live_node_num(m - 1);
                EG_unallocated_task_num(m) = num_task;
            end
            
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
        
        C_E_consume_node(count, t) = sum(C_WSN_E_consume_app);
        C_residual_E_node(count, t) = std_C_residual_E(m);
        C_live_node_node(count, t) = C_live_node_num(m);
        C_unallocated_task_node(count, t) = sum(C_unallocated_task_num);
        
        EG_E_consume_node(count, t) = sum(EG_WSN_E_consume_app);
        EG_residual_E_node(count, t) = std_EG_residual_E(m);
        EG_live_node_node(count, t) = EG_live_node_num(m);
        EG_unallocated_task_node(count, t) = sum(EG_unallocated_task_num);
        
        PSO_E_consume_node(count, t) = sum(PSO_WSN_E_consume_app);
        PSO_residual_E_node(count, t) = std_PSO_residual_E(m);
        PSO_live_node_node(count, t) = PSO_live_node_num(m);
        PSO_unallocated_task_node(count, t) = sum(PSO_unallocated_task_num);
        
        num_node_base = num_node;
        node_position_x_base = node_position_x;
        node_position_y_base = node_position_y;
        node_cluster_base = node_cluster;
        node_resource_base = node_resource;
    end
end

average_C_E_consume_node = sum(C_E_consume_node) ./ average_count;
average_C_residual_E_node = sum(C_residual_E_node) ./ average_count;
average_C_live_node_node = sum(C_live_node_node) ./ average_count;
average_C_unallocated_task_node = sum(C_unallocated_task_node) ./ average_count;
average_C_energy_consumption_per_task = average_C_E_consume_node ./ (num_app * num_task - average_C_unallocated_task_node);

average_EG_E_consume_node = sum(EG_E_consume_node) ./ average_count;
average_EG_residual_E_node = sum(EG_residual_E_node) ./ average_count;
average_EG_live_node_node = sum(EG_live_node_node) ./ average_count;
average_EG_unallocated_task_node = sum(EG_unallocated_task_node) ./ average_count;
average_EG_energy_consumption_per_task = average_EG_E_consume_node ./ (num_app * num_task - average_EG_unallocated_task_node);

average_PSO_E_consume_node = sum(PSO_E_consume_node) ./ average_count;
average_PSO_residual_E_node = sum(PSO_residual_E_node) ./ average_count;
average_PSO_live_node_node = sum(PSO_live_node_node) ./ average_count;
average_PSO_unallocated_task_node = sum(PSO_unallocated_task_node) ./ average_count;
average_PSO_energy_consumption_per_task = average_PSO_E_consume_node ./ (num_app * num_task - average_PSO_unallocated_task_node);

app_index = num_node_in_cluster_min : 2 : num_node_in_cluster_max;
figure(3);
plot(app_index, average_C_E_consume_node, 'r--', app_index, average_EG_E_consume_node, 'b-.', app_index, average_PSO_E_consume_node, 'g-');
xlabel('cluster size');
ylabel('energy consumption');
legend('TP-CAA-HH', 'G-CAA-HH', 'BPSO-CAA-H');

figure(5)
plot(app_index, average_C_unallocated_task_node, 'r--', app_index, average_EG_unallocated_task_node, 'b-.', app_index, average_PSO_unallocated_task_node, 'g-')
xlabel('cluster size');
ylabel('number of unallocated task');
legend('TP-CAA-HH', 'G-CAA-HH', 'BPSO-CAA-H');

figure(6)
plot(app_index, average_C_energy_consumption_per_task, 'r--', app_index, average_EG_energy_consumption_per_task, 'b-.', app_index, average_PSO_energy_consumption_per_task, 'g-')
xlabel('cluster size');
ylabel('energy consumption per task');
legend('TP-CAA-HH', 'G-CAA-HH', 'BPSO-CAA-H');

figure(8)
plot(app_index, average_C_live_node_node, 'r--', app_index, average_EG_live_node_node, 'b-.', app_index, average_PSO_live_node_node, 'g-')
xlabel('cluster size');
ylabel('live node number');
legend('TP-CAA-HH', 'G-CAA-HH', 'BPSO-CAA-H');

end