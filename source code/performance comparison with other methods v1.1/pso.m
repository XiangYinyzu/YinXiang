function [PSO_WSN_E_consume, PSO_task_allocated_node_result, unallocated_task_num] = pso(swarm_size, maxgen, w, w1, w2, c1, c2, L, E_elec, Efs, PSO_num_task, resource_type, dist_sink_node, T_cell, node_resource, task_resource, R_type)
%
% This function allocates tasks by PSO algorithm. 
%
% Record of revisions:
%   Data               Programmer            Description of language
%  ======            =============          =========================
% 21/12/2017          Xiang Yin                    Original
% 

unallocated_task_num = 0;
PSO_allocated_task_node = zeros(PSO_num_task, resource_type);
PSO_allocated_task_node_cell = cell(1, swarm_size);
for i = 1 : swarm_size
    PSO_node_resource = node_resource;
    for j = 1 : PSO_num_task
        for k = 1 : resource_type
            if task_resource(j, k) == 1
                PSO_node_index = find(PSO_node_resource == k);
                if length(PSO_node_index) == 0                       % if a task can not be allocated, this case occurs when some nodes have already been assigned with tasks, they become unavailable
                    unallocated_task_num = unallocated_task_num + 1;
                    %keyboard;
                    break;
                end
                temp = randi([1, length(PSO_node_index)]);
                PSO_allocated_task_node(j, k) = PSO_node_index(temp);         % assign task i to the first cluster that can complete it, this is an initial allocation
                PSO_node_resource(PSO_node_index(temp)) = 0;
            end  
        end
    end
    PSO_allocated_task_node_cell{i} = PSO_allocated_task_node;
    PSO_allocated_task_node = zeros(PSO_num_task, resource_type);
end

%��ʼ��
X = cell(1, swarm_size);
V = cell(1, swarm_size);
for i = 1 : swarm_size
    X{i} = PSO_allocated_task_node_cell{i};
    V{i} = v(PSO_num_task, resource_type); %�ٶ�
    %ÿ��λ�þ�����ٶȾ���һһ��Ӧ
end

%����fitnessֵ
%initial_fitness = fitness(Time, PSO_num_task_resource_type, PSO_allocated_task_node_cell, dist, task, num_node, P, L, PS, PC, X, B, Swarm_size);
initial_fitness = fitness(swarm_size, w1, w2, L, E_elec, Efs, PSO_num_task, dist_sink_node, PSO_allocated_task_node_cell);


%%����Ⱥ�Ż��㷨    
[best_fitness, best_index] = min(initial_fitness);%�ҳ���Ӧ��ֵ��С�ļ���λ��
Gbest = X{best_index};       %ȫ�����
Pbest = X;            %�������
fitness_Pbest = initial_fitness;     %���������Ӧ��ֵ
fitness_Gbest = best_fitness; %ȫ�������Ӧ��ֵ

%% ����Ѱ��
for i = 1 : maxgen
    i = i;
    fitness_value_updated = zeros(1, swarm_size);
    for j = 1 : swarm_size
        %�ٶȸ���
        V{j} = w * V{j} + c1 * rand * (Pbest{j} - X{j}) + c2 * rand * (Gbest - X{j});
        X{j} = X{j} + V{j};
        X{j} = round(X{j});%ȡ��
       
       for m = 1 : PSO_num_task
           for n = 1 : length(find(task_resource(m, : ) ~= 0))         %resource_type
               flag = ismember(X{j}(m, n), R_type{T_cell{m}(n)}); %�жϸ��º���Դ�����Ƿ񻹷���
               if flag == 0
                  X{j}(m, n) = Gbest(m, n);
               end
           end
       end
        
    %���º��������Ӧ��ֵ�ļ���
    %fitness_value_updated = fitness(Time, PSO_num_task_resource_type, node_allocation, dist, task, Num_node, P, L, PS, PC, X, B, swarm_size);
    fitness_value_updated = fitness(swarm_size, w1, w2, L, E_elec, Efs, PSO_num_task, dist_sink_node, PSO_allocated_task_node_cell);
    end        
  
    for k = 1 : swarm_size     
        if fitness_value_updated(k) < fitness_Pbest(k)
           Pbest{k} = X{k};
           fitness_Pbest(k) = fitness_value_updated(k);
        end
        if fitness_value_updated(k) < fitness_Gbest
           Gbest = X{k};
           fitness_Gbest = fitness_value_updated(k);
        end
    end
end

PSO_E_consumption_task_best = zeros(1, PSO_num_task);
for i = 1 : PSO_num_task
    PSO_task_allocated_node_best_index = find(Gbest(i, :) ~= 0);        % find out the nodes of allocation in best solution
    for j = 1 : length(PSO_task_allocated_node_best_index)
        PSO_E_consumption_task_best(i) =  PSO_E_consumption_task_best(i) + (L * E_elec + L * Efs * dist_sink_node(Gbest(i, PSO_task_allocated_node_best_index(j)))^2);
        %PSO_E_consumption_task_best(i) =  PSO_E_consumption_task_best(i) + (L * E_elec + L * Efs * dist_sink_node(PSO_task_allocated_node_best_index(j))^2);
    end
end
PSO_WSN_E_consume = sum(PSO_E_consumption_task_best);    % the energy consumption of the task allocatin result
PSO_task_allocated_node_result = Gbest;          % the final task allocatin result

end