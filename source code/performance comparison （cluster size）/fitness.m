function fitness_value = fitness(swarm_size, w1, w2, L, E_elec, Efs, num_task, dist_sink_node, PSO_allocated_task_node_cell)
%
% Record of revisions:
%   Data               Programmer            Description of language
%  ======            =============          =========================
% 21/12/2017          Xiang Yin                    Original
% 

PSO_STD_E_consumption = zeros(1, num_task);
fitness_value = zeros(1, swarm_size);
for i = 1 : swarm_size
    PSO_task_allocated_node_index = PSO_allocated_task_node_cell{i};        
    PSO_task_allocated_node_index(PSO_task_allocated_node_index == 0) = [];    
    PSO_E_consumption = zeros(1, length(PSO_task_allocated_node_index));
    for j = 1 : length(PSO_task_allocated_node_index)
        PSO_E_consumption(j) = (L * E_elec + L * Efs * dist_sink_node(PSO_task_allocated_node_index(j))^2);
    end
    PSO_E_consume_particle = sum(PSO_E_consumption);     % energy consumption for a particle
    PSO_STD_E_consumption_particle = sum(PSO_STD_E_consumption);
    fitness_value(i) = w1 * PSO_E_consume_particle + w2 * PSO_STD_E_consumption_particle;
end

end
