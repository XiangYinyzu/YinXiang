function y = v(task, task_resource_num)
Vmax = 2; 
Vmin = -2;

 for i = 1 : task
     for j = 1 : task_resource_num
         y(i,j) = (Vmax - Vmin) * rand + Vmin;
     end
 end
end