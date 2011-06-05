clear all; disp('busy');

volume(1) = 1;
target_volume = 2;

total_time = 1;
delta_t = 0.00001;

total_iterations = ceil(total_time/delta_t);

iteration = 1;

while iteration<total_iterations
    
    iteration = iteration +1;
    volume(iteration) = volume(iteration-1)*(1+delta_t*7*(1-volume(iteration-1)/target_volume));
    
end

% figure
time_range=linspace(0,total_time,total_iterations);
hold all
plot(time_range,volume)
    
    