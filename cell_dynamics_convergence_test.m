disp('busy');tic;close all; clear all;%profile on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_iterations = 2;
simulation_name = '';

grid_size = [7,7];
max_no_cells = 500;

delta_t = 50;
viscosity = 1;

%%%%%%%%%%%%%%%%%%%%%%%%% Initial configuration parameters %%%%%%%%%%%%%%%%%%%%%%%%%%

anneal_initial_configuration_logical = false;

compile_mex_functions_logical = false;
configuration_noise = 0.5;
file_to_load = 'convergence_test';

load_from_file_logical = true;
original_cells = 1:prod(grid_size);
% original_cells = [13 19 22 31 56 92 94 100];

% this doesn't seem to be working at the moment.
original_cells_colour = 'y';

if prod(grid_size) > max_no_cells
    error('Initial number of cells is greater than the maximum number allowed');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update position parameters %%%%%%%%%%%%%%%%%%%%%%%%%%

update_positions_logical = true;
update_positions_start = 0;

% 1 - forward euler 2 - 4th order runge-kutta
ode_solver_type = 1;

target_area_factor = 1.0;

% For ~50 cells

% initial_force_constant_magnitudes.area = 1e-0;
% initial_force_constant_magnitudes.deformation = 1e-4;
% initial_force_constant_magnitudes.elongation = 5e-7;
% initial_force_constant_magnitudes.perimeter = 1e-4;
% initial_force_constant_magnitudes.tension = 2e-4;

% For ~250 cells

initial_force_constant_magnitudes.area = 1e1;
initial_force_constant_magnitudes.deformation = 1e-4;
initial_force_constant_magnitudes.elongation = 5e-8;
initial_force_constant_magnitudes.perimeter = 1e-4;
initial_force_constant_magnitudes.tension = 2e-4;

boundary_force_constants.deformation = 1e-3;
boundary_force_constants.edge = 1e-1;

% initial_force_constant_magnitudes.area = 0;
% initial_force_constant_magnitudes.deformation = 0;
% initial_force_constant_magnitudes.elongation = 0;
% initial_force_constant_magnitudes.perimeter = 0;
% initial_force_constant_magnitudes.tension = 0;

% boundary_force_constants.deformation = 0;
% boundary_force_constants.edge = 0;

tension_anisotropy_factor = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Vertex rearrangement parameters %%%%%%%%%%%%%%%%%%%%%%%%%

T1_swaps_logical = true;
T1_swaps_start = 0;
T1_probability = 1.0;

threshold_T1_swaps_factor = 0.1;

% protection_time = 100;
protection_time = 0;

%%%%%%%%%%%%%%%%%%%%%%%%% Cell growth and mitosis parameters %%%%%%%%%%%%%%%%%%%%%%%%

cell_growth_logical = false;
cell_growth_Dpp_dependent = false;
mitosis_logical = false;

% growth speeds of medial (1) and lateral cells (2)
average_cell_growth_speeds(1) = 0.0005;
average_cell_growth_speeds(2) = 2*average_cell_growth_speeds(1);
% no_growth_time = 5000;
no_growth_time = 100;
% average_cell_growth_speeds = [5e-7 1e-6];
% medial_lateral_threshold_factor = 0.5;
medial_lateral_threshold_factor = 100;
% medial_lateral_threshold_factor = 0.0;

mitosis_start = 0;

% mitosis_angles_type can either be 'uniform', or a 1x2 vector specifying
% the mean and variance of a normal distribution
mitosis_angles_type = 'uniform';
% mitosis_angles_type = [0 0];

% mitosis_dependence can currently be either 'volume' or 'area'
mitosis_dependence = 'volume';

% determines whether mitosis takes place at a set volume (a certain
% fraction of the target volume) or stochastically. couldn't we just have a
% mitosis threshold and a probability? then could set this to one if we
% don't want it to be stochastic.
mitosis_random_logical = 0;

% determines a minimum threshold at which mitosis can occur. it is a
% fraction of the target parameter. it should be set to a number less than
% 1 if growth is logistic and mitosis is volume-dependent as a cell will
% never actually reach its target volume.
mitosis_threshold = 0.95;

target_volume_factor = 1.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cell death parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cell_death_logical = false;
cell_death_start = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEM parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FEM_solve_logical = false;

% degradation_constant = 0.0005;
degradation_constant = 0;
diffusion_speed = 0.00001;

% gradient type can be either 'x', 'y', or 'radial. this currently only
% works for the initial concentration. need to feed this into the
% FEM_solver to make it work for the source terms
gradient_type = 'radial';

% the initial concentration will be set to this value inside the source
% width. the exact nature depends on gradient_type above.
initial_concentration_magnitude = 0;

maximum_source_to_release = 0.12;

mesh_refinement_threshold_factor = 1.2;
% source_magnitude = 0.005;

source_magnitude = 0.001;
source_width = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Movie parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

movie_logical = 0;

% axis_values = [-1 1 -1 1];
axis_values = 'equal';
% axis_values_FEM = [-1 1 -1 1 -0.5 1.5];
axis_values_FEM = [-0.5 0.5 -0.5 0.5 -0.2 0.9];
% axis_values_FEM = 'equal';
extra_pause = 0.0;
% extra_pause = 0.1;
movie_name = simulation_name;
movie_period = 100;
movie_start = 0;
no_frames_for_statistical_plots = 100;
update_period_1 = 1e9;
update_period_2 = movie_period;
view_FEM_mesh = 0;
view_FEM_concentration = 1;
visualise_initial_configuration = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig_saves_logical = false;
fig_saves_name = simulation_name;

full_saves_logical = false;
full_saves_name = simulation_name;
full_saves_period = max(floor(max_iterations/3),1);

regular_tests_logical = true;

statistics_period = max(floor(max_iterations/1000),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compile mex functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

compile_mex_functions(compile_mex_functions_logical);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate initial configuration %%%%%%%%%%%%%%%%%%%%%%%%%%

[cell_growth_speeds_matrix,cells,FEM_elements,FEM_nodes,...
    refined_edge_matrix,vertices] = initial_configuration(anneal_initial_configuration_logical,...
    average_cell_growth_speeds,boundary_force_constants,configuration_noise,...
    gradient_type,grid_size,FEM_solve_logical,file_to_load,initial_concentration_magnitude,...
    initial_force_constant_magnitudes,load_from_file_logical,...
    max_no_cells,medial_lateral_threshold_factor,source_width);

% cell_growth_speeds(original_cells) = -0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save initial_save cells cell_growth_speeds_matrix FEM_elements FEM_nodes ...
    refined_edge_matrix vertices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialise simulations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

array_sizes = size(vertices.position,1);

[cells,fig_saves_location,figure_position,full_saves_location,mean_edge_length,...
    movie_location,stats,vertices] = initialise_simulation(array_sizes,...
    cells,fig_saves_logical,fig_saves_name,full_saves_logical,full_saves_name,...
    max_iterations,movie_logical,movie_name,statistics_period,target_area_factor,...
    max_no_cells,target_volume_factor,vertices);

if full_saves_logical
    
    full_saves_file_name = [full_saves_location,'initial_save'];
    save(full_saves_file_name);
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initial figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iteration = 0;

% axis_values = [min(node_positions(:,1))-0.5 max(node_positions(:,1))+0.5 ...
% 	min(node_positions(:,2))-0.5 max(node_positions(:,2))+0.5];

visualiser(cells,vertices,FEM_elements,FEM_nodes,axis_values,...
    axis_values_FEM,figure_position,false,iteration,...
    movie_logical,update_period_1,update_period_2,view_FEM_concentration,view_FEM_mesh,...
    visualise_initial_configuration)

if movie_logical == 2 && iteration > movie_start && ~rem(iteration,movie_period)
    
    M(1) = getframe(gcf);
    frame_counter = 2;
    
else
    
    frame_counter = 1;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vertices.previous_positions = repmat(vertices.position,1,10);
cells.previous_vertices = repmat({cells.vertices},1,10);

time = 0;
mitosis_counter = 0;
statistics_counter = 0;

no_mitosis_since_last_statistics = 0;
no_T1_swaps_since_last_statistics = 0;
no_deaths_since_last_statistics = 0;

total_source_released = 0;

while true
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    time = time + delta_t;
    iteration = iteration + 1;
    
    statistics_this_iteration_logical = false;
    
    if ~rem(iteration,statistics_period)
        statistics_this_iteration_logical = true;
        statistics_counter = statistics_counter+1;
    end
    
    vertices.previous_positions = [vertices.position vertices.previous_positions(:,1:18)];
    cells.previous_vertices = {cells.vertices,cells.previous_vertices{1:9}};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cells.target_area(strcmp(cells.state,'apoptotic')) = 0.99*cells.target_area(strcmp(cells.state,'apoptotic'));
    
    %     cells.force_constants.area(strcmp(cells.state,'apoptotic')) = 0.99*force_constants.area(cells.original_logical);
    %     cells.force_constants.deformation(strcmp(cells.state,'apoptotic')) = 0.99*force_constants.deformation(cells.original_logical);
    cells.force_constants.elongation(strcmp(cells.state,'apoptotic')) =...
        0.99^2*cells.force_constants.elongation(strcmp(cells.state,'apoptotic'));
    %     cells.force_constants.elongation(strcmp(cells.state,'apoptotic')) = 0;
    cells.force_constants.perimeter(strcmp(cells.state,'apoptotic')) =...
        1/0.99*cells.force_constants.perimeter(strcmp(cells.state,'apoptotic'));
    cells.force_constants.tension(strcmp(cells.state,'apoptotic')) =...
        1/0.99*cells.force_constants.tension(strcmp(cells.state,'apoptotic'));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    mesh_refinement_threshold = mesh_refinement_threshold_factor*mean_edge_length;
    
    [cells.area,cells.perimeter,cells.edge_lengths,current_cell_area_stats,...
        current_cell_perimeter_stats,current_shape_index_stats,...
        current_edge_length_stats,cells.shape_indices,long_edges] =...
        CalculateCellAreas(cells.vertices,vertices.position,...
        mesh_refinement_threshold);
    
    long_edges = long_edges(long_edges(:,1)>0,:);
    
    mean_cell_area = current_cell_area_stats(1);
    mean_edge_length = current_edge_length_stats(1);
    
    cells.state(cells.target_area-mean_cell_area>-1e-9&strcmp(cells.state,'baby')) = {'normal'};
    cells.target_area(strcmp(cells.state,'baby')) =...
        cells.target_area(strcmp(cells.state,'baby')) + 1/1000*mean_cell_area;
    
    cells = cell_volume_growth(array_sizes,cell_growth_logical,...
        cell_growth_Dpp_dependent,cells,delta_t,FEM_nodes.concentration,...
        no_growth_time,time);
    
    if statistics_this_iteration_logical
        
        stats.cell_area(statistics_counter,:) = current_cell_area_stats;
        stats.cell_perimeter(statistics_counter,:) = current_cell_perimeter_stats;
        stats.shape_index(statistics_counter,:) = current_shape_index_stats;
        stats.edge_length(statistics_counter,:) = current_edge_length_stats;
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEM Solve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if FEM_solve_logical
        
        FEM_nodes.position = UpdateFEMNodePositions(cells.vertices,vertices.position,cells.area,FEM_nodes.edge);
        
        [FEM_nodes.concentration,current_total_Dpp_stats,source_this_iteration,triangle_quality] = FEM_solver(...
            degradation_constant,delta_t,diffusion_speed,FEM_nodes.concentration,FEM_elements.nodes,...
            FEM_nodes.position,FEM_nodes.previous_position,source_magnitude,...
            source_width);
        
        FEM_nodes.previous_position = FEM_nodes.position;
        
        total_source_released = total_source_released+source_this_iteration;
        
        if statistics_this_iteration_logical
            
            stats.source_term(statistics_counter,:) = [total_source_released source_this_iteration];
            stats.total_Dpp(statistics_counter,:) = current_total_Dpp_stats;
            stats.triangle_quality(statistics_counter,:) = ...
                [mean(triangle_quality) max(triangle_quality) ...
                min(triangle_quality) std(triangle_quality)];
            
        end
        
        if total_source_released > maximum_source_to_release
            
            source_magnitude = 0;
            
        end
        
        [cells,FEM_elements,FEM_nodes,refined_edge_matrix] = remove_FEM_nodes_from_short_edges(cells,...
            FEM_elements,FEM_nodes,vertices,mesh_refinement_threshold,refined_edge_matrix);
        
        [cells,FEM_elements,FEM_nodes,refined_edge_matrix] = add_FEM_nodes_to_long_edges(...
            cells,FEM_elements,FEM_nodes,long_edges,refined_edge_matrix,vertices);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update Positions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if update_positions_logical && iteration > update_positions_start
        
        boundary_element =...
            CreateBoundaryElement(vertices.cells,cells.vertices,vertices.no_cells);
        
		[vertices,angle_deviations,current_total_force_stats,vertex_movements] =...
			ode_solver(cells,vertices,boundary_element,boundary_force_constants,cell_growth_logical,...
			cell_growth_Dpp_dependent,delta_t,FEM_nodes.concentration,mean_edge_length,...
			no_growth_time,ode_solver_type,tension_anisotropy_factor,time,viscosity);
        
        if test_for_nans(FEM_nodes.concentration,vertices.position,iteration,regular_tests_logical)
            break;
		end
        
		% should this be done inside ode_solver? might make things neater,
		% but possible less transparent what is going on
        if statistics_this_iteration_logical
            
            stats.total_boundary_nodes_stats(statistics_counter) =...
                length(boundary_element);
            
            stats.total_area_force(statistics_counter) = current_total_force_stats.area;
            stats.total_boundary_deformation_force(statistics_counter) =...
                current_total_force_stats.boundary_deformation;
            stats.total_boundary_edge_force(statistics_counter) =...
                current_total_force_stats.boundary_edge;
            stats.total_deformation_force(statistics_counter) =...
                current_total_force_stats.deformation;
            stats.total_elongation_force(statistics_counter) =...
                current_total_force_stats.elongation;
            stats.total_perimeter_force(statistics_counter) =...
                current_total_force_stats.perimeter;
            stats.total_tension_force(statistics_counter) = current_total_force_stats.tension;
            stats.total_force(statistics_counter) = current_total_force_stats.total;
            
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~rem(iteration,statistics_period)
        
        vertex_movement_magnitudes = sqrt(sum(vertex_movements.^2,2));
        
        [stats.angle_deviation(statistics_counter,:),stats.cell_height(...
            statistics_counter,:),stats.cell_height_to_area(statistics_counter,:),...
            stats.cell_volume(statistics_counter,:),stats.cells_per_node(...
            statistics_counter,:),stats.mesh_peclet_number(statistics_counter),...
            stats.node_movement(statistics_counter,:),stats.nodes_per_cell(...
            statistics_counter,:),stats.rosette(statistics_counter,:),...
            stats.total_no_cells(statistics_counter),stats.total_no_nodes(...
            statistics_counter),stats.xy_ratio(statistics_counter)] =...
            generate_statistics(angle_deviations,cells.area,cells.vertices,vertices.no_cells,...
            cells.volume,delta_t,diffusion_speed,FEM_solve_logical,mean_edge_length,vertex_movement_magnitudes,...
            vertices.position);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mitosis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    medial_lateral_threshold =...
        medial_lateral_threshold_factor*max(vertices.position(:,1));
    
    if mitosis_logical && iteration > mitosis_start
        
        
        [cells,vertices,...
            cell_growth_speeds_matrix,FEM_elements,FEM_nodes,...
            mitosis_location,refined_edge_matrix] = mitosis(cells,vertices,...
            cell_growth_speeds_matrix,delta_t,FEM_elements,FEM_nodes,...
            FEM_solve_logical,medial_lateral_threshold,mitosis_angles_type,...
            mitosis_dependence,mitosis_random_logical,mitosis_threshold,refined_edge_matrix,time);
        
        
        if mitosis_location ~= 0
            
            no_mitosis_since_last_statistics = no_mitosis_since_last_statistics+1;
            mitosis_counter = mitosis_counter+1;
            stats.mitosis_locations(mitosis_counter,:) = mitosis_location;
            %             save mitosis_test_save
            
        end
        
        if statistics_this_iteration_logical
            
            stats.no_mitosis(statistics_counter) = no_mitosis_since_last_statistics;
            no_mitosis_since_last_statistics = 0;
            
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% T1 Swaps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if T1_swaps_logical && iteration > T1_swaps_start
        
        threshold_T1_swaps = threshold_T1_swaps_factor*mean_edge_length;
        
        [cells.vertices,vertices.position,cells.FEM_elements,vertices.cells,vertices.no_cells,FEM_nodes.concentration,...
            FEM_elements.nodes,no_T1_swaps_this_iteration,FEM_nodes.previous_position,...
            vertices.time_created,FEM_nodes.edge] = T1Swaps(cells.vertices,vertices.position,cells.FEM_elements,...
            vertices.cells,vertices.no_cells,FEM_nodes.concentration,FEM_elements.nodes,FEM_solve_logical,...
            FEM_nodes.previous_position,protection_time,T1_probability,...
            threshold_T1_swaps,time,vertices.time_created,FEM_nodes.edge);
        
        no_T1_swaps_since_last_statistics = no_T1_swaps_since_last_statistics +...
            no_T1_swaps_this_iteration;
        
        %         if(no_T1_swaps_this_iteration)>0
        %             error('T1 swap');
        %         end
        
        if ~rem(iteration,statistics_period)
            
            stats.no_T1_swaps(statistics_counter) = no_T1_swaps_since_last_statistics;
            no_T1_swaps_since_last_statistics = 0;
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if cell_death_logical && iteration > cell_death_start
        
        threshold_area = 0.2*current_cell_area_stats(1);
        
        if rand<0.001
            random_cell = ceil(rand*length(cells.vertices));
			while cells.boundary_logical(random_cell)
				random_cell = ceil(rand*length(cells.vertices));
			end
            cells.state{random_cell} = 'apoptotic';
        end
        
        [cells,FEM_elements,FEM_nodes,no_deaths_this_iteration,vertices] = cell_death(...
            cells,FEM_elements,FEM_nodes,protection_time,threshold_area,time,vertices);
        
        no_deaths_since_last_statistics = no_deaths_since_last_statistics +...
            no_deaths_this_iteration;
        
        if ~rem(iteration,statistics_period)
            
            stats.no_deaths(statistics_counter) = no_deaths_since_last_statistics;
            no_deaths_since_last_statistics = 0;
            
        end
        
        %     if no_deaths_this_iteration > 0
        %         disp(['cell death at iteration ',num2str(iteration)])
        %         save death_save
        %         break;
        %     end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if test_for_nans(FEM_nodes.concentration,vertices.position,iteration,regular_tests_logical)
        break;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Movie %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 	axis_values = [min(node_positions(:,1))-0.5 max(node_positions(:,1))+0.5 ...
    % 		min(node_positions(:,2))-0.5 max(node_positions(:,2))+0.5];
    
    visualiser(cells,vertices,FEM_elements,FEM_nodes,axis_values,...
        axis_values_FEM,figure_position,false,iteration,...
        movie_logical,update_period_1,update_period_2,view_FEM_concentration,view_FEM_mesh,...
        visualise_initial_configuration);
    
    if movie_logical == 2 && iteration > movie_start && ~rem(iteration,movie_period)
        
        M(frame_counter) = getframe(gcf);
        frame_counter = frame_counter+1;
        
        if frame_counter > 10
            mpgwrite(M,jet,[movie_location,'update.mpg']);
            system(['cat ',movie_location,movie_name,'.mpg ',movie_location,...
                'update.mpg > ',movie_location,'temp.mpg']);
            system(['mv ',movie_location,'temp.mpg ',movie_location,...
                movie_name,'.mpg']);
            clear M
            frame_counter = 1;
        end
    end
    
    pause(extra_pause);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Full saves %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if full_saves_logical && ~rem(iteration,full_saves_period)
        
        full_saves_file_name = [full_saves_location,'iteration_',num2str(iteration)];
        save(full_saves_file_name);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if iteration == max_iterations
        disp('Maximum number of iterations reached');
        break;
    end
    
    if length(cells) == max_no_cells
        disp('Maximum number of cells reached')
        break;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save cells cell_growth_speeds_matrix FEM_elements FEM_nodes ...
    refined_edge_matrix vertices

if full_saves_logical
    
    full_saves_file_name = [full_saves_location,'final_save'];
    save(full_saves_file_name);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Final plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

visualiser(cells,vertices,FEM_elements,FEM_nodes,axis_values,...
    axis_values_FEM,figure_position,true,iteration,...
    movie_logical,update_period_1,update_period_2,view_FEM_concentration,view_FEM_mesh,...
    visualise_initial_configuration);

if movie_logical == 2
    M(frame_counter) = getframe(gcf);
    mpgwrite(M,jet,[movie_location,'update.mpg']);
    system(['cat ',movie_location,movie_name,'.mpg ',movie_location,...
        'update.mpg > ',movie_location,'temp.mpg']);
    system(['mv ',movie_location,'temp.mpg ',movie_location,...
        movie_name,'.mpg']);
end

statistical_plots(delta_t,fig_saves_location,fig_saves_logical,figure_position,...
    movie_location,movie_logical,movie_name,no_frames_for_statistical_plots,...
    statistics_counter,stats,time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display simulation info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_taken = toc;

disp(['Number of iterations: ',int2str(iteration)])
disp(['Time taken: ',num2str(round(time_taken)),' seconds'])

% profile viewer