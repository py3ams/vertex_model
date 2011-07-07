disp('busy');close all;clear all;tic;%profile on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

total_time = 1;

max_iterations = 1000;
no_refinements = 2;

% simulation_name = 'refinement_comparison/true_solution';
simulation_name = ['refinement_comparison/iterations_',num2str(max_iterations),...
   '_refinements_',num2str(no_refinements)];


grid_size = [10,10];
max_no_cells = 101;

delta_t = total_time/max_iterations;
viscosity = 0.01;

%%%%%%%%%%%%%%%%%%%%%%%%% Initial configuration parameters %%%%%%%%%%%%%%%%%%%%%%%%%%

anneal_initial_configuration_logical = false;

compile_mex_functions_logical = false;
configuration_noise = 0.5;
% can be 'square', 'random', or 'hexagonal'
configuration_type = 'hexagonal';

load_from_file_logical = true;
load_FEM_from_file_logical = false;
file_to_load = 'Saves/refinement_comparison/true_solution/initial_save';

% to set the colour of the original cells to be different in figures and
% movies, need to edit figure_loop.m. otherwise would have to pass a variable
% through the visualiser.
original_cells = 1:prod(grid_size);
% original_cells = [13 19 22 31 56 92 94 100];

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

initial_force_constant_magnitudes.area = 1e0;
initial_force_constant_magnitudes.deformation = 5e-4;
initial_force_constant_magnitudes.elongation = 1e-5;
initial_force_constant_magnitudes.perimeter = 1e-4;
initial_force_constant_magnitudes.tension = 2e-4;

% For ~250 cells

% initial_force_constant_magnitudes.area = 1e1;
% initial_force_constant_magnitudes.deformation = 1e-4;
% initial_force_constant_magnitudes.elongation = 1e-7;
% initial_force_constant_magnitudes.perimeter = 1e-4;
% initial_force_constant_magnitudes.tension = 2e-4;

boundary_force_constants.deformation = 5e-3;
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

T1_swaps_logical = false;
T1_swaps_start = 0;
T1_probability = 1.0;

threshold_T1_swaps_factor = 0.05;

% protection_time = 100;
protection_time = 0;

%%%%%%%%%%%%%%%%%%%%%%%%% Cell growth and mitosis parameters %%%%%%%%%%%%%%%%%%%%%%%%

cell_growth_logical = false;
cell_growth_start = 0;
cell_growth_concentration_dependent = false;
mitosis_logical = false;

% growth speeds of medial (1) and lateral cells (2)
average_cell_growth_speed(1) = 0.5;
average_cell_growth_speed(2) = 2*average_cell_growth_speed(1);

% n.b. this parameter is multiplied by the total internal chemical in each
% cell, i.e. cells.internal_chemical*cells.area. need to make sure the
% orders of magnitude are right
lambda = 5000;

% no_growth_time = 5000;
no_growth_time = 0;
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
% this is only used if mitosos_dependence is set to 'none';
mitosis_period = 0.2;

% determines whether mitosis takes place at a set volume (a certain
% fraction of the target volume) or stochastically. couldn't we just have a
% mitosis threshold and a probability? then could set this to one if we
% don't want it to be stochastic.
mitosis_random_logical = 1;

% determines a minimum threshold at which mitosis can occur. it is a
% fraction of the target parameter. it should be set to a number less than
% 1 if growth is logistic and mitosis is volume-dependent as a cell will
% never actually reach its target volume.
mitosis_threshold = 0.9;

% determines the fraction of the initial maximum cell volume that the
% target volume is set to. if less than one a load of divisions will likely
% occur at the start of a simulation.
target_volume_factor = 1.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cell death parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cell_death_logical = false;
cell_death_start = 0;

% sets the area threshold below which cells can die, as a fraction of the mean area.
% cells must be 3-sided for this to occur. most 3-sided cells will be well below
% this threshold anyway.
cell_death_area_threshold_factor = 0.2;
apoptosis_concentration_dependent = 0;
apoptosis_concentration_threshold = 0.00025;
% apoptosis_concentration_threshold = 1e6;
% apoptosis_no_put_pc = apoptosis no per unit time per cell 
apoptosis_no_put_pc_below_threshold = 0.2;
apoptosis_no_put_pc_above_threshold = 1e-3;
% apoptosis_baseline_prob_per_unit_time is for each cell, so the number of
% apoptotic cells will be this number multiplied by the number of cells.
% this only applied if apoptosis_concentration_depdendent is set to zero,
% even if FEM_solve_logical is false
apoptosis_no_baseline_put_pc = 0.025;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEM parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FEM_solve_logical = true;

% mesh_refinement_threshold_factor = 1.2;
mesh_refinement_threshold_factor = 10;
no_chemicals = 1;
chemical_to_view = 1;

% degradation_constant(1) = 0.0005;
degradation_constant = [0.00000 0.00002];
diffusion_speed = [0.1 0.00002];
% diffusion_speed(1) = 0;

% gradient type can be either 1 - in the x direction (with peak at x = 0), 2 -
% in the y direction with peak at y = 0, 3 - radial, or 4 - from the left-hand edge.
gradient_type = [4 3];

% the initial concentration will be set to this value inside the source
% width. the exact nature depends on gradient_type above.
initial_concentration_magnitude = [0.0 0.1];
initial_concentration_magnitude(1) = 0.1;

maximum_source_to_release = [0.1 0.1];
maximum_source_to_release(1) = 100;

% need to be careful with these numbers. the way that cell volume growth is
% set up at the moment requires concentration values to be of the order 1
% for them to have a suitable effect on growth. the trade-offs between
% source magnitude and degradation etc are therefore important.
source_magnitude = [0.001 0.002];
source_magnitude(1) = 0;
source_width = [0.15 0.1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Movie parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

movie_logical = 0;

axis_values = 0.8*[-1 1 -1 1];
% axis_values = [-1 2 -1.5 1.5];
% axis_values = 'equal';
% axis_values_FEM = [-1 1 -1 1 -0.5 1.5];
axis_values_FEM = [axis_values -0.05 0.15];
% axis_values_FEM = 'equal';
extra_pause = 0.0;
% extra_pause = 0.1;
include_statistical_plots_in_movie = false;
linewidth_cells = 5;
linewidth_elements = 1;
movie_name = simulation_name;
movie_start = 0;
no_frames_for_statistical_plots = 100;
update_period = 1;
view_FEM_mesh = 1;
view_FEM_concentration = 1;
view_initial_config = 0;
view_iteration_number = 0;
view_number_cells = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig_saves_logical = false;
fig_saves_name = simulation_name;

full_saves_logical = true;
full_saves_name = simulation_name;
% full_saves_period = max(floor(max_iterations/3),1);
full_saves_period = 1;

regular_tests_logical = false;

statistics_period = max(floor(max_iterations/1000),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compile mex functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

compile_mex_functions(compile_mex_functions_logical);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate initial configuration %%%%%%%%%%%%%%%%%%%%%%%%%%

[cell_growth_speeds_matrix,cells,FEM_elements,FEM_nodes,...
	refined_edge_matrix,vertices] = initial_configuration(anneal_initial_configuration_logical,...
	average_cell_growth_speed,boundary_force_constants,configuration_noise,configuration_type,...
	gradient_type,grid_size,FEM_solve_logical,file_to_load,initial_concentration_magnitude,...
	initial_force_constant_magnitudes,load_FEM_from_file_logical,load_from_file_logical,...
	max_no_cells,medial_lateral_threshold_factor,no_chemicals,no_refinements,source_width);

% cell_growth_speeds(original_cells) = -0.1;

save('initial_save')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialise simulations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

array_sizes = size(vertices.position,1);

[cells,fig_saves_location,figure_position,full_saves_location,mean_cell_area,mean_edge_length,...
	movie_location,stats,vertices] = initialise_simulation(array_sizes,...
	cells,fig_saves_logical,fig_saves_name,full_saves_logical,full_saves_name,...
	max_iterations,movie_logical,movie_name,statistics_period,target_area_factor,...
	max_no_cells,target_volume_factor,vertices);

mesh_refinement_threshold = mesh_refinement_threshold_factor*mean_edge_length;

if full_saves_logical
	
	full_saves_file_name = [full_saves_location,'initial_save'];
	save(full_saves_file_name);
	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initial figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iteration = 0;

% axis_values = [min(node_positions(:,1))-0.5 max(node_positions(:,1))+0.5 ...
% 	min(node_positions(:,2))-0.5 max(node_positions(:,2))+0.5];

visualiser(cells,vertices,FEM_elements,FEM_nodes,axis_values,...
    axis_values_FEM,chemical_to_view,false,iteration,...
    linewidth_cells,linewidth_elements,movie_logical,update_period,...
    view_FEM_concentration,view_FEM_mesh,view_initial_config,...
    view_iteration_number,view_number_cells);

if movie_logical == 2 && iteration > movie_start && ~rem(iteration,update_period)
	
	M(1) = getframe(gcf);
	frame_counter = 2;
	
else
	
	frame_counter = 1;
	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vertices.previous_positions = repmat(vertices.position,1,10);
cells.previous_vertices = repmat({cells.vertices},1,10);

time = 0;
death_counter = 0;
mitosis_counter = 0;
stats.counter = 0;

stats.no_deaths_since_last_statistics = 0;
stats.no_mitosis_since_last_statistics = 0;
no_T1_swaps_since_last_statistics = 0;
no_deaths_since_last_statistics = 0;

dead_chemical = 0;
total_source_released = 0;
total_ingestion = 0;

baseline_target_volume = cells.target_volume(1);

time_taken_to_start_of_loop = toc;
tic;

while true
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	time = time + delta_t;
	iteration = iteration + 1;
	
	stats.this_iteration_logical = false;
	
	if ~rem(iteration,statistics_period)
		stats.this_iteration_logical = true;
		stats.counter = stats.counter+1;
	end
	
	vertices.previous_positions = [vertices.position vertices.previous_positions(:,1:18)];
	cells.previous_vertices = {cells.vertices,cells.previous_vertices{1:9}};
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mitosis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
	medial_lateral_threshold =...
		medial_lateral_threshold_factor*max(vertices.position(:,1));
	
	% need edge_lengths to do mitosis so can't do on first iteration
	if mitosis_logical && iteration > min(mitosis_start,2)
		
		[cells,vertices,cell_growth_speeds_matrix,FEM_elements,FEM_nodes,...
			mitosis_counter,refined_edge_matrix,stats] = mitosis(cells,vertices,...
			cell_growth_speeds_matrix,delta_t,FEM_elements,FEM_nodes,...
			FEM_solve_logical,medial_lateral_threshold,mitosis_angles_type,...
			mitosis_counter,mitosis_dependence,mitosis_period,mitosis_random_logical,...
			mitosis_threshold,refined_edge_matrix,stats,time);
		
		%             save mitosis_test_save
		
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% T1 Swaps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% 	test_refined_edge_matrix
	
	if T1_swaps_logical && iteration > T1_swaps_start
		
		threshold_T1_swaps = threshold_T1_swaps_factor*mean_edge_length;
		
		%         error('stop')
		
		[cells.vertices,vertices.position,cells.FEM_elements,vertices.cells,vertices.no_cells,FEM_nodes.concentration,...
			FEM_elements.nodes,no_T1_swaps_this_iteration,FEM_nodes.previous_position,...
			vertices.time_created,FEM_nodes.edge,refined_edge_matrix_edits] =...
			T1Swaps(cells.vertices,vertices.position,cells.FEM_elements,...
			vertices.cells,vertices.no_cells,FEM_nodes.concentration,FEM_elements.nodes,FEM_solve_logical,...
			FEM_nodes.previous_position,protection_time,T1_probability,...
			threshold_T1_swaps,time,vertices.time_created,FEM_nodes.edge);
		
		%          error('stop')
		
		%         if(no_T1_swaps_this_iteration)>0
		%             error('T1 swap');
		%         end
		
		if sum(sum(refined_edge_matrix_edits))
			
			refined_edge_matrix_edits =...
				refined_edge_matrix_edits(refined_edge_matrix_edits(:,1)>0,:);
			
			refined_edge_matrix(sub2ind(size(refined_edge_matrix),...
				refined_edge_matrix_edits(:,1),refined_edge_matrix_edits(:,2))) =...
				refined_edge_matrix_edits(:,3);
			
			refined_edge_matrix(sub2ind(size(refined_edge_matrix),...
				refined_edge_matrix_edits(:,2),refined_edge_matrix_edits(:,1))) =...
				refined_edge_matrix_edits(:,3);
			
		end
		
		% 		test_refined_edge_matrix
		
		no_T1_swaps_since_last_statistics = no_T1_swaps_since_last_statistics +...
			no_T1_swaps_this_iteration;
		
		if stats.this_iteration_logical
			
			stats.no_T1_swaps(stats.counter) = no_T1_swaps_since_last_statistics;
			no_T1_swaps_since_last_statistics = 0;
			
		end
		
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cell death %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if cell_death_logical && iteration > cell_death_start
		
		% this function chooses which cells become apoptotic, and alters
		% the force constants etc on cells that are apoptotic
		[cells,dead_chemical,stats] = set_apoptotic_cells(cells,...
			apoptosis_no_baseline_put_pc,apoptosis_concentration_dependent,...
			apoptosis_concentration_threshold,apoptosis_no_put_pc_above_threshold,...
			apoptosis_no_put_pc_below_threshold,dead_chemical,delta_t,...
			FEM_solve_logical,stats);
		
		cell_death_area_threshold = cell_death_area_threshold_factor*mean_cell_area;
		
        [cells,death_counter,FEM_elements,FEM_nodes,no_deaths_this_iteration,...
            refined_edge_matrix,stats,vertices] = cell_death(cells,death_counter,...
            FEM_elements,FEM_nodes,FEM_solve_logical,protection_time,refined_edge_matrix,...
            stats,cell_death_area_threshold,time,vertices);		
		
		%     if no_deaths_this_iteration > 0
		%         disp(['cell death at iteration ',num2str(iteration)])
		%         save death_save
		%         break;
		%     end
		
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cell volume growth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if cell_growth_logical && iteration > cell_growth_start
		
		% 	concentration_before = FEM_nodes.concentration;
		
        [cells,FEM_elements,FEM_nodes,stats,vertices] = cell_volume_growth(...
            cell_growth_concentration_dependent,cells,delta_t,FEM_elements,...
            FEM_nodes,lambda,no_growth_time,stats,time,vertices);
		
		% 	concentration_after = FEM_nodes.concentration;
		
		%     total_chemical_internalised = total_chemical_internalised+chemical_internalised_stats(5);
		
		% 		test_refined_edge_matrix
		
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update Positions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if update_positions_logical && iteration > update_positions_start
		
		% make sure these values are up-to-date for ode_solver. all other
		% cell properties are edited directly within functions such as
		% mitosis/T1 swaps. we could make sure these are edited to so this
		% line would not be necessary.
		[cells.area,cells.perimeter,cells.edge_lengths] =...
			CalculateCellAreas(cells.vertices,vertices.position);
		
		boundary_element =...
			CreateBoundaryElement(vertices.cells,cells.vertices,vertices.no_cells);
		
		[vertices,angle_deviations,stats,vertex_movements] =...
			ode_solver(cells,vertices,boundary_element,boundary_force_constants,cell_growth_logical,...
			cell_growth_concentration_dependent,delta_t,FEM_nodes.concentration,mean_edge_length,...
			no_growth_time,ode_solver_type,stats,tension_anisotropy_factor,time,viscosity);
		
		if test_for_nans(FEM_nodes.concentration,vertices.position,iteration,regular_tests_logical)
			break;
		end

	end
	
	%%%%%%%%%%%%%%%% Update cell variables and refined_edge_matrix %%%%%%%%%%%%%%%%%%
	
	% re-calculate cell area, perimeter, and edge lengths, for use next
	% time around the loop. can also find all the necessary values for
	% stats now as well. we do this before FEM_solver so we can use
	% long_edges to change the refined_edge_matrix
	[cells.area,cells.perimeter,cells.edge_lengths,current_cell_area_stats,...
		current_cell_perimeter_stats,current_shape_index_stats,...
		current_edge_length_stats,cells.shape_indices,long_edges] =...
		CalculateCellAreas(cells.vertices,vertices.position,...
		mesh_refinement_threshold);
	
	mean_cell_area = current_cell_area_stats(1);
	mean_edge_length = current_edge_length_stats(1);
	
	mesh_refinement_threshold = mesh_refinement_threshold_factor*mean_edge_length;
	
	long_edges = long_edges(long_edges(:,1)>0,:);
	
% 	% check if target areas of baby cells have reached the mean cell area. don't be
% 	% tempted to change this to checking the absolute difference between
% 	% cells.target_area and mean_cell_area as it is possible for the target
% 	% area to go slightly above the mean cell area by more than 1e-9 which
% 	% then won't be caught by the logical statement
% 	cells.state(cells.state==2&cells.target_area-mean_cell_area>-1e-9) = 1;
    cells.state(cells.state==2&abs(time-cells.time_of_last_division)>0.1) = 1;

% 	
% 	cells.target_area(cells.state==2) =...
% 		cells.target_area(cells.state==2) + 1/1000*mean_cell_area;
	
%     cells.target_area =...
%         mean_cell_area*(1+2*cells.internal_chemical_quantity/max(cells.internal_chemical_quantity));

%     cells.target_area = 100*cells.internal_chemical_quantity;

	if stats.this_iteration_logical
		
		stats.cell_area(stats.counter,:) = current_cell_area_stats;
		stats.cell_perimeter(stats.counter,:) = current_cell_perimeter_stats;
		stats.shape_index(stats.counter,:) = current_shape_index_stats;
		stats.edge_length(stats.counter,:) = current_edge_length_stats;
		
	end
	
	% this seems to be the most logical place for this to go
	if FEM_solve_logical
		
% 		[cells,FEM_elements,FEM_nodes,refined_edge_matrix] = remove_FEM_nodes_from_short_edges(cells,...
% 			FEM_elements,FEM_nodes,vertices,mesh_refinement_threshold,refined_edge_matrix);
		
		% 		test_refined_edge_matrix
		
% 		[cells,FEM_elements,FEM_nodes,refined_edge_matrix] = add_FEM_nodes_to_long_edges(...
% 			cells,FEM_elements,FEM_nodes,long_edges,refined_edge_matrix,vertices);
		
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEM Solve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if FEM_solve_logical
		
		FEM_nodes.position = UpdateFEMNodePositions(cells.vertices,...
			vertices.position,cells.area,FEM_nodes.edge);
		
        cells = set_source_and_ingestion_rates(cells,FEM_nodes,gradient_type,...
            no_chemicals,source_magnitude,source_width);
        
		[cells,FEM_nodes,M,source_magnitude,stats,total_ingestion,total_source_released] = ...
			FEM_solver(cells,degradation_constant,delta_t,diffusion_speed,FEM_elements,...
			FEM_nodes,gradient_type,maximum_source_to_release,no_chemicals,source_magnitude,...
			source_width,refined_edge_matrix,stats,total_ingestion,total_source_released);
		
		FEM_nodes.previous_position = FEM_nodes.position;
		        
%         cells.target_volume = baseline_target_volume*(1+1000*cells.internal_chemical_quantity);
%         cells.target_volume = 1e-3*cells.internal_chemical_quantity/max(cells.internal_chemical_quantity);

        
	end
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	stats = generate_statistics(cells,vertices,angle_deviations,delta_t,...
		diffusion_speed,FEM_nodes,FEM_solve_logical,mean_edge_length,...
		stats,vertex_movements);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if test_for_nans(FEM_nodes.concentration,vertices.position,iteration,regular_tests_logical)
		break;
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Movie %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% 	axis_values = [min(node_positions(:,1))-0.5 max(node_positions(:,1))+0.5 ...
	% 		min(node_positions(:,2))-0.5 max(node_positions(:,2))+0.5];
	
    visualiser(cells,vertices,FEM_elements,FEM_nodes,axis_values,...
        axis_values_FEM,chemical_to_view,false,iteration,...
        linewidth_cells,linewidth_elements,movie_logical,update_period,...
        view_FEM_concentration,view_FEM_mesh,view_initial_config,...
        view_iteration_number,view_number_cells);
	
	if movie_logical == 2 && iteration > movie_start && ~rem(iteration,update_period)
		
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
% 		eval(['save ',full_saves_file_name,' FEM_nodes FEM_elements M']);
        save(full_saves_file_name)
		
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

save('final_save');

time_taken_in_main_loop = toc;

if full_saves_logical
	
	full_saves_file_name = [full_saves_location,'final_save'];
	save(full_saves_file_name);
	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Final plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

visualiser(cells,vertices,FEM_elements,FEM_nodes,axis_values,...
    axis_values_FEM,chemical_to_view,true,iteration,...
    linewidth_cells,linewidth_elements,movie_logical,update_period,...
    view_FEM_concentration,view_FEM_mesh,view_initial_config,...
    view_iteration_number,view_number_cells);

if movie_logical == 2
	M(frame_counter) = getframe(gcf);
	mpgwrite(M,jet,[movie_location,'update.mpg']);
	system(['cat ',movie_location,movie_name,'.mpg ',movie_location,...
		'update.mpg > ',movie_location,'temp.mpg']);
	system(['mv ',movie_location,'temp.mpg ',movie_location,...
		movie_name,'.mpg']);
end

statistical_plots(delta_t,fig_saves_location,fig_saves_logical,figure_position,...
    include_statistical_plots_in_movie,movie_location,movie_logical,...
    movie_name,no_frames_for_statistical_plots,stats.counter,stats,time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display simulation info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['Number of iterations: ',int2str(iteration)])
disp(['Time taken: ',num2str(round(time_taken_in_main_loop+time_taken_to_start_of_loop)),' seconds'])

% profile viewer