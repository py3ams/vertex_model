disp('busy');close all;clear all;tic;%profile on

% the first test solution will have this number of iterations and no refinements. the
% number of refinements will increase by one each time (equivalent to halving the
% spatial step) and the number of iterations doubles, halving the time step
iterations_in_first_test_solution = 500;
no_test_solutions = 5;

iterations_in_true_solution = 8000;
no_refinements_in_true_solution = 4;

load_initial_cells_from_file_logical = true;
file_to_load = 'Saves/refinement_comparison_movie/initial_save.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

root_directory = 'Saves/refinement_comparison_realtime/';
simulation_name = ['true_iterations_',num2str(iterations_in_true_solution),...
    '_refinements_',num2str(no_refinements_in_true_solution),'/'];

configuration_noise = 0.5;
configuration_type = 'hexagonal';
grid_size = [10,10];
viscosity = 1;

target_area_factor = 1.0;

initial_force_constant_magnitudes.area = 2000;
initial_force_constant_magnitudes.deformation = 1e-1;
initial_force_constant_magnitudes.elongation = 1e-3;
initial_force_constant_magnitudes.perimeter = 5e-2;
initial_force_constant_magnitudes.tension = 1e-1;

boundary_force_constants.deformation = 1e-1;
boundary_force_constants.edge = 5e1;

diffusion_speed = 0.1;
initial_concentration_magnitude = 0.1;
no_chemicals = 1;
source_width = 0.15;
gradient_type = 4;

tension_anisotropy_factor = 0.0;

movie_logical = false;
chemical_to_view = 1;
view_FEM_concentration = true;
view_FEM_mesh = true;
view_initial_config = false;
update_period = 100;
axis_values = 0.6*[-1 1 -1 1];
axis_values_FEM = [axis_values 0 0.1];
view_iteration_number = true;
view_number_cells = true;
linewidth_cells = 2;
linewidth_elements = 1;
solution_to_view = 1;

full_saves_logical = true;
full_saves_period = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_no_cells = prod(grid_size);
no_vertices = 4*max_no_cells;

if load_initial_cells_from_file_logical
   
   loaded_vars = load(file_to_load);
   cells.vertices = loaded_vars.cells.vertices;
   vertices.position = loaded_vars.vertices.position;
   
else

   [cells.vertices,vertices.position] =...
      initial_cell_mesh(no_vertices,configuration_noise,configuration_type,grid_size);

end

no_cells = length(cells.vertices);
cells.original_logical = true(no_cells,1);
cells.state = ones(no_cells,1);

[vertices.cells,vertices.no_cells] = CreateCellStore(cells.vertices,no_vertices);

cells.area = CalculateCellAreas(cells.vertices,vertices.position);
mean_volume = mean(sqrt(cells.area.^3));

cells.volume = max(mean_volume*(1+0.5*randn(no_cells,1)),...
    0.2*mean_volume*ones(no_cells,1));

cells.target_area = target_area_factor*mean(cells.area)*ones(no_cells,1);

cells.force_constants.area =...
    initial_force_constant_magnitudes.area*ones(no_cells,1);
cells.force_constants.deformation =...
    initial_force_constant_magnitudes.deformation*ones(no_cells,1);
cells.force_constants.elongation =...
    initial_force_constant_magnitudes.elongation*ones(no_cells,1);
cells.force_constants.perimeter =...
    initial_force_constant_magnitudes.perimeter*ones(no_cells,1);
cells.force_constants.tension =...
    initial_force_constant_magnitudes.tension*ones(no_cells,1);

no_cells = length(cells.vertices);

full_saves_location = [root_directory,simulation_name];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[true_solution_FEM_elements,true_solution_FEM_nodes] =...
    create_FEM_mesh(cells,vertices,no_refinements_in_true_solution,true);

true_solution_FEM_nodes.previous_position = true_solution_FEM_nodes.position;

true_solution_FEM_nodes.concentration = initialise_concentration(...
    vertices.no_cells,true_solution_FEM_nodes,gradient_type,...
    initial_concentration_magnitude,no_cells,no_chemicals,source_width);

no_nodes_in_true_solution = length(true_solution_FEM_nodes.position);
true_solution_delta_t = 1/iterations_in_true_solution;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we initialise these oddly so the first two lines of the following loop work
no_refinements_in_current_test_solution = -1;
iterations_in_current_test_solution = 1/2*iterations_in_first_test_solution;

FEM_elements = cell(no_test_solutions,1);
FEM_nodes = cell(no_test_solutions,1);
no_nodes_in_test_solution = zeros(no_test_solutions,1);
FEM_node_concentration_projection = cell(no_test_solutions,1);
node_differences = cell(no_test_solutions,1);
test_solution_delta_t = zeros(no_test_solutions,1);
solution_error_squared = zeros(no_test_solutions,1);
solution_error = zeros(no_test_solutions,1);
error_per_iteration = cell(no_test_solutions,1);
test_solution_period = zeros(no_test_solutions,1);
test_solution_no_refinements = zeros(no_test_solutions,1);
test_solution_iterations = zeros(no_test_solutions,1);
test_iteration = zeros(no_test_solutions,1);

for current_solution = 1:no_test_solutions
    
    no_refinements_in_current_test_solution = no_refinements_in_current_test_solution+1;
    iterations_in_current_test_solution = 2*iterations_in_current_test_solution;

    test_solution_no_refinements(current_solution) = no_refinements_in_current_test_solution;
    test_solution_iterations(current_solution) = iterations_in_current_test_solution;
    
    [FEM_elements{current_solution},FEM_nodes{current_solution}] =...
        create_FEM_mesh(cells,vertices,no_refinements_in_current_test_solution,true);
    
    FEM_nodes{current_solution}.previous_position = FEM_nodes{current_solution}.position;
    
    no_nodes_in_test_solution(current_solution) = length(FEM_nodes{current_solution}.position);
    
    FEM_nodes{current_solution}.concentration = initialise_concentration(...
        vertices.no_cells,FEM_nodes{current_solution},gradient_type,...
        initial_concentration_magnitude,no_cells,no_chemicals,source_width);
    
    FEM_node_concentration_projection{current_solution} =...
        find_concentration_projection(no_cells,no_nodes_in_test_solution(...
        current_solution),no_nodes_in_true_solution,FEM_nodes{...
        current_solution}.concentration,true_solution_FEM_nodes.edge);
    
    node_differences{current_solution} =...
        abs(FEM_node_concentration_projection{current_solution}-true_solution_FEM_nodes.concentration);

    assert(sum(sum(node_differences{current_solution}))==0)
    
    test_solution_delta_t(current_solution) = 1/iterations_in_current_test_solution;
    test_solution_period(current_solution) = iterations_in_true_solution/iterations_in_current_test_solution;
    
    error_per_iteration{current_solution} = zeros(iterations_in_true_solution,1);
    
end

if full_saves_logical
    if exist(full_saves_location,'dir')
        yesno = input(['Are you sure you want to overwrite directory ',...
            full_saves_location,'? '],'s');
        if strcmp(yesno,'yes') || strcmp(yesno,'y')
            delete([full_saves_location,'*']);
        else
            error(['Could not delete directory ',full_saves_location]);
        end
    else
        mkdir(full_saves_location);
    end
    save([full_saves_location,'initial_save'])
end

if movie_logical
    
    if strcmp(solution_to_view,'true')

        visualiser(cells,vertices,true_solution_FEM_elements,true_solution_FEM_nodes,...
            axis_values,axis_values_FEM,chemical_to_view,false,0,...
            linewidth_cells,linewidth_elements,movie_logical,update_period,...
            view_FEM_concentration,view_FEM_mesh,view_initial_config,...
            view_iteration_number,view_number_cells)
        
    else
        
        visualiser(cells,vertices,FEM_elements{solution_to_view},FEM_nodes{solution_to_view},...
            axis_values,axis_values_FEM,chemical_to_view,false,0,...
            linewidth_cells,linewidth_elements,movie_logical,update_period,...
            view_FEM_concentration,view_FEM_mesh,view_initial_config,...
            view_iteration_number,view_number_cells);
        
    end
    
end

boundary_element =...
    CreateBoundaryElement(vertices.cells,cells.vertices,vertices.no_cells);

for true_iteration = 1:iterations_in_true_solution
      
    [cells.area,cells.perimeter,cells.edge_lengths,~,~,~,current_edge_length_stats] =...
        CalculateCellAreas(cells.vertices,vertices.position);
    
    mean_edge_length = current_edge_length_stats(1);
    
    vertices.position =...
        UpdatePos(cells.vertices,vertices.position,boundary_element,cells.area,...
        cells.perimeter,cells.volume,vertices.no_cells,true_solution_delta_t,...
        cells.edge_lengths,mean_edge_length,cells.target_area,...
        tension_anisotropy_factor,viscosity,cells.force_constants.area,...
        boundary_force_constants.deformation,boundary_force_constants.edge,...
        cells.force_constants.deformation,cells.force_constants.elongation,...
        cells.force_constants.perimeter,cells.force_constants.tension);
    
    cells.area = CalculateCellAreas(cells.vertices,vertices.position);
    
    true_solution_FEM_nodes.position = UpdateFEMNodePositions(cells.vertices,...
        vertices.position,cells.area,true_solution_FEM_nodes.edge);
    
    [true_solution_FEM_nodes,M,real_nodes_logical] = FEM_solver_stripped(...
        true_solution_delta_t,diffusion_speed,true_solution_FEM_elements,...
        true_solution_FEM_nodes);
    
    true_solution_FEM_nodes.previous_position = true_solution_FEM_nodes.position;
    
    for current_solution = 1:no_test_solutions
        
        % for a period of 4, for example, this will cause the test solution to update at iteration 4,8,12,..
        if rem(true_iteration,test_solution_period(current_solution))==0
            
            FEM_nodes{current_solution}.position = UpdateFEMNodePositions(cells.vertices,...
                vertices.position,cells.area,FEM_nodes{current_solution}.edge);
            
            assert(FEM_nodes{current_solution}.position(1,1) ==...
                true_solution_FEM_nodes.position(1,1));
            
            assert(FEM_nodes{current_solution}.position(end-no_cells+1,1) ==...
                true_solution_FEM_nodes.position(end-no_cells+1,1));
            
            FEM_nodes{current_solution} = FEM_solver_stripped(test_solution_delta_t(current_solution),...
                diffusion_speed,FEM_elements{current_solution},FEM_nodes{current_solution});
            
            FEM_nodes{current_solution}.previous_position = FEM_nodes{current_solution}.position;
            
            % takes the FEM node concentrations of the test solution and uses interpolation
            % to find values at the extra nodes that are in the true solution but not the
            % test solution
            FEM_node_concentration_projection{current_solution} =...
                find_concentration_projection(no_cells,no_nodes_in_test_solution(...
                current_solution),no_nodes_in_true_solution,FEM_nodes{...
                current_solution}.concentration,true_solution_FEM_nodes.edge);
            
        end
        
        node_differences{current_solution} = abs(FEM_node_concentration_projection{...
            current_solution}-true_solution_FEM_nodes.concentration);
   
        % it is not obvious why this expression gives an error measure. indeed, it is not
        % obvious why M comes into it at all. need to look at the maths to see why this is
        % the case. we are doing int_{Omega} |c - hat{c}|^2 dx, where hat{c} is the true
        % solution and c is the solution.
        error_per_iteration{current_solution}(true_iteration) = true_solution_delta_t*...
            node_differences{current_solution}(real_nodes_logical)'*M*...
            node_differences{current_solution}(real_nodes_logical);
        
        solution_error_squared(current_solution) = solution_error_squared(current_solution) +...
            error_per_iteration{current_solution}(true_iteration);
        
        
        
    end
    
    if full_saves_logical && ~rem(true_iteration,full_saves_period)
    
        save([full_saves_location,'iteration_',num2str(true_iteration)]);
        
    end
        
    if movie_logical
        
        if strcmp(solution_to_view,'true')
            
            visualiser(cells,vertices,true_solution_FEM_elements,true_solution_FEM_nodes,...
                axis_values,axis_values_FEM,chemical_to_view,false,true_iteration,...
                linewidth_cells,linewidth_elements,movie_logical,update_period,...
                view_FEM_concentration,view_FEM_mesh,view_initial_config,...
                view_iteration_number,view_number_cells)
            
        else
            
            visualiser(cells,vertices,FEM_elements{solution_to_view},FEM_nodes{solution_to_view},...
                axis_values,axis_values_FEM,chemical_to_view,false,true_iteration,...
                linewidth_cells,linewidth_elements,movie_logical,update_period,...
                view_FEM_concentration,view_FEM_mesh,view_initial_config,...
                view_iteration_number,view_number_cells);
            
        end
        
    end
   
end

toc

for current_solution = 1:no_test_solutions
    
    solution_error(current_solution) = sqrt(solution_error_squared(current_solution));
    
    disp(['Iterations : ',num2str(test_solution_iterations(current_solution)),...
        ', Refinements : ',num2str(test_solution_no_refinements(current_solution)),...
        ', Solution error : ',num2str(solution_error(current_solution))])
end

if full_saves_logical
    
    save([full_saves_location,'final_save']);
    eval(['save ',root_directory,simulation_name,'solution_error solution_error'])
    
end
