disp('busy');close all;clear all;tic;%profile on

iterations_in_first_test_solution = 500;
no_test_solutions = 4;

iterations_in_true_solution = 8000;
no_refinements_in_true_solution = 4;

configuration_noise = 0.5;
configuration_type = 'hexagonal';
grid_size = [10,10];
max_no_vertices = prod(grid_size);

initial_force_constant_magnitudes.area = 1e0;
initial_force_constant_magnitudes.deformation = 5e-4;
initial_force_constant_magnitudes.elongation = 1e-5;
initial_force_constant_magnitudes.perimeter = 1e-4;
initial_force_constant_magnitudes.tension = 2e-4;

initial_concentration_magnitude = 0.1;
no_chemicals = 1;
source_width = 0.15;

[cells.vertices,vertices.position] =...
    initial_cell_mesh(max_no_vertices,configuration_noise,configuration_type,grid_size);

[vertices.cells,vertices.no_cells] = CreateCellStore(cells.vertices,max_no_vertices);

boundary_vertices =...
    CreateBoundaryElement(vertices.cells,cells.vertices,vertices.no_cells);

cells.volume = max(mean_volume*(1+0.5*randn(no_cells,1)),...
    0.2*mean_volume*ones(no_cells,1));

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

% find no cells in true solution and check this is the same as the test
% solution. if it is not something has gone very wrong!
no_cells = length(cells.vertices);

[true_solution_FEM_elements,true_solution_FEM_nodes] =...
    create_FEM_mesh(cells,vertices,no_refinements_in_true_solution);

true_solution_FEM_nodes.previous_position = true_solution_FEM_nodes.position;

true_solution_FEM_nodes.concentration = initialise_concentration(...
    vertices.no_cells,true_solution_FEM_nodes,gradient_type,...
    initial_concentration_magnitude,no_cells,no_chemicals,source_width);

real_nodes_logical = true_solution_FEM_nodes.position(:,1)~=0|true_solution_FEM_nodes.position(:,2)~=0;
no_nodes_in_true_solution = length(true_solution_FEM_nodes.position);
true_solution_delta_t = 1/iterations_in_true_solution;

% we initialise these oddly so the first two lines of the following loop work
no_refinements_in_current_test_solution = -1;
iterations_in_current_test_solution = 1/2*iterations_in_first_test_solution;

FEM_elements = cell(no_test_solutions,1);
FEM_nodes = cell(no_test_solutions,1);
no_nodes_in_test_solution = zeros(no_test_solutions,1);
FEM_node_concentration_projection = cell(no_test_solutions,1);
node_differences = cell(no_test_solutions,1);
test_solution_delta_t = zeros(no_test_solutions,1);
solution_error = zeros(no_test_solutions,1);
errors_per_iteration = cell(iterations_in_true_solution,1);
test_solution_period = zeros(iterations_in_true_solution,1);
test_solution_no_refinements = zeros(iterations_in_true_solution,1);
test_solution_iterations = zeros(iterations_in_true_solution,1);
test_iteration = zeros(iterations_in_true_solution,1);

for current_solution = 1:no_test_solutions
    
    no_refinements_in_current_test_solution = no_refinements_in_current_test_solution+1;
    iterations_in_current_test_solution = 2*iterations_in_current_test_solution;

    test_solution_no_refinements(current_solution) = no_refinements_in_current_test_solution;
    test_solution_iterations(current_solution) = iterations_in_current_test_solution;
    
    [FEM_elements{current_solution},FEM_nodes{current_solution}] =...
        create_FEM_mesh(cells,vertices,no_refinements_in_current_test_solution);
    
    FEM_nodes{current_solution}.previous_position = FEM_nodes{current_solution.position};
    
    no_nodes_in_test_solution(current_solution) = length(FEM_nodes{current_solution}.position);
    
    FEM_nodes{current_solution}.concentration = initialise_concentration(...
        vertices.no_cells,FEM_nodes{current_solution},gradient_type,...
        initial_concentration_magnitude,no_cells,no_chemicals,source_width);
    
    FEM_node_concentration_projection{current_solution} =...
        find_concentration_projection(no_cells,no_nodes_in_test_solution(...
        current_solution),no_nodes_in_true_solution,FEM_nodes{...
        current_solution.concentration},true_solution_FEM_nodes.concentration);
    
    node_differences{current_solution} =...
        abs(FEM_node_concentration_projection{current_solution}-true_solution_FEM_nodes.concentration);

    assert(sum(sum(node_differences{current_solution}))==0)
    
    test_solution_delta_t(current_solution) = 1/iterations_in_current_test_solution;
    test_solution_period(current_solution) = iterations_in_true_solution/iterations_in_current_test_solution;
    
    errors_per_iteration{current_solution} = zeros(iterations_in_true_solution,1);
    
end



% loop over every iteration in the true solution
for true_iteration = 1:iterations_in_true_solution
      
    UPDATEPOS
   
   % for a period of 4, for example, this will cause the test solution to update at iteration 4,8,12,.. 
   if rem(true_iteration,test_solution_period)==0
      
      test_solution_vars.FEM_nodes.position = UpdateFEMNodePositions(true_solution_vars.cells.vertices,...
			true_solution_vars.vertices.position,true_solution_vars.cells.area,test_solution_vars.FEM_nodes.edge);
      
      assert(test_solution_vars.FEM_nodes.position(1,1) ==...
         true_solution_vars.FEM_nodes.position(1,1));
      assert(test_solution_vars.FEM_nodes.position(end-no_cells+1,1) ==...
         true_solution_vars.FEM_nodes.position(end-no_cells+1,1));
      
      test_solution_vars.FEM_nodes = FEM_solver(test_solution_delta_t,...
         true_solution_vars.diffusion_speed(1),test_solution_vars.FEM_elements,...
         test_solution_vars.FEM_nodes);

      test_solution_vars.FEM_nodes.previous_position = test_solution_vars.FEM_nodes.position;
      
      % takes the FEM node concentrations of the test solution and uses interpolation
      % to find values at the extra nodes that are in the true solution but not the
      % test solution
      FEM_node_concentration_projection =...
         find_concentration_projection(no_cells,no_nodes_in_test_solution,...
         no_nodes_in_true_solution,test_solution_vars,true_solution_vars);

   end
   
   node_differences = abs(FEM_node_concentration_projection-true_solution_vars.FEM_nodes.concentration);
   
   % it is not obvious why this expression gives an error measure. indeed, it is not
   % obvious why M comes into it at all. need to look at the maths to see why this is
   % the case. we are doing int_{Omega} |c - hat{c}|^2 dx, where hat{c} is the true
   % solution and c is the solution.
   error_per_iteration(true_iteration) = true_solution_delta_t*...
      node_differences(real_nodes_logical)'*true_solution_vars.M*node_differences(real_nodes_logical);

   solution_error = solution_error + error_per_iteration(true_iteration);
   
end

toc

eval(['save ',root_directory,simulation_name,'solution_error_new solution_error error_per_iteration'])

figure
plot(error_per_iteration)
disp(['Total solution error : ',num2str(solution_error)])

% profile viewer