disp('busy');close all;clear all;tic;%profile on

iterations_in_test_solution = 500;
no_refinements = 4;

iterations_in_true_solution = 4000;

root_directory = 'Saves/refinement_comparison/';
simulation_name = ['iterations_',num2str(iterations_in_test_solution),...
   '_refinements_',num2str(no_refinements),'/'];

true_solution_vars = load([root_directory,'true_solution/initial_save']);

% find no cells in true solution and check this is the same as the test
% solution. if it is not something has gone very wrong!
no_cells = length(true_solution_vars.cells.vertices);
no_nodes_in_true_solution = length(true_solution_vars.FEM_nodes.position);

[test_solution_vars.FEM_elements,test_solution_vars.FEM_nodes] =...
   create_FEM_mesh(true_solution_vars.cells,true_solution_vars.vertices,no_refinements);

test_solution_vars.FEM_nodes.previous_position = test_solution_vars.FEM_nodes.position;
no_nodes_in_test_solution = length(test_solution_vars.FEM_nodes.position);

test_solution_vars.FEM_nodes.concentration = initialise_concentration(...
   true_solution_vars.vertices.no_cells,test_solution_vars.FEM_nodes,...
   true_solution_vars.gradient_type,true_solution_vars.initial_concentration_magnitude,no_cells,...
   true_solution_vars.no_chemicals,true_solution_vars.source_width);

FEM_node_concentration_projection =...
   find_concentration_projection(no_cells,no_nodes_in_test_solution,...
   no_nodes_in_true_solution,test_solution_vars.FEM_nodes.concentration,...
   true_solution_vars.FEM_nodes.edge);

node_differences =...
   abs(FEM_node_concentration_projection-true_solution_vars.FEM_nodes.concentration);
assert(sum(sum(node_differences))==0)
   
true_solution_delta_t = true_solution_vars.delta_t;
test_solution_delta_t = 1/iterations_in_test_solution;

% initialise errors to 0
solution_error = 0;
error_per_iteration = zeros(iterations_in_true_solution,1);

% figure out the number of true iterations per test iteration. initialise a
% counter to 0 to cycle through test_iterations as the main counter goes
% through true iterations
test_solution_period = iterations_in_true_solution/iterations_in_test_solution;
test_iteration = 0;

real_nodes_logical = true_solution_vars.FEM_nodes.position(:,1)~=0|...
   true_solution_vars.FEM_nodes.position(:,2)~=0;

% loop over every iteration in the true solution
for true_iteration = 1:iterations_in_true_solution
      
   true_solution_vars = load([root_directory,'true_solution/iteration_',num2str(true_iteration)]);
   
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
         no_nodes_in_true_solution,test_solution_vars.FEM_nodes.concetration,...
         true_solution_vars.FEM_nodes.edge);

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