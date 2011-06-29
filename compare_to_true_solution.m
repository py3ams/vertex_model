disp('busy');close all;clear all;tic;

iterations_in_test_solution = 100;
no_refinements = 2;

simulation_name = ['iterations_',num2str(iterations_in_test_solution),...
   '_refinements_',num2str(no_refinements)];

iterations_in_true_solution = 100;

true_solution_initial_vars = load('Saves/test_true_solution/initial_save');
test_solution_initial_vars = load(['Saves/',simulation_name,'/initial_save']);

true_solution_delta_t = true_solution_initial_vars.delta_t;

% find no cells in true solution and check this is the same as the test
% solution. if it is not something has gone very wrong!
no_cells = length(true_solution_initial_vars.cells.vertices);
assert(length(test_solution_initial_vars.cells.vertices)==no_cells);

no_nodes_in_test_solution = length(test_solution_initial_vars.FEM_nodes.concentration);
no_nodes_in_true_solution = length(true_solution_initial_vars.FEM_nodes.concentration);

% might want to divide error by this to understand relative magnitude
total_concentration = CalculateTotalDpp(true_solution_initial_vars.FEM_nodes.concentration,...
   true_solution_initial_vars.FEM_elements.nodes,true_solution_initial_vars.FEM_nodes.position);

% find which FEM nodes in the true solution are actually being used. this shouldn't change
% through time as there are no rearrangements going off
real_nodes_logical = true_solution_initial_vars.FEM_nodes.position(:,1)~=0;

% initialise errors to 0
solution_error = 0;
error_per_iteration = zeros(iterations_in_true_solution,1);

% figure out the number of true iterations per test iteration. initialise a
% counter to 0 to cycle through test_iterations as the main counter goes
% through true iterations
test_solution_period = iterations_in_true_solution/iterations_in_test_solution;
test_iteration = 0;

% for the first few true iterations, we use the initial solution for the
% test. this gets updated to the solution at the first iteration of test
% when the correct number of true iterations have passed corresponding to one
% test iteration. could instead update mid-way through or at the start - are these
% equally valid?
test_solution_vars.FEM_nodes = test_solution_initial_vars.FEM_nodes;
true_solution_vars.FEM_nodes = true_solution_initial_vars.FEM_nodes;

FEM_node_concentration_projection =...
   find_concentration_projection(no_cells,no_nodes_in_test_solution,...
   no_nodes_in_true_solution,test_solution_vars,true_solution_vars);

% the following are needed to calculate M for the initial condition, as this is
% not calculated in cell_dynamics so will not be found in true_solution_vars
no_real_nodes = sum(real_nodes_logical);
FEM_nodes_index_in_real_nodes = zeros(no_real_nodes,1);
FEM_nodes_index_in_real_nodes(real_nodes_logical) = 1:no_real_nodes;
FEM_elements_stripped = true_solution_initial_vars.FEM_elements.nodes(true_solution_initial_vars.FEM_elements.nodes(:,1)>0,:);
FEM_elements_real_node_indices = FEM_nodes_index_in_real_nodes(FEM_elements_stripped);

initial_node_differences = abs(FEM_node_concentration_projection-true_solution_vars.FEM_nodes.concentration);
   
real_node_positions = true_solution_vars.FEM_nodes.position(real_nodes_logical,:);
[I,J,MV] = Stiff2DMonly(FEM_elements_real_node_indices,real_node_positions);
    
M_initial = sparse(I,J,MV);
   
initial_solution_error = true_solution_delta_t*...
    initial_node_differences(real_nodes_logical)'*...
    M_initial*initial_node_differences(real_nodes_logical);

% should be 0!
disp(['Initial solution error = ',num2str(initial_solution_error)])

% loop over every iteration in the true solution
for true_iteration = 1:iterations_in_true_solution
      
   true_solution_vars = load(['Saves/test_true_solution/iteration_',num2str(true_iteration)]);
   
   % for a period of 4, for example, this will cause the test solution to update at iteration 4,8,12,.. 
   if ~rem(true_iteration,test_solution_period)
      
      test_iteration = test_iteration+1;
      test_solution_vars = load(['Saves/',simulation_name,'/iteration_',num2str(test_iteration)]);
      
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

eval(['save Saves/',simulation_name,'/solution_error solution_error error_per_iteration'])

plot(error_per_iteration)
disp(['Total solution error : ',num2str(solution_error)])