disp('busy');close all;clear all;

iterations = 500;
no_refinements = 0;

iterations_in_true_solution = 4000;

simulation_name = ['iterations_',num2str(iterations),...
   '_refinements_',num2str(no_refinements)];

true_solution_initial_vars = load('Saves/true_solution/initial_save');
test_solution_initial_vars = load(['Saves/',simulation_name,'/initial_save']);
  
no_vertices = length(true_solution_initial_vars.vertices.position);
assert(length(test_solution_initial_vars.vertices.position)==no_vertices);

no_cells = length(true_solution_initial_vars.cells.vertices);
assert(length(test_solution_initial_vars.cells.vertices)==no_cells);

no_nodes_in_true_solution = length(true_solution_initial_vars.FEM_nodes.concentration);

FEM_node_concentration_projection = zeros(no_nodes_in_true_solution,1);

for iteration = 1:iterations_in_true_solution

true_solution_vars = load('Saves/true_solution/final_save');
test_solution_vars = load(['Saves/',simulation_name,'/final_save']);

FEM_node_concentration_projection(1:no_vertices) = test_solution_vars.FEM_nodes.concentration(1:no_vertices);
FEM_node_concentration_projection(end-no_cells:end) = test_solution_vars.FEM_nodes.concentration(end-no_cells:end);

for current_node = 1:no_nodes_in_true_solution
   
   if true_solution_vars.FEM_nodes.edge(current_node,1)>0
      
      FEM_node_concentration_projection(current_node) = 0.5*...
         (FEM_node_concentration_projection(true_solution_vars.FEM_nodes.edge(current_node,1))+...
         FEM_node_concentration_projection(true_solution_vars.FEM_nodes.edge(current_node,2)));

   end
   
end

node_differences = abs(FEM_node_concentration_projection-true_solution_vars.FEM_nodes.concentration);

real_nodes_logical = true_solution_vars.FEM_nodes.position(:,1)~=0;
real_node_positions = true_solution_vars.FEM_nodes.position(real_nodes_logical,:);

no_real_nodes = sum(real_nodes_logical);
FEM_nodes_index_in_real_nodes(real_nodes_logical) = 1:no_real_nodes;

FEM_elements_stripped = true_solution_vars.FEM_elements.nodes(true_solution_vars.FEM_elements.nodes(:,1)>0,:);
FEM_elements_real_node_indices = FEM_nodes_index_in_real_nodes(FEM_elements_stripped);

[I,J,MV] = Stiff2DMonly(FEM_elements_real_node_indices,...
    real_node_positions);

M = sparse(I,J,MV);

error = node_differences(real_nodes_logical)'*M*node_differences(real_nodes_logical);