function FEM_node_concentration_projection =...
   find_concentration_projection(no_cells,no_nodes_in_test_solution,...
   no_nodes_in_true_solution,test_solution_FEM_node_concentrations,...
   true_solution_FEM_node_edges)

% create a vector the same size as true FEM_nodes.concentration to
% store the projection/interpolation of the test solution
FEM_node_concentration_projection = zeros(no_nodes_in_true_solution,1);

% put the values from the test solution into the appropriate places
% in the projection. cell centres go at the end, everything else at
% the start.
FEM_node_concentration_projection(1:no_nodes_in_test_solution-no_cells) =...
   test_solution_FEM_node_concentrations(1:no_nodes_in_test_solution-no_cells);

FEM_node_concentration_projection(end-no_cells+1:end) =...
   test_solution_FEM_node_concentrations(end-no_cells+1:end);

% loop over nodes that are not part of the test solution and find
% their value as a linear combination of neighbouring nodes. it is
% important not to loop over all nodes, as values will be overwritten
% in cases where refinement has occurred at a lower level to the true
% solution.
% for current_node = no_nodes_in_test_solution-no_cells+1:no_nodes_in_true_solution-no_cells
%    
%    if true_solution_vars.FEM_nodes.edge(current_node,1)>0
%       
%       FEM_node_concentration_projection(current_node) = 0.5*...
%          (FEM_node_concentration_projection(true_solution_vars.FEM_nodes.edge(current_node,1))+...
%          FEM_node_concentration_projection(true_solution_vars.FEM_nodes.edge(current_node,2)));
%       
%    end
%    
% end

FEM_node_concentration_projection = FindConcentrationProjectionEdgeNodes(...
   FEM_node_concentration_projection,true_solution_FEM_node_edges,...
   no_cells,no_nodes_in_test_solution);

