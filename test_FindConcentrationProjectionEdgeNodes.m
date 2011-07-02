disp('busy');close all;clear all;

no_nodes_in_true_solution = 500;
no_nodes_in_test_solution = 200;
no_cells = 100;

true_solution_FEM_edge_nodes = zeros(no_nodes_in_true_solution,2);
true_solution_FEM_edge_nodes(200:299,1) = ceil(100*rand(100,1));
true_solution_FEM_edge_nodes(200:249,2) = ceil(100*rand(50,1));
true_solution_FEM_edge_nodes(250:299,2) = ceil(100*rand(50,1))+400;

FEM_node_concentration_projection1 = zeros(no_nodes_in_true_solution,1);

FEM_node_concentration_projection1(1:100) = rand(100,1);
FEM_node_concentration_projection1(401:500) = rand(100,1);

FEM_node_concentration_projection2 = FEM_node_concentration_projection1;

for current_node = no_nodes_in_test_solution-no_cells+1:no_nodes_in_true_solution-no_cells
   
   if true_solution_FEM_edge_nodes(current_node,1)>0
      
      FEM_node_concentration_projection1(current_node) = 0.5*...
         (FEM_node_concentration_projection1(true_solution_FEM_edge_nodes(current_node,1))+...
         FEM_node_concentration_projection1(true_solution_FEM_edge_nodes(current_node,2)));
      
   end
   
end

mex FindConcentrationProjectionEdgeNodes.cpp

FEM_node_concentration_projection2 = FindConcentrationProjectionEdgeNodes(...
   FEM_node_concentration_projection2,true_solution_FEM_edge_nodes,...
   no_cells,no_nodes_in_test_solution);

assert(sum(sum(FEM_node_concentration_projection2-FEM_node_concentration_projection1))==0)