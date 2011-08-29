function [FEM_nodes,M,real_nodes_logical] =...
    FEM_solver(delta_t,diffusion_speed,FEM_elements,FEM_nodes)

% find the nodes that are actually currrently in use.
real_nodes_logical = FEM_nodes.position(:,1)~=0;
no_real_nodes = sum(real_nodes_logical);

real_node_positions = FEM_nodes.position(real_nodes_logical,:);
real_previous_node_positions =...
    FEM_nodes.previous_position(real_nodes_logical,:);

% the length of FEM_nodes_index_in_real_nodes is the same as FEM_nodes.position, assuming 
% that the last entry in FEM_nodes.position is a real node. for each FEM node,
% this vector contains its location within real_nodes_logical. for the nodes that
% aren't currently in use, and therefore not in real_nodes_logical, there will be a 0.
% it shouldn't matter if reverse indices is exactly the same length as
% FEM_nodes.positions, i.e. if the last cell has died, as in that case the node won't
% be in FEM_elements, so when we use FEM_nodes_index_in_real_nodes there won't be a problem
FEM_nodes_index_in_real_nodes(real_nodes_logical) = 1:no_real_nodes;

% we only need to worry about the FEM_elements that actually exist (i.e. are
% non-zero)
FEM_elements_stripped = FEM_elements.nodes(FEM_elements.nodes(:,1)>0,:);
FEM_elements_real_node_indices = FEM_nodes_index_in_real_nodes(FEM_elements_stripped);

% create M matrix for previous_node_positions. this is not the same as just storing
% the M matrix from the previous iteration, as the mesh may have been updated by T1
% swaps and deaths etc.
[I_prev,J_prev,MV_prev] = Stiff2DMonly(FEM_elements_real_node_indices,...
    real_previous_node_positions);

M_prev = sparse(I_prev,J_prev,MV_prev);

% create the M, A, and W matrices for the current node positions. M_const is the M
% matrix using piecewise constant basis functions, whereas M uses piecewise linear
% functions
[I,J,AV,MV,WV] =...
    Stiff2D(delta_t,FEM_elements_real_node_indices,real_node_positions,...
    real_previous_node_positions);

A = sparse(I,J,AV);
M = sparse(I,J,MV);
W = sparse(I,J,WV);

COEFF_MAT = M+delta_t*(diffusion_speed*A+W); 
rhs = M_prev*FEM_nodes.concentration(real_nodes_logical);

FEM_nodes.concentration(real_nodes_logical) = COEFF_MAT\rhs;
    