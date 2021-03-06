function [cells,FEM_nodes,M,source_magnitude,stats,total_ingestion,total_source_released] = ...
   FEM_solver(cells,degradation_constant,delta_t,diffusion_speed,FEM_elements,...
   FEM_nodes,gradient_type,internal_chemical_uptake_type,maximum_source_to_release,no_chemicals,source_magnitude,...
   source_type,source_width,refined_edge_matrix,stats,total_ingestion,total_source_released)

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
[I,J,AV,MV,WV,triangle_quality] =...
   Stiff2D(delta_t,FEM_elements_real_node_indices,real_node_positions,...
   real_previous_node_positions);

A = sparse(I,J,AV);
M = sparse(I,J,MV);
W = sparse(I,J,WV);

source_this_iteration = zeros(no_chemicals,1);

previous_concentration = FEM_nodes.concentration;

cell_ingestion_functions = zeros(length(cells.vertices),no_chemicals);

for current_chemical = 1:no_chemicals
   
   COEFF_MAT = M+delta_t*(diffusion_speed(current_chemical)*A+W);
   % COEFF_MAT = M+delta_t*(diffusion_speed(current_chemical)*A);
   
   if source_type == 1
      
      source_functions = zeros(no_real_nodes,1);
      
      if gradient_type(current_chemical) == 1
         
         source_functions(abs(real_node_positions(:,1))<(0.5*source_width(current_chemical))) = source_magnitude(current_chemical);
         
      elseif gradient_type(current_chemical) == 2
         
         source_functions(abs(real_node_positions(:,2))<(0.5*source_width(current_chemical))) = source_magnitude(current_chemical);
         
      elseif gradient_type(current_chemical) == 3
         
         source_functions(sqrt((real_node_positions(:,1).^2)+(real_node_positions(:,2).^2))<(0.5*source_width(current_chemical))) = source_magnitude(current_chemical);
         
      elseif gradient_type(current_chemical) == 4
         
         min_x_position = min(real_node_positions(:,1));
         source_functions(abs(real_node_positions(:,1)-min_x_position)<(0.5*source_width(current_chemical))) = source_magnitude(current_chemical);
         
      end
            
      rhs = M_prev*FEM_nodes.concentration(real_nodes_logical,current_chemical)*...
         (1-delta_t*degradation_constant(current_chemical)) + delta_t*M*source_functions;
      
      % to find the total source released during this iteration, we must
      % intergrate over space and time. CalculateTotalDpp does the spatial
      % integration, then we approximate the time integration by multiplying by
      % delta_t. note that we integrate the source functions, not
      % M*source_functions, which is totally different. we are only able to
      % do this because of the way we have set up the source functions
      % currently.
      source_this_iteration(current_chemical) =...
         delta_t*CalculateTotalDpp(source_functions,FEM_elements_real_node_indices,real_node_positions);

      no_activated_source_functions = sum(source_functions>0);
      
   elseif source_type == 2
      
      source_term = CalculateSourceTerm(cells.FEM_elements,cells.source_rate(:,...
         current_chemical),cells.area,FEM_elements.nodes,FEM_nodes.previous_position,...
         FEM_nodes_index_in_real_nodes,no_real_nodes);
      
      source_this_iteration(current_chemical) = delta_t*sum(cells.source_rate(:,current_chemical));
      
      cell_internal_chemical_to_maximum_ratios = cells.internal_chemical_quantity(:,...
         current_chemical)./cells.maximum_internal_chemical_quantity(:,current_chemical);
      
      if internal_chemical_uptake_type==1
         cell_ingestion_functions(:,current_chemical) = cells.ingestion_rate(:,current_chemical);
      elseif internal_chemical_uptake_type==2
         cell_ingestion_functions(:,current_chemical) = cells.ingestion_rate(:,current_chemical).*...
            (1-cell_internal_chemical_to_maximum_ratios);
      end
      
      ingestion_term = CalculateIngestionTerm(cells.FEM_elements,cell_ingestion_functions(:,current_chemical),...
         FEM_elements.nodes,FEM_nodes.concentration(:,current_chemical),FEM_nodes.previous_position,...
         FEM_nodes_index_in_real_nodes,no_real_nodes);
      
      rhs = M_prev*FEM_nodes.concentration(real_nodes_logical,current_chemical)...
         + delta_t*source_term - delta_t*ingestion_term;
      
      no_activated_source_functions = sum(cells.source_rate(:,current_chemical)>0);
      
   end
   
   FEM_nodes.concentration(real_nodes_logical,current_chemical) = COEFF_MAT\rhs;
   
   total_source_released(current_chemical) =...
      total_source_released(current_chemical)+source_this_iteration(current_chemical);
   
   if total_source_released(current_chemical) > maximum_source_to_release(current_chemical)
      
      source_magnitude(current_chemical) = 0;
      
   end
   
end

cells = calculate_ingested_amount(cells,cell_ingestion_functions,delta_t,...
   FEM_elements,FEM_nodes,no_chemicals,previous_concentration);

if stats.this_iteration_logical
   
   stats.no_activated_source_functions(stats.counter) = no_activated_source_functions;
   
   real_FEM_concentrations = FEM_nodes.concentration(real_nodes_logical,1);
   stats.concentration_node_values(stats.counter,:) =...
      [mean(real_FEM_concentrations) max(real_FEM_concentrations) ...
      min(real_FEM_concentrations) std(real_FEM_concentrations)];
   
   stats.total_concentration(stats.counter,:) = CalculateTotalDpp(...
      FEM_nodes.concentration(:,1),FEM_elements_stripped,FEM_nodes.position);
      
   stats.chemical_source(stats.counter,:) = [total_source_released(1) source_this_iteration(1)];
   
   real_cells_logical = cells.state~=3&cells.state~=4;
   stats.internal_chemical_quantity(stats.counter,:) = [mean(cells.internal_chemical_quantity(real_cells_logical,1)) ...
      max(cells.internal_chemical_quantity(real_cells_logical,1)) min(cells.internal_chemical_quantity(real_cells_logical,1)) ...
      std(cells.internal_chemical_quantity(real_cells_logical,1)) sum(cells.internal_chemical_quantity(real_cells_logical,1))];
   
   stats.chemical_in_cell_1(stats.counter) = cells.internal_chemical_quantity(1,1);
   
   stats.triangle_quality(stats.counter,:) = ...
      [mean(triangle_quality) max(triangle_quality) ...
      min(triangle_quality) std(triangle_quality)];
   stats.no_refined_edges(stats.counter) = sum(sum(refined_edge_matrix))/2;
   
end


