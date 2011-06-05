function source_term = calculate_source_term(cells,current_chemical,FEM_elements,...
	FEM_nodes,FEM_nodes_index_in_real_nodes,no_real_nodes)

source_term = zeros(no_real_nodes,1);
no_cells = length(cells.vertices);

for current_cell = 1:no_cells
	
	cell_FEM_elements = cells.FEM_elements{current_cell};
	
	for current_cell_element_local = 1:length(cell_FEM_elements)
		
		current_cell_element_global = cell_FEM_elements(current_cell_element_local);
		current_cell_element_nodes =...
			FEM_elements.nodes(current_cell_element_global,:);
		
		current_cell_element_nodes_in_real_nodes =...
			FEM_nodes_index_in_real_nodes(current_cell_element_nodes);
		
		node_positions_current_element =...
			FEM_nodes.previous_position(current_cell_element_nodes,:);
		
		jac = [node_positions_current_element(3,1)-node_positions_current_element(1,1), ...
			node_positions_current_element(2,1)-node_positions_current_element(1,1); ...
			node_positions_current_element(3,2)-node_positions_current_element(1,2), ...
			node_positions_current_element(2,2)-node_positions_current_element(1,2)];
		
		det_jac = abs(det(jac));
		
		% we are updating 3 elements of the source term vector, corresponding to the
		% 3 basis functions that are non-zero in this triangle
		source_term(current_cell_element_nodes_in_real_nodes) =...
			source_term(current_cell_element_nodes_in_real_nodes)+...
			1/6*det_jac*cells.source_rate(current_cell,current_chemical)/...
			cells.area(current_cell);
		
	end
	
end



