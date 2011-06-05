function [cells,FEM_elements,FEM_nodes,refined_edge_matrix] = add_FEM_nodes_to_long_edges(...
	cells,FEM_elements,FEM_nodes,long_edges,refined_edge_matrix,vertices)

% it is much quicker to check all the long edges in one go like this and
% then only loop over those that are not in the matrix than it is to loop
% over them all and put a logical statement within the for loop. it is
% faster as it saves a lot of operations.
non_refined_long_edges = find(~refined_edge_matrix(sub2ind(size(...
	refined_edge_matrix),long_edges(:,1),long_edges(:,2))));

for current_non_refined_long_edge_local = 1:size(non_refined_long_edges)
	
	node_1 = long_edges(non_refined_long_edges(...
		current_non_refined_long_edge_local),1);
	
	node_2 = long_edges(non_refined_long_edges(...
		current_non_refined_long_edge_local),2);
	
	% 	this line appears to check whether the current edge is already in
	% 	FEM_nodes.edge. if it is not we have to add a new node to the
	% 	mid-point of the edge and update all the elements etc. i don't
	% 	understand why we need this line - surely if the edge is not in
	% 	refined_edge_matrix it should not be in FEM_nodes.edge, otherwise
	% 	test_refined_edge_matrix would throw up an error.
	if ~any((FEM_nodes.edge(:,1)==node_1&FEM_nodes.edge(:,2)==node_2)|...
			(FEM_nodes.edge(:,1)==node_2&FEM_nodes.edge(:,2)==node_1))
		
		no_cells = length(cells.vertices);
		FEM_index_zero_cell = length(FEM_nodes.previous_position)-no_cells;
		no_vertices = length(vertices.position);
		new_FEM_node = find(FEM_nodes.previous_position(no_vertices+1:FEM_index_zero_cell,1)==0,1)+no_vertices;
		
		if isempty(new_FEM_node)
			error('no space for new FEM edge node');
		end
		
		FEM_nodes.previous_position(new_FEM_node,:) =...
			mean([FEM_nodes.previous_position(node_1,:); ...
			FEM_nodes.previous_position(node_2,:)]);
		
		FEM_nodes.concentration(new_FEM_node) =...
			mean([FEM_nodes.concentration(node_1); ...
			FEM_nodes.concentration(node_2)]);
		
		FEM_nodes.edge(new_FEM_node,:) = [node_1 node_2];
		
		node_1_cells = vertices.cells(node_1,:);
		node_1_cells = node_1_cells(node_1_cells>0);
		
		node_2_cells = vertices.cells(node_2,:);
		node_2_cells = node_2_cells(node_2_cells>0);
		
		cells_sharing_the_edge = node_1_cells(ismember(node_1_cells,node_2_cells));
		
		reusable_FEM_elements = find(FEM_elements.nodes(:,1)==0,2);
		no_reusable_FEM_elements = length(reusable_FEM_elements);
		
		no_FEM_elements = length(FEM_elements.nodes);
		
		if no_reusable_FEM_elements < 2
			
			if no_reusable_FEM_elements < 1
				
				elements_to_use = no_FEM_elements+1:no_FEM_elements+2;
				FEM_elements.nodes(elements_to_use,:) = zeros(2,3);
				
			else
				
				FEM_elements.nodes(no_FEM_elements+1,:) = zeros(1,3);
				elements_to_use = [reusable_FEM_elements no_FEM_elements+1];
				
			end
			
		else
			
			elements_to_use = reusable_FEM_elements;
			
		end
		
		for current_cell_local = 1:length(cells_sharing_the_edge)
			
			current_cell_global = cells_sharing_the_edge(current_cell_local);
			current_cell_elements = cells.FEM_elements{current_cell_global};
			
			for current_element_local = 1:length(current_cell_elements)
				
				current_element_global = current_cell_elements(current_element_local);
				
				if FEM_elements.nodes(current_element_global,1) == node_1 &&...
						FEM_elements.nodes(current_element_global,2) == node_2
					
					FEM_elements.nodes(current_element_global,2) = new_FEM_node;
					
					FEM_elements.nodes(elements_to_use(1),:) =...
						[new_FEM_node node_2 FEM_index_zero_cell+current_cell_global];
					
					cells.FEM_elements{current_cell_global} =...
						[current_cell_elements(1:current_element_local) ...
						elements_to_use(1) current_cell_elements(...
						current_element_local+1:length(current_cell_elements))];
					
					break;
					
				elseif FEM_elements.nodes(current_element_global,1) == node_2 &&...
						FEM_elements.nodes(current_element_global,2) == node_1
					
					FEM_elements.nodes(current_element_global,2) = new_FEM_node;
					
					FEM_elements.nodes(elements_to_use(2),:) =...
						[new_FEM_node node_1 FEM_index_zero_cell+current_cell_global];
					
					cells.FEM_elements{current_cell_global} =...
						[current_cell_elements(1:current_element_local) ...
						elements_to_use(2) current_cell_elements(...
						current_element_local+1:length(current_cell_elements))];
					
					break;
					
				end
			end
		end
	end
	
	refined_edge_matrix(node_1,node_2) = 1;
	refined_edge_matrix(node_2,node_1) = 1;
	
end