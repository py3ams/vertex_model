function [cells,FEM_elements,FEM_nodes,refined_edge_matrix] =...
	remove_FEM_node_from_specific_edge(cells,edge,FEM_elements,...
	FEM_nodes,refined_edge_matrix,vertices)

node_1 = edge(1);
node_2 = edge(2);

edge_node_index =...
	find((FEM_nodes.edge(:,1)==node_1&FEM_nodes.edge(:,2)==node_2)|...
	(FEM_nodes.edge(:,1)==node_2&FEM_nodes.edge(:,2)==node_1));

if ~isempty(edge_node_index)
	
	node_1_cells = vertices.cells(node_1,:);
	node_1_cells = node_1_cells(node_1_cells>0);
	
	node_2_cells = vertices.cells(node_2,:);
	node_2_cells = node_2_cells(node_2_cells>0);
	
	cells_sharing_the_edge = node_1_cells(ismember(node_1_cells,node_2_cells));
	
	for current_cell_local = 1:length(cells_sharing_the_edge)
		
		current_cell_global = cells_sharing_the_edge(current_cell_local);
		current_cell_elements = cells.FEM_elements{current_cell_global};
		
		for current_element_local = 1:length(current_cell_elements)
			
			current_element_global = current_cell_elements(current_element_local);
			
			if ismember(node_1,FEM_elements.nodes(current_element_global,:)) &&...
					ismember(edge_node_index,FEM_elements.nodes(current_element_global,:))
				
				FEM_elements.nodes(current_cell_elements(current_element_local),...
					FEM_elements.nodes(current_cell_elements(current_element_local),:)==edge_node_index) = node_2;
				
			elseif ismember(node_2,FEM_elements.nodes(current_element_global,:)) &&...
					ismember(edge_node_index,FEM_elements.nodes(current_element_global,:))
				
				FEM_elements.nodes(current_cell_elements(current_element_local),:) = 0;
				element_to_remove_local = current_element_local;
				
			end
			
		end
		
		current_cell_elements(element_to_remove_local) = [];
		cells.FEM_elements{current_cell_global} = current_cell_elements;
		%             FEM_elements.nodes(current_cell_elements,:)
		
	end
	
	FEM_nodes.previous_position(edge_node_index,:) = 0;
	FEM_nodes.concentration(edge_node_index) = 0;
	FEM_nodes.edge(edge_node_index,:) = 0;
	
	refined_edge_matrix(node_1,node_2) = 0;
	refined_edge_matrix(node_2,node_1) = 0;
	
end

