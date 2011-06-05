function [cells,FEM_elements,FEM_nodes,refined_edge_matrix] =...
	remove_FEM_nodes_from_short_edges(cells,FEM_elements,FEM_nodes,vertices,...
	mesh_refinement_threshold,refined_edge_matrix)

FEM_nodes_refined_edge_logical = FEM_nodes.edge(:,1)>0;
FEM_nodes_refined_edges = find(FEM_nodes_refined_edge_logical);

list_of_refined_edges = FEM_nodes.edge(FEM_nodes_refined_edge_logical,:);

length_refined_edges = sqrt((FEM_nodes.previous_position(...
	list_of_refined_edges(:,1),2)-FEM_nodes.previous_position(...
	list_of_refined_edges(:,2),2)).^2+(FEM_nodes.previous_position(...
	list_of_refined_edges(:,1),1)-FEM_nodes.previous_position(...
	list_of_refined_edges(:,2),1)).^2);

edges_from_list_below_threshold_logical = length_refined_edges<mesh_refinement_threshold;

index_edges_below_threshold_in_list_of_refined_edges =...
	find(edges_from_list_below_threshold_logical);

% edges_below_threshold = refined_edges(edges_below_threshold_logical,:);

for current_edge_node_local =...
		1:length(index_edges_below_threshold_in_list_of_refined_edges)
	
	current_node = FEM_nodes_refined_edges(...
		index_edges_below_threshold_in_list_of_refined_edges(current_edge_node_local));
			
	%     current_node = edge_nodes(current_edge_node_local);
	
	node_1 = FEM_nodes.edge(current_node,1);
	node_2 = FEM_nodes.edge(current_node,2);
	
	% 	% this line is pretty slow! it is a bit silly doing this with matlab
	% 	% code.
	%     length_current_edge =...
	%         norm(FEM_nodes.previous_position(FEM_nodes.edge(edge_nodes(current_edge_node_local),1),:)-...
	%         FEM_nodes.previous_position(FEM_nodes.edge(edge_nodes(current_edge_node_local),2),:));
	%
	% %     length_current_edge
	%
	%     if length_current_edge < mesh_refinement_threshold
	
	%         disp('yes')
	
	node_1_cells = vertices.cells(node_1,:);
	node_1_cells = node_1_cells(node_1_cells>0);
	
	node_2_cells = vertices.cells(node_2,:);
	node_2_cells = node_2_cells(node_2_cells>0);
	
	cells_sharing_the_edge = node_1_cells(ismember(node_1_cells,node_2_cells));
	
	for current_cell_local = 1:length(cells_sharing_the_edge)
		
		current_cell_global = cells_sharing_the_edge(current_cell_local);
		current_cell_elements = cells.FEM_elements{current_cell_global};
		
		%             FEM_elements.nodes(current_cell_elements,:)
		
		for current_element_local = 1:length(current_cell_elements)
			
			current_element_global = current_cell_elements(current_element_local);
			
			if ismember(node_1,FEM_elements.nodes(current_element_global,:)) &&...
					ismember(current_node,FEM_elements.nodes(current_element_global,:))
				
				FEM_elements.nodes(current_cell_elements(current_element_local),...
					FEM_elements.nodes(current_cell_elements(current_element_local),:)==current_node) = node_2;
				
			elseif ismember(node_2,FEM_elements.nodes(current_element_global,:)) &&...
					ismember(current_node,FEM_elements.nodes(current_element_global,:))
				
				FEM_elements.nodes(current_cell_elements(current_element_local),:) = 0;
				element_to_remove_local = current_element_local;
				
			end
			
		end
		
		current_cell_elements(element_to_remove_local) = [];
		cells.FEM_elements{current_cell_global} = current_cell_elements;
		%             FEM_elements.nodes(current_cell_elements,:)
		
	end
	
	FEM_nodes.previous_position(current_node,:) = 0;
	FEM_nodes.concentration(current_node) = 0;
	FEM_nodes.edge(current_node,:) = 0;
	
	refined_edge_matrix(node_1,node_2) = 0;
	refined_edge_matrix(node_2,node_1) = 0;
	
	%     end
end

