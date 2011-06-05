function distance_between_nodes =...
	calculate_distance_between_vertices(reference_node_positions,other_node_positions)

no_nodes = size(reference_node_positions,1);

difference_in_x_components =...
	other_node_positions(:,1) - reference_node_positions(:,1);

difference_in_y_components =...
	other_node_positions(:,2) - reference_node_positions(:,2);

distance_between_nodes = zeros(no_nodes,1);
for i = 1:no_nodes
	
	distance_between_nodes(i) =...
		norm([difference_in_x_components(i),difference_in_y_components(i)],2);
end