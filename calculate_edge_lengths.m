function edge_lengths = calculate_edge_lengths(cell_node_positions,edges_local)

if nargin == 1
    no_cell_nodes = size(cell_node_positions,1);
    edges_local = [(1:no_cell_nodes)' [2:no_cell_nodes 1]'];
end

edge_lengths = calculate_distance_between_nodes(cell_node_positions(...
	edges_local(:,1),:),cell_node_positions(edges_local(:,2),:));
