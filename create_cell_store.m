function [cell_store,cells_per_node] = create_cell_store(array_sizes,cells)

% find number of cells
no_cells = sum(~cellfun('isempty',cells));

% initialise arrays
cells_per_node = zeros(array_sizes,1);
cell_store = zeros(array_sizes,10);

% loop over all cells
for current_cell = 1:no_cells
	
	% find cell nodes associated with current cell
	cell_nodes = cells{current_cell};
	no_cell_nodes = length(cell_nodes);
	
	% loop over cell nodes
	for i = 1:no_cell_nodes
		
		% update cells_per_node and cell_store with the current node
		current_node = cell_nodes(i);
		cells_per_node(current_node) = cells_per_node(current_node) + 1;
		cell_store(current_node,cells_per_node(current_node)) = current_cell;
		
	end
end