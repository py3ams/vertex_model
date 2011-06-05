function node_positions = update_positions(array_sizes,cells_fixed_logical,...
	cells,delta_t,iteration,k_elas,node_positions,...
	nodes_fixed_logical,pressure_constant,target_area,time,...
	time_of_last_division,viscosity)

% initialise variables
forces = zeros(array_sizes,2);
previous_no_cell_nodes = 0;

% find number of cells
no_cells = sum(~cellfun('isempty',cells));

% loop over all cells
for current_cell = 1:no_cells
	
% 	if ~cells_fixed_logical(current_cell)
		
		% find nodes associated with current cell
		cell_nodes = cells{current_cell};
		no_cell_nodes = length(cell_nodes);
		
		% find node positions in 3D
		cell_node_positions = node_positions(cell_nodes,:);
		
		cell_centre = mean(cell_node_positions);
		
		% calculate cell area
		cell_area =...
			polyarea(cell_node_positions(:,1),cell_node_positions(:,2));
				
		% calculate edge lengths
		cell_edge_lengths =...
			calculate_edge_lengths(cell_node_positions);
		
		% calculate perimeter
		cell_perimeter = sum(cell_edge_lengths);
		
		% find tension and pressure normals. tension_unit_normals is
		% 2 x 2 x no_cell_nodes, whereas pressure_unit_normals is no_cell_nodes x 2.
		[tension_unit_normals,pressure_unit_normals] =...
			calculate_normals(cell_node_positions,no_cell_nodes);
		
		% create the reshape matrix to reshap variables into 2 x 2 x no_cell_nodes
		% matrices, so they can be multplied by tension unit normals.
		if no_cell_nodes ~= previous_no_cell_nodes
			
			clear reshape_matrix
			
			reshape_matrix(1,:,:) = [1:no_cell_nodes; 1:no_cell_nodes];
			reshape_matrix(2,:,:) = [[no_cell_nodes 1:no_cell_nodes-1];...
				[no_cell_nodes 1:no_cell_nodes-1]];
			
		end
		
		% reshape edge_lengths, migration_tension_factors, and tension_anisotropy
		edge_lengths_matrix = cell_edge_lengths(reshape_matrix);
		
		% calculate tension force acting at each node (a no_cell_nodes x 2 matrix)
		tension_force = cell_perimeter*k_elas*squeeze(sum(edge_lengths_matrix.*...
			tension_unit_normals))';
				
		pressure_force =...
			(1 + growth_function(cell_area,cell_centre,target_area,time,...
			time_of_last_division(current_cell)))*...
			pressure_unit_normals*pressure_constant/cell_area;
		
% 		if current_cell == 36
% 			display(pressure_force)
% 		end
		
		% calculate total force acting at each cell node, and add it to the total
		% force vector
		forces(cell_nodes,:) =...
			forces(cell_nodes,:) + pressure_force + tension_force;
		
		% store previous number of cell nodes (so we don't have to reformulate the
		% reshape matrix if it hasn't changed
		previous_no_cell_nodes = no_cell_nodes;
		
% 	end
	
end

node_positions(~nodes_fixed_logical,:) =...
	node_positions(~nodes_fixed_logical,:) +...
	delta_t/viscosity*forces(~nodes_fixed_logical,:);

% node_positions = node_positions + delta_t/viscosity*forces;

global total_force
total_force(iteration) = sum(sqrt(sum(forces.^2,2)));
