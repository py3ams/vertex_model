function [tension_unit_normals,pressure_unit_normals,internal_angles] =...
	calculate_normals(cell_node_positions,no_cell_nodes)
	
if nargin == 1
	no_cell_nodes = length(cell_node_positions);
end

% initialise output variables
concave_logical = false(no_cell_nodes,1);
internal_angles = zeros(no_cell_nodes,1);
tension_unit_normals = zeros(2,2,no_cell_nodes);
pressure_unit_normals = zeros(no_cell_nodes,2);

% loop over all cellnodes
for current_node_local = 1:no_cell_nodes
	
	current_node_position = cell_node_positions(current_node_local,:);
	
	% find clockwise and anti-clockwise node indices
	if current_node_local == 1
		
		clockwise_node_local = no_cell_nodes;
		anti_clockwise_node_local = current_node_local + 1;
		
	elseif current_node_local == no_cell_nodes
		
		clockwise_node_local = current_node_local - 1;
		anti_clockwise_node_local = 1;
		
	else
		
		clockwise_node_local = current_node_local - 1;
		anti_clockwise_node_local = current_node_local + 1;
		
	end
	
	anti_clockwise_node_position =...
		cell_node_positions(anti_clockwise_node_local,:);
	
	clockwise_node_position =...
		cell_node_positions(clockwise_node_local,:);
	
	unit_vector_to_anti_clockwise_node =...
		calculate_unit_vector(current_node_position,...
		anti_clockwise_node_position);
	
	unit_vector_to_clockwise_node =...
		calculate_unit_vector(current_node_position,...
		clockwise_node_position);
		
	% take cross product of vectors to anti-clockwise and clockwise nodes, and
	% use this to decide if the node is concave
	cross_product = cross([unit_vector_to_anti_clockwise_node 0],...
		[unit_vector_to_clockwise_node 0]);

	if cross_product(3) < 0
		concave_logical(current_node_local) = true;
	end	  
			  
	internal_angles(current_node_local) = acos(dot(...
		unit_vector_to_anti_clockwise_node,unit_vector_to_clockwise_node));
	
	% store the two unit vectors to neighbouring nodes as the tension
	% normals
	current_tension_normals =...
		[unit_vector_to_anti_clockwise_node; unit_vector_to_clockwise_node];
	
	sum_current_tension_normals = sum(current_tension_normals);
	
	current_pressure_normal =...
		-sum_current_tension_normals/sqrt(sum(sum_current_tension_normals.^2));
		
	% store the pressure normals
	pressure_unit_normals(current_node_local,:) = current_pressure_normal;
	tension_unit_normals(:,:,current_node_local) = current_tension_normals;
	
end

internal_angles(concave_logical) = 2*pi - internal_angles(concave_logical);

pressure_unit_normals(concave_logical,:) =...
	-pressure_unit_normals(concave_logical,:);