function node_positions = apply_boundary_force(boundary_constant,...
	boundary_element,delta_t,node_positions,viscosity)

no_boundary_nodes = length(boundary_element);
boundary_element_positions = node_positions(boundary_element,:);

boundary_force = zeros(no_boundary_nodes,2);

ideal_external_angle = pi+2*pi/no_boundary_nodes;

for current_node=1:no_boundary_nodes
	
	anti_clockwise_node = mod(current_node,no_boundary_nodes)+1;
	clockwise_node = mod(current_node-2,no_boundary_nodes)+1;
	
	current_node_position = boundary_element_positions(current_node,:);
	anti_clockwise_node_position = boundary_element_positions(anti_clockwise_node,:);
	clockwise_node_position = boundary_element_positions(clockwise_node,:);
	
	unit_vector_to_anti_clockwise_node =...
		calculate_unit_vector(current_node_position,...
		anti_clockwise_node_position);
	
	unit_vector_to_clockwise_node =...
		calculate_unit_vector(current_node_position,...
		clockwise_node_position);
	
	current_tension_normals =...
		[unit_vector_to_anti_clockwise_node; unit_vector_to_clockwise_node];
	
	sum_current_tension_normals = sum(current_tension_normals);
	
	internal_angle = acos(dot(...
		unit_vector_to_anti_clockwise_node,unit_vector_to_clockwise_node));
		
	if unit_vector_to_anti_clockwise_node(1)*unit_vector_to_clockwise_node(2)-...
			unit_vector_to_clockwise_node(1)*unit_vector_to_anti_clockwise_node(2) < 0

		external_angle = internal_angle;
		current_pressure_normal =...
			sum_current_tension_normals/sqrt(sum(sum_current_tension_normals.^2));
	
	else
		
		external_angle = 2*pi-internal_angle;
		current_pressure_normal =...
			-sum_current_tension_normals/sqrt(sum(sum_current_tension_normals.^2));
		
	end	  
	
	boundary_force(current_node,:) =...
		boundary_constant.*(ideal_external_angle-external_angle)^3.*...
		current_pressure_normal;
	
end

node_positions(boundary_element,:) =...
	node_positions(boundary_element,:) + delta_t/viscosity*boundary_force;
