function concentration = initialise_concentration(cells_per_vertex,...
	FEM_node_positions,gradient_type,initial_concentration_magnitude,no_cells,...
	no_chemicals,source_width)

no_FEM_nodes = length(FEM_node_positions);
concentration = zeros(no_FEM_nodes,1);

for current_chemical = 1:no_chemicals
	
	if gradient_type(current_chemical) == 3
		
		% first loop is over nodes that are also cell vertices
		for current_node_local = 1:length(cells_per_vertex)
			
			if cells_per_vertex(current_node_local) > 0 && ...
					sqrt((FEM_node_positions(current_node_local,1).^2)+...
					(FEM_node_positions(current_node_local,2).^2))<(0.5*source_width(current_chemical))
				
				concentration(current_node_local,current_chemical) = initial_concentration_magnitude(current_chemical);
				
			end
			
		end
		
		% second loop is over nodes at cell centroids
		for current_node_local = (no_FEM_nodes-no_cells):no_FEM_nodes
			
			if sqrt((FEM_node_positions(current_node_local,1).^2)+...
					(FEM_node_positions(current_node_local,2).^2))<(0.5*source_width(current_chemical))
				
				concentration(current_node_local,current_chemical) = initial_concentration_magnitude(current_chemical);
				
				
			end
		end
		
	elseif gradient_type(current_chemical) == 1
		
		for current_node_local = 1:length(cells_per_vertex)
			
			if cells_per_vertex(current_node_local) > 0 && ...
					abs(FEM_node_positions(current_node_local,1)) < (0.5*source_width(current_chemical))
				
				concentration(current_node_local,current_chemical) = initial_concentration_magnitude(current_chemical);
				
			end
			
		end
		
		for current_node_local = (no_FEM_nodes-no_cells):no_FEM_nodes
			
			if abs(FEM_node_positions(current_node_local,1)) < (0.5*source_width(current_chemical))
				
				concentration(current_node_local,current_chemical) = initial_concentration_magnitude(current_chemical);
				
			end
			
		end
		
	elseif gradient_type(current_chemical) == 2
		
		for current_node_local = 1:length(cells_per_vertex)
			
			if cells_per_vertex(current_node_local) > 0 && ...
					abs(FEM_node_positions(current_node_local,2)) < (0.5*source_width(current_chemical))
				
				concentration(current_node_local,current_chemical) = initial_concentration_magnitude(current_chemical);
				
			end
			
		end
		
		for current_node_local = (no_FEM_nodes-no_cells):no_FEM_nodes
			
			if abs(FEM_node_positions(current_node_local,2)) < (0.5*source_width(current_chemical))
				
				concentration(current_node_local,current_chemical) = initial_concentration_magnitude(current_chemical);
				
			end
			
		end
		
	elseif gradient_type(current_chemical) == 4
		
		min_x_position = min(FEM_node_positions(:,1));
		
		for current_node_local = 1:length(cells_per_vertex)
			
			if cells_per_vertex(current_node_local) > 0 && ...
					abs(FEM_node_positions(current_node_local,1)-min_x_position) < source_width(current_chemical)
				
				concentration(current_node_local,current_chemical) = initial_concentration_magnitude(current_chemical);
				
			end
			
		end
		
		for current_node_local = (no_FEM_nodes-no_cells):no_FEM_nodes
			
			if abs(FEM_node_positions(current_node_local,1)-min_x_position) < source_width(current_chemical)
				
				concentration(current_node_local,current_chemical) = initial_concentration_magnitude(current_chemical);
				
			end
			
		end
		
	else
		
		error('Invalid gradient_type');
		
	end
	
end

