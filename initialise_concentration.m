function concentration = initialise_concentration(cells_per_vertex,...
	FEM_nodes,gradient_type,initial_concentration_magnitude,no_cells,...
	no_chemicals,source_width)

no_FEM_nodes = length(FEM_nodes.position);
concentration = zeros(no_FEM_nodes,1);

for current_chemical = 1:no_chemicals
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%% gradient type 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if gradient_type(current_chemical) == 3
		
		% first loop is over nodes that are also cell vertices
		for current_node_local = 1:length(cells_per_vertex)
			
			if cells_per_vertex(current_node_local) > 0 && ...
					sqrt((FEM_nodes.position(current_node_local,1).^2)+...
					(FEM_nodes.position(current_node_local,2).^2))<(0.5*source_width(current_chemical))
				
				concentration(current_node_local,current_chemical) = initial_concentration_magnitude(current_chemical);
				
			end
			
		end
		
		% second loop is over nodes on edges between vertices
		for current_node_local = length(cells_per_vertex)+1:(no_FEM_nodes-no_cells-1)
			
			if FEM_nodes.edge(current_node_local,1)>0 && ...
					sqrt((FEM_nodes.position(current_node_local,1).^2)+...
					(FEM_nodes.position(current_node_local,2).^2))<(0.5*source_width(current_chemical))
				
				concentration(current_node_local,current_chemical) = initial_concentration_magnitude(current_chemical);
				
			end
			
		end
		
		% second loop is over nodes at cell centroids
		for current_node_local = (no_FEM_nodes-no_cells):no_FEM_nodes
			
			if sqrt((FEM_nodes.position(current_node_local,1).^2)+...
					(FEM_nodes.position(current_node_local,2).^2))<(0.5*source_width(current_chemical))
				
				concentration(current_node_local,current_chemical) = initial_concentration_magnitude(current_chemical);
				
				
			end
		end
		
	%%%%%%%%%%%%%%%%%%%%%%%%%%%% gradient type 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
	elseif gradient_type(current_chemical) == 1
		
		for current_node_local = 1:length(cells_per_vertex)
			
			if cells_per_vertex(current_node_local) > 0 && ...
					abs(FEM_nodes.position(current_node_local,1)) < (0.5*source_width(current_chemical))
				
				concentration(current_node_local,current_chemical) = initial_concentration_magnitude(current_chemical);
				
			end
			
		end
		
		for current_node_local = length(cells_per_vertex)+1:(no_FEM_nodes-no_cells-1)
			
			if FEM_nodes.edge(current_node_local,1)>0 && ...
					abs(FEM_nodes.position(current_node_local,1)) < (0.5*source_width(current_chemical))
				
				concentration(current_node_local,current_chemical) = initial_concentration_magnitude(current_chemical);
				
			end
			
		end
		
		for current_node_local = (no_FEM_nodes-no_cells):no_FEM_nodes
			
			if abs(FEM_nodes.position(current_node_local,1)) < (0.5*source_width(current_chemical))
				
				concentration(current_node_local,current_chemical) = initial_concentration_magnitude(current_chemical);
				
			end
			
		end
		
	%%%%%%%%%%%%%%%%%%%%%%%%%%%% gradient type 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
		
	elseif gradient_type(current_chemical) == 2
		
		for current_node_local = 1:length(cells_per_vertex)
			
			if cells_per_vertex(current_node_local) > 0 && ...
					abs(FEM_nodes.position(current_node_local,2)) < (0.5*source_width(current_chemical))
				
				concentration(current_node_local,current_chemical) = initial_concentration_magnitude(current_chemical);
				
			end
			
		end
		
		for current_t_node_local = length(cells_per_vertex)+1:(no_FEM_nodes-no_cells-1)
			
			if FEM_nodes.edge(current_node_local,1)>0 && ...
					abs(FEM_nodes.position(current_node_local,2)) < (0.5*source_width(current_chemical))
				
				concentration(current_node_local,current_chemical) = initial_concentration_magnitude(current_chemical);
				
			end
			
		end
		
		for current_node_local = (no_FEM_nodes-no_cells):no_FEM_nodes
			
			if abs(FEM_nodes.position(current_node_local,2)) < (0.5*source_width(current_chemical))
				
				concentration(current_node_local,current_chemical) = initial_concentration_magnitude(current_chemical);
				
			end
			
		end
		
	%%%%%%%%%%%%%%%%%%%%%%%%%%%% gradient type 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
		
	elseif gradient_type(current_chemical) == 4
		
		min_x_position = min(FEM_nodes.position(:,1));
		
		for current_node_local = 1:length(cells_per_vertex)
			
			if cells_per_vertex(current_node_local) > 0 && ...
					abs(FEM_nodes.position(current_node_local,1)-min_x_position) < source_width(current_chemical)
				
				concentration(current_node_local,current_chemical) = initial_concentration_magnitude(current_chemical);
				
			end
			
		end
		
		for current_node_local = length(cells_per_vertex)+1:(no_FEM_nodes-no_cells-1)
			
			if FEM_nodes.edge(current_node_local,1)>0 && ...
					abs(FEM_nodes.position(current_node_local,1)-min_x_position) < source_width(current_chemical)
				
				concentration(current_node_local,current_chemical) = initial_concentration_magnitude(current_chemical);
				
			end
			
		end
		
		for current_node_local = (no_FEM_nodes-no_cells):no_FEM_nodes
			
			if abs(FEM_nodes.position(current_node_local,1)-min_x_position) < source_width(current_chemical)
				
				concentration(current_node_local,current_chemical) = initial_concentration_magnitude(current_chemical);
				
			end
			
		end
		
	else
		
		error('Invalid gradient_type');
		
	end
	
end

