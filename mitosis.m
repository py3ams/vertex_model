function [cells,vertices,cell_growth_speeds_matrix,FEM_elements,FEM_nodes,...
	mitosis_counter,refined_edge_matrix,stats] = mitosis(cells,vertices,...
	cell_growth_speeds_matrix,delta_t,FEM_elements,FEM_nodes,...
	FEM_solve_logical,medial_lateral_threshold,mitosis_angles_type,...
	mitosis_counter,mitosis_dependence,mitosis_period,mitosis_random_logical,...
	mitosis_threshold,refined_edge_matrix,stats,time)

mitosis_location = 0;

if strcmp(mitosis_dependence,'area')
	
	division_probabilities = division_probability(cells.area,...
		cells.vertices,delta_t,mitosis_random_logical,mitosis_threshold,...
		target_area);
	
elseif strcmp(mitosis_dependence,'volume')
	
	division_probabilities = division_probability(cells.volume,...
		cells.vertices,delta_t,mitosis_random_logical,...
		mitosis_threshold,cells.target_volume);
	
elseif strcmp(mitosis_dependence,'none')
	
	division_probabilities = zeros(length(cells.vertices),1);
	
	% using the rem function doesn't work so well because as soon as time
	% is actually slightly lower than its true value (i.e 0.9999999 instead
	% of 1) the remainder becomes large
	if abs(round(time/mitosis_period)-time/mitosis_period)<1e-6
	
% 	if rem(time,mitosis_period)<1e-12
		
		% all cells will be included in cells_to_divide. the code will then
		% pick one of them randomly to actually become cell_to_divide.
		division_probabilities = ones(length(cells.vertices),1);
		
% 		% this may result in more than one cell being chosen to have its
% 		% division probability set to one, if cells have the same volume.
% 		% however the code later on ensures that only one is chosen to
% 		% divide at this iteration.
% 		[~,max_volume_indices] = max(cells.volume);
% 		division_probabilities(max_volume_indices) = 1;        
% 
	end
	
end

no_cells = length(cells.vertices);

% cells_to_divide = find((rand(no_cells,1) < division_probabilities)&~cells.mitosis_flag);
% this way cells can divide more than once, but only once they are no longer growing,
% and only if they are not apoptotic or dead.
cells_to_divide = find((rand(no_cells,1) < division_probabilities)&cells.state==1);

% cells_to_divide = 25;

no_cells_to_divide = length(cells_to_divide);

if no_cells_to_divide > 0
	
	% pick a random cell from those that have been chosen to divide. we only
	% allow one cell per iteration to divide. if mitosis_random_logical
	% is 0 then cells_to_divide should be of length one in most
	% circumstances, apart from possibly at the beginning of a simulation.
	% it is possible, though unlikely, that two cells go over the mitosis
	% threshold at the same iteration. in the case of stochastic mitosis it
	% is quite possible for cells_to_divide to contain more than one cell.
	cell_to_divide = cells_to_divide(ceil(rand*no_cells_to_divide));
	
%     cell_to_divide = 25;
    
	cell_vertices = cells.vertices{cell_to_divide};
	no_cell_vertices = length(cell_vertices);
	cell_vertex_positions = vertices.position(cell_vertices,:);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% disp(cell_vertex_positions);
	
	% figure
	% patchAS(cell_vertex_positions)
	% hold on
	% axis equal
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	new_cell_1 = cell_to_divide;
	new_cell_2 = no_cells+1;
	
	unused_vertices = find(vertices.no_cells == 0,2);
	new_vertex_1 = unused_vertices(1);
	new_vertex_2 = unused_vertices(2);
	
	cell_to_divide_centroid =...
		CalculateCentroid(cell_vertex_positions,cells.area(cell_to_divide));
	
	%     mitosis_location(1,1) = cell_to_divide_centroid(1)/max(vertices.position(:,1));
	%     mitosis_location(1,2) = cell_to_divide_centroid(2)/max(vertices.position(:,2));
	
	mitosis_location(1,1) = cell_to_divide_centroid(1);
	mitosis_location(1,2) = cell_to_divide_centroid(2);
	
	if strcmp(mitosis_angles_type,'uniform')
		angle_of_mitosis = pi*(rand-0.5);
	else
		angle_of_mitosis = mitosis_angles_type(1) +...
			max(mitosis_angles_type(2),1e-9)*randn;
	end
	
	
	% find gradient and y-intercept of the mitosis line, which passes
	% through the cell centroid and is at the prescribed angle. the line
	% takes the form y=mx+c.
	mitosis_line_gradient = tan(angle_of_mitosis);
	mitosis_line_intercept =...
		cell_to_divide_centroid(2)-mitosis_line_gradient*cell_to_divide_centroid(1);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% x = linspace(min(cell_vertex_positions(:,1)),max(cell_vertex_positions(:,1)),100);
	% y = mitosis_line_gradient*x + mitosis_line_intercept;
	% plot(x,y,'b','linewidth',3)
	% plot(cell_to_divide_centroid(1),cell_to_divide_centroid(2),'go','linewidth',3)
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	intersection_counter = 0;
	
	% loop over cell vertices to find which edges intersect with the mitosis
	% line
	for current_vertex_local = 1:no_cell_vertices
		
		clockwise_vertex_local = mod(current_vertex_local,no_cell_vertices)+1;
		
		current_vertex_position = cell_vertex_positions(current_vertex_local,:);
		clockwise_vertex_position = cell_vertex_positions(clockwise_vertex_local,:);
		
		% find gradient and y-intercept of the current edge
		current_edge_gradient =...
			(clockwise_vertex_position(2)-current_vertex_position(2))/...
			(clockwise_vertex_position(1)-current_vertex_position(1));
		
		current_edge_intercept =...
			current_vertex_position(2)-current_edge_gradient*current_vertex_position(1);
		
		% find the x-value of the intersection between the current edge and
		% the mitosis line. if this intersection takes place between the
		% min and max x-values of the edge, the mitosis line intersects the
		% edge
		x_intersection =...
			(current_edge_intercept-mitosis_line_intercept)/...
			(mitosis_line_gradient-current_edge_gradient);
		
		%         x_intersection
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		% 	x = linspace(current_vertex_position(1),clockwise_vertex_position(1),100);
		% 	y = current_edge_gradient*x + current_edge_intercept;
		%
		% 	plot(x,y,'y','linewidth',3)
		
		% 	y_intersection = current_edge_gradient*x_intersection + current_edge_intercept;
		% 	y_intersection2 = mitosis_line_gradient*x_intersection + mitosis_line_intercept;
		%
		% 	plot(x_intersection,y_intersection,'ko','linewidth',3)
		% 	plot(x_intersection,y_intersection2,'co','linewidth',3)
		
		% 	disp(x_intersection)
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		x_max = max(current_vertex_position(1),clockwise_vertex_position(1));
		x_min = min(current_vertex_position(1),clockwise_vertex_position(1));
		
		if x_intersection > x_min && x_intersection < x_max
			
			intersection_counter = intersection_counter+1;
			
			if intersection_counter == 1
				
				intersection_vertex_1_local = current_vertex_local;
				intersection_vertex_1_global = cell_vertices(current_vertex_local);
				
				clockwise_intersection_vertex_1_local = clockwise_vertex_local;
				clockwise_intersection_vertex_1_global = cell_vertices(clockwise_vertex_local);
				
				new_vertex_1_position =...
					[x_intersection mitosis_line_gradient*x_intersection+mitosis_line_intercept];
				
				distance_intersection_vertex_1_new_vertex_1 =...
					calculate_distance_between_vertices(current_vertex_position,new_vertex_1_position);
				
			elseif intersection_counter == 2
				
				intersection_vertex_2_local = current_vertex_local;
				intersection_vertex_2_global = cell_vertices(current_vertex_local);
				
				clockwise_intersection_vertex_2_local = clockwise_vertex_local;
				clockwise_intersection_vertex_2_global = cell_vertices(clockwise_vertex_local);
				
				new_vertex_2_position =...
					[x_intersection mitosis_line_gradient*x_intersection+mitosis_line_intercept];
				
				distance_intersection_vertex_2_new_vertex_2 =...
					calculate_distance_between_vertices(current_vertex_position,new_vertex_2_position);
				
				break;
				
			end
		end
	end
	
	if intersection_counter < 2
		
		error(['warning: cannot find intercepts for mitosis of cell ',...
			num2str(cell_to_divide),' at angle ',num2str(angle_of_mitosis)]);
		
	else
		
		vertices.position(new_vertex_1,:) = new_vertex_1_position;
		vertices.position(new_vertex_2,:) = new_vertex_2_position;
		
		vertices.time_created(new_vertex_1) = time;
		vertices.time_created(new_vertex_2) = time;
		
		% figure out which vertices go in each of the two new cells. we
		% know that intersection_vertex_1 comes before clockwise_intersection_vertex_1,
		% which comes before intersection_vertex_2. the only complication
		% is clockwise_intersection_vertex_2, which could potentially be
		% the first cell vertex, if intersection_vertex_2 is the last cell
		% vertex. the rest is just bookkeeping
		new_cell_1_vertices_global =...
			[new_vertex_1 cell_vertices(clockwise_intersection_vertex_1_local:...
			intersection_vertex_2_local) new_vertex_2];
		
		if intersection_vertex_2_local == no_cell_vertices
			new_cell_vertices_2_global =...
				[new_vertex_2 cell_vertices(1:intersection_vertex_1_local) new_vertex_1];
		else
			new_cell_vertices_2_global =...
				[new_vertex_2 cell_vertices(clockwise_intersection_vertex_2_local:end) ...
				cell_vertices(1:intersection_vertex_1_local) new_vertex_1];
		end
		
		cells.vertices{new_cell_1} = new_cell_1_vertices_global;
		cells.vertices{new_cell_2} = new_cell_vertices_2_global;
		
		cells.target_volume(new_cell_2) = cells.target_volume(new_cell_1);
		
		cells.original_logical(new_cell_2) = 0;
		cells.state(new_cell_1) = 2;
		cells.state(new_cell_2) = 2;
		
		cell_area_new_cells =...
			CalculateCellAreas(cells.vertices([new_cell_1,new_cell_2]),vertices.position);
		
		cells.area(new_cell_1) = cell_area_new_cells(1);
		cells.area(new_cell_2) = cell_area_new_cells(2);
		
		% target areas are initially set to be equal to the current area.
		% there should be code elsewhere that allows the target area of
		% 'baby' cells to increase over time
		cells.target_area(new_cell_1) = cell_area_new_cells(1);
		cells.target_area(new_cell_2) = cell_area_new_cells(2);
		
		cell_area_cell_to_divide = sum(cell_area_new_cells);
		cell_volume_cell_to_divide = cells.volume(cell_to_divide);
		
		% share the volume of the mother cell between the two daughter
		% cells based on the ratio of their areas.
		cells.volume(new_cell_1) =...
			cell_area_new_cells(1)/cell_area_cell_to_divide*cell_volume_cell_to_divide;
		cells.volume(new_cell_2) =...
			cell_area_new_cells(2)/cell_area_cell_to_divide*cell_volume_cell_to_divide;
		
		cells.force_constants.area(new_cell_2) = cells.force_constants.area(new_cell_1);
		cells.force_constants.deformation(new_cell_2) = cells.force_constants.deformation(new_cell_1);
		cells.force_constants.elongation(new_cell_2) = cells.force_constants.elongation(new_cell_1);
		cells.force_constants.perimeter(new_cell_2) = cells.force_constants.perimeter(new_cell_1);
		cells.force_constants.tension(new_cell_2) = cells.force_constants.tension(new_cell_1);
		
		% change time of last division for each of the new cells to the current time
		cells.time_of_last_division(new_cell_1) = time;
		cells.time_of_last_division(new_cell_2) = time;
		
		cells.mitosis_flag(new_cell_1) = 1;
		cells.mitosis_flag(new_cell_2) = 1;
		
		% centroid are used mainly in for FEM purposes but also to
		% determine whether cell is medial or lateral
		centroid_new_cell_1 = CalculateCentroid(vertices.position(...
			new_cell_1_vertices_global,:),cell_area_new_cells(1));
		
		centroid_new_cell_2 = CalculateCentroid(vertices.position(...
			new_cell_vertices_2_global,:),cell_area_new_cells(2));
		
		% this code is probably a bit obsolete - should consider getting
		% rid of cell_growth_speeds_matrix. it assigns the growths speeds
		% depending on whether they are medial or lateral. the matrix is
		% created before simulation starts.
		if abs(centroid_new_cell_1(1)) < medial_lateral_threshold
			
			new_medial_cell_growth_speed_index = find(cell_growth_speeds_matrix(:,1),1);
			
			cells.growth_speed(new_cell_1) =...
				cell_growth_speeds_matrix(new_medial_cell_growth_speed_index,1);
			
			cell_growth_speeds_matrix(new_medial_cell_growth_speed_index,1) = 0;
			
		else
			
			new_lateral_cell_growth_speed_index = find(cell_growth_speeds_matrix(:,2),1);
			
			cells.growth_speed(new_cell_1) =...
				cell_growth_speeds_matrix(new_lateral_cell_growth_speed_index,2);
			
			cell_growth_speeds_matrix(new_lateral_cell_growth_speed_index,2) = 0;
			
		end
		
		if abs(centroid_new_cell_2(1)) < medial_lateral_threshold
			
			new_medial_cell_growth_speed_index = find(cell_growth_speeds_matrix(:,1),1);
			
			cells.growth_speed(new_cell_2) =...
				cell_growth_speeds_matrix(new_medial_cell_growth_speed_index,1);
			
			cell_growth_speeds_matrix(new_medial_cell_growth_speed_index,1) = 0;
			
		else
			
			new_lateral_cell_growth_speed_index = find(cell_growth_speeds_matrix(:,2),1);
			
			cells.growth_speed(new_cell_2) =...
				cell_growth_speeds_matrix(new_lateral_cell_growth_speed_index,2);
			
			cell_growth_speeds_matrix(new_lateral_cell_growth_speed_index,2) = 0;
			
		end
		
		% we need to update the cells associated with the old cell vertices that are
		% now in new_cell_2. the vertices that are in new_cell_1 don't need
		% to be updated as they will already be associated with
		% cell_to_divide, which is the same as new_cell_1. we will then do the
		% new vertices, i.e. those that make up the edge of division,
		% separately, as they currently don't have any cells associated
		% with them.
		old_cell_vertices_now_in_cell_2 =...
			cell_vertices(ismember(cell_vertices,new_cell_1_vertices_global)==0);
		
		for temp_vertex_local = 1:length(old_cell_vertices_now_in_cell_2)
			
			temp_vertex_global = old_cell_vertices_now_in_cell_2(temp_vertex_local);
			cells_temp_vertex_global = vertices.cells(temp_vertex_global,:);
			cells_temp_vertex_global(cells_temp_vertex_global ==...
				cell_to_divide) = new_cell_2;
			vertices.cells(temp_vertex_global,:) = cells_temp_vertex_global;
			
		end
		
		vertices.no_cells(new_vertex_1) = 2;
		vertices.no_cells(new_vertex_2) = 2;
		
		vertices.cells(new_vertex_1,1:2) = [new_cell_1 new_cell_2];
		vertices.cells(new_vertex_2,1:2) = [new_cell_1 new_cell_2];
		
		% find out if either of the edges containing the new vertices are shared with other
		% cells. in most cases this will be true, unless the cell is on a
		% boundary
		[~,~,cell_with_same_edge_1] =...
			find_cells_containing_vertices(vertices.cells,cell_to_divide,...
			intersection_vertex_1_global,clockwise_intersection_vertex_1_global);
		
		[~,~,cell_with_same_edge_2] =...
			find_cells_containing_vertices(vertices.cells,cell_to_divide,...
			intersection_vertex_2_global,clockwise_intersection_vertex_2_global);
		
		% if the edges containing new vertices are shared with other cells, need to
		% add the new vertices into those cells at the correct position, and also add
		% cells to vertices.no_cells and vertices.cells for the new vertices
		if ~isempty(cell_with_same_edge_1)
			
			vertices.no_cells(new_vertex_1) = vertices.no_cells(new_vertex_1) + 1;
			vertices.cells(new_vertex_1,vertices.no_cells(new_vertex_1)) = cell_with_same_edge_1;
			
			cell_vertices_cell_with_same_edge_1 = cells.vertices{cell_with_same_edge_1};
			
			index_1 = find(cell_vertices_cell_with_same_edge_1 ==...
				clockwise_intersection_vertex_1_global);
			
			cell_vertices_cell_with_same_edge_1 =...
				[cell_vertices_cell_with_same_edge_1(1:index_1) new_vertex_1 ...
				cell_vertices_cell_with_same_edge_1(index_1+1:end)];
			
			cells.vertices{cell_with_same_edge_1} = cell_vertices_cell_with_same_edge_1;
			
		end
		
		if ~isempty(cell_with_same_edge_2)
			
			vertices.no_cells(new_vertex_2) = vertices.no_cells(new_vertex_2) + 1;
			vertices.cells(new_vertex_2,vertices.no_cells(new_vertex_2)) = cell_with_same_edge_2;
			
			cell_vertices_cell_with_same_edge_2 = cells.vertices{cell_with_same_edge_2};
			
			index_2 = find(cell_vertices_cell_with_same_edge_2 == ...
				clockwise_intersection_vertex_2_global);
			
			cell_vertices_cell_with_same_edge_2 =...
				[cell_vertices_cell_with_same_edge_2(1:index_2) new_vertex_2 ...
				cell_vertices_cell_with_same_edge_2(index_2+1:end)];
			
			cells.vertices{cell_with_same_edge_2} = cell_vertices_cell_with_same_edge_2;
			
		end
		
		% it is obviously important to do this step after we have set
		% vertices.no_cells etc for the new vertices.
		cells.boundary_logical(new_cell_1) = any(vertices.no_cells(cells.vertices{new_cell_1})<3);
		cells.boundary_logical(new_cell_2) = any(vertices.no_cells(cells.vertices{new_cell_2})<3);
		
		if FEM_solve_logical
			
			% this is the index within FEM nodes of the first cell-centre
			% node
			FEM_index_zero_cell = length(FEM_nodes.previous_position)-no_cells;
			
			% we must first remove any FEM_nodes from the intersection
			% edges. these are going to be replaced by new nodes corresponding
			% to the new vertices. it would be preferable to actually use the
			% concentration values of these edge FEM_nodes, if they exist, in
			% calculating the values at the new vertices.
			edge_1 = [intersection_vertex_1_global clockwise_intersection_vertex_1_global];
			
			[cells,FEM_elements,FEM_nodes,refined_edge_matrix] =...
				remove_FEM_node_from_specific_edge(cells,edge_1,FEM_elements,...
				FEM_nodes,refined_edge_matrix,vertices);
			
			edge_2 = [intersection_vertex_2_global clockwise_intersection_vertex_2_global];
			
			[cells,FEM_elements,FEM_nodes,refined_edge_matrix] =...
				remove_FEM_node_from_specific_edge(cells,edge_2,FEM_elements,...
				FEM_nodes,refined_edge_matrix,vertices);
			
			FEM_nodes.previous_position(new_vertex_1,:) = new_vertex_1_position;
			FEM_nodes.previous_position(new_vertex_2,:) = new_vertex_2_position;
			
			% find the values of concentration at the new vertices. this is a linear
			% interpolation between the values at the current vertices, based on
			% the ratio of the distance along the edge to the total edge length
			FEM_nodes.concentration(new_vertex_1) = FEM_nodes.concentration(intersection_vertex_1_global) +...
				distance_intersection_vertex_1_new_vertex_1/cells.edge_lengths{cell_to_divide}...
				(intersection_vertex_1_local)*(FEM_nodes.concentration(clockwise_intersection_vertex_1_global)-...
				FEM_nodes.concentration(intersection_vertex_1_global));
			
			FEM_nodes.concentration(new_vertex_2) = FEM_nodes.concentration(intersection_vertex_2_global) +...
				distance_intersection_vertex_2_new_vertex_2/cells.edge_lengths{cell_to_divide}...
				(intersection_vertex_2_local)*(FEM_nodes.concentration(clockwise_intersection_vertex_2_global)-...
				FEM_nodes.concentration(intersection_vertex_2_global));
			
			cells.internal_chemical_value(new_cell_2) = cells.internal_chemical_value(cell_to_divide);
			
% 			cells.ingestion_rate(new_cell_2) = cells.ingestion_rate(cell_to_divide);
			
			cells.internal_chemical_quantity(new_cell_1) =...
				cells.internal_chemical_value(new_cell_1)*cell_area_new_cells(1);
			cells.internal_chemical_quantity(new_cell_2) =...
				cells.internal_chemical_value(new_cell_2)*cell_area_new_cells(2);
			
			cells.maximum_internal_chemical_quantity(new_cell_2) = cells.maximum_internal_chemical_quantity(cell_to_divide);
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			% text(vertices.position(intersection_vertex_1_global,1),...
			% 	vertices.position(intersection_vertex_1_global,2),...
			% 	num2str(FEM_nodes.concentration(intersection_vertex_1_global)));
			%
			% text(vertices.position(clockwise_intersection_vertex_1_global,1),...
			% 	vertices.position(clockwise_intersection_vertex_1_global,2),...
			% 	num2str(FEM_nodes.concentration(clockwise_intersection_vertex_1_global)));
			%
			% text(new_vertex_1_position(1),new_vertex_1_position(2),num2str(FEM_nodes.concentration(new_vertex_1)));
			%
			% text(vertices.position(intersection_vertex_2_global,1),...
			% 	vertices.position(intersection_vertex_2_global,2),...
			% 	num2str(FEM_nodes.concentration(intersection_vertex_2_global)));
			%
			% text(vertices.position(clockwise_intersection_vertex_2_global,1),...
			% 	vertices.position(clockwise_intersection_vertex_2_global,2),...
			% 	num2str(FEM_nodes.concentration(clockwise_intersection_vertex_2_global)));
			%
			% text(new_vertex_2_position(1),new_vertex_2_position(2),num2str(FEM_nodes.concentration(new_vertex_2)));
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			% 		total_Dpp_in_cell_to_divide =...
			% 			CalculateTotalDpp(FEM_nodes.concentration,FEM_elements.nodes(cells.FEM_elements{cell_to_divide},:),...
			% 			previous_FEM_vertices.position);
			
			% this is where we create the temporary FEM mesh, described in
			% the ALE paper, using the new vertex nodes and the mother cell
			% centre, to calculate the total concentration in the sections
			% of the mother cell corresponding to the daughter cells. we
			% currently ignore edge nodes, this could be added
			% fairly simply.
			
			% new cell 1
			
			no_cell_vertices_new_cell_1 = length(new_cell_1_vertices_global);
			temp_FEM_elements_new_cell_1 = zeros(no_cell_vertices_new_cell_1,3);
			
			for temp_current_vertex_local=1:no_cell_vertices_new_cell_1
				
				temp_current_vertex_global =...
					new_cell_1_vertices_global(temp_current_vertex_local);
				
				temp_clockwise_vertex_global = new_cell_1_vertices_global(mod(...
					temp_current_vertex_local,no_cell_vertices_new_cell_1)+1);
				
				temp_FEM_elements_new_cell_1(temp_current_vertex_local,:) =...
					[temp_current_vertex_global temp_clockwise_vertex_global ...
					FEM_index_zero_cell+cell_to_divide];
				
			end
			
			total_Dpp_in_new_cell_1 = CalculateTotalDpp(...
				FEM_nodes.concentration,temp_FEM_elements_new_cell_1,...
				FEM_nodes.previous_position);
			
			% new cell 2
			
			no_cell_vertices_new_cell_2 = length(new_cell_vertices_2_global);
			temp_FEM_elements_new_cell_2 = zeros(no_cell_vertices_new_cell_2,3);
			
			for temp_current_vertex_local=1:no_cell_vertices_new_cell_2
				
				temp_current_vertex_global =...
					new_cell_vertices_2_global(temp_current_vertex_local);
				
				temp_clockwise_vertex_global = new_cell_vertices_2_global(mod(...
					temp_current_vertex_local,no_cell_vertices_new_cell_2)+1);
				
				temp_FEM_elements_new_cell_2(temp_current_vertex_local,:) =...
					[temp_current_vertex_global temp_clockwise_vertex_global ...
					FEM_index_zero_cell+cell_to_divide];
				
			end
			
			total_Dpp_in_new_cell_2 = CalculateTotalDpp(FEM_nodes.concentration,temp_FEM_elements_new_cell_2,...
				FEM_nodes.previous_position);
			
			% 		if abs(total_Dpp_in_cell_to_divide-...
			% 				(total_Dpp_in_new_cell_1+total_Dpp_in_new_cell_2))>1e-9
			% 			error('Dpp not conserved during mitosis');
			% 		end
			
% 			% we are going to create a new FEM_node, which is an edge node
% 			% at the mid-point on the newly created edge. we will use the
% 			% concentration of the mother cell centroid. the cell centroid
% 			% sits somewhere along this edge by definition of the edge, but
% 			% is not necessarily at it's mid-point.
% 			total_no_vertices = length(vertices.position);
% 			new_FEM_node = find(FEM_nodes.previous_position(total_no_vertices+1:FEM_index_zero_cell,1)==0,1)+total_no_vertices;
% 			
% 			if isempty(new_FEM_node)
% 				error('no space for new FEM edge node');
% 			end
% 			
% 			FEM_nodes.concentration(new_FEM_node) = FEM_nodes.concentration(FEM_index_zero_cell+cell_to_divide);
% 			
% 			% There seems to be something of a choice here. We could use the position
% 			% of the old centroid for the new FEM node. This will lie somewhere along
% 			% the mitosis line, but not neccessarily at the mid-point. The ALE can
% 			% then take care of the movement at the next iteration. We could
% 			% alternatively make an approximation and put the new node at the
% 			% mid-point with the concentration from the old centroid. Doing it this
% 			% way is consistent with what we do in T1 swaps.
% 			% it is really important to do this step before calculating the
% 			% concentration values at the centroids of the daughter cells, which
% 			% obviously depend on the position of this node. this must also
% 			% be done before we set the previous_position of the centroid
% 			% of new cell 1, as this is going to re-use the same spot in
% 			% the array as the cell_to_divide.
% 			
% 			FEM_nodes.previous_position(new_FEM_node,:) = FEM_nodes.previous_position(FEM_index_zero_cell+cell_to_divide,:);
% 			
% 			%             FEM_nodes.previous_position(new_FEM_node,:) = ...
% 			%                 0.5*(FEM_nodes.previous_position(new_vertex_1,:) +...
% 			%                 FEM_nodes.previous_position(new_vertex_2,:));
% 			
% 			FEM_nodes.edge(new_FEM_node,:) = [new_vertex_1 new_vertex_2];
% 			refined_edge_matrix(new_vertex_1,new_vertex_2) = 1;
% 			refined_edge_matrix(new_vertex_2,new_vertex_1) = 1;
			
			% it is important to do this step after setting the previous
			% position of the new_FEM_node, as this uses the
			% previous_position of cell_to_divide, which we are going to
			% overwrite with centroid_new_cell_1 (as new_cell_1 =
			% cell_to_divide)
			FEM_nodes.previous_position(FEM_index_zero_cell+new_cell_1,:) =...
				centroid_new_cell_1;
			
			FEM_nodes.previous_position(FEM_index_zero_cell+new_cell_2,:) =...
				centroid_new_cell_2;
			
			% it should already be true that edge for the FEM_node of
			% cell_to_divide/new_cell_1 is [0 0], as a centroid can't be an
			% edge node. for new_cell_2, we are simply adding an extra row
			% to the end of the array.
			FEM_nodes.edge(FEM_index_zero_cell+new_cell_1,:) = 0;
			FEM_nodes.edge(FEM_index_zero_cell+new_cell_2,:) = 0;
			
			% we are going to create variables to contain, in
			% clockwise order, all the FEM nodes of the new cells,
			% excluding the centroid nodes. in the most simple case this
			% will be exactly the same as the vertices of the new cells.
			% the code is similar to that used to create, for example,
			% new_cell_1_vertices_global.
			% however, when FEM edge nodes are present, these must also be
			% included. these variables are used to work out the
			% concentration at the central nodes, and also to create the new
			% list of FEM elements after mitosis.
			% wouldn't it be easier to just loop over the new cell vertices
			% directly and find any edge nodes? would have to make sure we
			% check the edge from the last vertex to the first as well.
			% could use a mod() function to ensure we cover this last edge.
			% we have already stored the old centroid node as a new edge node,
			% so this would be picked up during this loop. this would
			% simplify the code and remove the special cases before and after the
			% for loops. the current code uses cell_vertices of the mother
			% cell, even though we have already worked out which vertices
			% are in each of the daughter cells.
			
			% This code is really unneccessary and redundant - look at
			% changing it!
			
			%%% cell 1 %%%
			
			FEM_nodes_without_centre_new_cell_1_global = zeros(1,2*no_cell_vertices_new_cell_1);
			FEM_nodes_without_centre_new_cell_1_global(1:2) = [new_vertex_1 clockwise_intersection_vertex_1_global];
			
			no_FEM_nodes_found_cell_1 = 2;
			
			for temp_index = clockwise_intersection_vertex_1_local:intersection_vertex_2_local-1
				
				node_1 = cell_vertices(temp_index);
				node_2 = cell_vertices(temp_index+1);
				
				% find index where the current edge is in FEM_nodes.edge.
				% this could be empty.
				temp_edge_node = find((FEM_nodes.edge(:,1)==node_1&FEM_nodes.edge(:,2)==node_2)|...
					(FEM_nodes.edge(:,1)==node_2&FEM_nodes.edge(:,2)==node_1));
				
				if isempty(temp_edge_node)
					
					no_FEM_nodes_found_cell_1 = no_FEM_nodes_found_cell_1+1;
					FEM_nodes_without_centre_new_cell_1_global(no_FEM_nodes_found_cell_1) = node_2;
					
				else
					
					no_FEM_nodes_found_cell_1 = no_FEM_nodes_found_cell_1+1;
					FEM_nodes_without_centre_new_cell_1_global(no_FEM_nodes_found_cell_1) = temp_edge_node;
					
					no_FEM_nodes_found_cell_1 = no_FEM_nodes_found_cell_1+1;
					FEM_nodes_without_centre_new_cell_1_global(no_FEM_nodes_found_cell_1) = node_2;
					
				end
				
			end
			
			no_FEM_nodes_found_cell_1 = no_FEM_nodes_found_cell_1+1;
			FEM_nodes_without_centre_new_cell_1_global(no_FEM_nodes_found_cell_1) = new_vertex_2;
			
% 			no_FEM_nodes_found_cell_1 = no_FEM_nodes_found_cell_1+1;
% 			FEM_nodes_without_centre_new_cell_1_global(no_FEM_nodes_found_cell_1) = new_FEM_node;
			
			% remove the zeros
			FEM_nodes_without_centre_new_cell_1_global =...
				FEM_nodes_without_centre_new_cell_1_global(1:no_FEM_nodes_found_cell_1);
			
			%%% cell 2 %%%
			
			FEM_nodes_without_centre_new_cell_2_global = zeros(1,2*no_cell_vertices_new_cell_2);
			
			if intersection_vertex_2_local == no_cell_vertices
				
				FEM_nodes_without_centre_new_cell_2_global(1:2) = [new_vertex_2 cell_vertices(1)];
				no_FEM_nodes_found_cell_2 = 2;
				
				for temp_index = 1:intersection_vertex_1_local-1
					
					node_1 = cell_vertices(temp_index);
					node_2 = cell_vertices(temp_index+1);
					
					temp_edge_node = find((FEM_nodes.edge(:,1)==node_1&FEM_nodes.edge(:,2)==node_2)|...
						(FEM_nodes.edge(:,1)==node_2&FEM_nodes.edge(:,2)==node_1));
					
					if isempty(temp_edge_node)
						
						no_FEM_nodes_found_cell_2 = no_FEM_nodes_found_cell_2+1;
						FEM_nodes_without_centre_new_cell_2_global(no_FEM_nodes_found_cell_2) = node_2;
						
					else
						
						no_FEM_nodes_found_cell_2 = no_FEM_nodes_found_cell_2+1;
						FEM_nodes_without_centre_new_cell_2_global(no_FEM_nodes_found_cell_2) = temp_edge_node;
						
						no_FEM_nodes_found_cell_2 = no_FEM_nodes_found_cell_2+1;
						FEM_nodes_without_centre_new_cell_2_global(no_FEM_nodes_found_cell_2) = node_2;
						
					end
					
				end
				
				no_FEM_nodes_found_cell_2 = no_FEM_nodes_found_cell_2+1;
				FEM_nodes_without_centre_new_cell_2_global(no_FEM_nodes_found_cell_2) = new_vertex_1;
				
% 				no_FEM_nodes_found_cell_2 = no_FEM_nodes_found_cell_2+1;
% 				FEM_nodes_without_centre_new_cell_2_global(no_FEM_nodes_found_cell_2) = new_FEM_node;
				
			else
				
				FEM_nodes_without_centre_new_cell_2_global(1:2) = [new_vertex_2 clockwise_intersection_vertex_2_global];
				no_FEM_nodes_found_cell_2 = 2;
				
				for temp_index = clockwise_intersection_vertex_2_local:length(cell_vertices)-1
					
					%                     disp('here1')
					
					node_1 = cell_vertices(temp_index);
					node_2 = cell_vertices(temp_index+1);
					
					temp_edge_node = find((FEM_nodes.edge(:,1)==node_1&FEM_nodes.edge(:,2)==node_2)|...
						(FEM_nodes.edge(:,1)==node_2&FEM_nodes.edge(:,2)==node_1));
					
					if isempty(temp_edge_node)
						
						no_FEM_nodes_found_cell_2 = no_FEM_nodes_found_cell_2+1;
						FEM_nodes_without_centre_new_cell_2_global(no_FEM_nodes_found_cell_2) = node_2;
						
					else
						
						no_FEM_nodes_found_cell_2 = no_FEM_nodes_found_cell_2+1;
						FEM_nodes_without_centre_new_cell_2_global(no_FEM_nodes_found_cell_2) = temp_edge_node;
						
						no_FEM_nodes_found_cell_2 = no_FEM_nodes_found_cell_2+1;
						FEM_nodes_without_centre_new_cell_2_global(no_FEM_nodes_found_cell_2) = node_2;
						
					end
					
				end
				
				node_1 = cell_vertices(end);
				node_2 = cell_vertices(1);
				
				temp_edge_node = find((FEM_nodes.edge(:,1)==node_1&FEM_nodes.edge(:,2)==node_2)|...
					(FEM_nodes.edge(:,1)==node_2&FEM_nodes.edge(:,2)==node_1));
				
				if isempty(temp_edge_node)
					
					no_FEM_nodes_found_cell_2 = no_FEM_nodes_found_cell_2+1;
					FEM_nodes_without_centre_new_cell_2_global(no_FEM_nodes_found_cell_2) = node_2;
					
				else
					
					no_FEM_nodes_found_cell_2 = no_FEM_nodes_found_cell_2+1;
					FEM_nodes_without_centre_new_cell_2_global(no_FEM_nodes_found_cell_2) = temp_edge_node;
					
					no_FEM_nodes_found_cell_2 = no_FEM_nodes_found_cell_2+1;
					FEM_nodes_without_centre_new_cell_2_global(no_FEM_nodes_found_cell_2) = node_2;
					
				end
				
				for temp_index = 1:intersection_vertex_1_local-1
					
					%                     disp('here2')
					
					node_1 = cell_vertices(temp_index);
					node_2 = cell_vertices(temp_index+1);
					
					temp_edge_node = find((FEM_nodes.edge(:,1)==node_1&FEM_nodes.edge(:,2)==node_2)|...
						(FEM_nodes.edge(:,1)==node_2&FEM_nodes.edge(:,2)==node_1));
					
					if isempty(temp_edge_node)
						
						no_FEM_nodes_found_cell_2 = no_FEM_nodes_found_cell_2+1;
						FEM_nodes_without_centre_new_cell_2_global(no_FEM_nodes_found_cell_2) = node_2;
						
					else
						
						no_FEM_nodes_found_cell_2 = no_FEM_nodes_found_cell_2+1;
						FEM_nodes_without_centre_new_cell_2_global(no_FEM_nodes_found_cell_2) = temp_edge_node;
						
						no_FEM_nodes_found_cell_2 = no_FEM_nodes_found_cell_2+1;
						FEM_nodes_without_centre_new_cell_2_global(no_FEM_nodes_found_cell_2) = node_2;
						
					end
					
				end
				
				no_FEM_nodes_found_cell_2 = no_FEM_nodes_found_cell_2+1;
				FEM_nodes_without_centre_new_cell_2_global(no_FEM_nodes_found_cell_2) = new_vertex_1;
				
% 				no_FEM_nodes_found_cell_2 = no_FEM_nodes_found_cell_2+1;
% 				FEM_nodes_without_centre_new_cell_2_global(no_FEM_nodes_found_cell_2) = new_FEM_node;
				
			end
			
			FEM_nodes_without_centre_new_cell_2_global =...
				FEM_nodes_without_centre_new_cell_2_global(1:no_FEM_nodes_found_cell_2);
			
			% find the missing concentration values at the centroids
			FEM_nodes.concentration(FEM_index_zero_cell+new_cell_1) = FindCentroidDppValue(centroid_new_cell_1,...
				FEM_nodes_without_centre_new_cell_1_global,FEM_nodes.concentration,FEM_nodes.previous_position,total_Dpp_in_new_cell_1);
			
			FEM_nodes.concentration(FEM_index_zero_cell+new_cell_2) = FindCentroidDppValue(centroid_new_cell_2,...
				FEM_nodes_without_centre_new_cell_2_global,FEM_nodes.concentration,FEM_nodes.previous_position,total_Dpp_in_new_cell_2);
			
			% this is now a new section where we are going to construct all
			% the new elements that are necessitated by the division/
			
			no_elements_new_cell_1 = length(FEM_nodes_without_centre_new_cell_1_global);
			no_elements_new_cell_2 = length(FEM_nodes_without_centre_new_cell_2_global);
			
			reusable_FEM_elements = cells.FEM_elements{cell_to_divide};
			
			% there are also two extra elements in the cells sharing the
			% edges that get divided. not always - if the cell is on the
			% boundary this doesn't happen. we haven't yet altered the cell
			% vertices in these adjoining cells - worry about this later
			% on. we could search through the FEM elements to see if there
			% are any other unused elements, e.g. those from cells that
			% have died. should implement this.
			required_no_new_elements =...
				no_elements_new_cell_1+no_elements_new_cell_2-...
				length(reusable_FEM_elements);
			
			old_no_FEM_elements = length(FEM_elements.nodes);
			
			FEM_elements.nodes(old_no_FEM_elements+1:old_no_FEM_elements+required_no_new_elements,:) = 0;
			
			elements_to_use = [reusable_FEM_elements ...
				old_no_FEM_elements+1:old_no_FEM_elements+required_no_new_elements];
			
			% add the appropriate elements to cells.FEM_elements for each cell
			cells.FEM_elements{new_cell_1} = elements_to_use(1:no_elements_new_cell_1);
			
			cells.FEM_elements{new_cell_2} = ...
				elements_to_use(no_elements_new_cell_1+1:...
				no_elements_new_cell_1+no_elements_new_cell_2);
			
			% sort out the FEM_elements for the new cells by looping over
			% their FEM nodes and adding in the centroid node
			for temp_current_FEM_node_local = 1:no_elements_new_cell_1
				
				first_node_current_element =...
					FEM_nodes_without_centre_new_cell_1_global(temp_current_FEM_node_local);
				
				second_node_current_element =...
					FEM_nodes_without_centre_new_cell_1_global(mod(temp_current_FEM_node_local,no_elements_new_cell_1)+1);
				% it may not be immediately obvious why mod is used in this
				% way! but it works out correctly - for example type
				% mod((1:5)-1,4)+1. normally the clockwise node would be
				% temp_current_FEM_node_local+1.
				
				current_element_global = elements_to_use(temp_current_FEM_node_local);
				
				FEM_elements.nodes(current_element_global,:) =...
					[first_node_current_element second_node_current_element ...
					FEM_index_zero_cell+new_cell_1];
				
			end
			
			for temp_current_FEM_node_local = 1:no_elements_new_cell_2
				
				first_node_current_element =...
					FEM_nodes_without_centre_new_cell_2_global(temp_current_FEM_node_local);
				
				second_node_current_element =...
					FEM_nodes_without_centre_new_cell_2_global(mod(temp_current_FEM_node_local,no_elements_new_cell_2)+1);
				
				current_element_global = elements_to_use(no_elements_new_cell_1+temp_current_FEM_node_local);
				
				FEM_elements.nodes(current_element_global,:) =...
					[first_node_current_element second_node_current_element ...
					FEM_index_zero_cell+new_cell_2];
				
			end
			
			% we also need to split one of the FEM elements in the
			% cell_with_same_edge. this element is the one that contains
			% FEM nodes corresponding to the cell vertices along the edge
			% that has itself been split. this will require editing one
			% element and adding another. the choice of which element
			% the edited element becomes is arbitrary.
			if ~isempty(cell_with_same_edge_1)
				
				% again, we could instead look for unused FEM_elements here
				% rather than just tacking the new one on the end.
				new_FEM_element = length(FEM_elements.nodes)+1;
				
				% clockwise_intersection_vertex comes before intersection_vertex
				% in the FEM element because we are in the adjoining cell,
				% and any edge the is clockwise in one cell will be
				% anticlockwise in the adjacent cell.
				FEM_element_to_split = ...
					find((FEM_elements.nodes(:,1)==clockwise_intersection_vertex_1_global)&...
					(FEM_elements.nodes(:,2)==intersection_vertex_1_global));
				
				% edit this element
				FEM_elements.nodes(FEM_element_to_split,2) = new_vertex_1;
				
				FEM_elements.nodes(new_FEM_element,:) = ...
					[new_vertex_1 intersection_vertex_1_global ...
					FEM_index_zero_cell+cell_with_same_edge_1];
				
				% the final step is to edit the elements associated with
				% the cell to include the new element at the correct
				% location
				cell_with_same_edge_1_FEM_elements = cells.FEM_elements{cell_with_same_edge_1};
				
				index_element_to_split_in_cell_elements =...
					find(cell_with_same_edge_1_FEM_elements==FEM_element_to_split);
				
				cells.FEM_elements{cell_with_same_edge_1} =...
					[cell_with_same_edge_1_FEM_elements(1:index_element_to_split_in_cell_elements) new_FEM_element ...
					cell_with_same_edge_1_FEM_elements(index_element_to_split_in_cell_elements+1:end)];
				
			end
			
			if ~isempty(cell_with_same_edge_2)
				
				new_FEM_element = length(FEM_elements.nodes)+1;
				
				FEM_element_to_split = ...
					find((FEM_elements.nodes(:,1)==clockwise_intersection_vertex_2_global)&...
					(FEM_elements.nodes(:,2)==intersection_vertex_2_global));
				
				FEM_elements.nodes(FEM_element_to_split,2) = new_vertex_2;
				
				FEM_elements.nodes(new_FEM_element,:) = ...
					[new_vertex_2 intersection_vertex_2_global ...
					FEM_index_zero_cell+cell_with_same_edge_2];
				
				cell_with_same_edge_2_FEM_elements = cells.FEM_elements{cell_with_same_edge_2};
				
				index_element_to_split_in_cell_elements =...
					find(cell_with_same_edge_2_FEM_elements==FEM_element_to_split);
				
				cells.FEM_elements{cell_with_same_edge_2} =...
					[cell_with_same_edge_2_FEM_elements(1:index_element_to_split_in_cell_elements) new_FEM_element ...
					cell_with_same_edge_2_FEM_elements(index_element_to_split_in_cell_elements+1:end)];
				
			end
		end
	end
end

if mitosis_location ~= 0
	
	stats.no_mitosis_since_last_statistics = stats.no_mitosis_since_last_statistics+1;
	mitosis_counter = mitosis_counter+1;
	stats.mitosis_locations(mitosis_counter,:) = mitosis_location;
	
end

if stats.this_iteration_logical
	
	stats.no_mitosis(stats.counter) = stats.no_mitosis_since_last_statistics;
	stats.no_mitosis_since_last_statistics = 0;
	
end