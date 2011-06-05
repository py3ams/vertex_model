function stats = generate_statistics(cells,vertices,angle_deviations,delta_t,...
    diffusion_speed,FEM_nodes,FEM_solve_logical,mean_edge_length,...
    stats,vertex_movements)

if stats.this_iteration_logical
    
    statistics_counter = stats.counter;
    vertex_movement_magnitudes = sqrt(sum(vertex_movements.^2,2));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % we obviously don't want to include dead cells. apoptotic cells also present a
    % problem as their volume is set to 0 so cell height becomes NaN.
    cells_logical = cells.state==1|cells.state==2;
    no_alive_cells = sum(cells_logical);
    stats.total_no_cells(statistics_counter) = no_alive_cells;
    vertices_logical = vertices.no_cells > 0;
    stats.total_no_vertices(statistics_counter) = sum(vertices_logical);
    no_cells_per_vertex = vertices.no_cells(vertices_logical);
    stats.rosette(statistics_counter,:) = [sum(no_cells_per_vertex==4) ...
        sum(no_cells_per_vertex >= 5)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cell_heights = cells.volume(cells_logical)./cells.area(cells_logical);
    height_area_ratios = cell_heights./cells.area(cells_logical);
    vertices_per_cell = cellfun('length',cells.vertices);
    vertices_per_cell = vertices_per_cell(vertices_per_cell>0);
    vertex_movement_magnitudes = vertex_movement_magnitudes(vertices_logical);
    
    stat_values = MeanMaxMinStdTotal(angle_deviations,no_cells_per_vertex,cell_heights,...
        cells.volume(cells_logical),height_area_ratios,vertex_movement_magnitudes,vertices_per_cell,5);
    
    stats.angle_deviation(statistics_counter,:) = stat_values(1,1:4);
    stats.no_cells_per_vertex(statistics_counter,:) = stat_values(2,1:3);
    stats.cell_height(statistics_counter,:) = stat_values(3,1:4);
    stats.cell_volume(statistics_counter,:) = stat_values(4,1:4);
    stats.cell_height_to_area(statistics_counter,:) = stat_values(5,1:4);
    stats.vertex_movement(statistics_counter,:) = stat_values(6,1:2);
    stats.vertices_per_cell(statistics_counter,:) = stat_values(7,1:3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    vertex_positions = vertices.position;
    
    stats.xy_ratio(statistics_counter) = (max(vertex_positions(:,1))-min(...
        vertex_positions(:,1)))/(max(vertex_positions(:,2))-min(vertex_positions(:,2)));
    
    
    if FEM_solve_logical
        
        max_vertex_movement = stats.vertex_movement(statistics_counter,2);
        
        stats.mesh_peclet_number(statistics_counter) =...
            mean_edge_length*max_vertex_movement/(delta_t*min(diffusion_speed));
        
    else
        
        stats.mesh_peclet_number(statistics_counter) = 0;
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    no_FEM_nodes = size(FEM_nodes.concentration,1);
    
    % we might want these distributions even if FEM_solve_logical is false. in that
    % case would need to calculate the cell centres as they will not be stored in
    % FEM_node_positions. alternatively, find somewhere else where they are
    % calculated and add them to cells.
    if no_FEM_nodes>1
        % note that this is not the same as the number of alive cells that we use for
        % other stats
        length_cells = length(cells.vertices);
        cell_centre_positions = FEM_nodes.previous_position(no_FEM_nodes-length_cells+1:end,:);
        cell_centre_positions = cell_centre_positions(cells_logical,:);
        
        shifted_cell_centre_positions_x = cell_centre_positions(:,1)-min(cell_centre_positions(:,1));
        normalised_cell_centre_positions_x =...
            shifted_cell_centre_positions_x/max(shifted_cell_centre_positions_x);
        
        shifted_cell_centre_positions_y = cell_centre_positions(:,2)-min(cell_centre_positions(:,2));
        normalised_cell_centre_positions_y =...
            shifted_cell_centre_positions_y/max(shifted_cell_centre_positions_y);
        
%         stats.volume_distribution_x(statistics_counter,:) = polyfit(...
%             normalised_cell_centre_positions_x,cells.volume(cells_logical),1);
%         stats.volume_distribution_y(statistics_counter,:) = polyfit(...
%             normalised_cell_centre_positions_y,cells.volume(cells_logical),1);
%         
%         stats.area_distribution_x(statistics_counter,:) = polyfit(...
%             normalised_cell_centre_positions_x,cells.area(cells_logical),1);
%         stats.area_distribution_y(statistics_counter,:) = polyfit(...
%             normalised_cell_centre_positions_y,cells.area(cells_logical),1);
        
        % not really sure if it makes sense to do this
        normalised_cell_volumes = cells.volume(cells_logical)/...
            max(cells.volume(cells_logical));
        
        stats.normalised_volume_distribution_x(statistics_counter,:) = polyfit(...
            normalised_cell_centre_positions_x,normalised_cell_volumes,1);
        stats.normalised_volume_distribution_y(statistics_counter,:) = polyfit(...
            normalised_cell_centre_positions_y,normalised_cell_volumes,1);
        
        normalised_cell_areas = cells.area(cells_logical)/...
            max(cells.area(cells_logical));
        
        stats.normalised_area_distribution_x(statistics_counter,:) = polyfit(...
            normalised_cell_centre_positions_x,normalised_cell_areas,1);
        stats.normalised_area_distribution_y(statistics_counter,:) = polyfit(...
            normalised_cell_centre_positions_y,normalised_cell_areas,1);
    end
    
end
