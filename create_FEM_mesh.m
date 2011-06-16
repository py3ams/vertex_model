function [FEM_elements,FEM_nodes,cell_elements] =...
    create_FEM_mesh(cells,vertex_positions,no_refinements)

% size of the vectors to set aside to store node positions of the nodes
% associated with different levels of refinement. the total number of FEM_nodes after
% a refinement is up to four times the previous number (consider Euler's formula -
% there is a new node on each edge). therefore the cumulative number of FEM_nodes 
% set aside is equal to 4^refinement_level*initial_no_nodes. however, we do not take
% account of the fact that there are also nodes at cell centres, which we should.
% instead we are relying on the fact that we already leave plenty of space for new
% cells etc.
refinement_level_sizes = [3 13 47];

no_cells = length(cells);
cell_elements = cell(no_cells,1);

cell_areas = CalculateCellAreas(cells,vertex_positions);

FEM_elements.nodes = zeros(sum(cellfun('length',cells))*4^no_refinements,3);

total_refinement_size = sum(refinement_level_sizes(1:no_refinements))*length(vertex_positions);

FEM_nodes.position = ...
    [vertex_positions; zeros(total_refinement_size,2); zeros(no_cells,2)];

FEM_nodes.edge = zeros(size(FEM_nodes.position,1),2);

no_elements = 0;
cell_centre_node_index = size(FEM_nodes.position,1)-no_cells;

% there is an edge_node_counter for each level of refinement
edge_node_counters = length(vertex_positions)+([0 refinement_level_sizes(1:end-1)])*length(vertex_positions);

for current_cell = 1:no_cells
    
    cell_vertices = cells{current_cell};
    no_cell_vertices = length(cell_vertices);
    
    cell_vertex_positions = vertex_positions(cell_vertices,:);
    cell_centre_position = CalculateCentroid(cell_vertex_positions,cell_areas(current_cell));
    
    cell_centre_node_index = cell_centre_node_index+1;
    
    FEM_nodes.position(cell_centre_node_index,:) = cell_centre_position;
    
    % it works out easier to do this here rather than after no_elements has
    % changed at the end
    cell_elements{current_cell} = no_elements+1:no_elements+4^no_refinements*no_cell_vertices;
    
    for current_vertex_local = 1:length(cell_vertices)
        
        current_vertex_global = cell_vertices(current_vertex_local);
        clockwise_vertex_global = cell_vertices(mod(current_vertex_local,no_cell_vertices)+1);
        
        triangle_refinements = cell(no_refinements,1);
        triangle_refinements{1} = [current_vertex_global clockwise_vertex_global cell_centre_node_index];
        
        % 		if refinement_logical
        for refinement_level = 1:no_refinements
            
            no_triangles_at_current_level = 0;
            
            for current_outer_triangle_index = 1:size(triangle_refinements{refinement_level},1)
                
                current_triangle_nodes = triangle_refinements{refinement_level}(current_outer_triangle_index,:);
                
                edge_node_1_index =...
                    find((FEM_nodes.edge(:,1)==current_triangle_nodes(1)&FEM_nodes.edge(:,2)==current_triangle_nodes(2))|...
                    (FEM_nodes.edge(:,1)==current_triangle_nodes(2)&FEM_nodes.edge(:,2)==current_triangle_nodes(1)),1);
                
                if isempty(edge_node_1_index)
                    
                    edge_node_counters(refinement_level) = edge_node_counters(refinement_level)+1;
                    edge_node_1_index = edge_node_counters(refinement_level);
                    
                    FEM_nodes.position(edge_node_1_index,:) = 0.5*(FEM_nodes.position(current_triangle_nodes(1),:)+...
                        FEM_nodes.position(current_triangle_nodes(2),:));
                    
                    FEM_nodes.edge(edge_node_1_index,:) = [current_triangle_nodes(1) current_triangle_nodes(2)];
                    
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                edge_node_2_index =...
                    find((FEM_nodes.edge(:,1)==current_triangle_nodes(1)&FEM_nodes.edge(:,2)==current_triangle_nodes(3))|...
                    (FEM_nodes.edge(:,1)==current_triangle_nodes(3)&FEM_nodes.edge(:,2)==current_triangle_nodes(1)),1);
                
                if isempty(edge_node_2_index)
                    
                    edge_node_counters(refinement_level) = edge_node_counters(refinement_level)+1;
                    edge_node_2_index = edge_node_counters(refinement_level);
                    
                    FEM_nodes.position(edge_node_2_index,:) = 0.5*(FEM_nodes.position(current_triangle_nodes(1),:)+...
                        FEM_nodes.position(current_triangle_nodes(3),:));
                    
                    FEM_nodes.edge(edge_node_2_index,:) = [current_triangle_nodes(1) current_triangle_nodes(3)];
                    
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                edge_node_3_index =...
                    find((FEM_nodes.edge(:,1)==current_triangle_nodes(2)&FEM_nodes.edge(:,2)==current_triangle_nodes(3))|...
                    (FEM_nodes.edge(:,1)==current_triangle_nodes(3)&FEM_nodes.edge(:,2)==current_triangle_nodes(2)),1);
                
                if isempty(edge_node_3_index)
                    
                    edge_node_counters(refinement_level) = edge_node_counters(refinement_level)+1;
                    edge_node_3_index = edge_node_counters(refinement_level);
                    
                    FEM_nodes.position(edge_node_3_index,:) = 0.5*(FEM_nodes.position(current_triangle_nodes(2),:)+...
                        FEM_nodes.position(current_triangle_nodes(3),:));
                    
                    FEM_nodes.edge(edge_node_3_index,:) = [current_triangle_nodes(2) current_triangle_nodes(3)];
                    
                end
                
                no_triangles_at_current_level = no_triangles_at_current_level+1;
                triangle_refinements{refinement_level+1}(no_triangles_at_current_level,:) = ...
                    [current_triangle_nodes(1) edge_node_1_index edge_node_2_index];
                
                no_triangles_at_current_level = no_triangles_at_current_level+1;
                triangle_refinements{refinement_level+1}(no_triangles_at_current_level,:) = ...
                    [edge_node_1_index edge_node_3_index edge_node_2_index];
                
                no_triangles_at_current_level = no_triangles_at_current_level+1;
                triangle_refinements{refinement_level+1}(no_triangles_at_current_level,:) = ...
                    [edge_node_1_index current_triangle_nodes(2) edge_node_3_index];
                
                no_triangles_at_current_level = no_triangles_at_current_level+1;
                triangle_refinements{refinement_level+1}(no_triangles_at_current_level,:) = ...
                    [edge_node_2_index edge_node_3_index current_triangle_nodes(3)];
                
                
            end
            
        end
        
        if no_refinements == 0
            refinement_level = 0;
        end
        
        FEM_elements.nodes(no_elements+1:no_elements+4^no_refinements,:) = triangle_refinements{refinement_level+1};
        no_elements = no_elements+4^no_refinements;
        
    end
    
end

FEM_elements.nodes = FEM_elements.nodes(FEM_elements.nodes(:,1)>0,:);


