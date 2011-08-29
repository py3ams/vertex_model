function [FEM_elements,FEM_nodes,cell_elements] =...
   create_FEM_mesh(cells,vertices,no_refinements,no_mitosis_flag)

if nargin<4
   no_mitosis_flag = false;
end

% there are two distincts cases here. we can't have no_refinements > 0 if any
% rearrangements are going off. if no_refinements==0 and there are rearrangements
% (inc. mitosis/death), then we need to include space for the FEM edge nodes, which
% are created by default during mitosis. if no_refinements>0 these edge nodes are not
% required, but extra nodes for the refinements obviously are.

if no_refinements > 0 || no_mitosis_flag

    % size of the vectors to set aside to store node positions of the nodes
    % associated with different levels of refinement. these values are multiplied by
    % the size of the vector that stores vertex positions. the total number of FEM_nodes after
    % a refinement is up to four times the previous number (consider Euler's formula -
    % there is a new node on each edge). therefore the cumulative number of FEM_nodes
    % set aside is equal to 4^refinement_level*initial_no_nodes. we need to
    % remember to take into account the initial nodes that are cell centres,
    % which are why the numbers in this vector are bigger than
    % 4^refinement_level, as we multiply them by the initial_no_vertices, not
    % initial_no_nodes.
    refinement_level_sizes = [5 20 80 320];
    
    no_cells = length(cells.vertices);
    cell_elements = cell(no_cells,1);
    
    cell_areas = CalculateCellAreas(cells.vertices,vertices.position);
    
    % when there are no refinements the number of elements is equal to the
    % total number of cell vertices (i.e. counting each vertex multiple times
    % for the number of cells it is in). for each refinement the number of
    % elements increases by a factor of 4. this should be exactly right, there
    % should be no extra zeros on the end of FEM_elements.nodes at the end of
    % this function.
    FEM_elements.nodes = zeros(sum(cellfun('length',cells.vertices))*4^no_refinements,3);
    
    % find the total size of the vector needed to store the nodes created by
    % refinement
    total_refinement_size = sum(refinement_level_sizes(1:no_refinements))*length(vertices.position);
    
    % would it make more sense to arrange this as [vertices.position
    % zeros(max_no_cells,2) zeros(total_refinement_size,2)] and let the total
    % refinement size include cell centres up to the maximum number of cells?
    FEM_nodes.position = ...
        [vertices.position; zeros(total_refinement_size,2); zeros(no_cells,2)];
    
    FEM_nodes.edge = zeros(size(FEM_nodes.position,1),2);
    
    %index of the first cell centre within FEM_nodes.position. could also be
    %found by length(vertices.position)+total_refinement_size
    cell_centre_node_index = size(FEM_nodes.position,1)-no_cells;
    
    % finds the index for the first node associated with each level of refinement
    % (actually finds node before first node as we add 1 before we store a node)
    edge_node_indices = length(vertices.position)+([0 cumsum(refinement_level_sizes(...
        1:end-1))])*length(vertices.position);
    
    no_elements = 0;
    
    for current_cell = 1:no_cells
        
        cell_vertices = cells.vertices{current_cell};
        no_cell_vertices = length(cell_vertices);
        
        cell_vertex_positions = vertices.position(cell_vertices,:);
        cell_centre_position = CalculateCentroid(cell_vertex_positions,cell_areas(current_cell));
        
        cell_centre_node_index = cell_centre_node_index+1;
        
        FEM_nodes.position(cell_centre_node_index,:) = cell_centre_position;
        
        % it works out easier to do this here rather than after no_elements has
        % changed at the end. the number of elements added to the current cell
        % will be equal to 4^no_refinements*no_cell_vertices
        cell_elements{current_cell} = no_elements+1:no_elements+4^no_refinements*no_cell_vertices;
        
        for current_vertex_local = 1:length(cell_vertices)
            
            current_vertex_global = cell_vertices(current_vertex_local);
            clockwise_vertex_global = cell_vertices(mod(current_vertex_local,no_cell_vertices)+1);
            
            % we create this cell triangle_refinements each time round this loop.
            % it is associated with the triangle that would have been created in
            % the case of no_refinements, and the first entry in the cell is
            % equal to that big triangular element. for each refinement level we
            % create a new entry in the cell that contains the elements at the
            % current level. for must then loop over each of those elements at
            % the next level. the number of elements in the cell at each
            % refinement level is therefore equal to 4^no_refinements.
            triangle_refinements = cell(no_refinements+1,1);
            triangle_refinements{1} = [current_vertex_global clockwise_vertex_global cell_centre_node_index];
            
            for refinement_level = 1:no_refinements
                
                no_triangles_at_current_level = 0;
                
                % loop over each triangle at the previous refinement_level. n.b.
                % as we can't have an index 0 for no refinements, the indices are
                % shifted by 1 so when we index with refinement_level it is
                % actually the previous level!
                for current_outer_triangle_index = 1:size(triangle_refinements{refinement_level},1)
                    
                    current_triangle_nodes = triangle_refinements{refinement_level}(current_outer_triangle_index,:);
                    
                    % check with the first edge (arbitrarily made up of node 1 and
                    % node 2) has already been refined. if not, find the position
                    % of the new edge node and add it to FEM_nodes.position and
                    % FEM_nodes.edge.
                    edge_node_1_index =...
                        find((FEM_nodes.edge(:,1)==current_triangle_nodes(1)&FEM_nodes.edge(:,2)==current_triangle_nodes(2))|...
                        (FEM_nodes.edge(:,1)==current_triangle_nodes(2)&FEM_nodes.edge(:,2)==current_triangle_nodes(1)),1);
                    
                    if isempty(edge_node_1_index)
                        
                        edge_node_indices(refinement_level) = edge_node_indices(refinement_level)+1;
                        edge_node_1_index = edge_node_indices(refinement_level);
                        
                        FEM_nodes.position(edge_node_1_index,:) = 0.5*(FEM_nodes.position(current_triangle_nodes(1),:)+...
                            FEM_nodes.position(current_triangle_nodes(2),:));
                        
                        FEM_nodes.edge(edge_node_1_index,:) = [current_triangle_nodes(1) current_triangle_nodes(2)];
                        
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    % second edge of the big triangle is node 1 and node 3
                    edge_node_2_index =...
                        find((FEM_nodes.edge(:,1)==current_triangle_nodes(1)&FEM_nodes.edge(:,2)==current_triangle_nodes(3))|...
                        (FEM_nodes.edge(:,1)==current_triangle_nodes(3)&FEM_nodes.edge(:,2)==current_triangle_nodes(1)),1);
                    
                    if isempty(edge_node_2_index)
                        
                        edge_node_indices(refinement_level) = edge_node_indices(refinement_level)+1;
                        edge_node_2_index = edge_node_indices(refinement_level);
                        
                        FEM_nodes.position(edge_node_2_index,:) = 0.5*(FEM_nodes.position(current_triangle_nodes(1),:)+...
                            FEM_nodes.position(current_triangle_nodes(3),:));
                        
                        FEM_nodes.edge(edge_node_2_index,:) = [current_triangle_nodes(1) current_triangle_nodes(3)];
                        
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    % third edge is node 2 and node 3
                    edge_node_3_index =...
                        find((FEM_nodes.edge(:,1)==current_triangle_nodes(2)&FEM_nodes.edge(:,2)==current_triangle_nodes(3))|...
                        (FEM_nodes.edge(:,1)==current_triangle_nodes(3)&FEM_nodes.edge(:,2)==current_triangle_nodes(2)),1);
                    
                    if isempty(edge_node_3_index)
                        
                        edge_node_indices(refinement_level) = edge_node_indices(refinement_level)+1;
                        edge_node_3_index = edge_node_indices(refinement_level);
                        
                        FEM_nodes.position(edge_node_3_index,:) = 0.5*(FEM_nodes.position(current_triangle_nodes(2),:)+...
                            FEM_nodes.position(current_triangle_nodes(3),:));
                        
                        FEM_nodes.edge(edge_node_3_index,:) = [current_triangle_nodes(2) current_triangle_nodes(3)];
                        
                    end
                    
                    % add the four triangles to the triangle_refinements cell for
                    % this level (which is indexed at refinement_level+1}. can
                    % figure these out by drawing a triangle and refining it into 4
                    % smaller ones.
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
            
            % need this as the loop starts at refinement_level=1. could instead
            % loop from refinement_level =0 and change all the indices but this
            % way is simpler.
            if no_refinements == 0
                refinement_level = 0;
            end
            
            % add all the triangles at the correct refinement level to
            % FEM_elements.
            FEM_elements.nodes(no_elements+1:no_elements+4^no_refinements,:) = triangle_refinements{refinement_level+1};
            no_elements = no_elements+4^no_refinements;
            
        end
        
    end
    
else
    
    no_cells = length(cells.vertices);
    cell_elements = cell(no_cells,1);
    
    cell_areas = CalculateCellAreas(cells.vertices,vertices.position);
    
    FEM_elements.nodes = zeros(sum(cellfun('length',cells.vertices)),3);
    
    FEM_nodes.position = ...
        [vertices.position; zeros(3*length(vertices.position),2); zeros(no_cells,2)];
    
    FEM_nodes.edge = zeros(size(FEM_nodes.position,1),2);
        
    no_elements = 0;
    cell_centre_node_index = size(FEM_nodes.position,1)-no_cells;
    
    for current_cell = 1:no_cells
        
        cell_vertices = cells.vertices{current_cell};
        no_cell_vertices = length(cell_vertices);
        
        cell_vertex_positions = vertices.position(cell_vertices,:);
        cell_centre_position = CalculateCentroid(cell_vertex_positions,cell_areas(current_cell));
        
        cell_elements{current_cell} = no_elements+1:no_elements+no_cell_vertices;
        
        cell_centre_node_index = cell_centre_node_index+1;
        
        FEM_nodes.position(cell_centre_node_index,:) = cell_centre_position;
        
        for current_vertex_local = 1:length(cell_vertices)
            
            no_elements = no_elements+1;
            
            current_vertex_global = cell_vertices(current_vertex_local);
            clockwise_vertex_global = cell_vertices(mod(current_vertex_local,no_cell_vertices)+1);
            
            FEM_elements.nodes(no_elements,:) = [current_vertex_global clockwise_vertex_global cell_centre_node_index];
            
        end
        
    end
end

