function [FEM_elements,FEM_nodes,cell_elements] =...
    create_FEM_mesh(cells,vertex_positions)

refinement_logical = true;

no_cells = length(cells);
cell_elements = cell(no_cells,1);

cell_areas = CalculateCellAreas(cells,vertex_positions);

FEM_elements.nodes = zeros(sum(cellfun('length',cells)),3);

% the middle term is zeros we leave for possible mesh refinement. the number of edges
% in a triangular
FEM_nodes.position = ...
    [vertex_positions; zeros(3*length(vertex_positions),2); zeros(no_cells,2)];

FEM_nodes.edge = zeros(4*length(vertex_positions)+no_cells,2);

no_elements = 0;
cell_centre_node_index = 4*length(vertex_positions);
edge_node_index = length(vertex_positions);

for current_cell = 1:no_cells
    
    cell_vertices = cells{current_cell};
    no_cell_vertices = length(cell_vertices);
    
    cell_vertex_positions = vertex_positions(cell_vertices,:);
    cell_centre = CalculateCentroid(cell_vertex_positions,cell_areas(current_cell));
    
    cell_centre_node_index = cell_centre_node_index+1;
    
    FEM_nodes.position(cell_centre_node_index,:) = cell_centre;
    
    cell_elements{current_cell} = no_elements+1:no_elements+no_cell_vertices;
    
    for current_vertex_local = 1:length(cell_vertices)
        
        current_vertex_global = cell_vertices(current_vertex_local);
        clockwise_vertex_global = cell_vertices(mod(current_vertex_local,no_cell_vertices)+1);
        
        if refinement_logical
            
            % need to check first if the edge has already been done. to achieve this
            % need to record edges in FEM_nodes.edge once they are created. have got
            % the logic somewhere to find if an edge is already in this vector.
            edge_node_index = edge_node_index+1;
            
            edge_node_1_index = edge_node_index;
            
            edge_node_1_position = 0.5*(FEM_nodes.position(current_vertex_global,:)+...
                FEM_nodes.position(clockwise_vertex_global,:));           
            
            FEM_nodes.position(edge_node_index) = edge_node_1_position;
            
            %%%
            
            edge_node_index = edge_node_index+1;
            
            edge_node_2_index = edge_node_index;
            
            edge_node_2_position = 0.5*(FEM_nodes.position(current_vertex_global,:)+...
                cell_centre);

            FEM_nodes.position(edge_node_index) = edge_node_2_position;
            
            %%%
                       
            edge_node_index = edge_node_index+1;
            
            edge_node_3_index = edge_node_index;
            
            edge_node_3_position = 0.5*(FEM_nodes.position(clockwise_vertex_global,:)+...
                cell_centre);

            FEM_nodes.position(edge_node_index) = edge_node_3_position;
            
        end
        
        no_elements = no_elements+1;
        FEM_elements.nodes(no_elements,:) = [current_vertex_global edge_node_1_index edge_node_2_index];
        
        no_elements = no_elements+1;
        FEM_elements.nodes(no_elements,:) = [edge_node_1_index edge_node_3_index edge_node_2_index];
        
        no_elements = no_elements+1;
        FEM_elements.nodes(no_elements,:) = [edge_node_1_index clockwise_vertex_global edge_node_3_index];
        
        no_elements = no_elements+1;
        FEM_elements.nodes(no_elements,:) = [edge_node_2_index edge_node_3_index cell_centre_node_index];
        
        % need to all these elements to cell elements
        
    end
    
end
