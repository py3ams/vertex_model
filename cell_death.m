function [cells,death_counter,FEM_elements,FEM_nodes,no_deaths_this_iteration,...
   refined_edge_matrix,stats,vertices] = cell_death(cells,death_counter,...
   FEM_elements,FEM_nodes,FEM_solve_logical,protection_time,refined_edge_matrix,...
   stats,threshold_area,time,vertices)

length_cells = length(cells.vertices);

no_deaths_this_iteration = 0;

% it is much quicker to do these logical tests with a find command outside
% the loop rather than looping over all cells and doing logical tests for
% each one separately (even though that way the find command is not used).
% this way is also quicker than, for example, just finding the apoptotic
% cells outside the loop and checking their area and no_vertices inside the
% loop, though this may obviously depend on how frequently apoptotic cells
% occur in simulations. don't change this again unless you really know what
% you are doing!
potential_death_cells =...
   find((cells.area<threshold_area)&(cells.state==3)&...
   (cellfun('length',cells.vertices)==3));

for current_cell_local = 1:length(potential_death_cells)
   
   death_location = 0;
   current_cell = potential_death_cells(current_cell_local);
   cell_vertices = cells.vertices{current_cell};
   
   if all(time-vertices.time_created(cell_vertices) > protection_time)
      
      %         current_cell
      %         error('save variables now');
      
      no_deaths_this_iteration = no_deaths_this_iteration + 1;
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      %         disp(['Dying cell = ',num2str(current_cell)])
      %         vertices.position(cell_vertices,:)
      %         disp(['Vertex positions = ',num2str(vertices.position(cell_vertices,:))])
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      new_vertex_position = mean(vertices.position(cell_vertices,:));
      death_location = new_vertex_position;
      
      new_vertex = find(vertices.no_cells == 0,1);
      vertices.position(new_vertex,:) = new_vertex_position;
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      cells_with_same_edges = zeros(1,3);
      
      [~,~,cell_with_same_edges_1] =...
         find_cells_containing_vertices(vertices.cells,current_cell,...
         cell_vertices(1),cell_vertices(2));
      
      if ~isempty(cell_with_same_edges_1)
         
         cells_with_same_edges(1) = cell_with_same_edges_1;
         
         temp_cell_vertices = cells.vertices{cells_with_same_edges(1)};
         temp_cell_vertices(temp_cell_vertices==cell_vertices(1)) = [];
         temp_cell_vertices(temp_cell_vertices==cell_vertices(2)) = new_vertex;
         
         cells.vertices{cell_with_same_edges_1} = temp_cell_vertices;
         
         new_cell_area = CalculateCellAreas(cells.vertices(cell_with_same_edges_1),vertices.position);
         
         if FEM_solve_logical
            
            % should not need to update internal_chemical_quantity as
            % the whole point is that we are keeping it constant.
            cells.internal_chemical_value(cell_with_same_edges_1) =...
               cells.internal_chemical_value(cell_with_same_edges_1)*cells.area(cell_with_same_edges_1)/new_cell_area;
            
         end
         
         % important to do this after updating internal_chemical_value
         cells.area(cell_with_same_edges_1) = new_cell_area;
         
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      [~,~,cell_with_same_edges_2] =...
         find_cells_containing_vertices(vertices.cells,current_cell,...
         cell_vertices(2),cell_vertices(3));
      
      if ~isempty(cell_with_same_edges_2)
         
         cells_with_same_edges(2) = cell_with_same_edges_2;
         
         temp_cell_vertices = cells.vertices{cells_with_same_edges(2)};
         temp_cell_vertices(temp_cell_vertices==cell_vertices(2)) = [];
         temp_cell_vertices(temp_cell_vertices==cell_vertices(3)) = new_vertex;
         
         cells.vertices{cells_with_same_edges(2)} = temp_cell_vertices;
         
         new_cell_area = CalculateCellAreas(cells.vertices(...
            cell_with_same_edges_2),vertices.position);
         
         if FEM_solve_logical
            
            % should not need to update internal_chemical_quantity as
            % the whole point is that we are keeping it constant.
            cells.internal_chemical_value(cell_with_same_edges_2) =...
               cells.internal_chemical_value(cell_with_same_edges_2)*...
               cells.area(cell_with_same_edges_2)/new_cell_area;
            
         end
         
         % important to do this after updating internal_chemical_value
         cells.area(cell_with_same_edges_2) = new_cell_area;
         
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      [~,~,cell_with_same_edges_3] =...
         find_cells_containing_vertices(vertices.cells,current_cell,...
         cell_vertices(3),cell_vertices(1));
      
      if ~isempty(cell_with_same_edges_3)
         
         cells_with_same_edges(3) = cell_with_same_edges_3;
         
         temp_cell_vertices = cells.vertices{cells_with_same_edges(3)};
         temp_cell_vertices(temp_cell_vertices==cell_vertices(3)) = [];
         temp_cell_vertices(temp_cell_vertices==cell_vertices(1)) = new_vertex;
         
         cells.vertices{cells_with_same_edges(3)} = temp_cell_vertices;
         
         new_cell_area = CalculateCellAreas(cells.vertices(...
            cell_with_same_edges_3),vertices.position);
         
         if FEM_solve_logical
            
            % should not need to update internal_chemical_quantity as
            % the whole point is that we are keeping it constant.
            cells.internal_chemical_value(cell_with_same_edges_3) =...
               cells.internal_chemical_value(cell_with_same_edges_3)*...
               cells.area(cell_with_same_edges_3)/new_cell_area;
            
         end
         
         % important to do this after updating internal_chemical_value
         cells.area(cell_with_same_edges_3) = new_cell_area;
         
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      cells.state(current_cell) = 4;
      cells.vertices{current_cell} = [];
      cells.volume(current_cell) = 0;
      cells.area(current_cell) = 0;
      cells.growth_speed(current_cell) = 0;
      
      cells.force_constants.area(current_cell) = 0;
      cells.force_constants.deformation(current_cell) = 0;
      cells.force_constants.elongation(current_cell) = 0;
      cells.force_constants.perimeter(current_cell) = 0;
      cells.force_constants.tension(current_cell) = 0;
      
      cells.target_area(current_cell) = 0;
      cells.target_volume(current_cell) = 0;
      
      vertices.no_cells(cell_vertices) = 0;
      vertices.cells(cell_vertices,:) = 0;
      vertices.position(cell_vertices,:) = 0;
      
      no_cells_with_same_edges = sum(cells_with_same_edges>0);
      
      vertices.no_cells(new_vertex) = no_cells_with_same_edges;
      vertices.cells(new_vertex,1:no_cells_with_same_edges) =...
         cells_with_same_edges(cells_with_same_edges>0);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      if FEM_solve_logical
         
         cells.ingestion_rate(current_cell) = 0;
         
         cells.internal_chemical_value(current_cell) = 0;
         cells.internal_chemical_quantity(current_cell) = 0;
         
         new_node = new_vertex;
         
         length_FEM_nodes = length(FEM_nodes.previous_position);
         FEM_index_zero_cell = length_FEM_nodes-length_cells;
         
         FEM_nodes.previous_position(new_node,:) = FEM_nodes.previous_position(FEM_index_zero_cell+current_cell,:);
         FEM_nodes.previous_position(FEM_index_zero_cell+current_cell,:) = 0;
         FEM_nodes.previous_position(cell_vertices,:) = 0;
         
         FEM_nodes.concentration(new_node) = FEM_nodes.concentration(FEM_index_zero_cell+current_cell);
         FEM_nodes.concentration(FEM_index_zero_cell+current_cell) = 0;
         FEM_nodes.concentration(cell_vertices) = 0;
         
         % loop over the vertices of the cell to construct the 3 edges, and
         % remove any edge nodes from any of these edges. these are very
         % unlikely to exist, as they require long edges on cells with small
         % areas, but can happen.
         for current_vertex_local = 1:3
            
            current_edge = [cell_vertices(current_vertex_local) ...
               cell_vertices(mod(current_vertex_local,3)+1)];
            
            [cells,FEM_elements,FEM_nodes,refined_edge_matrix] =...
               remove_FEM_node_from_specific_edge(cells,current_edge,...
               FEM_elements,FEM_nodes,refined_edge_matrix,vertices);
            
         end
         
         % see if any of the cell vertices are members of a refined edge.
         % for each vertex there is only one edge which could be refined as
         % the two edges that are part of the cell should not be refined,
         % and each vertex has a maximum of three edges
         cell_vertices_refined_edge_logical = ismember(cell_vertices,FEM_nodes.edge);
         
         % loop over the cell vertices and replace any occurrences of them
         % in FEM_nodes.edge and refined_edge_matrix with the new node.
         % there is almost certainly a more elegant way of doing this.
         for current_vertex_local = 1:3
            
            if cell_vertices_refined_edge_logical(current_vertex_local)
               
               current_vertex_global = cell_vertices(current_vertex_local);
               
               index_vertex_in_FEM_nodes_edge = find(FEM_nodes.edge==current_vertex_global);
               
               if index_vertex_in_FEM_nodes_edge <= length_FEM_nodes
                  
                  index_refined_edge_partner = index_vertex_in_FEM_nodes_edge + length_FEM_nodes;
                  refined_edge_partner = FEM_nodes.edge(index_refined_edge_partner);
                  
               else
                  
                  index_refined_edge_partner = index_vertex_in_FEM_nodes_edge - length_FEM_nodes;
                  refined_edge_partner = FEM_nodes.edge(index_refined_edge_partner);
                  
               end
               
               % replace cell vertex with new node in both FEM_nodes.edge
               % and refined_edge_matrix
               FEM_nodes.edge(index_vertex_in_FEM_nodes_edge) = new_node;
               refined_edge_matrix(current_vertex_global,refined_edge_partner) = 0;
               refined_edge_matrix(refined_edge_partner,current_vertex_global) = 0;
               
               refined_edge_matrix(refined_edge_partner,new_node) = 1;
               refined_edge_matrix(new_node,refined_edge_partner) = 1;
               
            end
            
         end
         
         FEM_elements.nodes(cells.FEM_elements{current_cell},:) = 0;
         cells.FEM_elements{current_cell} = [];
         
         % loop over the 3 cells that share an edge each with the dying cell
         for cell_with_same_edge_local = 1:3
            
            cell_with_same_edge_global = cells_with_same_edges(cell_with_same_edge_local);
            
            if cell_with_same_edge_global>0
               
               cell_elements = cells.FEM_elements{cell_with_same_edge_global};
               
               for current_element_local = 1:length(cell_elements)
                  
                  current_element_global = cell_elements(current_element_local);
                  current_element_nodes = FEM_elements.nodes(current_element_global,:);
                  
                  element_overlaps_dying_cell_logical =...
                     ismember(current_element_nodes,cell_vertices);
                  
                  if sum(element_overlaps_dying_cell_logical)==2
                     
                     element_to_remove_local = current_element_local;
                     element_to_remove_global = current_element_global;
                     % can't do this in the loop as it will change size of
                     % cell_elements and mess things up
                     %                         FEM_elements.nodes(current_element_global,:) = 0;
                     %                         cell_elements(current_element_local) = [];
                     %                         cells.FEM_elements{cell_with_same_edge_global} = cell_elements;
                     
                  else
                     
                     current_element_nodes(element_overlaps_dying_cell_logical) = new_node;
                     FEM_elements.nodes(current_element_global,:) = current_element_nodes;
                     
                  end
                  
               end
               
               FEM_elements.nodes(element_to_remove_global,:) = 0;
               cell_elements(element_to_remove_local) = [];
               cells.FEM_elements{cell_with_same_edge_global} = cell_elements;
               
               %                 cells_with_same_edges
               
            end
         end
      end
   end
   
   if death_location ~= 0
      
      stats.no_deaths_since_last_statistics = stats.no_deaths_since_last_statistics+1;
      death_counter = death_counter+1;
      stats.death_locations(death_counter,:) = death_location;
      
   end
   
end

if stats.this_iteration_logical
   
   stats.no_deaths(stats.counter) = stats.no_deaths_since_last_statistics;
   stats.no_deaths_since_last_statistics = 0;
   
end


