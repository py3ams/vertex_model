function [cell_growth_speeds_matrix,cells,FEM_elements,FEM_nodes,refined_edge_matrix,...
   vertices] = initial_configuration(anneal_initial_configuration_logical,...
   average_cell_growth_speeds,boundary_force_constants,configuration_noise,configuration_type,...
   gradient_type,grid_size,FEM_solve_logical,file_to_load,initial_concentration_magnitude,...
   initial_force_constant_magnitudes,load_FEM_from_file_logical,load_from_file_logical,...
   max_no_cells,medial_lateral_threshold_factor,no_chemicals,no_refinements,source_width)

if load_from_file_logical
   
   %     error('load from file not currently available');
   
   variables = load(file_to_load);
   
   cells = variables.cells;
   vertices = variables.vertices;
   
   cell_growth_speeds_matrix = variables.cell_growth_speeds_matrix;
   
   % can load cells but not FEM if we want to
   if load_FEM_from_file_logical
      
      FEM_elements = variables.FEM_elements;
      FEM_nodes = variables.FEM_nodes;
      refined_edge_matrix = variables.refined_edge_matrix;
      
   end
   
   max_no_vertices = length(vertices.position);
   
else
   
   % cannot load FEM from file if cells are not also being loaded
   
   if load_FEM_from_file_logical
      disp('Warning: setting load_FEM_from_file_logical to false as load_from_file_logical is false');
      load_FEM_from_file_logical = false;
   end
   
   max_no_vertices = 4*max_no_cells;
   
   [cells.vertices,vertices.position] =...
      initial_cell_mesh(max_no_vertices,configuration_noise,configuration_type,grid_size);
   
   no_cells = length(cells.vertices);
   
   [vertices.cells,vertices.no_cells] = CreateCellStore(cells.vertices,max_no_vertices);
   
   cells.area = CalculateCellAreas(cells.vertices,vertices.position);
   
   cells.original_logical = true(no_cells,1);
   
   % 1 = normal cell; 2 = baby cell; 3 = apoptotic cell; 4 = dead cell;
   cells.state = ones(no_cells,1);
   
   % we are using a bit of a fix here to identify cells on the boundary,
   % by first using the CreateBoundaryElement function to find the
   % vertices on the boundary, then finding cells containing any of those
   % vertices. the cellfun is very slow, which is why it was removed from
   % the main loop in cell_dynamics.
   boundary_vertices =...
      CreateBoundaryElement(vertices.cells,cells.vertices,vertices.no_cells);
   
   cells.boundary_logical = cellfun(@(x)any(ismember(x,boundary_vertices)),cells.vertices);
   
   mean_volume = mean(sqrt(cells.area.^3));
   
   cells.volume = max(mean_volume*(1+0.5*randn(no_cells,1)),...
      0.2*mean_volume*ones(no_cells,1));
   
   %     cells.volume = 2*mean_volume*rand(no_cells,1);
   
   %     cells.volume = mean_volume*ones(no_cells,1);
   
   cell_growth_speeds_matrix =...
      [2*average_cell_growth_speeds(1)*rand(10*max_no_cells,1) ...
      2*average_cell_growth_speeds(2)*rand(10*max_no_cells,1)];
   
   % 	cell_growth_speeds_matrix =...
   % 		[average_cell_growth_speeds(1)*ones(10*max_no_cells,1) ...
   % 		average_cell_growth_speeds(2)*ones(10*max_no_cells,1)];
   
   medial_cells = abs(cellfun(@(x)mean(vertices.position(x,1)),cells.vertices))<...
      medial_lateral_threshold_factor*max(vertices.position(:,1));
   
   lateral_cells = ~medial_cells;
   
   no_medial_cells = sum(medial_cells);
   no_lateral_cells = no_cells - no_medial_cells;
   
   cells.growth_speed(medial_cells,1) =...
      cell_growth_speeds_matrix(1:no_medial_cells,1);
   
   cells.growth_speed(lateral_cells,1) =...
      cell_growth_speeds_matrix(1:no_lateral_cells,2);
   
   cell_growth_speeds_matrix(1:no_medial_cells,1) = 0;
   cell_growth_speeds_matrix(1:no_lateral_cells,2) = 0;
   
   cells.force_constants.area =...
      initial_force_constant_magnitudes.area*ones(no_cells,1);
   cells.force_constants.deformation =...
      initial_force_constant_magnitudes.deformation*ones(no_cells,1);
   cells.force_constants.elongation =...
      initial_force_constant_magnitudes.elongation*ones(no_cells,1);
   cells.force_constants.perimeter =...
      initial_force_constant_magnitudes.perimeter*ones(no_cells,1);
   cells.force_constants.tension =...
      initial_force_constant_magnitudes.tension*ones(no_cells,1);
   
   cells.mitosis_flag = zeros(no_cells,1);
   
   % these are set for real in set_source_and_ingestion.m
   cells.source_rate = zeros(no_cells,no_chemicals);
   cells.ingestion_rate = zeros(no_cells,no_chemicals);
   
   %     figure
   %     for i =1:no_cells
   %         hold on
   %         patchAS(vertices.position(cells.vertices{i},:),'r',2)
   %     end
   
   if anneal_initial_configuration_logical
      
      delta_t = 1;
      viscosity = 1;
      tension_anisotropy_factor = 0.0;
      
      boundary_element =...
         CreateBoundaryElement(vertices.cells,cells.vertices,vertices.no_cells);
      
      for iter = 1:1000
         
         [cells.area,cells.perimeter,cells.edge_lengths,cell_area_stats,~,~,...
            edge_length_stats] = CalculateCellAreas(cells.vertices,vertices.position);
         
         mean_edge_length = edge_length_stats(1);
         cells.target_area = cell_area_stats(1)*ones(no_cells,1);
         
         %             figure
         %             for i =1:no_cells
         %                 hold on
         %                 patchAS(vertices.position(cells.vertices{i},:),'r',2)
         %             end
         
         vertices.position =...
            UpdatePos(cells.vertices,vertices.position,boundary_element,cells.area,...
            cells.perimeter,cells.volume,vertices.no_cells,delta_t,cells.edge_lengths,...
            mean_edge_length,cells.target_area,tension_anisotropy_factor,viscosity,...
            cells.force_constants.area,boundary_force_constants.deformation,...
            boundary_force_constants.edge,cells.force_constants.deformation,...
            cells.force_constants.elongation,cells.force_constants.perimeter,...
            cells.force_constants.tension);
         
         
      end
      
   end
   
end

% if load_FEM_from_file_logical is true at this point we will already have
% FEM_elements and FEM_nodes loaded
if FEM_solve_logical && ~load_FEM_from_file_logical
   
   [FEM_elements,FEM_nodes,cells.FEM_elements] =...
      create_FEM_mesh(cells,vertices,no_refinements);
   
   %         figure('position',[100 100 800 800])
   %         trisurf(FEM_elements.nodes,FEM_nodes.position(:,1),FEM_nodes.position(:,2),...
   %             zeros(length(FEM_nodes.position(:,1)),1),'linewidth',1);
   %         axis off;grid off;view([0 90]);axis equal
   %         hold on;
   %         for current_cell = 1:length(cells.vertices)
   %             hold on
   %             patchAS(vertices.position(cells.vertices{current_cell},:),[50,200,50]/256,4)
   %         end
   
   FEM_nodes.concentration =...
      initialise_concentration(vertices.no_cells,FEM_nodes,...
      gradient_type,initial_concentration_magnitude,...
      length(cells.vertices),no_chemicals,source_width);
   
   no_cells = length(cells.vertices);
   
   cells.internal_chemical_value = zeros(no_cells,no_chemicals);
   cells.internal_chemical_quantity = zeros(no_cells,no_chemicals);
   cells.maximum_internal_chemical_quantity = 5e-4*ones(no_cells,no_chemicals);
   
end

% if FEM_solve_logical is false all these variables should be 0 even if
% load_FEM_from_file_logical is true (i.e. we are overwriting previous
% declarations. if load_FEM_from_file_logical and FEM_solve_logical are false
% this will be the first time these variables have been set.
if ~FEM_solve_logical
   
   FEM_elements.nodes = [0 0 0; 0 0 0];
   FEM_nodes.position = [0 0];
   FEM_nodes.concentration = 0;
   cells.FEM_elements = {0};
   FEM_nodes.edge = [0 0];
   cells.internal_chemical_value = 0;
   cells.internal_chemical_quantity = 0;
   % 	Brk = ones(length(node_positions)+no_cells,1);
   
end

% these variables are set in this way regardless of whether FEM_solve_logical
% is true.
if ~load_FEM_from_file_logical
   refined_edge_matrix = sparse(max_no_vertices,max_no_vertices);
   %     refined_edge_matrix = sparse(1,1);
   FEM_nodes.previous_position = FEM_nodes.position;
end



