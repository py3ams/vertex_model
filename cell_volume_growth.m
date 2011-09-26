function [cells,FEM_elements,FEM_nodes,stats,vertices] = cell_volume_growth(...
   cell_growth_concentration_dependent,cells,concentration_dependence_type,delta_t,FEM_elements,...
   FEM_nodes,growth_solver_type,lambda,no_growth_time,stats,time,vertices)

% find cells that are not apoptotic or dead and have had sufficient time since
% division. could we not include this third condition in baby cells (i.e.
% cehck if cells.state==2)
cells.growing_logical =...
   time-cells.time_of_last_division > no_growth_time & cells.state~=3 & ...
   cells.state~=4;

if cell_growth_concentration_dependent
   
   if numel(FEM_nodes.concentration)==1
      
      error('Cell growth cannot be concentration dependent if you are not solving for concentration!')
      
   else
      
      if concentration_dependence_type == 1
         
         no_cells = length(cells.vertices);
         FEM_cell_centre_concentrations = FEM_nodes.concentration(end-no_cells:end);
         
         cells.volume(cells.growing_logical) =...
            cells.volume(cells.growing_logical).*(1+delta_t*cells.growth_speed(...
            cells.growing_logical).*(1+lambda*FEM_cell_centre_concentrations(...
            cells.growing_logical)).*(1-cells.volume(cells.growing_logical)./...
            cells.target_volume(cells.growing_logical)));
         
      elseif concentration_dependence_type == 2
         
         cells.volume(cells.growing_logical) =...
            cells.volume(cells.growing_logical).*(1+delta_t*cells.growth_speed(...
            cells.growing_logical).*(1+lambda*cells.internal_chemical_quantity(...
            cells.growing_logical)).*(1-cells.volume(cells.growing_logical)./...
            cells.target_volume(cells.growing_logical)));
         
      end
      
      
   end
   
else
   
   growing_cells_target_volume = cells.target_volume(cells.growing_logical);
   growing_cells_initial_volume = cells.initial_volume(cells.growing_logical);
   growing_cells_volume = cells.volume(cells.growing_logical);
   growing_cells_speed = cells.growth_speed(cells.growing_logical);
   
   if growth_solver_type == 1
      
      cells.volume(cells.growing_logical) =...
         growing_cells_volume.*(1+delta_t*growing_cells_speed.*...
         (1-growing_cells_volume./growing_cells_target_volume));
      
   elseif growth_solver_type == 2
      
      cells.volume(cells.growing_logical) =...
         growing_cells_initial_volume.*growing_cells_target_volume./...
         (growing_cells_initial_volume+(growing_cells_target_volume-...
         growing_cells_initial_volume).*exp(-growing_cells_speed*time));
      
   else
      
      error('Invalid growth_solver_type');
      
   end
   
end


