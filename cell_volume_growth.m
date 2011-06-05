function [cells,FEM_elements,FEM_nodes,stats,vertices] = cell_volume_growth(...
    cell_growth_concentration_dependent,cells,delta_t,FEM_elements,...
    FEM_nodes,lambda,no_growth_time,stats,time,vertices)

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
		
		cells.volume(cells.growing_logical) =...
			cells.volume(cells.growing_logical).*(1+delta_t*cells.growth_speed(...
			cells.growing_logical).*(1+lambda*cells.internal_chemical_quantity(...
			cells.growing_logical)).*(1-cells.volume(cells.growing_logical)./...
			cells.target_volume(cells.growing_logical)));
            
        
    end
    
else
    
    cells.volume(cells.growing_logical) =...
        cells.volume(cells.growing_logical).*(1+delta_t*cells.growth_speed(...
        cells.growing_logical).*(1-cells.volume(cells.growing_logical)./...
        cells.target_volume(cells.growing_logical)));
    
end


