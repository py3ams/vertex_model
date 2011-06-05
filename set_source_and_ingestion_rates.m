function cells = set_source_and_ingestion_rates(cells,FEM_nodes,gradient_type,...
	no_chemicals,source_magnitude,source_width)

no_cells = length(cells.vertices);
cells.source_rate = zeros(no_cells,no_chemicals);
cells.ingestion_rate = zeros(no_cells,no_chemicals);
% apoptotic and dead cells cannot produce or ingest chemical
cells_logical = cells.state~=3 & cells.state~=4;
cell_centre_positions = FEM_nodes.previous_position(end-no_cells+1:end,:);

for current_chemical = 1:no_chemicals
	
	if gradient_type(current_chemical) == 1
		
		source_producing_logical = abs(cell_centre_positions(:,1))<(0.5*...
			source_width(current_chemical))&cells_logical;
		
	elseif gradient_type(current_chemical) == 2
		
		source_producing_logical = abs(cell_centre_positions(:,2))<(0.5*...
			source_width(current_chemical))&cells_logical;
		
	elseif gradient_type(current_chemical) == 3
		
		source_producing_logical = sqrt((cell_centre_positions(:,1).^2)+(...
			cell_centre_positions(:,2).^2))<(0.5*source_width(current_chemical))&...
            cells_logical;
		
	elseif gradient_type(current_chemical) == 4
		
		min_x_position = min(cell_centre_positions(:,1));
		source_producing_logical = abs(cell_centre_positions(:,1)-min_x_position)<...
			(0.5*source_width(current_chemical))&cells_logical;
		
	else
		
		error('Invalid gradient_type');
		
	end
	
% 	source_producing_logical = rand(no_cells,1)>0.9;
	
	cells.source_rate(source_producing_logical,current_chemical) =...
		source_magnitude(current_chemical);
    
    cells.ingestion_rate(cells_logical,current_chemical) = 0;

%     cells.ingestion_rate(~source_producing_logical,current_chemical) = 1;
	
%     no_FEM_nodes = size(FEM_nodes.position(:,1),1);
%     
%     cell_centre_positions = FEM_nodes.previous_position(no_FEM_nodes-no_cells+1:end,:);
%     cell_centre_positions = cell_centre_positions(cells_logical,:);
%     
%     shifted_cell_centre_positions_x = cell_centre_positions(:,1)-min(cell_centre_positions(:,1));
%     normalised_cell_centre_positions_x =...
%         shifted_cell_centre_positions_x/max(shifted_cell_centre_positions_x);
%     
%     cells.ingestion_rate(cells_logical,current_chemical) = 1e6*normalised_cell_centre_positions_x;
	
end

