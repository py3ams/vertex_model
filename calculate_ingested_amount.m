function cells = calculate_ingested_amount(cells,cell_ingestion_functions,...
	delta_t,FEM_elements,FEM_nodes,no_chemicals,previous_concentration)

% could also check cells.ingestion_rate>0? this would rely on ensuring that
% the ingestion rate is set to 0 when a cell becomes apoptotic. we should
% ensure that we are testing the right cell states here, AND that we set
% the ingestion rate to 0 when they move into that state. if we want, for
% example, apoptotic cells to still ingest, we should change the following
% logical statement, and change what happens to the ingestion rate in
% set_apoptotic_cells.
cells_ingesting_logical = cells.state~=4&cells.state~=3;

real_cell_areas =...
	CalculateCellAreas(cells.vertices(cells_ingesting_logical),FEM_nodes.position);
previous_real_cell_areas =...
	CalculateCellAreas(cells.vertices(cells_ingesting_logical),FEM_nodes.previous_position);

M = diag(real_cell_areas);

for current_chemical = 1:no_chemicals
	
	% find the amount of chemical in the same regions of space as each
	% cell. this is not the internal chemical.
	real_cell_concentrations = CalculateCellConcentrations(cells.FEM_elements(...
		cells_ingesting_logical),FEM_elements.nodes,previous_concentration(:,current_chemical),...
		FEM_nodes.previous_position);
	
	% need to look at the maths to see why the ingestion term is given by
	% this. see chemical_production_and_ingestion.tex.
	ingestion_term =...
		cell_ingestion_functions(cells_ingesting_logical,current_chemical).*real_cell_concentrations;
	
	rhs = previous_real_cell_areas.*cells.internal_chemical_value(...
		cells_ingesting_logical,current_chemical) + delta_t*ingestion_term;
	
	cells.internal_chemical_value(cells_ingesting_logical,current_chemical) = M\rhs;
	
	% the quantity is always given by value*cell_area. it is questionable
	% whether we need this as a separate property, but it is now used in
	% several places so best to keep it.
	cells.internal_chemical_quantity(cells_ingesting_logical,current_chemical) =...
		cells.internal_chemical_value(cells_ingesting_logical,current_chemical).*...
		real_cell_areas;
	
end