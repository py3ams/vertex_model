function [cells,dead_chemical,stats] = set_apoptotic_cells(cells,...
	apoptosis_baseline_prob_per_unit_time,apoptosis_concentration_dependent,...
	apoptosis_concentration_threshold,apoptosis_period,apoptosis_prob_above_threshold,...
	apoptosis_prob_below_threshold,apoptosis_type,dead_chemical,delta_t,FEM_solve_logical,stats,time)
	
length_cells = length(cells.vertices);

% 	if rand<0.001
% 		random_cell = ceil(rand*length_cells);
% 		while cells.boundary_logical(random_cell)
% 			random_cell = ceil(rand*length_cells);
% 		end
% 		cells.state(random_cell) = 3;
% 	end

if apoptosis_concentration_dependent

	probability_apoptosis(cells.internal_chemical_quantity<=apoptosis_concentration_threshold,1) = apoptosis_prob_below_threshold*delta_t;
	probability_apoptosis(cells.internal_chemical_quantity>apoptosis_concentration_threshold,1) = apoptosis_prob_above_threshold*delta_t;
	
	new_apoptotic_cells = rand(length_cells,1)<probability_apoptosis&~cells.boundary_logical&cells.state==1;
	
	
else
   
   if strcmp(apoptosis_type,'random')
	
      new_apoptotic_cells = rand(length_cells,1)<(apoptosis_baseline_prob_per_unit_time*delta_t)&...
         ~cells.boundary_logical&cells.state~=3;
      
   elseif strcmp(apoptosis_type,'regular')
      
      if abs(round(time/apoptosis_period)-time/apoptosis_period)<1e-6
         
         new_apoptotic_cells = find(~cells.boundary_logical&cells.state==1);
         new_apoptotic_cells = new_apoptotic_cells(ceil(rand*length(new_apoptotic_cells)));
         
      else
         
         new_apoptotic_cells = [];
         
      end
      
   end
	
end

cells.state(new_apoptotic_cells) = 3;
cells.volume(new_apoptotic_cells) = 0;

if FEM_solve_logical
	
	cells.ingestion_rate(new_apoptotic_cells) = 0;
	
	dead_chemical = dead_chemical + sum(cells.internal_chemical_quantity(...
		new_apoptotic_cells));
	
	cells.internal_chemical_value(new_apoptotic_cells) = 0;
	cells.internal_chemical_quantity(new_apoptotic_cells) = 0;
	
	if stats.this_iteration_logical
		
		stats.dead_chemical(stats.counter) = dead_chemical;
		
	end
	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
apoptotic_cells = cells.state==3;

% if FEM_solve_logical
cells.state(new_apoptotic_cells) = 3;
cells.volume(new_apoptotic_cells) = 0;

if FEM_solve_logical
	
	cells.ingestion_rate(new_apoptotic_cells) = 0;
	
	dead_chemical = dead_chemical + sum(cells.internal_chemical_quantity(...
		new_apoptotic_cells));
	
	cells.internal_chemical_value(new_apoptotic_cells) = 0;
	cells.internal_chemical_quantity(new_apoptotic_cells) = 0;
	
	if stats.this_iteration_logical
		
		stats.dead_chemical(stats.counter) = dead_chemical;
		
	end
	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
apoptotic_cells = cells.state==3;

% if FEM_solve_logical
%
% 	cells.internal_chemical_value(apoptotic_cells) =...
% 		0.99*cells.internal_chemical_value(apoptotic_cells);
%
% 	cells.internal_chemical_quantity(apoptotic_cells) =...
% 		0.99*cells.internal_chemical_quantity(apoptotic_cells);
%
% end

% alter the target area and force constants on apoptotic cells
cells.target_area(apoptotic_cells) = 0.999*cells.target_area(apoptotic_cells);
cells.volume(apoptotic_cells) = 0;

cells.force_constants.area(apoptotic_cells) = 0;
cells.force_constants.deformation(apoptotic_cells) = 0;
cells.force_constants.elongation(apoptotic_cells) = 0;
% cells.force_constants.elongation(apoptotic_cells) =...
% 	0.99^2*cells.force_constants.elongation(apoptotic_cells);
cells.force_constants.perimeter(apoptotic_cells) =...
	1/0.999*cells.force_constants.perimeter(apoptotic_cells);
cells.force_constants.tension(apoptotic_cells) =...
	1/0.999*cells.force_constants.tension(apoptotic_cells);
