function p = division_probability(cell_parameter,cell_vertices,delta_t,...
	division_period,mitosis_random_logical,mitosis_threshold,target_parameter)

if mitosis_random_logical

	p = (cell_parameter./target_parameter)*delta_t/division_period;
	
	% cells smaller than the threshold are unable to divide.
	p(cell_parameter < mitosis_threshold*target_parameter) = 0;
	
else
	
	p = cell_parameter > mitosis_threshold*target_parameter;
	
end

if sum(cellfun('isempty',cell_vertices)&p>0)>0
	disp('warning: an empty cell has a non-zero division probability')
end
