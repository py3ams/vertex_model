function p = division_probability(cell_parameter,cell_vertices,delta_t,...
	mitosis_random_logical,mitosis_threshold,target_parameter)

if mitosis_random_logical

	p = ((cell_parameter./target_parameter).^2)*delta_t/10;
	
	% cells smaller than the threshold are unable to divide.
	p(cell_parameter < mitosis_threshold*target_parameter) = 0;
	
else
	
	p = cell_parameter > mitosis_threshold*target_parameter;
	
end

if sum(cellfun('isempty',cell_vertices)&p>0)>0
	disp('warning: an empty cell has a non-zero division probability')
end
