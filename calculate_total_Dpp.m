function total_Dpp = calculate_total_Dpp(Dpp,FEM_elements,gauss_precision,...
	FEM_node_positions)

no_elements = size(FEM_elements,1);      % Number of elements

total_Dpp = 0;

for current_element = 1:no_elements
	
	nodes_global = FEM_elements(current_element,:);                       
	element_node_positions = FEM_node_positions(nodes_global,:);
	element_Dpp = Dpp(nodes_global);
	
	jac = [element_node_positions(2,1)-element_node_positions(1,1), ...
		element_node_positions(3,1)-element_node_positions(1,1); ...
		element_node_positions(2,2)-element_node_positions(1,2), ...
		element_node_positions(3,2)-element_node_positions(1,2)];
	
	det_jac = det(jac);
	
	[Gauss_pointsx Gauss_pointsy Gauss_weights] = Gauss(gauss_precision);

	for i = 1:gauss_precision
	
		X = Gauss_pointsx(i);
		Y = Gauss_pointsy(i);
		W = Gauss_weights(i);
		
		total_Dpp = total_Dpp + W*det_jac*(element_Dpp(1)*(1-X-Y) +...
			element_Dpp(2)*X + element_Dpp(3)*Y);
		
	end
	
end