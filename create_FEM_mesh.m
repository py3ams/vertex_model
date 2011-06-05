function [FEM_elements,FEM_nodes,cell_elements] =...
	create_FEM_mesh(cells,node_positions)

no_cells = length(cells);
cell_elements = cell(no_cells,1);

cell_areas = CalculateCellAreas(cells,node_positions);

FEM_elements.nodes = zeros(sum(cellfun('length',cells)),3);

% the middle term is zeros we leave for possible mesh refinement. the number of edges
% in a triangular 
FEM_nodes.position = ...
    [node_positions; zeros(3*length(node_positions),2); zeros(no_cells,2)];

FEM_nodes.edge = zeros(4*length(node_positions)+no_cells,2);

no_elements = 0;
new_node = 4*length(node_positions);

for i = 1:no_cells
	
	cell_nodes = cells{i};
	
	cell_node_positions = node_positions(cell_nodes,:);
	cell_centre = CalculateCentroid(cell_node_positions,cell_areas(i));
	
	cell_elements{i} = no_elements+1:no_elements+length(cell_nodes);
	
	new_node = new_node+1;
	
	FEM_nodes.position(new_node,:) = cell_centre;
	
	for j = 1:length(cell_nodes)-1
		
		no_elements = no_elements+1;
		FEM_elements.nodes(no_elements,:) = [cell_nodes(j) cell_nodes(j+1) new_node];
		
	end
	
	no_elements = no_elements+1;
	FEM_elements.nodes(no_elements,:) = [cell_nodes(end) cell_nodes(1) new_node];
	
end
	