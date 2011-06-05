function [cells_containing_node_1,cells_containing_node_2,cell_with_same_edge] =...
    find_cells_containing_nodes(cell_store,current_cell,node_1,node_2)
    
% look for other cells containing the first node
cell_store_node_1 = cell_store(node_1,:);
cells_containing_node_1 =...
	cell_store_node_1((cell_store_node_1>0)&(cell_store_node_1~=current_cell));     
        
% look for all other cells containing the second node
cell_store_node_2 = cell_store(node_2,:);
cells_containing_node_2 =...
	cell_store_node_2((cell_store_node_2>0)&(cell_store_node_2~=current_cell));
        
% look if any cells share both nodes
cell_with_same_edge =...
	cells_containing_node_1(ismember(...
	cells_containing_node_1,cells_containing_node_2));

if ~isempty(cell_with_same_edge)

	cells_containing_node_1(cells_containing_node_1 == cell_with_same_edge) = [];
	cells_containing_node_2(cells_containing_node_2 == cell_with_same_edge) = [];

end