function centre_of_mass = calculate_centroid(cell_node_positions,cell_area)
	
no_cell_nodes = size(cell_node_positions,1);

cx = 0;
cy = 0;

for i = 1:no_cell_nodes
       
    x_i = cell_node_positions(i,1);
    x_iplus1 = cell_node_positions(mod(i,no_cell_nodes)+1,1);
    y_i = cell_node_positions(i,2);
    y_iplus1 = cell_node_positions(mod(i,no_cell_nodes)+1,2);
    
    cx = cx + (x_i+x_iplus1)*(x_i*y_iplus1-x_iplus1*y_i);
    cy = cy + (y_i+y_iplus1)*(x_i*y_iplus1-x_iplus1*y_i);
    
end

cx = cx/(6*cell_area);
cy = cy/(6*cell_area);

centre_of_mass = -[cx cy];