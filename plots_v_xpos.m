cell_centre_indices =...
    length(FEM_nodes.concentration)-(length(cells.vertices)-1:-1:0);
cell_centre_positions_x = FEM_nodes.position(cell_centre_indices,1);
a = polyfit(cell_centre_positions_x,cells.volume,1);
figure;
plot(cell_centre_positions_x,cells.volume,'x');
hold on;
plot(cell_centre_positions_x,a(1)*cell_centre_positions_x+a(2))