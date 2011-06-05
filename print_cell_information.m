function print_cell_information(cell_nos,cells)

for i = 1:length(cell_nos)
    
    current_cell = cell_nos(i);

    disp(' ');
    disp(['cell number: ',num2str(current_cell)])
    disp(['state: ',num2str(cells.state(current_cell))])
    disp(['vertices: ',num2str(cells.vertices{current_cell})])
    disp(['area: ',num2str(cells.area(current_cell))])
    disp(['target area: ',num2str(cells.target_area(current_cell))])
    disp(['volume: ',num2str(cells.volume(current_cell))])
    disp(['target volume: ',num2str(cells.target_volume(current_cell))])
    disp(['perimeter: ',num2str(cells.perimeter(current_cell))])
    disp(['edge lengths: ',num2str([cells.edge_lengths{current_cell}]')])
    disp(['shape index: ',num2str(cells.shape_indices(current_cell))])
    disp(['growing logical: ',num2str(cells.growing_logical(current_cell))])
    disp(['growth speed: ',num2str(cells.growth_speed(current_cell))])
    disp(['time of last division: ',num2str(cells.time_of_last_division(current_cell))])
    disp(['original logical: ',num2str(cells.original_logical(current_cell))])
    disp(['boundary logical: ',num2str(cells.boundary_logical(current_cell))])
    disp(['FEM elements: ',num2str(cells.FEM_elements{current_cell})])
    
    disp(' ');
    
end
