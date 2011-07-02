function [cells,node_positions,voronoi_points] = initial_cell_mesh(...
	array_sizes,configuration_noise,configuration_type,grid_size)
	
if nargin < 4
   configuration_type = 'square';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% set up the tessellation points %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialise arrays
node_positions = zeros(array_sizes,2);

% we make the initial width 1 (-0.5 to 0.5) so all distances are normalized to this
% initial width
x_max = 0.5;
spacing = 1/grid_size(1);

% find y_max based on grid size and spacing - there is probably a better
% way of doing this
y_max = floor(grid_size(2)/2)*spacing-0.5*spacing*rem(grid_size(2)+1,2);

% cell centres we want are from -x_max+1/2*spacing to x_max-1/2*spacing -
% we add three extra cells either side
x_coords = linspace(-x_max-2.5*spacing,x_max+2.5*spacing,grid_size(1)+6);
y_coords = linspace(-y_max-2.5*spacing,y_max+2.5*spacing,grid_size(2)+6)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% sort the coordinates out and do the tessellation %%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% by default tessellation is square based

x_coords_matrix = repmat(x_coords,length(y_coords),1);
y_coords_matrix = repmat(y_coords,1,length(x_coords));

% add random noise to the inner matrix excluding the outer two rows and
% columns. it is important to have one layer outside the cells that has
% random noise, otherwise vertices get created at exact locations like
% 0,0.1 etc. also note grid_size(2) relates to number of rows as it is the
% number of cells in the y-axis, and grid_size(1) relates to columns
size_inner_matrix = size(x_coords_matrix(3:grid_size(2)+4,3:grid_size(1)+4));

x_coords_matrix(3:grid_size(2)+4,3:grid_size(1)+4) = x_coords_matrix(3:...
   grid_size(2)+4,3:grid_size(1)+4) + spacing*configuration_noise*(...
   rand(size_inner_matrix)-0.5);
   
y_coords_matrix(3:grid_size(2)+4,3:grid_size(1)+4) = y_coords_matrix(3:...
   grid_size(2)+4,3:grid_size(1)+4) + spacing*configuration_noise*(...
   rand(size_inner_matrix)-0.5);

% for the random configuration, completely randomise all points excluding
% the outer three rows and columns of the matrix. note that in this case
% there will be one layer of semi-randomised points outside the cells, then
% two completely regular layers.
if strcmp(configuration_type,'random')

   size_matrix_to_randomise = size(x_coords_matrix(4:grid_size(2)+3,4:grid_size(1)+3));
   
   x_coords_matrix(4:grid_size(2)+3,4:grid_size(1)+3) = rand(size_matrix_to_randomise)-0.5;
   y_coords_matrix(4:grid_size(2)+3,4:grid_size(1)+3) = rand(size_matrix_to_randomise)-0.5;
   
elseif strcmp(configuration_type,'hexagonal')
    
   y_coords_matrix(:,1:2:end) = y_coords_matrix(:,1:2:end)+0.5*spacing;

end

% need to create this vector as we output it from the function
voronoi_points = [x_coords_matrix(:) y_coords_matrix(:)];

[temp_node_positions, temp_cells] = voronoin(voronoi_points);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% figure out which cells and nodes we actually want %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

no_cells = grid_size(1)*grid_size(2);
cells = cell(no_cells,1);

% cells_fixed_logical = false(5*length(temp_cells),1);

% keep a record of the nodes we actually use
node_store = zeros(length(temp_node_positions),1);

no_nodes = 0;
cell_counter = 0;

% loop through the right temp_cells using clever indexing
for i = 4:grid_size(1)+3
	
	first_column_cell = (i-1)*(grid_size(2)+6)+3;
	
	for j = 1:grid_size(2)
		
		cell_counter = cell_counter + 1;
		cell_nodes = temp_cells{first_column_cell + j};
		cells{cell_counter} = cell_nodes;
		
		new_nodes = cell_nodes(~ismember(cell_nodes,node_store));
		node_store(no_nodes+1:no_nodes+length(new_nodes)) = new_nodes;
		no_nodes = no_nodes + length(new_nodes);
		
	end
end

% remove blank entries in node store
node_store = node_store(node_store>0);

% store positions of nodes used
node_positions(1:no_nodes,:) = temp_node_positions(node_store,:);

% loop over all cells to give cells correct node indices
for current_cell = 1:no_cells
	
	old_cell_nodes = cells{current_cell};
	cell_nodes = zeros(1,length(old_cell_nodes));
	
	% check cell nodes are orientated clockwise, if not then rearrange
	% them
	temp_cell_node_positions = [temp_node_positions(old_cell_nodes,:);...
		temp_node_positions(old_cell_nodes(1),:)] + 100;
	
	centre_of_mass_x_component = 0;
	
	for j = 1:length(old_cell_nodes)
		
		x_j = temp_cell_node_positions(j,1);
		x_jplus1 = temp_cell_node_positions(j+1,1);
		y_j = temp_cell_node_positions(j,2);
		y_jplus1 = temp_cell_node_positions(j+1,2);
		
		centre_of_mass_x_component = centre_of_mass_x_component +...
			(x_j+x_jplus1)*(x_j*y_jplus1-x_jplus1*y_j);
		
		cell_nodes(j) = find(node_store==old_cell_nodes(j));
		
	end
	
	if centre_of_mass_x_component > 0
		cell_nodes = flipdim(cell_nodes,2);
	end
	
	% store cell
	cells{current_cell} = cell_nodes;

end
