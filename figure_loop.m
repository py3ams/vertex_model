function figure_loop(cells,vertices,linewidth,normal_colour)

if nargin < 4
	normal_colour = 'r';
	if nargin < 3
		linewidth = 1;
	end
end

length_cells = length(cells.vertices);

patch_colour = repmat(normal_colour,length_cells,1); 
patch_colour(cells.original_logical,:) = normal_colour;
patch_colour(cells.state==3) = 'r';
patch_colour(cells.state==2) = 'r';

hold on
for current_cell =1:length_cells
    patchAS(vertices.position(cells.vertices{current_cell},:),patch_colour(current_cell),linewidth)
end

% cellfun(@(x)patchAS(vertices.position(x,:),patch_colour,linewidth),cells.vertices);

% if nargin > 2
% 	
% 	cellfun(@(x)patchAS(node_positions(x,:),original_cells_colour,linewidth),...
% 		cells(original_cells));
% 	
% 	cellfun(@(x)patchAS(node_positions(x,:),'r',linewidth),cells.vertices(cells.original_logical));
% 	
% else
% 	
% 	cellfun(@(x)patchAS(node_positions(x,:),'r',linewidth),cells);
% 	
% end
		