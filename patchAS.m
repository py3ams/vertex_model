function patchAS(cell_node_positions,colour,linewidth)

if nargin < 3
	linewidth = 1;
    if nargin < 2
        colour = 'r';
    end    
end

if size(cell_node_positions,2) == 2
	patch(cell_node_positions(:,1),cell_node_positions(:,2),colour,'linewidth',linewidth)
elseif size(cell_node_positions,2) == 3
	patch(cell_node_positions(:,1),cell_node_positions(:,2),...
		cell_node_positions(:,3),colour,'linewidth',linewidth)
end
