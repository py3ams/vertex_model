close all

c = ['b','r','y','k','c'];
for i = 1:no_cells
	cell_nodes = cells{i};
	for j = 1:length(cell_nodes)
		plot(node_positions(cell_nodes(j),1),node_positions(cell_nodes(j),2),[c(mod(d,5)+1),'o'])
		hold on
		axis([-0.2 0.6 0 0.6])
		pause(0.1)
	end
	pause(0.3)
	d = d+1;
end
