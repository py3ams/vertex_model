function [cell_area_stats,cell_height_stats,cell_height_to_area_stats...
	cell_perimeter_stats,cell_volume_stats,cells_per_node_stats,...
	edge_length_stats,nodes_per_cell_stats,rosette_stats,shape_index_stats,...
	total_no_cells_stats,total_no_nodes_stats,xy_ratio_stats] = final_statistics(...
	cells,cells_per_node,cell_volumes,node_positions)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

no_cells = length(cells);
total_no_cells_stats = no_cells;
nodes_logical = cells_per_node > 0;
total_no_nodes_stats = sum(nodes_logical);
cells_per_node = cells_per_node(nodes_logical);
rosette_stats = [sum(cells_per_node==4) sum(cells_per_node >= 5)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[cell_areas,~,~,cell_area_stats,cell_perimeter_stats,shape_index_stats,...
	edge_length_stats] = CalculateCellAreas(cells,node_positions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cell_heights = cell_volumes./cell_areas;
height_area_ratios = cell_heights./cell_areas;
nodes_per_cell = cellfun('length',cells);

stat_values = MeanMaxMinStdTotal(cells_per_node,cell_heights,...
	cell_volumes,height_area_ratios,nodes_per_cell,...
	5);

cells_per_node_stats = stat_values(1,1:3);
cell_height_stats = stat_values(2,1:4);
cell_volume_stats = stat_values(3,1:4);
cell_height_to_area_stats = stat_values(4,1:4);
nodes_per_cell_stats = stat_values(5,1:3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xy_ratio_stats = (max(node_positions(:,1))-min(node_positions(:,1)))/...
    (max(node_positions(:,2))-min(node_positions(:,2)));

