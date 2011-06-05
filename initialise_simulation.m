function [cells,fig_saves_location,figure_position,full_saves_location,...
	initial_mean_cell_area,initial_mean_edge_length,movie_location,stats,vertices] =...
    initialise_simulation(array_sizes,cells,fig_saves_logical,fig_saves_name,...
	full_saves_logical,full_saves_name,max_iterations,movie_logical,movie_name,...
	statistics_period,target_area_factor,target_no_cells,target_volume_factor,...
    vertices)
	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_statistics = floor(max_iterations/statistics_period);

no_cells = length(cells.vertices);

stats.total_no_cells = zeros(max_statistics,1);
stats.no_cells_per_vertex = zeros(max_statistics,3);
stats.total_boundary_vertices = zeros(max_statistics,1);

[cells.area,~,~,cell_area_stats,~,~,edge_length_stats] =...
    CalculateCellAreas(cells.vertices,vertices.position);

initial_mean_cell_area = cell_area_stats(1);
initial_mean_edge_length = edge_length_stats(1);

cells.target_area = target_area_factor*mean(cells.area)*ones(no_cells,1);
cells.target_volume = target_volume_factor*max(cells.volume)*ones(no_cells,1);

stats.total_no_vertices = zeros(max_statistics,1);
stats.vertices_per_cell = zeros(max_statistics,3);

vertices.time_created = -1e9*ones(array_sizes,1);
cells.time_of_last_division = -1e9*ones(no_cells,1);

stats.rosette = zeros(max_statistics,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stats.cell_area = zeros(max_statistics,5);
stats.cell_perimeter = zeros(max_statistics,4);
stats.shape_index = zeros(max_statistics,4);
stats.edge_length = zeros(max_statistics,3);

stats.cell_height = zeros(max_statistics,4);
stats.cell_height_to_area = zeros(max_statistics,4);
stats.cell_volume = zeros(max_statistics,4);

stats.angle_deviation = zeros(max_statistics,4);
stats.vertex_movement = zeros(max_statistics,2);

stats.volume_distribution_x = zeros(max_statistics,2);
stats.volume_distribution_y = zeros(max_statistics,2);

stats.area_distribution_x = zeros(max_statistics,2);
stats.area_distribution_y = zeros(max_statistics,2);

stats.normalised_volume_distribution_x = zeros(max_statistics,2);
stats.normalised_volume_distribution_y = zeros(max_statistics,2);

stats.normalised_area_distribution_x = zeros(max_statistics,2);
stats.normalised_area_distribution_y = zeros(max_statistics,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stats.death_locations = zeros(10*target_no_cells,2);
stats.mitosis_locations = zeros(10*target_no_cells,2);
stats.no_mitosis = zeros(max_statistics,1);
stats.no_T1_swaps = zeros(max_statistics,1);
stats.no_deaths = zeros(max_statistics,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stats.total_area_force = zeros(max_statistics,1);
stats.total_boundary_deformation_force = zeros(max_statistics,1);
stats.total_boundary_edge_force = zeros(max_statistics,1);
stats.total_deformation_force = zeros(max_statistics,1);
stats.total_elongation_force = zeros(max_statistics,1);
stats.total_perimeter_force = zeros(max_statistics,1);
stats.total_tension_force = zeros(max_statistics,1);
stats.total_force = zeros(max_statistics,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stats.no_refined_edges = zeros(max_statistics,1);
stats.mesh_peclet_number = zeros(max_statistics,1);
stats.node_speed = zeros(max_statistics,1);
stats.xy_ratio = zeros(max_statistics,1);
stats.triangle_quality = zeros(max_statistics,4);

stats.concentration_node_values = zeros(max_statistics,4);
stats.total_concentration = zeros(max_statistics,1);

% this records the total source released and the source released during
% individual iterations. the first value is not the cumsum of the second as
% we don't take statistics every iteration. the second value is a snapshot
% of the source being released per iteration, whereas the first value is an
% accurate record of the total released to that point. 
stats.chemical_source = zeros(max_statistics,2);

% this records mean etc per cell of internal chemical. this is very
% different from what is being recorded for the source as it is cumulative,
% i.e. it is the sum of all the chemical that has been ingested so far. it
% is not the amounts ingested in a given iteration (though these could be
% worked out).
stats.internal_chemical_quantity = zeros(max_statistics,5);
stats.chemical_in_cell_1 = zeros(max_statistics,1);
stats.dead_chemical = zeros(max_statistics,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

movie_location = 'Movies/';

if movie_logical == 2
	if exist([movie_location,movie_name,'.mpg'],'file')
		yesno = input(['Are you sure you want to overwrite file ',...
			movie_location,movie_name,'.mpg? '],'s');
		if strcmp(yesno,'yes') || strcmp(yesno,'y')
			delete([movie_location,movie_name,'.mpg']);
		else
			error(['Could not delete movie ',movie_name,'.mpg'])
		end
	end
	system(['touch ',movie_location,movie_name,'.mpg']);
end

screen_size = get(0,'screensize');
figure_position =...
	[(screen_size(3)-5/4*screen_size(4))/2 0 5/4*screen_size(4) screen_size(4)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

full_saves_location = ['Saves/',full_saves_name,'/'];

if full_saves_logical
	if exist(full_saves_location,'dir')
		yesno = input(['Are you sure you want to overwrite directory ',...
			full_saves_location,'? '],'s');
		if strcmp(yesno,'yes') || strcmp(yesno,'y')
			delete([full_saves_location,'*']);
		else
			error(['Could not delete directory ',full_saves_location]);
		end
	else
		mkdir(full_saves_location);
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig_saves_location = ['Figs/',fig_saves_name,'/'];

if fig_saves_logical
	if exist(fig_saves_location,'dir')
		yesno = input(['Are you sure you want to overwrite directory ',...
			fig_saves_location,'? '],'s');
		if strcmp(yesno,'yes') || strcmp(yesno,'y')
			delete([fig_saves_location,'*']);
		else
			error(['Could not delete directory ',fig_saves_location]);
		end
	else
		mkdir(fig_saves_location);
	end
end
