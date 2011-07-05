function statistical_plots(delta_t,fig_saves_location,fig_saves_logical,...
	figure_position,include_statistical_plots_in_movie,movie_location,...
	movie_logical,movie_name,no_frames_for_statistical_plots,...
   statistics_counter,stats,time)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_range = linspace(delta_t,time,statistics_counter);
frame_counter = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

total_no_cells = stats.total_no_cells(1:statistics_counter,:);
figure('outerposition',figure_position);
subplot(2,3,1)
stairs(time_range,total_no_cells);
xlabel('Time')
ylabel('Total number of cells')
axis([0 time min(total_no_cells)-2 max(total_no_cells)+2])

total_no_vertices = stats.total_no_vertices(1:statistics_counter,:);
subplot(2,3,2)
stairs(time_range,total_no_vertices);
xlabel('Time')
ylabel('Total number of vertices')
axis([0 time min(total_no_vertices)-2 max(total_no_vertices)+2])

no_cells_per_vertex = stats.no_cells_per_vertex(1:statistics_counter,:);
subplot(2,3,3)
plot(time_range,no_cells_per_vertex(:,1))
hold on
stairs(time_range,no_cells_per_vertex(:,2),'r')
stairs(time_range,no_cells_per_vertex(:,3),'r')
xlabel('Time')
ylabel('Mean no cells per vertex')
axis([0 time 0 max(no_cells_per_vertex(:,2))+2])

vertices_per_cell = stats.vertices_per_cell(1:statistics_counter,:);
subplot(2,3,4)
plot(time_range,vertices_per_cell(:,1))
hold on
stairs(time_range,vertices_per_cell(:,2),'r')
stairs(time_range,vertices_per_cell(:,3),'r')
xlabel('Time')
ylabel('Mean vertices per cell')
axis([0 time 0 max(vertices_per_cell(:,2))+2])

% total_Dpp = stats.total_Dpp(1:statistics_counter,:);
% subplot(2,3,5)
% plot(time_range,total_Dpp)
% xlabel('Time')
% ylabel('Total Dpp')
% if max(total_Dpp)>0
%     axis([0 time 0.9*min(total_Dpp) 1.1*max(total_Dpp)])
% else
%     axis_values = axis;
%     axis([axis_values(1) time axis_values(3) axis_values(4)])
% end
	
xy_ratio = stats.xy_ratio(1:statistics_counter,:);
subplot(2,3,5)
plot(time_range,xy_ratio)
xlabel('Time')
ylabel('XY ratio')
axis([0 time 0.9*min(xy_ratio) 1.1*max(xy_ratio)])

% mesh_peclet_number = stats.mesh_peclet_number(1:statistics_counter,:);
% subplot(2,3,6)
% plot(time_range,mesh_peclet_number)
% hold on
% xlabel('Time')
% ylabel('Maximum mesh peclet number')
% if max(mesh_peclet_number)>0
%     axis([0 time 0 1.1*max(mesh_peclet_number)])
% end

time_to_solve_FEM = stats.time_to_solve_FEM(1:statistics_counter,:);
subplot(2,3,6)
plot(total_no_cells,time_to_solve_FEM)
hold on
xlabel('No cells')
ylabel('Time to solve FEM')
% if max(mesh_peclet_number)>0
%     axis([0 time 0 1.1*max(mesh_peclet_number)])
% end

% add the figure to the movie if necessary
if movie_logical == 2
	M(frame_counter:frame_counter+no_frames_for_statistical_plots) = getframe(gcf);
	frame_counter = frame_counter+no_frames_for_statistical_plots+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('outerposition',figure_position);

axis_distribution_gradients = [0 time -0.2 0.2];

% min_gradient_value = min([min(stats.normalised_volume_distribution_x(:,1)) ...
%     min(stats.normalised_volume_distribution_y(:,1)) min(stats.normalised_area_distribution_x(:,1)) ...
%     min(stats.normalised_area_distribution_y(:,1))]);
% 
% max_gradient_value = max([max(stats.normalised_volume_distribution_x(:,1)) ...
%     max(stats.normalised_volume_distribution_y(:,1)) max(stats.normalised_area_distribution_x(:,1)) ...
%     max(stats.normalised_area_distribution_y(:,1))]);

% if min_gradient_value < 0
%     min_axis_value = 1.1*min_gradient_value;
% else
%     min_axis_value = 0.9*min_gradient_value;
% end
% 
% if max_gradient_value < 0
%     max_axis_value = 0.9*max_gradient_value;
% else
%     max_axis_value = 1.1*max_gradient_value;
% end

gradient_volume_distribution_x = stats.normalised_volume_distribution_x(:,1);
subplot(2,2,1)
plot(time_range,gradient_volume_distribution_x)
xlabel('Time')
ylabel('Volume distribution gradient (x)')
% axis(axis_distribution_gradients)
a = axis;
axis([a(1) time a(3) a(4)]) 

gradient_volume_distribution_y = stats.normalised_volume_distribution_y(:,1);
subplot(2,2,2)
plot(time_range,gradient_volume_distribution_y)
xlabel('Time')
ylabel('Volume distribution gradient (y)')
% axis(axis_distribution_gradients)
a = axis;
axis([a(1) time a(3) a(4)]) 

gradient_area_distribution_x = stats.normalised_area_distribution_x(:,1);
subplot(2,2,3)
plot(time_range,gradient_area_distribution_x)
xlabel('Time')
ylabel('Area distribution gradient (x)')
% axis(axis_distribution_gradients)
a = axis;
axis([a(1) time a(3) a(4)]) 

gradient_area_distribution_y = stats.normalised_area_distribution_y(:,1);
subplot(2,2,4)
plot(time_range,gradient_area_distribution_y)
xlabel('Time')
ylabel('Area distribution gradient (y)')
% axis(axis_distribution_gradients)
a = axis;
axis([a(1) time a(3) a(4)]) 

if movie_logical == 2
	M(frame_counter:frame_counter+no_frames_for_statistical_plots) = getframe(gcf);
	frame_counter = frame_counter+no_frames_for_statistical_plots+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('outerposition',figure_position);

cell_area = stats.cell_area(1:statistics_counter,:);
subplot(2,3,1)
plot(time_range,cell_area(:,1))
hold on
% plot(time_range,cell_area(:,2),'r')
% plot(time_range,cell_area(:,3),'r')
plot(time_range,cell_area(:,1)+2*cell_area(:,4),'r')
plot(time_range,cell_area(:,1)-2*cell_area(:,4),'r')
% plot(time_range,target_area*ones(size(time_range)),'g--')
xlabel('Time')
ylabel('Mean cell area')
axis([0 time min(cell_area(:,1)-2.5*cell_area(:,4)) ...
	max(cell_area(:,1)+2.5*cell_area(:,4))])

subplot(2,3,2)
plot(time_range,cell_area(:,5))
xlabel('Time')
ylabel('Total area')
axis([0 time 0.9*min(cell_area(:,5)) 1.1*max(cell_area(:,5))])

cell_height = stats.cell_height(1:statistics_counter,:);
subplot(2,3,3)
plot(time_range,cell_height(:,1))
hold on
% plot(time_range,cell_height(:,2),'r')
% plot(time_range,cell_height(:,3),'r')
plot(time_range,cell_height(:,1)+2*cell_height(:,4),'r')
plot(time_range,cell_height(:,1)-2*cell_height(:,4),'r')
xlabel('Time')
ylabel('Mean cell height')
axis([0 time min(cell_height(:,1)-2.5*cell_height(:,4)) ...
	max(cell_height(:,1)+2.5*cell_height(:,4))])

cell_height_to_area = stats.cell_height_to_area(1:statistics_counter,:);
subplot(2,3,4)
plot(time_range,cell_height_to_area(:,1))
hold on
% plot(time_range,cell_height_to_area(:,2),'r')
% plot(time_range,cell_height_to_area(:,3),'r')
plot(time_range,cell_height_to_area(:,1)+2*cell_height_to_area(:,4),'r')
plot(time_range,cell_height_to_area(:,1)-2*cell_height_to_area(:,4),'r')
xlabel('Time')
ylabel('Height area ratio')
axis([0 time ...
	min(cell_height_to_area(:,1)-2.5*cell_height_to_area(:,4)) ...
	max(cell_height_to_area(:,1)+2.5*cell_height_to_area(:,4))])

cell_perimeter = stats.cell_perimeter(1:statistics_counter,:);
subplot(2,3,5)
plot(time_range,cell_perimeter(:,1))
hold on
% plot(time_range,cell_perimeter(:,2),'r')
% plot(time_range,cell_perimeter(:,3),'r')
plot(time_range,cell_perimeter(:,1)+2*cell_perimeter(:,4),'r')
plot(time_range,cell_perimeter(:,1)-2*cell_perimeter(:,4),'r')
xlabel('Time')
ylabel('Mean cell perimeter')
axis([0 time min(cell_perimeter(:,1)-2.5*cell_perimeter(:,4)) ...
	max(cell_perimeter(:,1)+2.5*cell_perimeter(:,4))])

cell_volume = stats.cell_volume(1:statistics_counter,:);
subplot(2,3,6)
plot(time_range,cell_volume(:,1))
hold on
% plot(time_range,cell_volume(:,2),'r')
% plot(time_range,cell_volume(:,3),'r')
plot(time_range,cell_volume(:,1)+2*cell_volume(:,4),'r')
plot(time_range,cell_volume(:,1)-2*cell_volume(:,4),'r')
% plot(time_range,target_volume*ones(size(time_range)),'g--')
xlabel('Time')
ylabel('Mean cell volume')
if(max(cell_volume(:,4)>0))
	axis([0 time min(cell_volume(:,1)-2.5*cell_volume(:,4)) ...
		max(cell_volume(:,1)+2.5*cell_volume(:,4))])
end

if movie_logical == 2
	M(frame_counter:frame_counter+no_frames_for_statistical_plots) = getframe(gcf);
	frame_counter = frame_counter+no_frames_for_statistical_plots+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('outerposition',figure_position);

edge_length = stats.edge_length(1:statistics_counter,:);
subplot(2,3,1)
plot(time_range,edge_length(:,1))
hold on
plot(time_range,edge_length(:,2),'r')
plot(time_range,edge_length(:,3),'r')
xlabel('Time')
ylabel('Mean edge length')
axis([0 time 0 1.1*max(edge_length(:,2))])

shape_index = stats.shape_index(1:statistics_counter,:);
subplot(2,3,2)
plot(time_range,shape_index(:,1))
hold on
% plot(time_range,shape_index(:,2),'r')
% plot(time_range,shape_index(:,3),'r')
plot(time_range,shape_index(:,1)+2*shape_index(:,4),'r')
plot(time_range,shape_index(:,1)-2*shape_index(:,4),'r')
xlabel('Time')
ylabel('Mean shape index')
axis([0 time min(shape_index(:,1)-2.5*shape_index(:,4)) ...
	max(shape_index(:,1)+2.5*shape_index(:,4))])

angle_deviation = stats.angle_deviation(1:statistics_counter,:);
subplot(2,3,3)
plot(time_range,angle_deviation(:,1))
hold on
% plot(time_range,angle_deviation(:,2),'r')
% plot(time_range,angle_deviation(:,3),'r')
plot(time_range,angle_deviation(:,1)+2*angle_deviation(:,4),'r')
plot(time_range,angle_deviation(:,1)-2*angle_deviation(:,4),'r')
xlabel('Time')
ylabel('Angle deviation')
if max(angle_deviation(:,1))>0
	axis([0 time min(angle_deviation(:,1)-2.5*angle_deviation(:,4)) ...
		max(angle_deviation(:,1)+2.5*angle_deviation(:,4))])
end

vertex_movement = stats.vertex_movement(1:statistics_counter,:);
subplot(2,3,4)
plot(time_range,vertex_movement(:,1))
hold all
plot(time_range,vertex_movement(:,2),'r')
xlabel('Time')
ylabel('Absolute vertex movements')
axis([0 time 0 1.1*max(vertex_movement(:,2))])

subplot(2,3,5)
relative_mean_vertex_movement = vertex_movement(:,1)./edge_length(:,1);
plot(time_range,relative_mean_vertex_movement)
xlabel('Time')
ylabel('Relative mean vertex movement')
axis([0 time 0 1.1*max(relative_mean_vertex_movement)])

subplot(2,3,6)
relative_maximum_vertex_movement = vertex_movement(:,2)./edge_length(:,3);
plot(time_range,relative_maximum_vertex_movement)
xlabel('Time')
ylabel('Relative maximum vertex movement')
axis([0 time 0 1.1*max(relative_maximum_vertex_movement)])


if movie_logical == 2
	M(frame_counter:frame_counter+no_frames_for_statistical_plots) = getframe(gcf);
	frame_counter = frame_counter+no_frames_for_statistical_plots+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('outerposition',figure_position);

no_mitosis = stats.no_mitosis(1:statistics_counter,:);
subplot(2,3,1)
stairs(time_range,no_mitosis)
xlabel('Time')
ylabel('Number of cells experiencing mitosis')
axis([0 time 0 max(no_mitosis)+2])
	
% mitosis_locations = stats.mitosis_locations(stats.mitosis_locations(:,1)~=0,:);
% no_bars = 10;
% subplot(2,3,2)
% % [hist_data,bin_centres] = hist(mitosis_locations(:,1),no_bars);
% [hist_data,bin_centres] =...
%     hist(sqrt(mitosis_locations(:,1).^2+mitosis_locations(:,2).^2),no_bars);
% hist_data = hist_data/length(mitosis_locations)*100;
% bar(bin_centres,hist_data);
% xlabel('Mitosis location')

% mitosis_radii = sqrt(mitosis_locations(:,1).^2+mitosis_locations(:,2).^2);
% mitosis_radii_within_original_radius = mitosis_radii(mitosis_radii<0.5);
% no_bars = 10;
% subplot(2,3,3)
% [hist_data,bin_centres] =...
% 	hist(mitosis_radii_within_original_radius,no_bars);
% hist_data = hist_data./bin_centres;
% hist_data = hist_data./sum(hist_data)*100;
% % hist_data = hist_data/length(mitosis_locations)*100;
% % radii_bin_centres = linspace(0,0.5,11);
% bar(bin_centres,hist_data);
% xlabel('Normalised mitosis locations')
% axis([0 0.5 0 20])

mitosis_locations = stats.mitosis_locations(stats.mitosis_locations(:,1)~=0,:);
no_bars = 10;
subplot(2,3,2)
[hist_data,bin_centres] = hist(mitosis_locations(:,1),no_bars);
hist_data = hist_data/length(mitosis_locations)*100;
bar(bin_centres,hist_data);
xlabel('Mitosis location')
% axis([-0.5 0.5 0 20])

% mitosis_radii = sqrt(mitosis_locations(:,1).^2+mitosis_locations(:,2).^2);
% mitosis_within_original_radius = mitosis_locations(mitosis_radii<0.5,1);
% if sum(mitosis_within_original_radius)>1
%     no_bars = 10;
%     subplot(2,3,3)
%     [hist_data,bin_centres] =...
%         hist(mitosis_within_original_radius,no_bars);
%     hist_data = hist_data./(2*sqrt(0.5^2 - bin_centres.^2));
%     hist_data = hist_data./sum(hist_data)*100;
%     % hist_data = hist_data/length(mitosis_locations)*100;
%     % radii_bin_centres = linspace(0,0.5,11);
%     bar(bin_centres,hist_data);
%     xlabel('Normalised mitosis locations')
%     axis([-0.5 0.5 0 20])
% end

no_deaths = stats.no_deaths(1:statistics_counter,:);
subplot(2,3,3)
stairs(time_range,no_deaths)
xlabel('Time')
ylabel('Number of cell deaths')
axis([0 time 0 max(no_deaths)+2])

no_T1_swaps = stats.no_T1_swaps(1:statistics_counter,:);
subplot(2,3,4)
stairs(time_range,no_T1_swaps)
xlabel('Time')
ylabel('Number of T1 transitions')
axis([0 time 0 max(no_T1_swaps)+2])

triangle_quality_stats = stats.triangle_quality(1:statistics_counter,:);
subplot(2,3,5)
plot(time_range,triangle_quality_stats(:,1))
hold on
% plot(time_range,triangle_quality_stats(:,2),'r')
% plot(time_range,triangle_quality_stats(:,3),'r')
plot(time_range,triangle_quality_stats(:,1)+2*triangle_quality_stats(:,4),'r')
% plot(time_range,triangle_quality_stats(:,1)-2*triangle_quality_stats(:,4),'r')
xlabel('Time')
ylabel('Triangle Quality')
if max(triangle_quality_stats(:,1))>0
    axis([0 time 0 1.1*max(triangle_quality_stats(:,1)+2*triangle_quality_stats(:,4))])
end

death_locations = stats.death_locations(stats.death_locations(:,1)~=0,:);
no_bars = 10;
subplot(2,3,6)
[hist_data,bin_centres] = hist(death_locations(:,1),no_bars);
hist_data = hist_data/length(death_locations)*100;
bar(bin_centres,hist_data);
xlabel('Death location')

if movie_logical == 2
	M(frame_counter:frame_counter+no_frames_for_statistical_plots) = getframe(gcf);
	frame_counter = frame_counter+no_frames_for_statistical_plots+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('outerposition',figure_position);

total_concentration = stats.total_concentration(1:statistics_counter,:);
total_chemical_released = stats.chemical_source(1:statistics_counter,1);
total_internal_chemical = stats.internal_chemical_quantity(1:statistics_counter,5);
total_dead_chemical = stats.dead_chemical(1:statistics_counter);

initial_concentration_in_system = total_concentration(1);
initial_internal_chemical = stats.internal_chemical_quantity(1,5);

chemical_internalised = total_internal_chemical-initial_internal_chemical;

% dead chemical starts at 0 for a new simulation so don't need to subtract
% initial value
net_chemical = initial_concentration_in_system + total_chemical_released -...
	chemical_internalised - total_dead_chemical;

subplot(2,3,1)
plot(time_range,total_chemical_released)
hold all
plot(time_range,total_internal_chemical)
plot(time_range,total_dead_chemical)
plot(time_range,net_chemical)
xlabel('Time')
ylabel('Net chemical released')
legend('Source','Internal','Dead','Net','Location','Northwest')
% if max(net_chemical)>0
%     axis([0 time 0.9*min(net_chemical) 1.1*max(net_chemical)])
% else
    axis_values = axis;
    axis([axis_values(1) time axis_values(3) axis_values(4)])
% end
a = axis;

subplot(2,3,2)
plot(time_range,total_concentration)
xlabel('Time')
ylabel('Total concentration in system')
% if max(total_concentration)>0
%     axis([0 time 0.9*min(total_concentration) 1.1*max(total_concentration)])
% else
%     axis_values = axis;
%     axis([axis_values(1) time axis_values(3) axis_values(4)])
% end
axis(a)

% concentration_released_per_iteration = stats.chemical_source(1:statistics_counter,2);
% subplot(2,3,3)
% plot(time_range,concentration_released_per_iteration)
% xlabel('Time')
% ylabel('Current source term')
% if max(concentration_released_per_iteration)>0
%     axis([0 time 0.9*min(concentration_released_per_iteration) 1.1*max(concentration_released_per_iteration)])
% else
%     axis_values = axis;
%     axis([axis_values(1) time axis_values(3) axis_values(4)])
% end

mean_internal_chemical_quantity = stats.internal_chemical_quantity(1:statistics_counter,1);
max_internal_chemical_quantity = stats.internal_chemical_quantity(1:statistics_counter,2);
min_internal_chemical_quantity = stats.internal_chemical_quantity(1:statistics_counter,3);
std_internal_chemical_quantity = stats.internal_chemical_quantity(1:statistics_counter,4); 
subplot(2,3,3)
plot(time_range,mean_internal_chemical_quantity)
hold on
plot(time_range,max_internal_chemical_quantity,'r')
plot(time_range,min_internal_chemical_quantity,'r')
plot(time_range,mean_internal_chemical_quantity+2*std_internal_chemical_quantity,'g')
% plot(time_range,mean_internal_chemical_quantity-2*std_internal_chemical_quantity,'g')
xlabel('Time')
ylabel('Mean internal chemical quantity')
if max(min_internal_chemical_quantity)>0
    axis([0 time 0.9*min(min_internal_chemical_quantity) 1.1*max(max_internal_chemical_quantity)])
else
    axis_values = axis;
    axis([axis_values(1) time axis_values(3) axis_values(4)])
end

mean_concentration = stats.concentration_node_values(1:statistics_counter,1);
max_concentration = stats.concentration_node_values(1:statistics_counter,2);
min_concentration = stats.concentration_node_values(1:statistics_counter,3);
subplot(2,3,4)
plot(time_range,mean_concentration)
hold on
plot(time_range,max_concentration,'r')
plot(time_range,min_concentration,'r')
xlabel('Time')
ylabel('Mean node concentration value')
if max(min_concentration)>0
    axis([0 time 0.9*min(min_concentration) 1.1*max(max_concentration)])
else
    axis_values = axis;
    axis([axis_values(1) time axis_values(3) axis_values(4)])
end

no_refined_edges = stats.no_refined_edges(1:statistics_counter);
subplot(2,3,5)
plot(time_range,no_refined_edges)
xlabel('Time')
ylabel('Number of refined edges')
if max(no_refined_edges)>0
    axis([0 time 0 1.1*max(no_refined_edges)])
else
    axis_values = axis;
    axis([axis_values(1) time axis_values(3) axis_values(4)])
end

chemical_in_cell_1 = stats.chemical_in_cell_1(1:statistics_counter);
subplot(2,3,6)
plot(time_range,chemical_in_cell_1)
xlabel('Time')
ylabel('Chemical in cell 1')
if max(chemical_in_cell_1)>0
    axis([0 time 0 1.1*max(chemical_in_cell_1)])
else
    axis_values = axis;
    axis([axis_values(1) time axis_values(3) axis_values(4)])
end

% mean_chemical_ingestion_rate = stats.chemical_ingestion_rate(1:statistics_counter,1);
% max_chemical_ingestion_rate = stats.chemical_ingestion_rate(1:statistics_counter,2);
% min_chemical_ingestion_rate = stats.chemical_ingestion_rate(1:statistics_counter,3);
% subplot(2,3,6)
% plot(time_range,mean_chemical_ingestion_rate)
% hold on
% plot(time_range,max_chemical_ingestion_rate,'r')
% plot(time_range,min_chemical_ingestion_rate,'r')
% xlabel('Time')
% ylabel('Mean chemical ingested rate')
% if min(max_chemical_ingestion_rate)>0
%     axis([0 time 0.9*min(min_chemical_ingestion_rate) 1.1*max(max_chemical_ingestion_rate)])
% else
%     axis_values = axis;
%     axis([axis_values(1) time axis_values(3) axis_values(4)])
% end

if movie_logical == 2
	M(frame_counter:frame_counter+no_frames_for_statistical_plots) = getframe(gcf);
	frame_counter = frame_counter+no_frames_for_statistical_plots+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('outerposition',figure_position);

total_area_force = stats.total_area_force(1:statistics_counter,:);
total_boundary_deformation_force =...
	stats.total_boundary_deformation_force(1:statistics_counter,:);
total_boundary_edge_force = stats.total_boundary_edge_force(1:statistics_counter,:);
total_deformation_force = stats.total_deformation_force(1:statistics_counter,:);
total_elongation_force = stats.total_elongation_force(1:statistics_counter,:);
total_perimeter_force = stats.total_perimeter_force(1:statistics_counter,:);
total_tension_force = stats.total_tension_force(1:statistics_counter,:);
total_force = stats.total_force(1:statistics_counter,:);
total_boundary_vertices = stats.total_boundary_vertices(1:statistics_counter,:);

plot(time_range,total_area_force./total_no_vertices)
hold all
plot(time_range,total_boundary_deformation_force./total_boundary_vertices)
plot(time_range,total_boundary_edge_force./total_boundary_vertices)
plot(time_range,total_deformation_force./total_no_vertices)
plot(time_range,total_elongation_force./total_no_vertices)
plot(time_range,total_perimeter_force./total_no_vertices)
plot(time_range,total_tension_force./total_no_vertices)
plot(time_range,total_force./total_no_vertices)
axis_values = axis;
axis([axis_values(1) time axis_values(3) axis_values(4)])
xlabel('Time')
ylabel('Average force per node')
legend('Area','Boundary Deformation','Boundary Edge','Deformation','Elongation',...
	'Perimeter','Tension','Total')

if movie_logical == 2 && include_statistical_plots_in_movie
	M(frame_counter:frame_counter+no_frames_for_statistical_plots) = getframe(gcf);
	mpgwrite(M,jet,[movie_location,'update.mpg']);
	system(['cat ',movie_location,movie_name,'.mpg ',movie_location,...
		'update.mpg > ',movie_location,'temp.mpg']);
	system(['mv ',movie_location,'temp.mpg ',movie_location,movie_name,'.mpg']);
	system(['rm ',movie_location,'update.mpg']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fig_saves_logical
	
	figurecount = length(findobj('Type','figure'));
	
	for currentfig = 1:figurecount
		saveas(currentfig,[fig_saves_location,'Fig',num2str(currentfig)]);
	end
	
end
