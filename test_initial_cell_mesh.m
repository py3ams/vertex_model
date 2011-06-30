disp('busy');clear all; close all;

array_sizes = 1000;
configuration_noise = 0.5;
grid_size = [10,10];
tessellation_type = 'hexagonal';

[cell_vertices,vertex_positions,voronoi_points] =...
	initial_cell_mesh(array_sizes,configuration_noise,grid_size,tessellation_type);

figure('position',[200 200 500 500],'color','white')
% we put the axes just inside the figure so there is enough space for the
% box and it can be seen clearly within a matlab figure
axes('position',[0.01 0.01 0.98 0.98],'linewidth',2)

for current_cell = 1:length(cell_vertices)
    hold on
    patchAS(vertex_positions(cell_vertices{current_cell},:),'r',2)
end

plot(voronoi_points(:,1),voronoi_points(:,2),'bo','linewidth',3,'markersize',3)

set(gca,'xtick',[])     
set(gca,'ytick',[])
set(gca,'ticklength',[0 0])
box on

axis(0.6*[-1 1 -1 1])

figure('position',[200 200 500 500],'color','white')
axes('position',[0.01 0.01 0.98 0.98],'linewidth',2)

for current_cell = 1:length(cell_vertices)
    hold on
    patchAS(vertex_positions(cell_vertices{current_cell},:),'r',2)
end

set(gca,'xtick',[])     
set(gca,'ytick',[])
set(gca,'ticklength',[0 0])
box on

axis(0.55*[-1 1 -1 1])