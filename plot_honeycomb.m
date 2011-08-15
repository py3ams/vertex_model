disp('busy');clear all; %close all;

array_sizes = 1000;
configuration_noise = 0.0;
configuration_type = 'hexagonal';
grid_size = [10,12];
boxon = false;

[cell_vertices,vertex_positions,voronoi_points] =...
	initial_cell_mesh(array_sizes,configuration_noise,configuration_type,grid_size);

zero_entries_in_position_matrx = find(vertex_positions(:,1)==0);
% check that the first time we see a zero, the next entry is also zero. in
% other words there are no zeros in the middle of all the vertex positions.
assert(vertex_positions(zero_entries_in_position_matrx(1)+1,1)==0);

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
if boxon
    box on
else
    axis off
end
    
axis(1*[-1 1 -1 1])

figure('position',[75 45 1180 779],'color','white')
axes('position',[0.01 0.01 0.98 0.98],'linewidth',2)

for current_cell = 1:length(cell_vertices)
    hold on
    patchAS(vertex_positions(cell_vertices{current_cell},:),'w',10)
end

set(gca,'xtick',[])     
set(gca,'ytick',[])
set(gca,'ticklength',[0 0])
if boxon
    box on
else
    axis off
end

axis([-0.4460 0.4657 -0.0099 0.4344])