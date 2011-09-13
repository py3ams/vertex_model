disp('busy');clear all; close all;

array_sizes = 1000;
configuration_noise = 0.0;
configuration_type = 'random';
grid_size = [10,10];
linewidth=3;

[cell_vertices,vertex_positions,voronoi_points] =...
	initial_cell_mesh(array_sizes,configuration_noise,configuration_type,grid_size);

zero_entries_in_position_matrx = find(vertex_positions(:,1)==0);
% check that the first time we see a zero, the next entry is also zero. in
% other words there are no zeros in the middle of all the vertex positions.
assert(vertex_positions(zero_entries_in_position_matrx(1)+1,1)==0);

figure('position',[100 100 500 500],'color','white','PaperPositionMode','auto')
axes('position',[0.01 0.01 0.98 0.98],'linewidth',2,'xcolor','white','ycolor','w',...
    'zcolor','w','ticklength',[0 0],'xtick',[],'ytick',[])

for current_cell = 1:length(cell_vertices)
    hold on
    patchAS(vertex_positions(cell_vertices{current_cell},:),'r',linewidth)
end

plot(voronoi_points(:,1),voronoi_points(:,2),'bo','linewidth',linewidth,'markersize',3)

axis(0.7*[-1 1 -0.995 1.015])
% axis(0.7*[-1 1 -1 1])

figure('position',[100 100 500 500],'color','white','PaperPositionMode','auto')
axes('position',[0.01 0.01 0.98 0.98],'linewidth',2,'xcolor','white','ycolor','w',...
    'zcolor','w','ticklength',[0 0],'xtick',[],'ytick',[])

for current_cell = 1:length(cell_vertices)
    hold on
    patchAS(vertex_positions(cell_vertices{current_cell},:),'r',linewidth)
end

axis(0.55*[-1 1 -0.95 1.05])