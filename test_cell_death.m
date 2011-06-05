disp('busy');clear all; close all;

before = load('death_save');

[after.cells,death_counter,after.FEM_elements,after.FEM_nodes,no_deaths_this_iteration,...
    refined_edge_matrix,stats,after.vertices] = cell_death(before.cells,before.death_counter,...
    before.FEM_elements,before.FEM_nodes,before.FEM_solve_logical,before.protection_time,before.refined_edge_matrix,...
    before.stats,before.cell_death_area_threshold,before.time,before.vertices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axis_values = [-0.05 0.1 -0.2 -0.1];
axis_values_FEM = [axis_values 0 0.1];

% if you are using a new save (because mitosis has changed and no longer works with
% this save, need to plot view all cells and figure out where cell 25 is.
figure('units','normalized','outerposition',[0 0.25 1 0.5])
subplot(1,2,1)
figure_loop(before.cells,before.vertices)
subplot(1,2,2)
figure_loop(after.cells,after.vertices)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[100 100 300 300],'color','white')
axes('position',[0 0 1 1])
figure_loop(before.cells,before.vertices)
axis(axis_values)

figure('position',[100 100 300 300],'color','white')
axes('position',[0 0 1 1])
figure_loop(before.cells,before.vertices)
axis(axis_values)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[100 100 300 300],'color','white')
axes('position',[0 0 1 1])
% cellfun(@(x)patch(before.vertices.position(x,1),before.vertices.position(x,2),'r','linewidth',2.5,'FaceAlpha',0),before.cells.vertices);
% hold on;
trisurf(before.FEM_elements.nodes(before.FEM_elements.nodes(:,1)>0,:),before.FEM_nodes.position(:,1),...
    before.FEM_nodes.position(:,2),zeros(length(before.FEM_nodes.position),1))
% trisurf(FEM_element_nodes,before.FEM_nodes.previous_position(:,1),...
%     before.FEM_nodes.previous_position(:,2),before.FEM_nodes.concentration)
axis(axis_values_FEM)
view(0,90)

figure('position',[100 100 300 300],'color','white')
axes('position',[0 0 1 1])
% cellfun(@(x)patch(after.vertices.position(x,1),after.vertices.position(x,2),'r','linewidth',2.5,'FaceAlpha',0),after.cells.vertices);
% hold on;
trisurf(after.FEM_elements.nodes(after.FEM_elements.nodes(:,1)>0,:),after.FEM_nodes.position(:,1),...
    after.FEM_nodes.position(:,2),zeros(length(after.FEM_nodes.position),1))
% trisurf(FEM_element_nodes,after.FEM_nodes.previous_position(:,1),...
%     after.FEM_nodes.previous_position(:,2),after.FEM_nodes.concentration)
axis(axis_values_FEM)
view(0,90)