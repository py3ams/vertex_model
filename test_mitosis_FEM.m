disp('busy');clear all; close all;

% n.b. need to set cells to divide to 25 in mitosis.m for this to work

before = load('mitosis_save');
no_cells = length(before.cells.vertices);

mitosis_dependence = 'none';
mitosis_period = 0.001;

[after.cells,after.vertices,after.cell_growth_speeds_matrix,after.FEM_elements,after.FEM_nodes,...
    after.mitosis_counter,after.refined_edge_matrix,after.stats] = mitosis(before.cells,before.vertices,...
    before.cell_growth_speeds_matrix,before.delta_t,before.FEM_elements,before.FEM_nodes,...
    before.FEM_solve_logical,before.medial_lateral_threshold,before.mitosis_angles_type,...
    before.mitosis_counter,mitosis_dependence,mitosis_period,before.mitosis_random_logical,...
    before.mitosis_threshold,before.refined_edge_matrix,before.stats,before.time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axis_values = [-0.35 -0.1 -0.17 0.08];
axis_values_FEM = [axis_values 0 0.1];

% if you are using a new save (because mitosis has changed and no longer works with
% this save, need to plot view all cells and figure out where cell 25 is.
% figure
% subplot(1,2,1)
% figure_loop(before.cells,before.vertices)
% subplot(1,2,2)
% figure_loop(after.cells,after.vertices)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[100 100 300 300],'color','white')
axes('position',[0 0 1 1])
figure_loop(before.cells,before.vertices)
axis(axis_values)

figure('position',[100 100 300 300],'color','white')
axes('position',[0 0 1 1])
figure_loop(after.cells,after.vertices)
axis(axis_values)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[100 100 300 300],'color','white')
axes('position',[0 0 1 1])
trisurf(before.FEM_elements.nodes(before.FEM_elements.nodes(:,1)>0,:),before.FEM_nodes.previous_position(:,1),...
    before.FEM_nodes.previous_position(:,2),zeros(length(before.FEM_nodes.previous_position),1),'linewidth',2.5)
hold on;
cellfun(@(x)patch(before.vertices.position(x,1),before.vertices.position(x,2),'r','linewidth',5,'FaceAlpha',0),before.cells.vertices);
axis([axis_values 0 0.1])
axis off
view(0,90)

figure('position',[100 100 300 300],'color','white')
axes('position',[0 0 1 1])
trisurf(after.FEM_elements.nodes(after.FEM_elements.nodes(:,1)>0,:),after.FEM_nodes.previous_position(:,1),...
    after.FEM_nodes.previous_position(:,2),zeros(length(after.FEM_nodes.previous_position),1),'linewidth',2.5)
hold on;
cellfun(@(x)patch(after.vertices.position(x,1),after.vertices.position(x,2),'r','linewidth',5,'FaceAlpha',0),after.cells.vertices);
axis([axis_values 0 0.1])
axis off
view(0,90)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[100 100 300 300])
axes('position',[0 0 1 1])

trisurf(before.FEM_elements.nodes(before.cells.FEM_elements{25,:},:),before.FEM_nodes.previous_position(:,1),...
    before.FEM_nodes.previous_position(:,2),before.FEM_nodes.concentration)
axis(axis_values_FEM)
grid off
view(0,10)

% FEM_element_nodes = after.FEM_elements.nodes(after.FEM_elements.nodes(:,1)>0,:);
FEM_element_nodes = [after.FEM_elements.nodes(after.cells.FEM_elements{25,:},:);
    after.FEM_elements.nodes(after.cells.FEM_elements{no_cells+1,:},:)];

figure('position',[100 100 300 300],'color','white')
axes('position',[0 0 1 1])
trisurf(FEM_element_nodes,after.FEM_nodes.previous_position(:,1),...
    after.FEM_nodes.previous_position(:,2),after.FEM_nodes.concentration)
axis(axis_values_FEM)
grid off
view(0,10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long

disp(['Total concentration in cell before division = ',num2str(CalculateTotalDpp(...
    before.FEM_nodes.concentration,before.FEM_elements.nodes(...
    before.cells.FEM_elements{25},:),before.FEM_nodes.previous_position))]);

disp(['Total concentration in cells after division = ',num2str(CalculateTotalDpp(...
    after.FEM_nodes.concentration,after.FEM_elements.nodes(...
    after.cells.FEM_elements{25},:),after.FEM_nodes.previous_position)+...
    CalculateTotalDpp(after.FEM_nodes.concentration,after.FEM_elements.nodes(...
    after.cells.FEM_elements{52},:),after.FEM_nodes.previous_position))]);