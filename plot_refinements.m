disp('busy')

clear FEM_elements FEM_nodes cells.FEM_elements

no_refinements = 3;

[FEM_elements,FEM_nodes,cells.FEM_elements] =...
    create_FEM_mesh(cells.vertices,vertices.position,no_refinements);

FEM_element_nodes = FEM_elements.nodes(FEM_elements.nodes(:,1)>0,:);

figure('position',[100 100 800 800])
trisurf(FEM_element_nodes,FEM_nodes.position(:,1),FEM_nodes.position(:,2),...
    zeros(length(FEM_nodes.position(:,1)),1),'linewidth',2);
axis off;grid off;view([0 90]);axis equal
hold on;
for current_cell = 1:length(cells.vertices)
    hold on
    patchAS(vertices.position(cells.vertices{current_cell},:),[50,200,50]/256,5)
end
axis([-0.2254 -0.0688 -0.3391 -0.1744])