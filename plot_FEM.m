function plot_FEM(cells,FEM_elements,FEM_nodes,axis_values_FEM,chemical_to_view,...
	linewidth_cells,linewidth_elements,vertices,view_FEM_concentration)

set_caxis_logical = true;
set_colormap_logical = true;

if view_FEM_concentration
	FEM_angle = [0 60];
else
	FEM_angle = [0 90];
end
shading_style = 'faceted';

caxis_vals = [0 axis_values_FEM(6)];
green1 = [50,180,50]/255;
white = [255,255,0]/255;
colormap_val = [linspace(green1(1),white(1),300)' ...
    linspace(green1(2),white(2),300)' linspace(green1(3),white(3),300)'];

FEM_element_nodes = FEM_elements.nodes(FEM_elements.nodes(:,1)>0,:);

if view_FEM_concentration
	
   trisurf(FEM_element_nodes,FEM_nodes.previous_position(:,1),FEM_nodes.previous_position(:,2),...
		FEM_nodes.concentration(:,chemical_to_view),'linewidth',linewidth_elements)
	% 			coloraxis = caxis;
	% 			coloraxis = [0 0.1];
	axis off;grid off;axis(axis_values_FEM);view(FEM_angle);shading(shading_style);
   
   if set_caxis_logical
      caxis(caxis_vals);
   end
   
   if set_colormap_logical
      colormap(colormap_val);
   end
      
else
	
	trisurf(FEM_element_nodes,FEM_nodes.previous_position(:,1),FEM_nodes.previous_position(:,2),...
		zeros(length(FEM_nodes.previous_position(:,1)),1),'linewidth',linewidth_elements)
	axis off;grid off;axis(axis_values_FEM);view(FEM_angle);shading(shading_style);
      
   if set_caxis_logical
      caxis(caxis_vals);
   end
   
   if set_colormap_logical
      colormap(colormap_val);
   end
   
	hold on;
	for current_cell=1:length(cells.vertices)
		hold on
		patchAS(vertices.position(cells.vertices{current_cell},:),[50,200,50]/256,linewidth_cells)
	end
	
end

drawnow;