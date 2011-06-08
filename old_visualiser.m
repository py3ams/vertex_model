function old_visualiser(cells,vertices,FEM_elements,FEM_nodes,axis_values,...
	axis_values_FEM,figure_position,final_iteration_logical,iteration,...
	movie_logical,update_period_1,update_period_2,view_FEM_concentration,...
	view_FEM_mesh,view_iteration_number,visualise_initial_configuration)

% this code is really confused now! need to go through and sort out what i mean by
% left_axes_1, right_axes_2 etc, particularly in cases of more than one chemical
% being present.

chemicals_to_view = [1];

if view_FEM_concentration
	FEM_angle = [0 50];
else
	FEM_angle = [0 90];
	%     axis_values_FEM = 'equal';
end
linewidth = 3;
shading_style = 'faceted';

global cell_figure_handle right_text_axes left_axes_1 left_axes_2 ...
	right_axes_1 right_axes_2 left_text_axes

FEM_element_nodes = FEM_elements.nodes(FEM_elements.nodes(:,1)>0,:);

% find the number of cells
no_cells = sum(~cellfun('isempty',cells.vertices));


% set everything up during the first iteration
if ~iteration
	
	cell_figure_handle = figure('outerposition',figure_position,'color','white');
	
	if visualise_initial_configuration
		
		left_text_axes = axes('position',[0 0.9 0.5 0.05]);
		axis([0 1 0 1])
		axis off
		
		if view_iteration_number
			figure_text = [' Iteration ',num2str(iteration),'; Number of cells = ',...
				int2str(no_cells),' '];
		else
			figure_text = ['Number of cells = ',int2str(no_cells),' '];
		end
		
		%         figure_text = [' Number of cells = ',...
		%             int2str(no_cells),' '];
		text(0.5,0.7,figure_text,'FontName','FixedWidth','FontSize',14,...
			'FontWeight','bold','HorizontalAlignment','center','EdgeColor',...
			'black','LineWidth',2)
		
		right_text_axes = axes('position',[0.5 0.9 0.5 0.05]);
		axis([0 1 0 1])
		axis off
		
		if view_iteration_number
			figure_text = [' Iteration ',num2str(iteration),'; Number of cells = ',...
				int2str(no_cells),' '];
		else
			figure_text = [' Number of cells = ',int2str(no_cells),' '];
		end
		%         figure_text = [' Number of cells = ',...
		%             int2str(no_cells),' '];
		text(0.5,0.7,figure_text,'FontName','FixedWidth','FontSize',14,...
			'FontWeight','bold','HorizontalAlignment','center','EdgeColor',...
			'black','LineWidth',2)
		
		if view_FEM_mesh
			
			left_axes_1 = axes('position',[0.01 0.45 0.48 0.45]);
			figure_loop(cells,vertices,linewidth)
			axis(axis_values)
			axis off;
			
			left_axes_2 = axes('position',[0.01 0 0.48 0.45]);
			figure_loop(cells,vertices,linewidth)
			axis(axis_values)
			axis off;
			
			right_axes_1 = axes('position',[0.51 0.45 0.48 0.44]);
			trisurf(FEM_element_nodes,FEM_nodes.previous_position(:,1),FEM_nodes.previous_position(:,2),FEM_nodes.concentration(:,1))
			axis off;grid off;axis(axis_values_FEM);view(FEM_angle);shading(shading_style);
			
			right_axes_2 = axes('position',[0.51 0.01 0.48 0.44]);
			trisurf(FEM_element_nodes,FEM_nodes.previous_position(:,1),FEM_nodes.previous_position(:,2),FEM_nodes.concentration(:,1))
			axis off;grid off;axis(axis_values_FEM);view(FEM_angle);shading(shading_style);
			
			%           right_axes_3 = axes('position',[0.75 0 0.24 0.45]);
			% 			trisurf(FEM_elements,FEM_node_positions(:,1),FEM_node_positions(:,2),Brk)
			% 			grid off;
			% 			shading(shading_style)
			% 			axis(axis_values_FEM);
			% 			axis off;
			% 			view(FEM_angle);
			
		else
			
			left_axes_1 = axes('position',[0.01 0 0.48 0.9]);
			figure_loop(cells,vertices,linewidth)
			axis(axis_values)
			axis off;
			
			left_axes_2 = axes('position',[0.51 0 0.48 0.9]);
			figure_loop(cells,vertices,linewidth)
			axis(axis_values);
			axis off;
			
		end
		
	else
		
		right_text_axes = axes('position',[0 0.9 1 0.09]);
		axis([0 1 0 1])
		axis off
		
		if view_iteration_number
			figure_text = [' Iteration ',num2str(iteration),'; Number of cells = ',...
				int2str(no_cells),' '];
		else
			figure_text = [' Number of cells = ',int2str(no_cells),' '];
		end
		%         figure_text = [' Number of cells = ',...
		%             int2str(no_cells),' '];
		text(0.5,0.7,figure_text,'FontName','FixedWidth','FontSize',14,...
			'FontWeight','bold','HorizontalAlignment','center','EdgeColor',...
			'black','LineWidth',2)
		
		if view_FEM_mesh
			
			% to make this axes square have to remember that figure itself
			% is in 5/4 ratio (horizontal/vertices), so have to make
			% vertices axis 5/4 bigger than the horizontal axis in relative
			% terms
			left_axes_1 = axes('position',[0.01 0.5*(0.9-1.25*0.49) 0.49 1.25*0.49]);
			figure_loop(cells,vertices,linewidth)
			axis(axis_values)
			axis off;
			
			if ~view_FEM_concentration
				
				right_axes_1 = axes('position',[0.51 0.5*(0.9-1.25*0.49) 0.49 1.25*0.49]);
				trisurf(FEM_element_nodes,FEM_nodes.previous_position(:,1),...
					FEM_nodes.previous_position(:,2),zeros(length(FEM_nodes.previous_position),1),'linewidth',1.5)
				hold on;
				cellfun(@(x)patch(vertices.position(x,1),vertices.position(x,2),[255,30,30]/255,'linewidth',5,'FaceAlpha',0),cells.vertices);
				axis off;grid off;axis(axis_values_FEM);view(FEM_angle);
				
			else
				
				if length(chemicals_to_view) == 1
					
					right_axes_1 = axes('position',[0.51 0.5*(0.9-1.25*0.49) 0.49 1.25*0.49]);
					trisurf(FEM_element_nodes,FEM_nodes.previous_position(:,1),FEM_nodes.previous_position(:,2),...
						FEM_nodes.concentration(:,chemicals_to_view(1)));
					axis off;grid off;axis(axis_values_FEM);view(FEM_angle);shading(shading_style);
					
				elseif length(chemicals_to_view) == 2
					
					right_axes_1 = axes('position',[0.51 0.01 0.49 0.44]);
					trisurf(FEM_element_nodes,FEM_nodes.previous_position(:,1),FEM_nodes.previous_position(:,2),...
						FEM_nodes.concentration(:,chemicals_to_view(1)));
					axis off;grid off;axis(axis_values_FEM);view(FEM_angle);shading(shading_style);
					
					right_axes_2 = axes('position',[0.51 0.45 0.49 0.44]);
					trisurf(FEM_element_nodes,FEM_nodes.previous_position(:,1),FEM_nodes.previous_position(:,2),...
						FEM_nodes.concentration(:,chemicals_to_view(2)));
					axis off;grid off;axis(axis_values_FEM);view(FEM_angle);shading(shading_style);
					
				else
					error('Not enough/Too many chemicals to view')
				end
			end
			
			%             right_axes_3 = axes('position',[0.51 0 0.48 0.49]);
			%             trisurf(FEM_elements,FEM_node_positions(:,1),FEM_node_positions(:,2),Brk)
			%             grid off;
			%             shading(shading_style)
			%             axis(axis_values_FEM);
			%             axis off;
			%             view(FEM_angle);
			
			%             axes('position',[0.17 0.42 0.01 0.01]);
			%             axis([0 1 0 1]);
			%             axis off
			%             set(gca,'FontSize',20);
			%             text(0.5,0.5,'Dpp','FontName','FixedWidth','FontSize',18,...
			%                 'FontWeight','bold','HorizontalAlignment','center')
			
			%             axes('position',[0.87 0.42 0.01 0.01]);
			%             axis([0 1 0 1]);
			%             axis off
			%             %             set(gca,'FontSize',20);
			%             text(0.5,0.5,'Brk','FontName','FixedWidth','FontSize',18,...
			%                 'FontWeight','bold','HorizontalAlignment','center')
			
		else
			
			left_axes_1 = axes('position',[0.01 0 0.98 0.9]);
			figure_loop(cells,vertices,linewidth)
			axis(axis_values)
			axis off;
			
		end
		
	end
	
else
	
	if movie_logical && visualise_initial_configuration && ~rem(iteration,update_period_2)
		
		figure(cell_figure_handle);
		
		if view_iteration_number
			figure_text = [' Iteration ',num2str(iteration),'; Number of cells = ',...
				int2str(no_cells),' '];
		else
			figure_text = [' Number of cells = ',int2str(no_cells),' '];
		end
		%         figure_text = [' Number of cells = ',...
		%             int2str(no_cells),' '];
		
		set(cell_figure_handle,'CurrentAxes',left_text_axes);cla;
		text(0.5,0.7,figure_text,'FontName','FixedWidth','FontSize',14,...
			'FontWeight','bold','HorizontalAlignment','center','EdgeColor',...
			'black','LineWidth',2)
		
		set(cell_figure_handle,'CurrentAxes',left_axes_2);cla;
		figure_loop(cells,vertices,linewidth)
		axis(axis_values)
		
		if view_FEM_mesh
			
			set(cell_figure_handle,'CurrentAxes',right_axes_2);cla;
			
			if ~view_FEM_concentration
				trisurf(FEM_element_nodes,FEM_nodes.previous_position(:,1),...
					FEM_nodes.previous_position(:,2),zeros(length(FEM_nodes.previous_position),1),'linewidth',1.5)
				hold on;
				cellfun(@(x)patch(vertices.position(x,1),vertices.position(x,2),[255,30,30]/255,'linewidth',5,'FaceAlpha',0),cells.vertices);
				axis off;grid off;axis(axis_values_FEM);view(FEM_angle);
			else
				trisurf(FEM_element_nodes,FEM_nodes.previous_position(:,1),FEM_nodes.previous_position(:,2),FEM_nodes.concentration(:,1));
				axis off;grid off;axis(axis_values_FEM);view(FEM_angle);shading(shading_style);
			end
			
			
		end
		
		drawnow;
		
	end
	
	if (movie_logical && ~rem(iteration,update_period_1)) || final_iteration_logical
		
		figure(cell_figure_handle);
		if view_iteration_number
			figure_text = [' Iteration ',num2str(iteration),'; Number of cells = ',...
				int2str(no_cells),' '];
		else
			figure_text = [' Number of cells = ',int2str(no_cells),' '];
		end
		%         figure_text = [' Number of cells = ',...
		%             int2str(no_cells),' '];
		
		set(cell_figure_handle,'CurrentAxes',right_text_axes);cla;
		text(0.5,0.7,figure_text,'FontName','FixedWidth','FontSize',14,...
			'FontWeight','bold','HorizontalAlignment','center','EdgeColor',...
			'black','LineWidth',2)
		
		set(cell_figure_handle,'CurrentAxes',left_axes_1);cla;
		figure_loop(cells,vertices,linewidth)
		axis(axis_values)
		
		if view_FEM_mesh
			
			if ~view_FEM_concentration
				set(cell_figure_handle,'CurrentAxes',right_axes_1);cla;
				trisurf(FEM_element_nodes,FEM_nodes.previous_position(:,1),...
					FEM_nodes.previous_position(:,2),zeros(length(FEM_nodes.previous_position),1),'linewidth',1.5)
				hold on;
				cellfun(@(x)patch(vertices.position(x,1),vertices.position(x,2),[255,30,30]/255,'linewidth',5,'FaceAlpha',0),cells.vertices);
				
			else
				set(cell_figure_handle,'CurrentAxes',right_axes_1);cla;
				trisurf(FEM_element_nodes,FEM_nodes.previous_position(:,1),FEM_nodes.previous_position(:,2),...
					FEM_nodes.concentration(:,chemicals_to_view(1)));
				axis off;grid off;view(FEM_angle);axis(axis_values_FEM);
				if length(chemicals_to_view) == 2
					set(cell_figure_handle,'CurrentAxes',right_axes_2);cla;
					trisurf(FEM_element_nodes,FEM_nodes.previous_position(:,1),FEM_nodes.previous_position(:,2),...
						FEM_nodes.concentration(:,chemicals_to_view(2)));
					axis off;grid off;view(FEM_angle);axis(axis_values_FEM);
				end
			end
			
			
			%             set(cell_figure_handle,'CurrentAxes',right_axes_3);cla;
			%             trisurf(FEM_elements,FEM_node_positions(:,1),FEM_node_positions(:,2),Brk);
			%             grid off;
			%             shading(shading_style);
			%             axis(axis_values_FEM);
			%             axis off;
			%             view(FEM_angle);
			%
		end
		
		drawnow;
		
	end
	
end

