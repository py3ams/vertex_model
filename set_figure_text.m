function set_figure_text(cells,iteration,view_iteration_number,view_number_cells)

no_cells = sum(~cellfun('isempty',cells.vertices));

figure_text = [];
if view_iteration_number
	figure_text = [' Iteration ',num2str(iteration),' '];
end

if view_number_cells
	figure_text = [figure_text,' Number of cells = ',int2str(no_cells),' '];
end

text(0.5,0.7,figure_text,'FontName','FixedWidth','FontSize',14,...
	'FontWeight','bold','HorizontalAlignment','center','EdgeColor',...
	'black','LineWidth',2)