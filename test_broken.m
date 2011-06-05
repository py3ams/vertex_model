function test_broken(cell_store,cells,Dpp,node_positions,iteration,...
	function_name,regular_tests_logical)

global cell_figure_handle

array_sizes = size(node_positions,1);

if regular_tests_logical
	
	% display error message if any cells have lenghts less than 3 but greater than
	% 0, or if any node positions are NaN.
	if any(abs(cellfun('length',cells)-1.5)<1)
		
		setWindowState(cell_figure_handle,'minimize')
		error(['**** Cells have disappeared in ',function_name,' at iteration ',...
			int2str(iteration),' ****'])
		
	elseif any(any(isnan(node_positions)))
		
		setWindowState(cell_figure_handle,'minimize')
		error(['**** node_position NaNs in ',function_name,' at iteration ',...
			int2str(iteration),' ****'])
		
	elseif any(any(isnan(Dpp)))
		
		setWindowState(cell_figure_handle,'minimize')
		error(['**** Dpp NaNs in ',function_name,' at iteration ',...
			int2str(iteration),' ****'])
		
	else
		
		[dummy_cell_store] =...
			CreateCellStore(cells,array_sizes);
		
		logical_equal = sum(cell_store,2)==sum(dummy_cell_store,2);
		
		if ~all(logical_equal)
			error(['**** cell_store broken in ',function_name,' at iteration ',...
				int2str(iteration),' ****'])
		end
		
	end
	
end