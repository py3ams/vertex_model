function nans_logical =...
	test_for_nans(Dpp,node_positions,iteration,regular_tests_logical)

nans_logical = false;

if regular_tests_logical
	
	if any(any(isnan(node_positions)))
		
		disp(['**** node_position NaNs in at iteration ',int2str(iteration),' ****'])
		nans_logical = true;
		
	elseif any(any(isnan(Dpp)))
		
		disp(['**** Dpp NaNs in at iteration ',int2str(iteration),' ****'])
		nans_logical = true;
		
	end
	
end