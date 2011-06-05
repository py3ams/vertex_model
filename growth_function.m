function G = growth_function(cell_area,cell_centre,target_area,time,...
	time_of_last_division)
	
growth_speed = 10*exp(-abs(cell_centre(1)));

G = growth_speed*(time-time_of_last_division)*cell_area/target_area;