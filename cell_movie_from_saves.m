disp('busy');clear all;close all;

% this is better than using movie_name or simulation_name as they get overwritten
% when we load up the saved sims!
folder_name = 'area_force_only';
saved_iterations = [1:1:1000];
% saved_iterations = [1000 10000];

save_movie_logical = true;
axis_size = 1.1;

linewidth = 5;
frame_counter = 1;

if save_movie_logical
    if exist(['Movies/',folder_name,'.mpg'],'file')
        yesno = input(['Are you sure you want to overwrite file ',...
            folder_name,'.mpg? '],'s');
        if strcmp(yesno,'yes') || strcmp(yesno,'y')
            delete(['Movies/',folder_name,'.mpg']);
        else
            error(['Could not delete movie ',folder_name,'.mpg'])
        end
    end
    system(['touch Movies/',folder_name,'.mpg']);
end

load(['Saves/',folder_name,'/initial_save.mat'])

figure('position',[100 100 500 500],'PaperPositionMode','auto','color','white')
axes('position',[0 0 1 1])
figure_loop(cells,vertices,linewidth);
% only at view(90,0) does ellipsoid actually fill the axes due to matlab's rendering, so 
% we make the y axis a bit more than twice the x
axis(axis_size*[-0.5 0.5 -0.5 0.5]);
axis off;
drawnow;

if save_movie_logical
    
    M(frame_counter) = getframe(gcf);
    frame_counter = frame_counter+1;
    
end

for unused_variable = 1:length(saved_iterations)
    
    if ~exist(['Saves/',folder_name,'/iteration_',num2str(saved_iterations(unused_variable)),'.mat'],'file')
        error(['No save for iteration_',num2str(saved_iterations(unused_variable)),'.mat'])
    end
    
    load(['Saves/',folder_name,'/iteration_',num2str(saved_iterations(unused_variable)),'.mat']);
    
    cla;
    figure_loop(cells,vertices,linewidth);
    axis(axis_size*[-0.5 0.5 -0.5 0.5]);
    axis off;
    drawnow;
    
    if save_movie_logical
        
        M(frame_counter) = getframe(gcf);
        frame_counter = frame_counter+1;
        
        if frame_counter > 100 || unused_variable == length(saved_iterations)
            mpgwrite(M,jet,[movie_location,'update.mpg']);
            system(['cat ',movie_location,folder_name,'.mpg ',movie_location,...
                'update.mpg > ',movie_location,'temp.mpg']);
            system(['mv ',movie_location,'temp.mpg ',movie_location,...
                folder_name,'.mpg']);
            clear M
            frame_counter = 1;
        end
        
    end
    
end